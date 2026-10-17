#' Align two or more cell datasets
#'
#' Given two or more cell feature datasets, the function will try to align the spatial coordinates of samples shared within each dataset.
#'
#'
#' @param N_cores Integer. Number of cores to parallelize your computation.
#' @param Cell_feature_list A named list containing the datasets to be aligned. Each dataset must have been adequately formatted using the [Data_arrange_function()].
#' @param Iterations The number of iterations of the ICP algorithm (see details) will run before stopping.
#' @param Phenotypes_included (OPTIONAL) If the datasets contain information regarding cell phenotype labels, a character vector indicating the cell phenotypes that will be used in the aligning process.
#' @param Apply_DBSCAN_denoise A logical value indicating if density-based filtering should be performed before alignment.
#' @inheritParams Clustering_Phenotyper
#' @param Concave_Hull A logical value indicating if alignment should be performed on the bare points, or based on the concave hull computed from the points.
#' @inheritParams Image_size_calculator
#' @param Polygon_sampling_density A numeric value indicating the sampling density of the concave hull polygon. Higher values result in lower density sampling.
#' @param Perform_final_ICP_refinement A logical value indicating if final ICP refinement taking into account all the points in the provided datasets should be performed.
#'
#' @details
#' The [Cell_aligner_function()] assumes that there is spatial correlation between the samples provided. The function also assumes that spatial patterns contained within each
#' sample can be used as landmarks to enable spatial alignment.
#'
#' The function will first identify samples that are shared between the datasets provided according to the Subject_Names. The samples present in the first item of the list
#' will be used as spatial reference. This means that the spatial X Y coordinates of the samples present in the second list item will be aligned to match the sample in the first
#' list item. Afterwards, the samples of the third list item will be matched with the one in the second list item and so on. Any missing sample will be skipped by the algorithm.
#'
#' The user can reduce the spatial noise and computational burden by selecting only certain cell phenotype labels in the alignment process, by performing DBSCAN filtering and by
#' aligning the concave hull polygon in stead of the whole point pattern.
#'
#' The algorithm implemented in the [Cell_aligner_function()] will first align the centroids of the reference and moving point datasets. Then the Morpho::icpmat function from the
#' Morpho package is executed to perform fine alignment. Afterwards, image flipping (either X or Y) followed by rotation (1 to 360 degrees) is performed, re-aligned using the Morpho::icpmat function
#' and model error tested by computing average distance to closest neighbor to the reference points. The alignment with the smallest error is selected and the alignment model is recorded to be
#' used to align the whole dataset if required.
#'
#'
#' @returns A list containing two elements. The first is the Cell_feature_list containing the aligned X Y coordinates, the second is a tibble summarizing the alignment
#' process for every sample included and the error of the model.
#'
#'
#' @examples
#' \dontrun{
#'  #Get the original dataset
#' Original_dataset <- CSM::CSM_Phenotypecell_test
#'
#'#Generate a randomly rotated version
#'Rotated_dataset <-
#'  purrr::map_dfr(unique(Original_dataset$Subject_Names), function(Subject){
#'    #obtain each individual sample
#'  Individual_sample <- Original_dataset %>% dplyr::filter(Subject_Names == Subject)
#'  #Compute a random degree
#'  Random_degrees <- sample(15:345, size = 1)
#'  Random_radians <- Random_degrees * pi / 180
#'
#'  Individual_sample[c("X", "Y")] <- Rotate_points(as.matrix(Individual_sample[c("X", "Y")]), theta = Random_radians, non_negative = TRUE)
#'
#'  return(Individual_sample)
#'})
#'#Generate a randomly flipped and translated
#'Flipped_Translated_dataset <-
#'purrr::map_dfr(unique(Original_dataset$Subject_Names), function(Subject){
#'  #obtain each individual sample
#'  Individual_sample <- Original_dataset %>% dplyr::filter(Subject_Names == Subject)
#'  #Compute a random flipping
#'  Random_flipping <- sample(c("x", "y"), size = 1)
#'  Individual_sample[c("X", "Y")] <- Flip_xy(as.matrix(Individual_sample[c("X", "Y")]), mode = Random_flipping)
#'
#'  #Compute a random translation
#'  Random_X_translation <- sample(seq(from = 25, to = 100, by = 25), size = 1)
#'  Random_Y_translation <- sample(seq(from = 25, to = 100, by = 25), size = 1)
#'  Individual_sample$X <- Individual_sample$X + Random_X_translation
#'  Individual_sample$Y <- Individual_sample$Y + Random_Y_translation
#'
#'  return(Individual_sample)
#'})
#'
#'#Place it al into a named list
#'DATASETS_LIST <-
#' list(Reference = Original_dataset,
#'      Rotated = Rotated_dataset,
#'      Flipped_translated = Flipped_Translated_dataset)
#'
#'
#' #Run the alignment protocol
#'Aligned <-
#' Cell_aligner_function(N_cores = 1,
#'                        Cell_feature_list = DATASETS_LIST,
#'                        Iterations = 5,
#'                        Phenotypes_included = NULL,
#'
#'                        Apply_DBSCAN_denoise = FALSE,
#'
#'                        Concave_Hull = TRUE,
#'                        Hull_ratio = 0.1,
#'                        Polygon_sampling_density = 5,
#'
#'                        Perform_final_ICP_refinement = FALSE
#'                        )
#'
#'
#'
#' #Plot the results
#'plot(
#'   purrr::map_dfr(1:length(TEST$Cell_feature_list), function(List_item){
#'   Item <- TEST$Cell_feature_list[[List_item]] %>% dplyr::mutate(Modification = names(TEST$Cell_feature_list)[List_item])
#' }) %>%
#'   ggplot(aes(X, Y, color = Modification)) +
#'   facet_wrap(~Subject_Names) +
#' geom_point(size = 3)
#' )
#'}
#'
#' @import patchwork
#'
#' @export

Cell_aligner_function <-
  function(
    N_cores = 1,
    Cell_feature_list,
    Iterations = 5,

    Phenotypes_included = NULL,

    Apply_DBSCAN_denoise = FALSE,
    Min_cell_no = NULL,
    Distance_radius = NULL,

    Concave_Hull = FALSE,
    Hull_ratio = 0.3,
    Polygon_sampling_density = 10,

    Perform_final_ICP_refinement = FALSE
  ){

    #######Check suggested packages#######
    {
      if(!requireNamespace("Rvcg", quietly = FALSE)) stop(
        paste0("Rvcg CRAN package is required to execute the function. Please install using the following code: ",
               expression(install.packages("Rvcg")))
      )

      if(!requireNamespace("Morpho", quietly = FALSE)) stop(
        paste0("Morpho CRAN package is required to execute the function. Please install using the following code: ",
               expression(install.packages("Morpho")))
      )

      #If concave hull required check sf
      if(Concave_Hull){
        if(!requireNamespace("sf", quietly = FALSE)) stop(
          paste0("sf CRAN package is required to execute the function. Please install using the following code: ",
                 expression(install.packages("sf")))
        )
      }
      #If DBSCAN filtering is required
      if(Apply_DBSCAN_denoise){
        if(!requireNamespace("dbscan", quietly = FALSE)) stop(
          paste0("dbscan CRAN package is required to execute the function. Please install using the following code: ",
                 expression(install.packages("dbscan")))
        )
      }
    }

    #######Argument check general arguments#######
    #What to do on exit
    on.exit({
      future::plan("future::sequential")
      gc()
      })

    #Check cores
    if(!all(N_cores >= 1 & N_cores%%1 == 0)) stop("N_cores must be an integer value > 0")

    #Check iterations
    if(!all(Iterations >= 1 & Iterations%%1 == 0)) stop("Iterations must be an integer value > 0")

    #Check Perform_final_ICP_refinement
    if(!all(is.logical(Perform_final_ICP_refinement), length(Perform_final_ICP_refinement))) stop("Perform_final_ICP_refinement must be a single logical value")

    #Check hull-related arguments
    if(!is.logical(Concave_Hull)) stop("Concave_Hull must be a logical value")
    if(Concave_Hull){
      if(!all(is.numeric(Hull_ratio), Hull_ratio >= 0, Hull_ratio <= 1)) stop("Hull_ratio must be a numeric value between 0 and 1")
      if(!all(is.numeric(Polygon_sampling_density), Polygon_sampling_density >= 1)) stop("Polygon_sampling_density must be a numeric value >= 1")
    }

    #Check DBSCAN related arguments
    if(Apply_DBSCAN_denoise){
      #Check other arguments
      if(!all(is.numeric(Min_cell_no), Min_cell_no%%1 == 0, Min_cell_no > 1)) stop("Min_cell_no must be an integer value > 1")
      if(!all(is.numeric(Distance_radius), Distance_radius > 0)) stop("Distance_radius must be a numeric value > 0")
    }

    #Check that the Cell_feature_list is a list, with non repeated names
    if(!is.list(Cell_feature_list)) stop("Cell_feature_list must be a named list")
    if(length(Cell_feature_list) < 2) stop("Cell_feature_list must be of length >= 2")
    if(is.null(names(Cell_feature_list)) | "" %in% names(Cell_feature_list)) stop("Cell_feature_list must be a named list")
    if(any(duplicated(names(Cell_feature_list)))) stop("Names of Cell_feature_list must be unique")

    #Check that all the cell_feature list contains tibbles with Cell_no, X, Y, Subject_Names
    purrr::walk(1:length(Cell_feature_list), function(List_item){
      if(!identical(names(Cell_feature_list[[List_item]])[1:4],  c("Cell_no", "X", "Y", "Subject_Names"))) { #Check if Data is correctly formatted
        stop(paste0("Issue with list item n", List_item, " ", names(Cell_feature_list)[List_item], ":",
                    "\n",
                    "Data is not adequately formatted."))
      }
    })

    #Check phenotypes
    if(!is.null(Phenotypes_included)){
      #Check that phenotypes are contained within each item provided if necessary
      purrr::walk(1:length(Cell_feature_list), function(List_item){
        if(!"Phenotype" %in% names(Cell_feature_list[[List_item]])) stop(paste0("Issue with list item n", List_item, " ", names(Cell_feature_list)[List_item], ":",
                                                                                    "\n",
                                                                                    "A column named Phenotype containing cell phenotype labels has not been identified."))

        if(sum(Cell_feature_list[[List_item]]$Phenotype %in% Phenotypes_included) < 1) stop(paste0("Issue with list item n", List_item, " ", names(Cell_feature_list)[List_item], ":",
                                                                                                   "\n",
                                                                                                   "No ", stringr::str_c(Phenotypes_included, collapse = ", "), " cell types present in the sample."))
      })

      #Remove samples without phenotype cells
      Final_Cell_feature_list <- purrr::map(Cell_feature_list, function(List_item){
        return(List_item %>% dplyr::filter(Phenotype %in% Phenotypes_included))
      })
      names(Final_Cell_feature_list) <- names(Cell_feature_list)

      #Compute for every list item how many samples have been removed, If any removed print a warning message
      purrr::walk(1:length(Final_Cell_feature_list), function(List_item){
        Final_subjects <- unique(Final_Cell_feature_list[[List_item]]$Subject_Names)
        Original_subjects <- unique(Cell_feature_list[[List_item]]$Subject_Names)

        if(!all(Original_subjects %in% Final_subjects)){
          Removed_subjects <- Original_subjects[!Original_subjects %in% Final_subjects]

          cat(paste0("\n",
                     "The following subjects have been removed from ", names(Final_Cell_feature_list)[List_item], ". ",
                     "No ", stringr::str_c(Phenotypes_included, collapse = ", "), " cell types present in sample:" )
              )
          Colored_print(paste0("\n",
                               str_c(Removed_subjects, collapse = "\n")),
                        color = "red")

        }
      })
    }

    else{
      Final_Cell_feature_list <- Cell_feature_list
    }

    #Paste a message indicating the order of the alignment
    cat(paste("\nThe following order will be attempted to perform alignment:\n"))
    Main_ref_name <- stringr::str_c(names(Final_Cell_feature_list)[1], " (Main spatial reference)")
    Alignment_order <- stringr::str_c(c(Main_ref_name, names(Final_Cell_feature_list)[-1]), collapse = "\n\u21D1\n")
    Centered_print(Alignment_order)

    ####BUILD THE LOOK UP TABLE####
    #Find the dataset within the list with the largest amount of samples
    List_index <- which.max(purrr::map_dbl(Final_Cell_feature_list, ~length(unique(.$Subject_Names))))
    Reference_sample_name <- unique(Final_Cell_feature_list[[List_index]]$Subject_Names)

    #Iterate across the Final_Cell_feature_list to find the closest match
    Look_up_tibble <-
      purrr::map_dfc(1:length(Final_Cell_feature_list), function(List_item_index){
        List_item <- Final_Cell_feature_list[[List_item_index]]
        Names_in_item <- unique(List_item$Subject_Names)

        #Generate the vector
        Closest_vector <- vector(mode = "character", length = length(Reference_sample_name))
        Closest_value <- vector(mode = "integer", length = length(Reference_sample_name))

        #For every item in reference_sample_name find the closest match
        for(i in 1:length(Reference_sample_name)){
          Target_name <- Reference_sample_name[i]
          Distance_vector <- adist(Names_in_item, Target_name, fixed = TRUE)

          Closest_value[i] <- min(Distance_vector)
          Closest_vector[i] <- Names_in_item[which.min(Distance_vector)]
        }

        Final_tibble <- tibble::tibble(Target = Reference_sample_name,  Name = Closest_vector, Dist = Closest_value)
        Final_tibble <- Final_tibble %>%
          dplyr::group_by(Name) %>%
          dplyr::mutate(
            Name = dplyr::if_else(
              dplyr::row_number() == which.min(Dist),
              Name,
              NA_character_
            )
          ) %>%
          dplyr::ungroup()

        #Get only the names required
        Final_tibble <- Final_tibble["Name"]
        names(Final_tibble) <- names(Final_Cell_feature_list)[List_item_index]
        return(Final_tibble)
      })

    #Create a list with the samples
    Samples_list <-
      purrr::map(names(Look_up_tibble), function(List_item_name){
        Item <- Final_Cell_feature_list[[List_item_name]]
        Samples <- na.omit(Look_up_tibble[[List_item_name]])
        Final_list <- purrr::map(Samples, function(sample) Item %>% dplyr::filter(Subject_Names == sample) %>% dplyr::select(1:4))
        names(Final_list) <- Samples
        return(Final_list)
      })
    names(Samples_list) <- names(Look_up_tibble)

    #####DBSCAN FILTERING#####
    if(Apply_DBSCAN_denoise){
      cat("\n Running a test of DBSCAN denoising")

      Test_completed <- "NO" #prepare to enter the loop
      while(Test_completed == "NO"){
        #Choose a random row from the look up tibble
        Random_row <- sample(1:nrow(Look_up_tibble), size = 1)
        #Choose the random column to initate the alignment process
        Random_column <- sample(1:ncol(Look_up_tibble), size = 1)
        #Get the dataset info
        Dataset_list_name <- names(Look_up_tibble)[Random_column]
        Dataset_sample_name <- Look_up_tibble[[Random_row, Random_column]]
        #If is NA then repeat the process
        if(is.na(Dataset_sample_name)) next

        #Get the dataset
        Dataset <- Samples_list[[Dataset_list_name]][[Dataset_sample_name]]

        #Proceed with algorithm
        DB_results <- dbscan::dbscan(Dataset[c("X", "Y")], eps = Distance_radius, minPts = Min_cell_no, borderPoints = FALSE)

        #Plot the result
        if(length(DB_results$cluster) <= 100000){
          plot(
            tibble::tibble(X = Dataset[["X"]], Y = Dataset[["Y"]], Cluster = DB_results$cluster) %>%
              dplyr::mutate(Type = case_when(Cluster == 0 ~ "Noise",
                                             TRUE ~ "Approved")) %>%
              ggplot(aes(x = X, y = Y, color = Type)) + geom_point(size = 0.8) +
              scale_color_manual("", values = c(Approved = "black", Noise = "grey")) +
              cowplot::theme_cowplot()+
              ggtitle(paste0(Dataset_list_name, "\n", Dataset_sample_name)) +
              theme(panel.grid = element_blank())
          )
        }

        #Subsample plot for large dataset
        if(length(DB_results$cluster) > 100000){
          message(">100K observations to generate plots. A random subset containing 10% of the dataset will be selected.")
          plot(
            tibble::tibble(X = Dataset[["X"]], Y = Dataset[["Y"]], Cluster = DB_results$cluster) %>%
              dplyr::mutate(Cluster = case_when(Cluster == 0 ~ "Noise",
                                                TRUE ~ "Approved")) %>%
              dplyr::slice_sample(n = 100000) %>%
              ggplot(aes(x = X, y = Y, color = Type)) + geom_point(size = 1.5) +
              scale_color_manual("", values = c(Approved = "black", Noise = "grey"))+
              cowplot::theme_cowplot() +
              ggtitle(paste0(Dataset_list_name, "\n", Dataset_sample_name)) +
              theme(panel.grid = element_blank())
          )
        }

        #Generate a menu to proceed with compuation
        answer <- menu(c("Proceed", "Abort", "Test again"), title = "Check parameters provided and sample image generated. Should the analysis proceed?")
        #If user decides to stop then abort function and return stop message
        if(answer == 2) stop("The function has been stopped. Please tune parameters and try again")
        if(answer == 1) Test_completed <- "Proceed"
      }

      #If proceed then modify the samples list
      cat("\n Running DBSCAN denoising")
      Samples_list <-
        purrr::map(Samples_list, function(List_item){
          #Obtain the sample names list
          Sample_names <- names(List_item)
          #Run the DBSCAN within each specific sample
          Subsample_list <- purrr::map(List_item, function(Sample_item){
            DB_results <- dbscan::dbscan(Sample_item[c("X", "Y")], eps = Distance_radius, minPts = Min_cell_no, borderPoints = FALSE)
            return(Sample_item %>% dplyr::filter(DB_results$cluster != 0))
          })
          #Modify names
          names(Subsample_list) <- Sample_names
          return(Subsample_list)
      })
      names(Samples_list) <- names(Look_up_tibble)
    }

    #Turn this sample list into a shared object
    Samples_list <- mori::share(Samples_list)

    #####RUN AN EXAMPLE AND PROCEED IF EVERYTHING IS OK####
    cat("\n Running a test of alignment parameters")
    Test_completed <- "NO" #prepare to enter the loop
    while(Test_completed == "NO"){
      #Choose a random row from the look up tibble
      Random_row <- sample(1:nrow(Look_up_tibble), size = 1)
      #Choose the random column to initate the alignment process
      Random_column <- sample(1:(ncol(Look_up_tibble)-1), size = 1)

      #Get the alignment dataset info
      Target_dataset_list_name <- names(Look_up_tibble)[Random_column]
      Moving_dataset_list_name <- names(Look_up_tibble)[Random_column + 1]
      Target_dataset_sample_name <- Look_up_tibble[[Random_row, Random_column]]
      Moving_dataset_sample_name <- Look_up_tibble[[Random_row, Random_column + 1]]

      #If any is NA then repeat the process
      if(anyNA(c(Target_dataset_sample_name, Moving_dataset_sample_name))) next

      #Obtain the datasets
      Target_dataset <- Samples_list[[Target_dataset_list_name]][[Target_dataset_sample_name]]
      Moving_dataset <- Samples_list[[Moving_dataset_list_name]][[Moving_dataset_sample_name]]

      #Transform the Target and moving datasets to concave hull sampled polygon lines if required
      if(Concave_Hull){

        #First compute the concave hull of both datasets and sample the concave hull polygon
        Target_dataset <- sf::st_as_sf(Target_dataset, coords = c("X", "Y")) #Generate sf object
        Target_dataset <- sf::st_union(Target_dataset) #Turn it into a geometry object
        Target_dataset <- sf::st_concave_hull(Target_dataset, ratio = Hull_ratio) #Obtain the concave hull
        Target_dataset <- sf::st_boundary(Target_dataset) #Turn it into a line
        Target_dataset <- sf::st_line_sample(Target_dataset, density = 1/Polygon_sampling_density) #Sample the line according to the desired density
        Target_dataset <- sf::st_coordinates(sf::st_cast(Target_dataset, "POINT")) #Obtain the coordinates of the sampled points

        Moving_dataset <- sf::st_as_sf(Moving_dataset, coords = c("X", "Y"))
        Moving_dataset <- sf::st_union(Moving_dataset)
        Moving_dataset <- sf::st_concave_hull(Moving_dataset, ratio = Hull_ratio)
        Moving_dataset <- sf::st_boundary(Moving_dataset)
        Moving_dataset <- sf::st_line_sample(Moving_dataset, density = 1/Polygon_sampling_density)
        Moving_dataset <- sf::st_coordinates(sf::st_cast(Moving_dataset, "POINT"))

      }

      #If not turn it into a matrix
      else{
        Target_dataset <- as.matrix(Target_dataset %>% dplyr::select(X, Y))
        Moving_dataset <- as.matrix(Moving_dataset %>% dplyr::select(X, Y))
      }

      #Run the ICP aligner
      Alignment_results <- ICP_cell_aligner(Moving_coords = Moving_dataset, Target_coords = Target_dataset, Iterations = Iterations)

      #Prepare to plot some results
      X_max <- max(c(Target_dataset[,1], Moving_dataset[,1]))
      Y_max <- max(c(Target_dataset[,2], Moving_dataset[,2]))

      Target_plot <-
        Target_dataset %>%
        ggplot(aes(X, Y)) +
        geom_point(color = "blue", size = 1.5) +
        coord_cartesian(xlim = c(0, X_max), ylim = c(0, Y_max)) +
        theme_minimal()+
        ggtitle("Target", paste0(Target_dataset_list_name, "\n", Target_dataset_sample_name)) +
        theme(panel.grid = element_blank())

      Moving_plot <-
        Moving_dataset %>%
        ggplot(aes(X, Y)) +
        geom_point(color = "red", size = 1) +
        coord_cartesian(xlim = c(0, X_max), ylim = c(0, Y_max)) +
        theme_minimal()+
        ggtitle("Moving", paste0(Moving_dataset_list_name, "\n", Moving_dataset_sample_name)) +
        theme(panel.grid = element_blank(),
              aspect.ratio = Y_max/X_max)

      Corrected_plot <-
        ggplot() +
        geom_point(aes(X, Y), color = "blue", size = 1.5, data = Target_dataset) +
        geom_point(aes(X, Y), color = "red", size = 1, alpha = 0.5, data = Alignment_results$Coordinates)+
        coord_cartesian(xlim = c(0, X_max), ylim = c(0, Y_max)) +
        theme_minimal()+
        ggtitle("Alignment result") +
        theme(panel.grid = element_blank(),
              aspect.ratio = Y_max/X_max)

      Model_results_plot <-
        ggplot() +
        annotate("text", label = paste("Flipping action", Alignment_results$Model$Flip_action, sep = ": "), x = 0, y = 6) +
        annotate("text", label = paste("Rotation (degrees)", round(Alignment_results$Model$Angle_degrees, digits = 2), sep = ": "), x = 0, y = 5) +
        annotate("text", label = paste("X shift", round(Alignment_results$Model$Translation[[1]], digits = 2), sep = ": "), x = 0, y = 4) +
        annotate("text", label = paste("Y shift", round(Alignment_results$Model$Translation[[2]], digits = 2), sep = ": "), x = 0, y = 3) +

        annotate("text", label = paste("Error", round(Alignment_results$Error, digits = 2), sep = ": "), x = 0, y = 1) +

        theme_minimal() +
        theme(panel.grid = element_blank(),
              axis.text = element_blank(),
              axis.title = element_blank())

      plot(Target_plot + Moving_plot + Corrected_plot + Model_results_plot + plot_layout(ncol = 2, nrow = 2))

      #Generate a menu to proceed with compuation
      answer <- menu(c("Proceed", "Abort", "Test again"), title = "Check parameters provided and sample image generated. Should the analysis proceed?")
      #If user decides to stop then abort function and return stop message
      if(answer == 2) stop("The function has been stopped. Please tune parameters and try again")
      if(answer == 1) Test_completed <- "Proceed"
    }


    ######PERFORM ALIGNMENT#####
    future::plan("future::multisession", workers = N_cores)
    options(future.globals.maxSize = Inf, future.rng.onMisuse = "ignore")
    furrr::furrr_options(scheduling = Inf)

    #Execute the function
    #Generate a tibble containing the alignment results (the new coordinates and the model)
    cat("\n Performing alignment")
    Alignment_result_tibble <-
      furrr::future_map_dfr(1:nrow(Look_up_tibble), function(Look_up_tibble_index){

        #Generate the sample list
        Sample_names <- na.omit(unlist(Look_up_tibble[Look_up_tibble_index, ]))
        Sample_list <- purrr::map(1:length(Sample_names), function(Index){
          List_name <- names(Sample_names)[Index]
          Sample_name <- Sample_names[Index]

          #Obtain the individual sample
          Sample <- Samples_list[[List_name]][[Sample_name]]

          #If Hull based then return the hull points
          if(Concave_Hull){
            #First compute the concave hull of both datasets and sample the concave hull polygon
            Sample <- sf::st_as_sf(Sample, coords = c("X", "Y")) #Generate sf object
            Sample <- sf::st_union(Sample) #Turn it into a geometry object
            Sample <- sf::st_concave_hull(Sample, ratio = Hull_ratio) #Obtain the concave hull
            Sample <- sf::st_boundary(Sample) #Turn it into a line
            Sample <- sf::st_line_sample(Sample, density = 1/Polygon_sampling_density) #Sample the line according to the desired density
            Sample <- sf::st_coordinates(sf::st_cast(Sample, "POINT")) #Obtain the coordinates of the sampled points
          }
          #If not turn it into a matrix
          else{
            Sample <- as.matrix(Sample %>% dplyr::select(X, Y))
          }

          #Return the final productr
          return(Sample)
        })

        #The names will be the ones of the list element where they come from, not the sample name
        names(Sample_list) <- names(Sample_names)

        Model_list <- list()
        Error_list <- list()

        #This sample list will be modified in place by the for loop
        for(i in 1:(length(Sample_list)-1)){
          #Define Target and Moving
          Taget <- Sample_list[[i]]
          Moving <- Sample_list[[i+1]]

          #Perform ICP alignment
          Alignment_results <- ICP_cell_aligner(Moving_coords = Moving, Target_coords = Taget, Iterations = Iterations)

          #Transform the Sample_list item with the new coords
          Sample_list[[i+1]] <- Alignment_results$Coordinates

          #Add the model and erro results to the list
          Model_list[[i]] <- Alignment_results$Model
          Error_list[[i]] <- Alignment_results$Error
        }

        #Remove the first element of the Sample_list
        Sample_list <- Sample_list[-1]

        #Generate the tibble with the results
        Reference_tibble <- tibble::tibble(Subject_Names = Sample_names[1],
                                           List_Names = names(Sample_names)[1],
                                           New_coords = NA,
                                           Model = NA,
                                           Error = NA,
                                           Reference_Name = NA_character_,
                                           Reference_List = NA_character_)

        Cases_tibble <- tibble::tibble(Subject_Names = Sample_names[-1],
                                       List_Names = names(Sample_names)[-1],
                                       New_coords = Sample_list,
                                       Model = Model_list,
                                       Error = unlist(Error_list),
                                       Reference_Name = Sample_names[-length(Sample_names)],
                                       Reference_List = names(Sample_names)[-length(Sample_names)])

        Result <- dplyr::bind_rows(Reference_tibble, Cases_tibble)

        return(Result)

      }, .progress = TRUE)
    #Remove the sample list, the Final_Cell_feature_list and run gc
    rm(Samples_list)
    rm(Final_Cell_feature_list)
    gc()

    ####SUBSTITUTE THE SPATIAL COORDS####
    cat("\n Substituting spatial coords")
    #We will need to redo the Cell_feature_list
    Cell_feature_list <- mori::share(Cell_feature_list)
    List_names <- names(Cell_feature_list)

    Cell_feature_list <-
      furrr::future_map(1:length(Cell_feature_list), function(List_item_index){
        #Get the list item name and the actual list item
        Current_list_name <- names(Cell_feature_list)[List_item_index]
        Current_list_item <- Cell_feature_list[[List_item_index]]

        #Iterate over the subject names of the list items
        Final_list <-
          purrr::map(unique(Current_list_item$Subject_Names), function(Sample_name){
            #Get the current sample and the model results
            Current_sample <- Current_list_item %>% dplyr::filter(Subject_Names == Sample_name)
            Alignment_sample_result <- Alignment_result_tibble %>% dplyr::filter(Subject_Names == Sample_name, List_Names == Current_list_name)

            #If the alignment result has not been performed
            if(nrow(Alignment_sample_result) == 0) return(list(Sample = Current_sample,
                                                               Translation_correction = NA,
                                                               New_error = NA))

            #If the current sample is the absolute reference then return it as it is
            if(is.na(Alignment_sample_result$Reference_Name)) return(list(Sample = Current_sample,
                                                                          Translation_correction = NA,
                                                                          New_error = NA))

            #If the alignment has been performed with the whole cells (No phenotypes, no DBSCAN, no Concave hull), just substitute the X Y coords
            if(all(is.null(Phenotypes_included), !Apply_DBSCAN_denoise, !Concave_Hull)){
              Current_sample[(c("X", "Y"))] <- Alignment_sample_result[["New_coords"]]
              Translation_correction <- NA
              New_error <- NA
            }

            #Else apply the models
            else{
              #Model
              Model <- Alignment_sample_result[["Model"]]
              Final_coords <- as.matrix(Current_sample[(c("X", "Y"))])

              #Perform Flip
              Final_coords <- Flip_xy(Final_coords, mode = Model[[1]]$Flip_action)
              #Perform rotation
              Final_coords <- Rotate_points(Final_coords, theta = Model[[1]]$Angle_radians)
              #Perform translation
              Final_coords[,1] <- Final_coords[,1] + Model[[1]]$Translation[1]
              Final_coords[,2] <- Final_coords[,2] + Model[[1]]$Translation[2]

              #The initial translation correction is NA
              Translation_correction <- NA
              New_error <- NA

              #If final ICP required
              if(Perform_final_ICP_refinement){
                #Get the target coords
                Target_list_name <- Alignment_sample_result$Reference_List
                Target_sample_name <- Alignment_sample_result$Reference_Name
                Target_coords <- as.matrix(Cell_feature_list[[Target_list_name]] %>% dplyr::filter(Subject_Names == Target_sample_name) %>% dplyr::select(X, Y))

                #Perform the final ICP
                ICP_final_correction <- Morpho::icpmat(x = Final_coords, y = Target_coords, iterations = Iterations, type = "rigid")
                colnames(ICP_final_correction) <- c("X", "Y")

                #Obtain the final modification (trnaslation and final error)
                Added_X_translation <- as.numeric(min(ICP_final_correction[,1] - Final_coords[,1]))
                Added_Y_translation <- as.numeric(min(ICP_final_correction[,2] - Final_coords[,2]))
                Translation_correction <- c(X_shift = Added_X_translation, Y_shift = Added_Y_translation)
                Final_error_sense <- sqrt(sum(Rvcg::vcgKDtree(Target_coords, ICP_final_correction, k=1)$distance^2))
                Final_error_antisense <- sqrt(sum(Rvcg::vcgKDtree(ICP_final_correction, Target_coords, k=1)$distance^2))
                New_error <- max(c(Final_error_sense, Final_error_antisense))

                #Substitute the original value
                Final_coords <- ICP_final_correction
              }

              #Assign the coordinates to the current sample
              Current_sample[(c("X", "Y"))] <- Final_coords
            }

            #Return the current sample
            return(
              list(Sample = Current_sample,
                   Translation_correction = Translation_correction,
                   New_error = New_error)
            )
          })

        #Add names to the list
        names(Final_list) <- unique(Current_list_item$Subject_Names)

        #Return the final tibble with all the Samples
        return(Final_list)
      }, .progress = TRUE)
    future::plan("future::sequential")
    gc()

    #Change the names
    names(Cell_feature_list) <- List_names

    #If no final ICP performed then return the Cell feature list merging the samples
    if(!Perform_final_ICP_refinement){
      #Obtain the cell feature list
      Cell_feature_list <-
        purrr::map(1:length(Cell_feature_list), function(List_item_index){
          Current_sample <- Cell_feature_list[[List_item_index]]
          return(purrr::map_dfr(Current_sample, ~.[["Sample"]]))
        })

      #Check that all X and Y coordinates are above 0
      Min_x <- min(purrr::map_dbl(Cell_feature_list, function(List_item) min(List_item$X)))
      if(Min_x < 0){
        Cell_feature_list <- purrr::map(Cell_feature_list, function(List_item) List_item %>% dplyr::mutate(X = X + abs(Min_x)))
        cat(paste0("\n", "The minimum X value of some of the cells after alignment was below zero",
                   "\n", "Adding ", abs(Min_x), " to all X values"))
      }

      Min_y <- min(purrr::map_dbl(Cell_feature_list, function(List_item) min(List_item$Y)))
      if(Min_y < 0){
        Cell_feature_list <- purrr::map(Cell_feature_list, function(List_item) List_item %>% dplyr::mutate(Y = Y + abs(Min_y)))
        cat(paste0("\n", "The minimum Y value of some of the cells after alignment was below zero",
                   "\n", "Adding ", abs(Min_y), " to all X values"))
      }

      names(Cell_feature_list) <- List_names

      #Return the final elements
      return(list(Cell_feature_list = Cell_feature_list,
                  Alignment_results = Alignment_result_tibble)
      )
    }

    else{
      #Substitute the model and error for every item of the look up table
      for(i in 1:nrow(Alignment_result_tibble)){
        #If the item belongs to an absolute reference sample don´t do anything
        if(is.na(Alignment_result_tibble$Reference_Name[i])) next

        #Else find the new error and x and Y shift
        List_name <- Alignment_result_tibble$List_Names[i]
        Sample_name <- Alignment_result_tibble$Subject_Names[i]
        New_error <- Cell_feature_list[[List_name]][[Sample_name]][["New_error"]]
        New_X_shift <- Cell_feature_list[[List_name]][[Sample_name]][["Translation_correction"]][1]
        New_Y_shift <- Cell_feature_list[[List_name]][[Sample_name]][["Translation_correction"]][2]

        #Subtitute the new values in the Alignment_result_tibble
        Alignment_result_tibble$Error[i] <- New_error
        Alignment_result_tibble$Model[i][[1]][["Translation"]][1] <- Alignment_result_tibble$Model[i][[1]][["Translation"]][1] + New_X_shift
        Alignment_result_tibble$Model[i][[1]][["Translation"]][2] <- Alignment_result_tibble$Model[i][[1]][["Translation"]][2] + New_Y_shift
      }

      #Obtain the cell feature list
      Cell_feature_list <-
        purrr::map(1:length(Cell_feature_list), function(List_item_index){
          Current_sample <- Cell_feature_list[[List_item_index]]
          return(purrr::map_dfr(Current_sample, ~.[["Sample"]]))
        })
      #Check that all X and Y coordinates are above 0
      Min_x <- min(purrr::map_dbl(Cell_feature_list, function(List_item) min(List_item$X)))
      if(Min_x < 0){
        Cell_feature_list <- purrr::map(Cell_feature_list, function(List_item) List_item %>% dplyr::mutate(X = X + abs(Min_x)))
        cat(paste0("\n", "The minimum X value of some of the cells after alignment was below zero",
                   "\n", "Adding ", abs(Min_x), " to all X values"))
      }

      Min_y <- min(purrr::map_dbl(Cell_feature_list, function(List_item) min(List_item$Y)))
      if(Min_y < 0){
        Cell_feature_list <- purrr::map(Cell_feature_list, function(List_item) List_item %>% dplyr::mutate(Y = Y + abs(Min_y)))
        cat(paste0("\n", "The minimum Y value of some of the cells after alignment was below zero",
                   "\n", "Adding ", abs(Min_y), " to all X values"))
      }


      names(Cell_feature_list) <- List_names

      #Return the final elements
      return(list(Cell_feature_list = Cell_feature_list,
                  Alignment_results = Alignment_result_tibble)
      )
    }
  }
