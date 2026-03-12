#' Calculate the distance from a cell to the closest positive pixel
#'
#' Given a thresholded image (either binary or multilevel) and the cell feature information matrix, the function will compute the distance from cell to the closest non-zero value pixel.
#' This can be used to analyze the spatial interaction of cells with extra-cellular elements in the image. Images can be thresholded using the [Pixel_Threshold_calculator()] function.
#'
#'
#' @param N_cores Integer. Number of cores to parallelize your computation.
#' @param Directory Character string specifying the path to the folder containing thresholded images.
#' @param Image_rotate (OPTIONAL) A integer value indicating the degrees of rotation of the image. If NULL no rotation is applied.
#' @param Image_x_flip A logical value indicating if X image flip should be performed.
#' @param Image_y_flip A logical value indicating if X image flip should be performed.
#' @param DATA A dataframe or tibble containing cell feature data.
#' @param Phenotypes_included A character vector indicating the phenotype labels that will be included in the analysis.
#' @param Pixel_distance_ratio (OPTIONAL) A numeric value indicating the ratio between pixel size / cell coordinates. Use this argument when cell coordinates have been transformed to a distance unit like microns.
#'
#' @returns A tibble containing Cell feature information along with the distance to the closest positive pixel.
#'
#' @seealso [Pixel_Threshold_calculator()]
#'
#' @examples
#' \dontrun{
#' Cell_to_pixel_distance_calculator(
#'     N_cores = 1,
#'     Directory = "Path/To/Thresholded_images",
#'     Image_rotate = 90,
#'     Image_x_flip = FALSE,
#'     Image_y_flip = TRUE,
#'     DATA = CSM_Phenotypecell_test,
#'     Phenotypes_included = unique(CSM_Phenotypecell_test$Phenotype),
#'     Pixel_distance_ratio = NULL
#'  )
#'}
#'
#' @export

Cell_to_pixel_distance_calculator <-
    function(N_cores = 1,
             Directory,
             Image_rotate = NULL,
             Image_x_flip = FALSE,
             Image_y_flip = FALSE,
             DATA,
             Phenotypes_included,
             Pixel_distance_ratio = NULL){

      #######Check suggested packages#######
      {
        if(!requireNamespace("magick", quietly = FALSE)) stop(
          paste0("magick CRAN package is required to execute the function. Please install using the following code: ",
                 expression(install.packages("magick")))
        )
        if(!requireNamespace("EBImage", quietly = TRUE)) stop(
          paste0("EBImage Bioconductor package is required to execute the function. Please install using the following code: ",
                 expression({
                   if (!require("BiocManager", quietly = TRUE))
                     install.packages("BiocManager")

                   BiocManager::install("EBImage")
                 })
          )
        )
        if(!requireNamespace("rtree", quietly = FALSE)) stop(
          paste0("rtree GitHub package is required to execute the function. Please install using the following code: ",
                 expression(remotes::install_github("akoyabio/rtree")))
        )
      }


      #######Specify that on exit return to single core and run gc#######
      on.exit({
        future::plan("future::sequential")
        gc()
      })

      #######Argument check general arguments#######
      Argument_checker <- c(N_cores_OK = (N_cores >= 1 & N_cores%%1 == 0),
                            Empty_directory_OK = length(dir(Directory)) >= 1,
                            Image_rotate_OK = if(!is.null(Image_rotate)) {
                              all(is.numeric(Image_rotate), Image_rotate >= 0, Image_rotate <= 360, Image_rotate%%1 == 0)
                            } else(TRUE),
                            Image_x_flip_OK = is.logical(Image_x_flip),
                            Image_y_flip_OK = is.logical(Image_y_flip),
                            DATA_OK = all(identical(names(DATA)[c(1:4)], c("Cell_no", "X", "Y", "Subject_Names")), "Phenotype" %in% names(DATA)),
                            Phenotypes_included_OK = all(Phenotypes_included %in% DATA[["Phenotype"]]),
                            Pixel_distance_ratio_OK = if(!is.null(Pixel_distance_ratio)) {
                              all(is.numeric(Pixel_distance_ratio), Pixel_distance_ratio > 0)
                            } else(TRUE)
      )

      Stop_messages <- c(N_cores_OK = "N_cores must be an integer value > 0",
                         Empty_directory_OK = "No files found at the directory provided. Please check out the path.",
                         Image_rotate_OK = "Image_rotate must be either NULL or a integer value between 0 and 360",
                         Image_x_flip_OK = "Image_x_flip must be a logical value",
                         Image_y_flip_OK = "Image_y_flip must be a logical value",
                         DATA_OK = "DATA must be adequately formatted and must contain a column containing Phenotype information",
                         Phenotypes_included_OK = paste0("Phenotypes_included must be any of the following: ", stringr::str_c(DATA[["Phenotype"]], collapse = ", ")),
                         Pixel_distance_ratio_OK = "Pixel_distance_ratio must be either NULL or a numeric value > 0"
      )
      #Check arguments and stop if necessary
      if(!all(Argument_checker)){
        stop(cat(Stop_messages[!Argument_checker],
                 fill = sum(!Argument_checker)))
      }

      #Check specifically that directory contains adequate files, and that these are present in Subject_Names of data
      Full_image_names <- dir(Directory, full.names = TRUE)
      Image_names <- dir(Directory, full.names = FALSE)

      #Select only the thresholded and not the tissue masks
      Full_image_names_selected <- Full_image_names[stringr::str_detect(Image_names, "Thresholded")]
      Image_names_selected <- Image_names[stringr::str_detect(Image_names, "Thresholded")]
      Image_names_selected_list <- stringr::str_split(Image_names_selected, "_")


      #Check that all have been processed using the Pixel_Threshold_calculator
      if(!all(purrr::map_lgl(Image_names_selected_list, ~.[[1]] == "Processed"))){
        Problematic_images <- Image_names_selected[!map_lgl(Image_names_selected_list, ~.[[1]] == "Processed")]
        stop(paste0("According to image names, the following images have not been processed using the Pixel_Threshold_calculator: ", stringr::str_c(Problematic_images, collapse = ", ")))
      }
      if(!all(purrr::map_lgl(Image_names_selected_list, ~length(.) == 4))){
        Problematic_images <- Image_names_selected[!map_lgl(Image_names_selected_list, ~length(.) == 4)]
        stop(paste0("According to image names, the following images have not been processed using the Pixel_Threshold_calculator: ",
                    stringr::str_c(Problematic_images, collapse = ", ")))
      }
      #Check that all images in the directory are unique
      if(length(unique(purrr::map_chr(Image_names_selected_list, ~.[[2]]))) != length(Image_names_selected_list)){
        Problematic_images <-purrr::map_chr(Image_names_selected_list, ~.[[2]])[duplicated(purrr::map_chr(Image_names_selected_list, ~.[[2]]))]
        stop(paste0("The following image appear to be duplicated in the directory: ",
                    stringr::str_c(Problematic_images, collapse = ", "),
                    ". Images must be unique."))
      }
      #Check that all images in the directory are from the same target
      if(length(unique(purrr::map_chr(Image_names_selected_list, ~.[[3]]))) != 1){
        Number_markers <- sort(table(purrr::map_chr(Image_names_selected_list, ~.[[3]])), decreasing = TRUE)
        stop(paste0("Besides ", names(Number_markers)[1], " the directory contains images using the following markers: ",
                    stringr::str_c(names(Number_markers)[-1], collapse = ", "),
                    ". A single marker type must be included in the directory"))
      }
      #Check that all images in the directory have been thresholded homogenously
      if(length(unique(purrr::map_chr(Image_names_selected_list, ~.[[4]]))) != 1){
        Number_threshold_strategies <- sort(table(purrr::map_chr(Image_names_selected_list, ~.[[4]])), decreasing = TRUE)
        stop(paste0("Besides ", names(Number_threshold_strategies)[1], " the directory contains images thresholded using the following strategy: ",
                    stringr::str_c(names(Number_markers)[-1], collapse = ", "),
                    ". A single threshold strategy must be selected."))
      }
      #Check that image names are present in DATA Subject_Names calculating the closest
      Subject_names_in_images <- purrr::map_chr(Image_names_selected_list, ~.[[2]])
      Subject_names_in_data <- unique(DATA$Subject_Names)
      
      #######Look up table generation#######
      #Generate a look up tibble with the image name, the closest name in data$Subject_Names and the image path
      #Evaluate which subject name in data is the closest one to the subject name in images
      Closest_name_vector <-purrr::map_dbl(Subject_names_in_images, function(Image_name){
        Distance_vector <- adist(Image_name, Subject_names_in_data, fixed = TRUE, ignore.case = TRUE)
        which.min(Distance_vector)
      })
      Names_tibble <- tibble(Images_names = Subject_names_in_images,
                             Subject_names_in_data = Subject_names_in_data[Closest_name_vector])
      Names_tibble$Identical <- apply(Names_tibble, MARGIN = 1, function(Row) identical(Row[[1]], Row[[2]]))

      #If Subject_Names in data are duplicated then stop the computing
      if(sum(duplicated(Names_tibble$Subject_names_in_data)) > 0){
        Duplicated_names_list <-
          purrr::map(Names_tibble$Subject_names_in_data[duplicated(Names_tibble$Subject_names_in_data)],
                     function(Duplicated_names) unname(unlist(Names_tibble %>% dplyr::filter(Subject_names_in_data == Duplicated_names) %>% dplyr::select(Images_names))))
        names(Duplicated_names_list) <- stringr::str_c(Names_tibble$Subject_names_in_data[duplicated(Names_tibble$Subject_names_in_data)], " in DATA_matched by: ")

        print(Duplicated_names_list)
        stop("Multiple image names are matched with same Subject_Names in data. Please check image names in directory.")
      }

      #If match is not exact print a message
      if(sum(!Names_tibble$Identical) > 0){
        message("Subject names in the following images do not exactly. Approximate match will be used")
        print(Names_tibble %>% dplyr::filter(!Identical))
      }

      #Remove Subject_Names in data not present in Names_tibble
      if(sum(!unique(DATA$Subject_Names) %in% Names_tibble$Subject_names_in_data) > 0){
        Absent_Subject_Names <- unique(DATA$Subject_Names)[!unique(DATA$Subject_Names) %in% Names_tibble$Subject_names_in_data]
        message("The following Subject_Names in DATA are not present in image names: ",
                stringr::str_c(Absent_Subject_Names, collapse = ", "),
                ". They will be removed before analysis")
        DATA <- DATA %>% dplyr::filter(Subject_Names %in% Names_tibble$Subject_names_in_data)
      }

      #Generate the final Names_tibble including image path
      Names_tibble$Image_URL <- Full_image_names_selected

      #Also note the target being measured (used to name distance column in final result)
      Target_being_measured <- Image_names_selected_list[[1]][[3]]

      #Check that thresholded images have positive pixels and check that images have adequate number of target cells
      print("Removing samples without target cells")
      #If no phenotypes present in image remove images
      Images_with_cells <-
        purrr::map_lgl(Names_tibble$Subject_names_in_data, function(Image_name_in_subject){
          DATA <- DATA %>% dplyr::filter(Subject_Names == Image_name_in_subject)
          DATA <- DATA %>% dplyr::filter(Phenotype %in% Phenotypes_included)
          nrow(DATA) > 0
        })
      if(!all(Images_with_cells)){
        message(paste0("The following image do not have target cells. They will be removed: ",
                       stringr::str_c(Names_tibble$Subject_names_in_data[!Images_with_cells], collapse = ", "))
        )
        Names_tibble <- Names_tibble[Images_with_cells, ]
      }
      if(!nrow(Names_tibble) > 0) stop("No images with adequate number of cells.")
      
      #######Check that images have positive pixels#######
      print("Removing images without non-cero value pixels")

      #Check that images have positive pixels
      #Will iterate for every image in the names tibble
      future::plan("future::multisession", workers = N_cores)
      options(future.globals.maxSize = Inf, future.rng.onMisuse = "ignore")
      furrr::furrr_options(scheduling = Inf)
      Images_with_pixels <-
        furrr::future_map_lgl(Names_tibble$Image_URL, function(Image){
          Image <- magick::image_read(Image) %>% magick::as_EBImage()
          sum(Image) > 0
        })
      #Return to single core
      future::plan("future::sequential")
      gc()
      if(!all(Images_with_pixels)){
        Colored_print(paste0("\n", "The following image/s do not have positive pixels. They will be removed: ",
                       stringr::str_c("\n", Names_tibble$Subject_names_in_data[!Images_with_pixels], collapse = ", ")),
                      color = "red"
                      )
        Names_tibble <- Names_tibble[Images_with_pixels, ]
      }
      if(!nrow(Names_tibble) > 0) stop("No images with adequate number of pixels and/or cells.")
      
      
      #######RANDOM TEST#######
      #Run a random test with the image with the lowest cell counts
      #Find the names of all the samples
      Sample_names <- (DATA %>% dplyr::filter(Subject_Names %in% unique(Names_tibble$Subject_names_in_data)) %>%
                            dplyr::count(Subject_Names) %>% dplyr::arrange(n))[[1]]
      
      Answer_value <- 3
      while(Answer_value == 3){
        Selected_sample <- sample(Sample_names, size = 1)
        #Obtain the points in the sample
        DATA_test <- DATA %>% dplyr::filter(Subject_Names == Selected_sample)
        #Obtain the path to the image
        Image_path <- Names_tibble$Image_URL[which(Names_tibble$Subject_names_in_data == Selected_sample)]
        
        #get the image
        Image <- magick::image_read(as.character(Image_path))
        
        #Execute rotation if required
        if(!is.null(Image_rotate)) Image <- Image %>% magick::image_rotate(degrees = Image_rotate)
        #Execute x or y flip if required
        if(Image_x_flip) Image <- Image %>% magick::image_flop()
        if(Image_y_flip) Image <- Image %>% magick::image_flip() #Here we will flip the image as required by user as annotation_raster will flip it for plotting purposes
        
        #obtain the max width and height
        X_max <- magick::image_info(Image)$width
        Y_max <- magick::image_info(Image)$height
        
        #plot both results
        print(paste0("Generating a sample overlay image using ", Selected_sample))
        plot(
          ggplot() +
            annotation_raster(Image, xmin = 1, xmax = X_max, ymin = 0, ymax = Y_max) +
            geom_point(aes(x = X, y = Y), color = "red", data = DATA_test) +
            theme_minimal() +
            guides(color = "none") +
            scale_x_continuous("", labels = NULL, limits = c(1, X_max)) +
            scale_y_continuous("", labels = NULL, limits = c(1, Y_max)) +
            ggtitle(Selected_sample) +
            theme(panel.grid = element_blank(),
                  panel.background = element_rect(fill = "black"),
                  plot.title = element_text(hjust = 0.5))
        )
        
        #Generate a menu to proceed with compuation
        answer <- menu(c("Proceed", "Abort", "Test again"), title = "Check parameters provided and sample image generated. Should the analysis proceed?")
        #If user decides to stop then abort function and return stop message
        if(answer == 2) stop("The function has been stopped. Please tune parameters and try again")
        if(answer == 1) Answer_value <- "Proceed"
      }
      

      #If OK then run the final analysis
      print("Running distance to positive pixel computation")
      
      
      #######COMPUTE THE ACTUAL DISTANCE BETWEEN CELLS AND PIXELS#######

      #Will iterate for every image in the names tibble
      future::plan("future::multisession", workers = N_cores)
      options(future.globals.maxSize = Inf, future.rng.onMisuse = "ignore")
      furrr::furrr_options(scheduling = Inf)

      RESULTS <- suppressMessages(
        furrr::future_map_dfr(seq_along(1:nrow(Names_tibble)), function(Index_image){
          #Import DATA
          DATA_image <- DATA %>% dplyr::filter(Subject_Names == Names_tibble[["Subject_names_in_data"]][[Index_image]])

          #Remove Phenotypes not included
          DATA_image <- DATA_image %>% dplyr::filter(Phenotype %in% Phenotypes_included)

          #Import Image
          Image <- magick::image_read(Names_tibble$Image_URL[[Index_image]])
          
          #Execute rotation if required
          if(!is.null(Image_rotate)) Image <- Image %>% magick::image_rotate(degrees = Image_rotate)
          #Execute x or y flip if required
          if(Image_x_flip) Image <- Image %>% magick::image_flop()
          if(!Image_y_flip) Image <- Image %>% magick::image_flip() #Since the function is prepared to work with native Y axis keep it as it is if the user decides to flip it
          
          #Turn to an EBImage object
          Image <- magick::as_EBImage(Image)
          #Colapse it to a tibble
          Image_tibble <- as_tibble(expand.grid(X = 1:dim(Image)[[1]], Y = 1:dim(Image)[[2]]))
          Image_tibble$Value <- as.vector(Image)
          
          #remove the image and run gc
          rm(Image)
          gc()
          
          #Keep non-0 values 
          Image_tibble <- Image_tibble %>% dplyr::filter(Value != 0)
          #If the pixel distance ratio is provided modify the values of the positive pixel values
          if(!is.null(Pixel_distance_ratio)) Image_mage_tibble <- Image_mage_tibble %>%dplyr::mutate(X = X*Pixel_distance_ratio, Y = Y*Pixel_distance_ratio)

          #Cells will be origin, pixels will be targets.
          COO_info <- cbind(DATA_image[["X"]], DATA_image[["Y"]])

          #First calculate min distance to closest neighbors for binary thresholded images
          if(length(unique(Image_tibble$Value)) == 1){
            #Calculate the targets
            Targets <- rtree::RTree(cbind(Image_tibble[["X"]], Image_tibble[["Y"]]))

            #Calculate the closest pixel
            Index <- rtree::knn(Targets, COO_info, 1L)
            Closest_pixel_tibble <- Image_tibble[unlist(Index), c(1,2)]
            names(Closest_pixel_tibble) <- c("Target_X", "Target_Y")

            Closest_pixel_tibble <- suppressMessages(dplyr::bind_cols(COO_info, Closest_pixel_tibble))
            names(Closest_pixel_tibble)[c(1:2)] <- c("X", "Y")
            Closest_pixel_tibble <- Closest_pixel_tibble %>% dplyr::mutate(Dist_X = (X - Target_X)^2,
                                                                           Dist_Y = (Y - Target_Y)^2,
                                                                           DIST = sqrt(Dist_X + Dist_Y))
            DATA_image$DIST <- Closest_pixel_tibble$DIST
            names(DATA_image)[ncol(DATA_image)] <- stringr::str_c(Target_being_measured, "_DIST", collapse = "_")
            return(DATA_image)
          }

          #Then calculate min distance for every value of multithresholded images
          if(length(unique(Image_tibble$Value)) > 1){
            #Obtain the unique values
            Unique_values <- unique(Image_tibble$Value)

            #Calculate a tibble containing distances to every unique value
            Closest_pixel_tibble <-
              purrr::map_dfc(seq_along(1:length(Unique_values)), function(Indivudal_target){
                Image_tibble <- Image_tibble %>% dplyr::filter(Value == Unique_values[Indivudal_target])

                Targets <- rtree::RTree(cbind(Image_tibble[["X"]], Image_tibble[["Y"]]))

                #Calculate the closest pixel
                Index <- rtree::knn(Targets, COO_info, 1L)
                Closest_pixel_tibble <- Image_tibble[unlist(Index), c(1,2)]
                names(Closest_pixel_tibble) <- c("Target_X", "Target_Y")

                Closest_pixel_tibble <- suppressMessages(bind_cols(COO_info, Closest_pixel_tibble))
                names(Closest_pixel_tibble)[c(1:2)] <- c("X", "Y")
                Closest_pixel_tibble <- Closest_pixel_tibble %>% dplyr::mutate(Dist_X = (X - Target_X)^2,
                                                                               Dist_Y = (Y - Target_Y)^2,
                                                                               DIST = sqrt(Dist_X + Dist_Y))
                Closest_pixel_tibble <- Closest_pixel_tibble %>% dplyr::select(DIST)
                names(Closest_pixel_tibble) <- stringr::str_c(Target_being_measured, "_DIST_", round(Unique_values[Indivudal_target], digits = 3), collapse = "")
                return(Closest_pixel_tibble)
              })
            DATA_image <-dplyr::bind_cols(DATA_image, Closest_pixel_tibble)
            return(DATA_image)
          }
        }, .progress = TRUE)
      )

      #Return to single core
      future::plan("future::sequential")
      gc()

      return(RESULTS)

    }

