#' Computes morphology metrics from tiled image cell data
#'
#' Given a list containing cell data from tiled images (generated using the [Image_tiling_processing_function()] or the [Tiled_Image_Clustering_function()]) the function
#' identifies discrete structures and calculates a variety of morphology metrics (see details).
#'
#'
#' @param N_cores Integer. Number of cores to parallelize your computation.
#'
#' @param Tiled_Images A list containing tiled images obtained using [Image_tiling_processing_function()] or the [Tiled_Image_Clustering_function()].
#' @param Variable_name A character value indicating the column name containing the cell labels.
#' @param Feature_to_analyze A character value indicating the cell type or neighborhood type that will be morphologically analyzed
#'
#' @param Min_tile_number A integer value specifying the minimum number of tiles any morphological structure must contain to be considered real.
#' @param Smooth_amount A numeric value specifying the amount of Gaussian smoothing to be applied to the tile information (higher values yield larger and more smoothened structures).
#' @param Tolerance_value A numeric value specifying the minimum height of the structure between its highest point (seed) and the point where it contacts another object (MAY BE NULL).
#' @param Neighborhood_distance A integer value indicating the radius of the neighborhood in pixels for the detection of neighboring objects. Higher value smooths out small objects.
#'
#' @param Minimum_cell_no_per_tile (If Tiled_Images are generated using the [Image_tiling_processing_function()]) A integer value indicating the minimum number of cells of interest a tile must contain to be considered positive (default is 1).
#'
#' @details
#' The function will identify discrete structures using the watershed algorithm implemented in the simpleSeg Bioconductor package. To this end, the image tiles will be assigned a
#' positive (1) or negative (0) value according to the presence or absence of the desired cell feature. Afterwards, these positive/negative tile pattern is transformed
#' into a pseudo-image and fed to the watershed algorithm to compute the structure mask image. This image is then analyzed to obtain the morphological features of the structures idenfified
#' using the cytomapper Bioconductor package.
#'
#' @returns A list containing a tibble with the structure morphology data and the structure mask.
#'
#' @seealso [Image_tiling_processing_function()], [Tiled_Image_Clustering_function()], [Structure_morphometry_analyzer()], [Structure_morphometry_graph_maker()]
#'
#' @examples
#' \dontrun{
#' #Divide cells into discrete tiles
#'Tiled_Images <-
#' Image_tiling_processing_function(
#'  N_cores = 1,
#'  DATA = CSM_Phenotypecell_test,
#'  Tile_width = 25,
#'  Tile_height = 25,
#'  Variables_to_keep = "Phenotype")
#'
#' #Try to identify morphological structures for OTHER cells according to tile distribution
#'Structures <-
#' Structure_morphometry_calculator(N_cores = 1,
#'
#'                                 Tiled_Images = Tiled_Images,
#'                                 Variable_name = "Phenotype",
#'                                 Feature_to_analyze = "OTHER",
#'
#'                                 Min_tile_number = 1,
#'                                 Smooth_amount = 0.01,
#'                                 Tolerance_value = 1,
#'                                 Neighborhood_distance = 1,
#'                                 Minimum_cell_no_per_tile = 2)
#'
#' #Analyze the results
#' Structure_morphometry_analyzer(DATA = Structures)
#'
#' #Graph the results
#' Structure_morphometry_graph_maker(Tiled_Images = Tiled_Images,
#'                                   Variable_name = "Phenotype",
#'                                   DATA_Morphometry = Structures,
#'                                   Color_by = "Structure_ID",
#'                                   Image_name = "ABCM22001_B09_MiniCrop.tif")
#'}
#'
#'
#' @export


Structure_morphometry_calculator <-
  function(N_cores = 1,

           Tiled_Images,
           Variable_name,
           Feature_to_analyze,


           Min_tile_number = 4,
           Smooth_amount = 0.0000001,
           Tolerance_value = 0.0000001,
           Neighborhood_distance = 1,

           Minimum_cell_no_per_tile = 1
  ){

    ########Check suggested packages########
    {
      if(!requireNamespace("EBImage", quietly = TRUE)) stop(
        paste0("EBImage Bioconductor package is required to execute the function. Please install using the following code: ",
               expression({
                 if (!require("BiocManager", quietly = TRUE))
                   install.packages("BiocManager")

                 BiocManager::install("EBImage")
               })
        )
      )
      if(!requireNamespace("simpleSeg", quietly = TRUE)) stop(
        paste0("simpleSeg Bioconductor package is required to execute the function. Please install using the following code: ",
               expression({
                 if (!require("BiocManager", quietly = TRUE))
                   install.packages("BiocManager")

                 BiocManager::install("simpleSeg")
               })
        )
      )
      if(!requireNamespace("S4Vectors", quietly = TRUE)) stop(
        paste0("S4Vectors Bioconductor package is required to execute the function. Please install using the following code: ",
               expression({
                 if (!require("BiocManager", quietly = TRUE))
                   install.packages("BiocManager")

                 BiocManager::install("S4Vectors")
               })
        )
      )
      if(!requireNamespace("cytomapper", quietly = TRUE)) stop(
        paste0("cytomapper Bioconductor package is required to execute the function. Please install using the following code: ",
               expression({
                 if (!require("BiocManager", quietly = TRUE))
                   install.packages("BiocManager")

                 BiocManager::install("cytomapper")
               })
        )
      )
    }

    ########Check arguments########
    if(!all(N_cores >= 1 & N_cores%%1 == 0)) stop("N_cores must be an integer value > 0")

    if(!all(length(Variable_name) == 1, is.character(Variable_name))) stop("Variable_name should be a single character value")
    if(!all(length(Feature_to_analyze) == 1, is.character(Feature_to_analyze))) stop("Feature_to_analyze should be a single character value")

    #Check arguments by generating a argument check vector and message vector
    Argument_checker <- c(Min_tile_number_OK = all(is.numeric(Min_tile_number), Min_tile_number%%1 == 0, Min_tile_number >= 1),
                          Smooth_amount_OK = all(is.numeric(Smooth_amount), Smooth_amount >= 0),
                          Tolerance_value_OK = if(!is.null(Tolerance_value)) {
                            all(is.numeric(Tolerance_value), Tolerance_value >= 0)
                          } else(TRUE),
                          Neighborhood_distance_OK = all(Neighborhood_distance%%1 == 0, Neighborhood_distance >= 0)
    )


    Stop_messages <- c(Min_tile_number_OK = "Min_tile_number must be a positive integer value larger than 0",
                       Smooth_amount_OK = "Smooth amount must be a positive numeric value",
                       Tolerance_value_OK = "Tolerance value must be NULL or a positive numeric value",
                       Neighborhood_distance_OK = "Neighborhood distance must be a positive integer value"
    )

    #Check arguments and stop if necessary
    if(!all(Argument_checker)){
      stop(cat(Stop_messages[!Argument_checker],
               fill = sum(!Argument_checker)))
    }

    ########Check datasets provided########
    #We will check if the dataset provided is a list
    if(!is.list(Tiled_Images)) stop("Tiled_Images should be a list created with the Image_tiling_processing_function or Tiled_Image_Clustering_function")

    List_content_tibble <- unique(purrr::map_lgl(names(Tiled_Images), ~tibble::is_tibble(Tiled_Images[[.]])))
    if(length(List_content_tibble) > 1) stop("Tiled_Images should be a list created with the Image_tiling_processing_function or Tiled_Image_Clustering_function")

    #If is tibble then they are derived from Tiled_Image_Clustering_function. Check specific variables
    if(List_content_tibble){
      #If variable name not found stop and print an error
      if(!all(purrr::map_lgl(Tiled_Images, ~(Variable_name %in% names(.))))) stop(paste0(Variable_name, " not found in the names of one or more of Tiled_Images elements"))

      #If no target found in any image then stop the computation
      if(!Feature_to_analyze %in% unique(unlist(purrr::map(Tiled_Images, ~.[[Variable_name]])))) stop(paste0(Feature_to_analyze,
                                                                                                             " not found in any of the Tiled_Images provided"))
      #Else proceed
      Tiled_Images <- purrr::map(Tiled_Images, function(Image){
        #Select the required columns
        Image <- Image %>%
          dplyr::select(Subject_Names, tile_X_centroid, tile_Y_centroid, tile_id, tile_xmin, tile_xmax, tile_ymin, tile_ymax, dplyr::all_of(Variable_name))

        names(Image)[9] <- "Target_column"

        Image <- Image %>% dplyr::mutate(Target_column = dplyr::case_when(Target_column == Feature_to_analyze ~ 1,
                                                                          TRUE ~ 0))
      })

      #Apply a correction as tiles without clustering have been removed from the analysis
      Names <- names(Tiled_Images)#Preserve names
      Tiled_Images <-
        purrr::map(1:length(Tiled_Images), function(Image_index){
          Image <- Tiled_Images[[Image_index]]

          #Compute all the possible tiles
          Full_tile_tibble <- tidyr::expand_grid(tile_xmin = unique(Image$tile_xmin), tile_ymin = unique(Image$tile_ymin))
          #Compute the tile width and height
          Tile_width <- unique(Image$tile_xmax - Image$tile_xmin)
          Tile_height <- unique(Image$tile_ymax - Image$tile_ymin)

          Full_tile_tibble <- Full_tile_tibble %>% dplyr::mutate(tile_xmax = tile_xmin + Tile_width,
                                                                 tile_ymax = tile_ymin + Tile_height)
          #Join with existing data
          Full_tile_tibble <-
            dplyr::left_join(Full_tile_tibble, Image, dplyr::join_by(tile_xmin, tile_xmax, tile_ymin, tile_ymax))

          Full_tile_tibble$tile_id[is.na(Full_tile_tibble$tile_id)] <- "ZGenerated_tile"
          Full_tile_tibble$tile_X_centroid[is.na(Full_tile_tibble$tile_X_centroid)] <- Full_tile_tibble$tile_xmin[is.na(Full_tile_tibble$tile_X_centroid)] + (Tile_width/2)
          Full_tile_tibble$tile_Y_centroid[is.na(Full_tile_tibble$tile_Y_centroid)] <- Full_tile_tibble$tile_ymin[is.na(Full_tile_tibble$tile_Y_centroid)] + (Tile_height/2)
          Full_tile_tibble$Subject_Names <- names(Tiled_Images)[[Image_index]]
          Full_tile_tibble$Target_column[is.na(Full_tile_tibble$Target_column)] <- 0


          Full_tile_tibble <- Full_tile_tibble %>% dplyr::select(Subject_Names, tile_X_centroid, tile_Y_centroid, tile_id, tile_xmin, tile_xmax, tile_ymin, tile_ymax, Target_column)
          return(Full_tile_tibble)
        })

      names(Tiled_Images) <- Names

      #Print that Minimum_cell_no_per_tile will be ignored
      print("Tiled_Images are derived from Tiled_Image_Clustering_function. Minimum_cell_no_per_tile will be ignored.")
    }

    #If the list elements are not tibbles (are lists) then they are derived from Image_tiling_processing_function and required modification
    if(!List_content_tibble){
      #Check that Variable_name is present in all images
      if(!all(purrr::map_lgl(Tiled_Images, ~(Variable_name %in% names(.[[2]]))))) stop(paste0(Variable_name, " not found in the names of one or more of Tiled_images elements"))

      #If no target found in any image then stop the computation
      #If no target found in any image then stop the computation
      if(!Feature_to_analyze %in% unique(unlist(purrr::map(Tiled_Images, ~.[[2]][[Variable_name]])))) stop(paste0(Feature_to_analyze,
                                                                                                                  " not found in any of the Tiled_images provided"))

      #If everything is OK then proceed by reducing the Tiled_Images to a list with elements of length 1
      Names <- names(Tiled_Images)#Preserve names
      Tiled_Images <- purrr::map(1:length(Tiled_Images), function(Image_index){
        Interim <- Tiled_Images[[Image_index]][[2]] %>% dplyr::rename("Target" = dplyr::all_of(Variable_name))

        Count_per_tile <- Interim %>% dplyr::filter(Target == Feature_to_analyze) %>% dplyr::count(tile_id)
        Tile_info <- left_join(Tiled_Images[[Image_index]][[1]], Count_per_tile, by = "tile_id")
        Tile_info <- Tile_info %>% mutate(Target_column = case_when(n >= Minimum_cell_no_per_tile ~ 1,
                                                                    TRUE ~ 0))
        Tile_info$Subject_Names <- names(Tiled_Images)[Image_index]

        Tile_info <- Tile_info %>%
          dplyr::select(Subject_Names, tile_X_centroid, tile_Y_centroid, tile_id, tile_xmin, tile_xmax, tile_ymin, tile_ymax, Target_column)

      })
      names(Tiled_Images) <- Names
    }

    ########Check that tiles are squared and remove any image without positive tiles########
    #Compute if tiles are squared
    Is_squared_tiles <-
      purrr::map_lgl(Tiled_Images, function(Image){
        Tile_width <- Image$tile_xmax[1] - Image$tile_xmin[1]
        Tile_height <- Image$tile_ymax[1] - Image$tile_ymin[1]

        return(near(Tile_width, Tile_height, tol = 2)) #sometimes tiles are minimally different even if created with equal size and height
      })

    #If any non-squared stop the computation
    if(!all(Is_squared_tiles)){
      Problematic_tiles <- names(Is_squared_tiles)[!Is_squared_tiles]

      stop(paste0("The tiles of the following images are not squared: ", stringr::str_c(Problematic_tiles, collapse = ", "), ". Morphometric analysis can only work with square tiles"))
    }

    #Compute if any image has no positive tiles
    Devoid_tiles <-
      purrr::map_lgl(Tiled_Images, function(Image){
        sum(Image$Target_column) == 0
      })
    if(any(Devoid_tiles)){
      Images_without_tiles <- names(Tiled_Images)[Devoid_tiles]

      Colored_print(paste0("\nThe following images do not contain any evaluable tiles and will be romved from the analysis: ",
                           "\n",
                           stringr::str_c(Images_without_tiles, collapse = "\n")),
                    color = "red")

      #Remove the actual images
      Tiled_Images <- Tiled_Images[!Devoid_tiles]
    }

    ########RUN TEST########
    Test_OK <- 3
    while(Test_OK == 3){

      #Enter another loop to adequately compute the segmentation mask
      Segmentation_mask_ok <- "NOT_OK"
      Images_for_test_remaining <- names(Tiled_Images)

      while(Segmentation_mask_ok == "NOT_OK"){
        #If no images to test remaining stop the function
        if(length(Images_for_test_remaining) == 0) stop("No images are adequate for morphometry computation given the provided parameters. Try modifying tiling strategy or morphometry parameters.")

        #Generate a test with a random sample
        Random_image_name <- sample(Images_for_test_remaining, size = 1)
        Example_image_tiles <- Tiled_Images[[Random_image_name]]

        #Calculate the number of rows and columns
        x_length <- length(unique(Example_image_tiles$tile_xmin))#Calculte the X pixel length
        y_length <- length(unique(Example_image_tiles$tile_ymin))#Calculate the Y pixel length

        #Arrange the image according to the X centroid
        Interim <- Example_image_tiles %>% dplyr::arrange(tile_xmin, desc(tile_ymin))#arrange the tibble in an adequate format

        #Compute the pseudo-image (a matrix)
        Pseudo_image <- matrix(Interim[["Target_column"]], nrow = y_length, ncol = x_length)#Calculate the matrix
        colnames(Pseudo_image) <- unique(Interim$tile_xmin)
        rownames(Pseudo_image) <- unique(Interim$tile_ymin)

        #Turn the pseudo-image to an EBImage and Cytomapper object
        Pseudo_image <- EBImage::Image(Pseudo_image) %>% EBImage::flip()
        Image <- cytomapper::CytoImageList(Pseudo_image)
        cytomapper::channelNames(Image) <- "COUNTS" #Modify the channel name
        S4Vectors::mcols(Image)$imageID <- "Pseudo_Image"#Modify name

        #Compute the structure mask and the object morphometry in a secure manner
        Seg_results <- try(simpleSeg::simpleSeg(Image,
                                                nucleus = "COUNTS",
                                                cellBody = "none",
                                                sizeSelection = Min_tile_number,
                                                smooth = Smooth_amount,
                                                transform = NULL,
                                                watershed = "intensity",
                                                tolerance = Tolerance_value,
                                                ext = Neighborhood_distance,
                                                discSize = 1,
                                                tissue = NULL,
                                                pca = FALSE,
                                                cores = 1)
        )


        try(
          #Change the ID of the mask
          S4Vectors::mcols(Seg_results)$imageID <- "Pseudo_Image"
        )

        #Compute the mask features
        Position_morphology_mean <- try(
          cytomapper::measureObjects(mask = Seg_results,
                                     image = Image,
                                     img_id = "imageID",
                                     feature_types = c("basic", "shape", "moment"),
                                     shape_feature = c("area", "perimeter", "radius.mean", "radius.sd", "radius.max", "radius.min"),
                                     moment_feature = c("cx", "cy", "eccentricity", "majoraxis"),
                                     basic_feature = "mean")
        )

        Morphology_tibble <- try(
          #Generate a tibble with the mask features
          tibble::as_tibble(cbind(as_tibble(SummarizedExperiment::colData(Position_morphology_mean)),
                                  as_tibble(t(SummarizedExperiment::assays(Position_morphology_mean)[[1]])))) %>% dplyr::select(-objectNum)
        )
        try(names(Morphology_tibble)[-c(1:12)] <- stringr::str_c(names(Morphology_tibble)[-c(1:12)], "_AVERAGE"))


        #Check errors and proceed if required
        if(berryFunctions::is.error(Seg_results)){
          #Print message and remove image
          print(paste0("Unable to compute structure mask for ", Random_image_name, ". Trying another random sample."))
          Images_for_test_remaining <- Images_for_test_remaining[-which(Images_for_test_remaining == Random_image_name)]

        } else if(max(Seg_results[[1]]) == 0){
          #Print message and remove image
          print(paste0("No structures found in ", Random_image_name, ". Trying another random sample."))
          Images_for_test_remaining <- Images_for_test_remaining[-which(Images_for_test_remaining == Random_image_name)]

        } else if(berryFunctions::is.error(Morphology_tibble)){
          #Print message and remove image
          print(paste0("No structures found in ", Random_image_name, ". Trying another random sample."))
          Images_for_test_remaining <- Images_for_test_remaining[-which(Images_for_test_remaining == Random_image_name)]
        } else Segmentation_mask_ok <- "PROCEED"
      }

      #Modify the results to match the original dataset
      #The dataset is usually reversed and flipped
      X_centroids <- Morphology_tibble$m.cx
      Y_centroids <- Morphology_tibble$m.cy
      Morphology_tibble <- Morphology_tibble %>% dplyr::mutate(m.cx = Y_centroids, m.cy = X_centroids)
      Morphology_tibble <- Morphology_tibble %>%
        mutate(X_final = -m.cx,
               Y_final = -m.cy)

      #Now map the X Y points to the original cell coordinates
      Morphology_tibble$X_final <- Morphology_tibble$X_final + dim(Pseudo_image)[[2]]
      Morphology_tibble$X_final <- Morphology_tibble$X_final/dim(Pseudo_image)[[2]]
      Morphology_tibble$X_final <- min(Example_image_tiles$tile_xmin) +
        (max(Example_image_tiles$tile_xmax) - min(Example_image_tiles$tile_xmin)) * Morphology_tibble$X_final

      Morphology_tibble$Y_final <- Morphology_tibble$Y_final + dim(Pseudo_image)[[1]]
      Morphology_tibble$Y_final <- Morphology_tibble$Y_final/dim(Pseudo_image)[[1]]
      Morphology_tibble$Y_final <- min(Example_image_tiles$tile_ymin) +
        (max(Example_image_tiles$tile_ymax) - min(Example_image_tiles$tile_ymin)) * Morphology_tibble$Y_final

      #Generate the plots
      Mask_plot <-
        ggplot() +
        annotation_raster(EBImage::colorLabels(Seg_results[[1]]) %>% EBImage::rotate(angle = -90) %>% EBImage::flip() %>% EBImage::flop(),
                          xmin = min(Example_image_tiles$tile_xmin), xmax = max(Example_image_tiles$tile_xmax),
                          ymin = min(Example_image_tiles$tile_ymin), ymax = max(Example_image_tiles$tile_ymax)) +


        geom_point(aes(x = X_final, y =  Y_final), shape = 21, fill = "white", color = "black", size = 3, data = Morphology_tibble)+
        scale_x_continuous(limits = c(min(Example_image_tiles$tile_xmin), max(Example_image_tiles$tile_xmax))) +
        scale_y_continuous(limits = c(min(Example_image_tiles$tile_ymin), max(Example_image_tiles$tile_ymax))) +
        ggtitle("Structure mask", Random_image_name) +
        theme(aspect.ratio =
                (max(Example_image_tiles$tile_xmax) - min(Example_image_tiles$tile_xmin)) /
                (max(Example_image_tiles$tile_ymax) - min(Example_image_tiles$tile_ymin))
        )


      Tiles_plot <-
        Example_image_tiles %>%
        ggplot() +
        geom_rect(aes(xmin = tile_xmin, xmax = tile_xmax, ymin = tile_ymin, ymax = tile_ymax, fill = as.factor(Target_column))) +
        scale_fill_manual("", labels = c("Negative", "Positive"), values = c("black", "white")) +
        geom_point(aes(x = X_final, y =  Y_final), shape = 21,  fill = "red", color = "white", size = 3, data = Morphology_tibble) +
        ggtitle("Tile pattern", Random_image_name) +
        theme(aspect.ratio =
                (max(Example_image_tiles$tile_xmax) - min(Example_image_tiles$tile_xmin)) /
                (max(Example_image_tiles$tile_ymax) - min(Example_image_tiles$tile_ymin))
        )

      plot(Tiles_plot + Mask_plot + patchwork::plot_layout(nrow = 1, ncol = 2))


      Test_OK <- menu(choices = c("Proceed", "Abort", "Test again"), title = "Check the structures identified by the algorithm")
      if(Test_OK == 2) stop("Analysis aborted")
    }

    ########Perform the actual computation########
    #Set the on exit statement
    on.exit({
      future::plan("future::sequential")
      gc()
    })

    #Set the cores
    future::plan("future::multisession", workers = N_cores)
    options(future.globals.maxSize = Inf, future.rng.onMisuse = "ignore")
    furrr::furrr_options(scheduling = Inf)

    #Execute the algorithm
    FINAL_RESULT_LIST <-
      furrr::future_map(Tiled_Images, function(Image){

        Example_image_tiles <- Image

        #Compute the tile size
        Tile_size <- unique(Image$tile_xmax - Image$tile_xmin)

        #Calculate the number of rows and columns
        x_length <- length(unique(Example_image_tiles$tile_xmin))#Calculte the X pixel length
        y_length <- length(unique(Example_image_tiles$tile_ymin))#Calculate the Y pixel length

        #Arrange the image according to the X centroid
        Interim <- Example_image_tiles %>% dplyr::arrange(tile_xmin, desc(tile_ymin))#arrange the tibble in an adequate format

        #Compute the pseudo-image (a matrix)
        Pseudo_image <- matrix(Interim[["Target_column"]], nrow = y_length, ncol = x_length)#Calculate the matrix
        colnames(Pseudo_image) <- unique(Interim$tile_xmin)
        rownames(Pseudo_image) <- unique(Interim$tile_ymin)

        #Turn the pseudo-image to an EBImage and Cytomapper object
        Pseudo_image <- EBImage::Image(Pseudo_image) %>% EBImage::flip()
        Image <- cytomapper::CytoImageList(Pseudo_image)
        cytomapper::channelNames(Image) <- "COUNTS" #Modify the channel name
        S4Vectors::mcols(Image)$imageID <- "Pseudo_Image"#Modify name

        #Compute the structure mask and the object morphometry in a secure manner
        Seg_results <- try(simpleSeg::simpleSeg(Image,
                                                nucleus = "COUNTS",
                                                cellBody = "none",
                                                sizeSelection = Min_tile_number,
                                                smooth = Smooth_amount,
                                                transform = NULL,
                                                watershed = "intensity",
                                                tolerance = Tolerance_value,
                                                ext = Neighborhood_distance,
                                                discSize = 1,
                                                tissue = NULL,
                                                pca = FALSE,
                                                cores = 1)
        )


        try(
          #Change the ID of the mask
          S4Vectors::mcols(Seg_results)$imageID <- "Pseudo_Image"
        )

        #Compute the mask features
        Position_morphology_mean <- try(
          cytomapper::measureObjects(mask = Seg_results,
                                     image = Image,
                                     img_id = "imageID",
                                     feature_types = c("basic", "shape", "moment"),
                                     shape_feature = c("area", "perimeter", "radius.mean", "radius.sd", "radius.max", "radius.min"),
                                     moment_feature = c("cx", "cy", "eccentricity", "majoraxis", "theta"),
                                     basic_feature = "mean")
        )

        Morphology_tibble <- try(
          #Generate a tibble with the mask features
          tibble::as_tibble(cbind(as_tibble(SummarizedExperiment::colData(Position_morphology_mean)),
                                  as_tibble(t(SummarizedExperiment::assays(Position_morphology_mean)[[1]])))) %>% dplyr::select(-objectNum)
        )
        try(names(Morphology_tibble)[-c(1:12)] <- stringr::str_c(names(Morphology_tibble)[-c(1:12)], "_AVERAGE"))


        #Check errors and proceed if required
        if(berryFunctions::is.error(Seg_results)){
          #return message
          return("Unable to perform morphometry analysis")

        } else if(max(Seg_results[[1]]) == 0){
          #return message
          return("Unable to perform morphometry analysis")

        } else if(berryFunctions::is.error(Morphology_tibble)){
          #return message
          return("Unable to perform morphometry analysis")
        } else {

          #Modify the results to match the original dataset
          #The dataset is usually reversed and flipped
          X_centroids <- Morphology_tibble$m.cx
          Y_centroids <- Morphology_tibble$m.cy
          Morphology_tibble <- Morphology_tibble %>% dplyr::mutate(m.cx = Y_centroids, m.cy = X_centroids)
          Morphology_tibble <- Morphology_tibble %>%
            mutate(X_final = -m.cx,
                   Y_final = -m.cy)

          #Now map the X Y points to the original cell coordinates
          Morphology_tibble$X_final <- Morphology_tibble$X_final + dim(Pseudo_image)[[2]]
          Morphology_tibble$X_final <- Morphology_tibble$X_final/dim(Pseudo_image)[[2]]
          Morphology_tibble$X_final <- min(Example_image_tiles$tile_xmin) +
            (max(Example_image_tiles$tile_xmax) - min(Example_image_tiles$tile_xmin)) * Morphology_tibble$X_final

          Morphology_tibble$Y_final <- Morphology_tibble$Y_final + dim(Pseudo_image)[[1]]
          Morphology_tibble$Y_final <- Morphology_tibble$Y_final/dim(Pseudo_image)[[1]]
          Morphology_tibble$Y_final <- min(Example_image_tiles$tile_ymin) +
            (max(Example_image_tiles$tile_ymax) - min(Example_image_tiles$tile_ymin)) * Morphology_tibble$Y_final

          #Select desired variables
          Morphology_tibble <- Morphology_tibble %>% dplyr::select(object_id, X_final, Y_final,
                                                                   s.area, s.perimeter, s.radius.mean, s.radius.sd, s.radius.min, s.radius.max,
                                                                   m.majoraxis, m.eccentricity, m.theta_AVERAGE)
          #Change names
          names(Morphology_tibble)[1:3] <- c("Structure_ID", "X_centroid", "Y_centroid")

          #Modify variables to pixel/distance metric and not tile based
          Morphology_tibble <- Morphology_tibble %>%
            dplyr::mutate(Cell_analyzed = Feature_to_analyze,
                          s.area = s.area * (Tile_size^2), #To account for tile area
                          s.perimeter = s.perimeter * Tile_size,
                          s.radius.mean = s.radius.mean * Tile_size,
                          s.radius.sd = s.radius.sd * Tile_size,
                          s.radius.min = s.radius.min * Tile_size,
                          s.radius.max = s.radius.max * Tile_size,
                          m.majoraxis = m.majoraxis * Tile_size)
          return(list(Morphology_tibble = Morphology_tibble,
                      Structure_mask = Seg_results[[1]] %>% EBImage::rotate(angle = -90) %>% EBImage::flip() %>% EBImage::flop()
                      )
                 )
        }
      },
      .progress = TRUE)

    future::plan("future::sequential")
    gc()

    ####REMOVE CONFLICTIVE IMAGES and RETURN result####
    #Remove images without adequate analysis
    Conflictive_images <-
      purrr::map_lgl(FINAL_RESULT_LIST, function(Image){

        if(length(Image[[1]]) == 1){ #If length of first element is 1 then test if the message
          if(Image[[1]] == "Unable to perform morphometry analysis") return(TRUE) # Test if the message is the expected one

        } else return(FALSE) #If both conditions are not met then return FALSE
      })

    #If there is any image without a result then remove it before exiting the function
    if(any(Conflictive_images)){
      Conflictive_names <- names(FINAL_RESULT_LIST)[Conflictive_images]

      Colored_print(paste0("\nThe following images could not be analyzed by the function and will be removed from the final result: ",
                           "\n",
                           stringr::str_c(Conflictive_names, collapse = "\n")),
                    color = "red")

      FINAL_RESULT_LIST <- FINAL_RESULT_LIST[!Conflictive_images]

    }

    #Return the final result
    return(FINAL_RESULT_LIST)
  }



