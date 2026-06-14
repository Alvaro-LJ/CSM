#' Imports large images to the environment
#'
#' The function is aimed to import large images or image subsets to the environment for graphication purposes (generate images that can be easily plotted). 
#'
#' @param N_cores Integer. Number of cores to parallelize your computation.
#' @param Image_directory Character specifying the path to the image file that needs to be imported.
#' @param Log10_pixel_output A single numeric value larger than 1. It indicates the size (in pixels) of the processed image. For example, if set to 6 (the default) the output image will have a size of 10^6 pixels.
#' @param Percentage_downsize (OPTIONAL) Alternatively, the user can provide the downsize factor as the percentage of original pixel number. Used only if Log10_pixel_output is NULL.
#' @param X_crop_coordinates (OPTIONAL) If only a subset of the image needs to be imported, a numeric vector of length = 2, indicating the min and max X coordinates to perform image subsetting.
#' @param Y_crop_coordinates (OPTIONAL) If only a subset of the image needs to be imported, a numeric vector of length = 2, indicating the min and max Y coordinates to perform image subsetting.
#' @param Ordered_Channels (OPTIONAL) If image has multiple channel, the user can provide channel names. These names will be used with Channels_to_keep (see below) to filter the resulting image.
#' @param Channels_to_keep (OPTIONAL) If image has multiple channel, the user can provide channel names. These names will be used with Channels_to_keep to filter the resulting image.
#' @param Save_processed_images A logical value indicating if the resulting image should be written in disk.
#' @param Output_Directory If Save_processed_images is TRUE, the path to the folder where the resulting image should be written.
#'
#'
#' @returns A list containing the resulting image, the native image pixel dimensions, the processed image pixel dimensions and the subsetting coordinates (if required).
#'
#' @details
#' The function is able to read tif files. It can also work with pyramidal tif images. The `Smart_image_importer()` is designed to import any type of image irrespective of it's size. In order to 
#' work with large images, the user needs to set a reasonable pixel size limit. Large tif images will be imported in tiles, downsized and then reconstructed to give rise to the output image. For pyramidal
#' tif files the resolution layer demonstrating the closest pixel size to the desired pixel-limit will be returned. If images are not exceedingly large (less than 1.5Gb), they will be imported and directly 
#' downsized. 
#' 
#' Tile importing is powered by the RBioFormats Bioconductor package that parses the Bioformats Java library to R. Sometimes even if the user#' has installed the RBioFormats package from 
#' bioconductor the function will still return an error. If this happens, the issue can be solved checking the installed Java version in the computer and manually installing the CRAN rJava package.
#' If the function returns an error due to RAM saturation, sometimes it can be solved by allowing Java to use an absurd amount of RAM (for example rJava::.jinit(parameters="-Xmx500g")).
#'
#'
#' @examples
#' \dontrun{
#' #Create temporary input directory------------------------------
#' Input_Dir <- tempfile(pattern = "tempdir1_Input")
#' dir.create(Input_Dir, recursive = TRUE)
#'
#' #Save images in Input directory
#' purrr::map(1:2,
#' function(Image){
#'    EBImage::writeImage(CSM_MiniMultiTiff_test[[Image]],
#'    file.path(Input_Dir, names(CSM_MiniMultiTiff_test)[Image]))
#' })
#' 
#' #reduce the pixel content of the image
#' Smart_image_importer(Image_directory = dir(Input_Dir, full.names = TRUE)[1], 
#'                      Log10_pixel_output = 5
#' ) 
#'
#'}
#' @export


Smart_image_importer <- 
  function(
    N_cores = 1,
    Image_directory,
    Log10_pixel_output = 6,
    Percentage_downsize = NULL,
    X_crop_coordinates = NULL,
    Y_crop_coordinates = NULL,
    Ordered_Channels = NULL,
    Channels_to_keep = NULL,
    Save_processed_images = FALSE,
    Output_Directory = NULL
  ){
    
    #####WHAT TO DO ON EXIT#####
    on.exit({
      future::plan("future::sequential")
      gc()
    }
    )
    
    #####CHECK SUGGESTED PACKAGES####
    #Check installation of suggested packages
    {
      if(!requireNamespace("magick", quietly = FALSE)) stop(
        paste0("magick CRAN package is required to execute the function. Please install using the following code: ",
               expression(install.packages("magick"))
               )
        )
      if(!requireNamespace("RBioFormats", quietly = TRUE)) stop(
        paste0("RBioFormats Bioconductor package is required to execute the function. Please install using the following code: ",
               expression({
                 if (!require("BiocManager", quietly = TRUE))
                   install.packages("BiocManager")
                 
                 BiocManager::install("RBioFormats")
               })
        )
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
    }
    
    #####ARGUMENT CHECK and basic data obtention#####
    #Check the cores
    if(!all(N_cores >= 1 & N_cores%%1 == 0)) stop("N_cores must be an integer value > 0")
    
    #Image directory
    Image_name <- stringr::str_extract(Image_directory, "[^/]+$")
    Image_extension <- stringr::str_extract(Image_name, "[^.]+$")
  
    #Get the image type
    Image_type <- "Unkwown"
    if(stringr::str_detect(Image_extension, stringr::regex("tiff?", ignore_case = TRUE))) Image_type <- "TIFF"
    if(Image_type == "Unkwown") stop(paste0(Image_name, " has an unkwown file extension. Only tif or OME-tiff or pyramidal tiff are supported"))
    
    #Get the Image size and metadata
    Image_size <- file.info(Image_directory)$size / 1024^3 #Get the image size in Gb
    Image_metadata <- RBioFormats::read.metadata(Image_directory, filter.metadata = TRUE)

    #Check if is pyramidal TIFF with the following check
    if(!is.null(unlist(purrr::map(Image_metadata@.Data, function(Metadata_item) names(Metadata_item$coreMetadata))))) Image_type <- "PYRAMIDAL"
    
    #Check if image is RGB
    Is_RGB <- FALSE
    try({
      #If image type is common tiff
      if(Image_type == "TIFF"){
        if(Image_metadata@.Data[[2]]$PhotometricInterpretation == "RGB") Is_RGB <- TRUE
      }
      #If image type is pyramidal tiff
      if(Image_type == "PYRAMIDAL"){
        if(Image_metadata@.Data[[1]][[2]]$PhotometricInterpretation == "RGB") Is_RGB <- TRUE
      }
    }, silent = TRUE)
      
    #If is RGB and channels have been substracted then stop the function
    if(Is_RGB & any(!is.null(Ordered_Channels), !is.null(Channels_to_keep))) stop("If image is RGB channels can not be subsetted")
    
    #Log10_pixel_output or Percentage_downsize
    #Check LOG10_pixel_output if required
    if(!is.null(Log10_pixel_output)){
      if(!all(length(Log10_pixel_output) == 1, is.numeric(Log10_pixel_output), Log10_pixel_output > 1)) stop("Log10_pixel_output must be a single numeric value higher than 1")
    }
    #If both LOG10 and percentage provided print a message
    if(all(!is.null(Log10_pixel_output), !is.null(Percentage_downsize))) message("\nBoth Log10_pixel_output and Percentage_downsize provided. Percentage_downsize will be ignored.")
    
    #If only Percentage_downsize provided compute LOG10_pixel_output
    if(all(is.null(Log10_pixel_output), !is.null(Percentage_downsize))){
      
      if(!all(length(Percentage_downsize) == 1, is.numeric(Percentage_downsize), Percentage_downsize > 0, Percentage_downsize <= 1)) stop("Percentage_downsize must be a numeric value between 0 and 1")
      
      #If pyramidal
      if(Image_type == "PYRAMIDAL"){
        #Obtain the info from the layer with the highest resolution of get the total pixel number\
        Max_Layer_data <- Image_metadata@.Data[[1]]$coreMetadata
        Total_pixel_number <- Max_Layer_data$sizeX * Max_Layer_data$sizeY
      }
      #If only TIFF
      else Total_pixel_number <- Image_metadata$coreMetadata$sizeX * Image_metadata$coreMetadata$sizeY
      
      #Compute the LOG_10_pixel_output
      Log10_pixel_output <- round(log10(Total_pixel_number * Percentage_downsize), digits = 3)
      
    }
    
    #Crop coordinates
    if(any(!is.null(X_crop_coordinates), !is.null(Y_crop_coordinates))){
      if(!all(!is.null(X_crop_coordinates), !is.null(Y_crop_coordinates))) stop("Both X_crop_coordinates and Y_crop_coordinates must be provided")
      if(!all(length(X_crop_coordinates) == 2, length(Y_crop_coordinates) == 2)) stop("Two coordinates must be provided for both X_crop_coordinates and Y_crop_coordinates")
      if(!all(is.numeric(X_crop_coordinates), is.numeric(Y_crop_coordinates))) stop("X_crop_coordinates and Y_crop_coordinates must be numeric vectors")
      if(!all(c(X_crop_coordinates, Y_crop_coordinates) > 0)) stop("X_crop_coordinates and Y_crop_coordinates must be higher than 0")
      if(!all(c(X_crop_coordinates, Y_crop_coordinates)%%1 == 0)) stop("X_crop_coordinates and Y_crop_coordinates must be numeric integer values")
      if(X_crop_coordinates[2] <= X_crop_coordinates[1]) stop("The first pixel coordinate of X_crop_coordinates must be smaller than the second pixel coordinate")
      if(Y_crop_coordinates[2] <= Y_crop_coordinates[1]) stop("The first pixel coordinate of Y_crop_coordinates must be smaller than the second pixel coordinate")
    }
    
    #Channels (if provided)
    if(any(!is.null(Ordered_Channels), !is.null(Channels_to_keep)) & !Is_RGB){
      
      #Check that they are both provided
      if(!all(!is.null(Ordered_Channels), !is.null(Channels_to_keep))) stop("Ordered_Channels and Channels_to_keep must be provided")
      
      #Check that ordered channels and channels to keep are unique
      if(length(Ordered_Channels) != length(unique(Ordered_Channels))) stop("Ordered_Channels must contain non-duplicated character values")
      if(length(Channels_to_keep) != length(unique(Channels_to_keep))) stop("Channels_to_keep must contain non-duplicated character values")
      
      #Check that ordered channels have the same length as channels in image
      if(Image_type == "TIFF") Number_channels <-  max(c(Image_metadata$coreMetadata$sizeZ, Image_metadata$coreMetadata$sizeC, Image_metadata$coreMetadata$imageCount))
      if(Image_type == "PYRAMIDAL") Number_channels <-  max(c(Image_metadata@.Data[[1]]$coreMetadata$sizeZ , Image_metadata@.Data[[1]]$coreMetadata$sizeC, Image_metadata@.Data[[1]]$coreMetadata$imageCount))
      if(length(Ordered_Channels) != Number_channels) stop(paste0("Mismatch in the number of Ordered_Channels provided. ", 
                                                                  length(Ordered_Channels), " channel names provided but ",
                                                                  Number_channels, " identified in image"))
      
      #Check that Channels_to_keep is contained in Ordered_Channels
      if(!all(Channels_to_keep %in% Ordered_Channels)){
        Problematic_channels <- Channels_to_keep[!Channels_to_keep %in% Ordered_Channels]
        stop(paste0("The following Channels_to_keep are not present in Ordered_channels:\n", stringr::str_c(Problematic_channels, collapse = "\n")))
      }
      
      #Generate the channel_to_keep_index
      Channels_to_keep_index <- match(Channels_to_keep, Ordered_Channels)
    }
    
    #Save processed image and directory provided
    if(!all(length(Save_processed_images) == 1, is.logical(Save_processed_images))) stop("Save_processed_images must be a logical value")
    if(Save_processed_images){
      if(!dir.exists(Output_Directory)) stop("Directory provided does not work. Please review")
    }
    
    
    ####Pyramidal (OME) TIFF####
    if(Image_type == "PYRAMIDAL"){
      #Get the number of layers and the resolution for each layer
      #Compute the number of total pixels per layer and the approx size of each layer
      Layer_information <- 
        purrr::map_dfr(1:length(Image_metadata@.Data), function(Resolution_layer_index){
          Layer_data <- Image_metadata@.Data[[Resolution_layer_index]]$coreMetadata
          
          tibble::tibble(resolutionLevel = Layer_data$resolutionLevel,
                         Image_Xmax = Layer_data$sizeX,
                         Image_Ymax = Layer_data$sizeY,
                         Total_pixel_count = Layer_data$sizeX*Layer_data$sizeY,
                         Total_channels = max(Layer_data$sizeC, Layer_data$sizeZ, Layer_data$imageCount)
          )
        })
      Layer_information$Layer_Ratio <- Layer_information$Total_pixel_count / max(Layer_information$Total_pixel_count)
      Layer_information$Layer_approx_size <- Image_size * (Layer_information$Total_pixel_count / sum(Layer_information$Total_pixel_count))
      
      #Get the basic data
      Image_Xmax <- Layer_information$Image_Xmax[1]
      Image_Ymax <- Layer_information$Image_Ymax[1]
      
      
      #COMPLETE IMAGE, NO SUBSETTING
      if(all(is.null(X_crop_coordinates), is.null(Y_crop_coordinates))){
        
        #To increase JAVA RAM use
        rJava::.jinit(parameters="-Xmx200g")
        
        #Compute which layer is closer to the required pixel amount
        Selected_layer <- Layer_information$resolutionLevel[which.min(abs(Layer_information$Total_pixel_count - 10^Log10_pixel_output))]
        Candidate_layer_data <- Layer_information %>% dplyr::filter(resolutionLevel >= Selected_layer)
        
        Uncomplete_processing <- TRUE #Set the status of the image
        #Generate a while loop, the function will try to return the closest possible layer
        while(Uncomplete_processing){
          #Select the resolution level
          Resolution_layer_being_tested <- Candidate_layer_data$resolutionLevel[1]
          
          cat(paste0("\nAttempting to obtain resolution level ", Resolution_layer_being_tested))
          #Attempt to obtain image in that resolution level
          try({
            Final_image <- RBioFormats::read.image(Image_directory,
                                                   proprietary.metadata = FALSE,
                                                   read.metadata = TRUE,
                                                   normalize = TRUE, 
                                                   resolution = Resolution_layer_being_tested
                                                   )
            #Change the color mode
            if(!Is_RGB) EBImage::colorMode(Final_image) <- "Grayscale"
            if(Is_RGB) EBImage::colorMode(Final_image) <- "Color"
            
            #Subset channels if required
            if(any(!is.null(Ordered_Channels), !is.null(Channels_to_keep))) Final_image <- Final_image[,,Channels_to_keep_index]
            #Obtain the reduction factor (the layer ratio compared to the maximum layer info)
            Reduction_factor <- Candidate_layer_data$Layer_Ratio[1]
            
            Final_image <- magick::image_read(Final_image)
            Uncomplete_processing <- FALSE
          })
          
          Candidate_layer_data <- Candidate_layer_data[-1,]
          if(nrow(Candidate_layer_data) == 0) stop("All pyramidal tif layers have been attempted to be extracted without success")
        }
      }
      
      #SUBSETTING REQUIRED
      if(any(!is.null(X_crop_coordinates), !is.null(Y_crop_coordinates))){
        
        #If subset coordinates are out of bounds stop the computation
        if(any(X_crop_coordinates[2] > Image_Xmax, Y_crop_coordinates[2] > Image_Ymax)) stop(paste0("Max X_crop_coordinates and y_crop_coordinates allowed are ", Image_Xmax, " and ", Image_Ymax, ", respectively"))
        
        #To increase JAVA RAM use
        rJava::.jinit(parameters="-Xmx200g")
        
        #Compute the pixel size for the subset in all of the layers
        Absolute_Subset_pixel_count <- (X_crop_coordinates[2] - X_crop_coordinates[1]) * (Y_crop_coordinates[2] - Y_crop_coordinates[1])
        Layer_information$Subset_pixel_count <- ceiling(Layer_information$Layer_Ratio * Absolute_Subset_pixel_count)
        
        #Compute which layer is closer to the required pixel amount
        Selected_layer <- Layer_information$resolutionLevel[which.min(abs(Layer_information$Subset_pixel_count - 10^Log10_pixel_output))]
        Candidate_layer_data <- Layer_information %>% dplyr::filter(resolutionLevel >= Selected_layer)
        Candidate_layer_data$X_crop_coordinates_min <- ceiling(X_crop_coordinates[1] * sqrt(Candidate_layer_data$Layer_Ratio)) #Sqroot to account for the two dimensions
        Candidate_layer_data$X_crop_coordinates_max <- ceiling(X_crop_coordinates[2] * sqrt(Candidate_layer_data$Layer_Ratio))
        Candidate_layer_data$Y_crop_coordinates_min <- ceiling(Y_crop_coordinates[1] * sqrt(Candidate_layer_data$Layer_Ratio))
        Candidate_layer_data$Y_crop_coordinates_max <- ceiling(Y_crop_coordinates[2] * sqrt(Candidate_layer_data$Layer_Ratio))
        
        Uncomplete_processing <- TRUE #Set the status of the image
        #Generate a while loop, the function will try to return the closest possible layer
        while(Uncomplete_processing){
          #Select the resolution level
          Resolution_layer_being_tested <- Candidate_layer_data$resolutionLevel[1]
          
          cat(paste0("\nAttempting subsetting from resolution level ", Resolution_layer_being_tested))
          #Attempt to obtain image in that resolution level
          try({
            Final_image <- RBioFormats::read.image(Image_directory,
                                                   proprietary.metadata = FALSE,
                                                   read.metadata = TRUE,
                                                   normalize = TRUE, 
                                                   resolution = Resolution_layer_being_tested,
                                                   subset = list(X = c(Candidate_layer_data$X_crop_coordinates_min[1]:Candidate_layer_data$X_crop_coordinates_max[1]),
                                                                 Y = c(Candidate_layer_data$Y_crop_coordinates_min[1]:Candidate_layer_data$Y_crop_coordinates_max[1]))
            )
            #Change color mode as required
            if(!Is_RGB) EBImage::colorMode(Final_image) <- "Grayscale"
            if(Is_RGB) EBImage::colorMode(Final_image) <- "Color"
            
            #Subset channels if required
            if(any(!is.null(Ordered_Channels), !is.null(Channels_to_keep))) Final_image <- Final_image[,,Channels_to_keep_index]
            
            #Obtain the reduction factor (the subset pixel count in the selected layer compared to the absolute subset pixel count)
            Reduction_factor <- Candidate_layer_data$Subset_pixel_count[1] / Absolute_Subset_pixel_count
            
            Final_image <- magick::image_read(Final_image)
            Uncomplete_processing <- FALSE
          })
          
          Candidate_layer_data <- Candidate_layer_data[-1,]
          if(nrow(Candidate_layer_data) == 0) stop("All pyramidal tif layers have been attempted to be subsetted without success")
        }
          
      }
    }
    
    ####TIFF####
    if(Image_type == "TIFF"){
      #Basic image metadata obtained for non_piramidal TIFF
      Image_Xmax <- Image_metadata$coreMetadata$sizeX
      Image_Ymax <- Image_metadata$coreMetadata$sizeY
      
      #Compute the reduction factor
      Reduction_factor <- (10^Log10_pixel_output) / (Image_Xmax * Image_Ymax)
      #If reduction factor is equal or larger than 1 and there is no need to Crop, then stop the function because it makes no sense
      if(Reduction_factor >= 1 & is.null(X_crop_coordinates)){
        message("Image already has lower pixel count than indicated. No image resizing will be attempted")
        Reduction_factor <- 1
      }
      
      #If reduction factor is equal or larger than 1 and the user wants to crop, set the reduction factor to 1 ()
      if(Reduction_factor >= 1 & !is.null(X_crop_coordinates)) Reduction_factor <- 1
      
      #COMPLETE IMAGE, NO SUBSETTING
      if(all(is.null(X_crop_coordinates), is.null(Y_crop_coordinates))){
        
        Uncomplete_processing <- TRUE #Set the status of the image
        
        #If the image size is smaller than 1.5 Gb attempt to reduce it using magick (It can deal with up to 3Gb according to benchmark)
        if(Image_size <= 1.5){
          cat("\nAttempting direct image resizing...")
          try({
            Final_image <- magick::image_read(Image_directory)
            
            #Keep desired channels only
            if(any(!is.null(Ordered_Channels), !is.null(Channels_to_keep))) Final_image <- Final_image[Channels_to_keep_index]
            
            #Apply reduction
            Final_width <- ceiling(Image_Xmax * sqrt(Reduction_factor))
            Final_height <- ceiling(Image_Ymax * sqrt(Reduction_factor))
            Final_image <- magick::image_resize(Final_image, magick::geometry_size_pixels(width = Final_width, height = Final_height))
            Uncomplete_processing <- FALSE #Turn to FALSE
          })
        }
        
        #If Image_size is larger or not completely processed (failed computation)
        if(Image_size > 1.5 | Uncomplete_processing){
          cat("\nAttempting resizing by tiling...")
          #Compute tile size
          N_rows_columns <- ceiling(sqrt(Image_size/0.5)) #Each tile is aimed to be 0.5Gb of size (OPTIMAL BASED ON BENCHMARK)
          
          #Compute the tile coordinates
          #Define the X and Y positions to visit
          X_Positions <- seq(from = 1, to = Image_Xmax, by = ceiling(Image_Xmax/N_rows_columns))
          Y_Positions <- seq(from = 1, to = Image_Ymax, by = ceiling(Image_Ymax/N_rows_columns))
          
          #If the Max X and Y is equal to the Image MAX then remove the final X or Y position
          if(X_Positions[length(X_Positions)] == Image_Xmax) X_Positions <- X_Positions[-length(X_Positions)]
          if(Y_Positions[length(Y_Positions)] == Image_Ymax) Y_Positions <- Y_Positions[-length(Y_Positions)]
          
          Tile_Position_tibble <- tidyr::expand_grid(X_Positions, Y_Positions)
          names(Tile_Position_tibble) <- c("X_Position_Min", "Y_Position_Min")
          Tile_Position_tibble$X_Position_Max <- Tile_Position_tibble$X_Position_Min + ceiling(Image_Xmax/N_rows_columns)-1
          Tile_Position_tibble$Y_Position_Max <- Tile_Position_tibble$Y_Position_Min + ceiling(Image_Ymax/N_rows_columns)-1
          Tile_Position_tibble$Tile_ID <- stringr::str_c("Tile", 1:nrow(Tile_Position_tibble))
          #If the final X or Y position is larger than the Xmax or Ymax of the image make the ammendments
          Tile_Position_tibble$X_Position_Max[Tile_Position_tibble$X_Position_Max > Image_Xmax] <- Image_Xmax
          Tile_Position_tibble$Y_Position_Max[Tile_Position_tibble$Y_Position_Max > Image_Ymax] <- Image_Ymax
          
          #Generate the column and the row index
          Col_index_tibble <- tibble::tibble(X_Position_Min = unique(Tile_Position_tibble$X_Position_Min),
                                             Column_index = stringr::str_c("Col_", 1:length(unique(Tile_Position_tibble$X_Position_Min))))
          
          Row_index_tibble <- tibble::tibble(Y_Position_Min = unique(Tile_Position_tibble$Y_Position_Min),
                                             Row_index = stringr::str_c("Row_", 1:length(unique(Tile_Position_tibble$Y_Position_Min))))
          #Bind to the tile_Position_tibble
          Tile_Position_tibble <- dplyr::left_join(Tile_Position_tibble, Col_index_tibble, by = "X_Position_Min") 
          Tile_Position_tibble <- dplyr::left_join(Tile_Position_tibble, Row_index_tibble, by = "Y_Position_Min") 
          
          
          #To increase JAVA RAM use
          rJava::.jinit(parameters="-Xmx200g") 
          #Compute all the tiles
          cat("\nTiling the image")
          
          #Plan the multicore session
          future::plan("future::multisession", workers = N_cores)
          options(future.globals.maxSize = Inf, future.rng.onMisuse = "ignore")
          furrr::furrr_options(scheduling = Inf)
          
          Tiles <- 
            furrr::future_map(1:nrow(Tile_Position_tibble), function(Image_index){
              rJava::.jinit(parameters="-Xmx200g") #To increase JAVA RAM use
              Tile_x_min <- Tile_Position_tibble$X_Position_Min[[Image_index]]
              Tile_x_max <- Tile_Position_tibble$X_Position_Max[[Image_index]]
              
              Tile_y_min <- Tile_Position_tibble$Y_Position_Min[[Image_index]]
              Tile_y_max <- Tile_Position_tibble$Y_Position_Max[[Image_index]]
              
              library(RBioFormats)
              
              Tile_image <- RBioFormats::read.image(Image_directory,
                                                    proprietary.metadata = FALSE,
                                                    read.metadata = TRUE,
                                                    normalize = TRUE,
                                                    subset = list(X = c(Tile_x_min:Tile_x_max),
                                                                  Y = c(Tile_y_min:Tile_y_max))
              )
              
              #Change color mode as required
              if(!Is_RGB) EBImage::colorMode(Tile_image) <- "Grayscale"
              if(Is_RGB) EBImage::colorMode(Tile_image) <- "Color"
              
              #Keep desired channels only
              if(any(!is.null(Ordered_Channels), !is.null(Channels_to_keep))) Tile_image <- Tile_image[,,Channels_to_keep_index]
              
              #Reduce the amount required
              Tile_image <- EBImage::resize(Tile_image, w = ceiling(dim(Tile_image)[1]*sqrt(Reduction_factor)), h = ceiling(dim(Tile_image)[2]*sqrt(Reduction_factor)))#Sqroot to account for double square reduction
              
              return(Tile_image)
              
            },
            .progress = TRUE)
          future::plan("future::sequential")
          gc()
          
          
          #Build the columns
          cat("\nMerging tiles")
          Cols <- purrr::map(unique(Tile_Position_tibble$Column_index), function(Current_column_index){
            Tile_index <- Tile_Position_tibble$Column_index == Current_column_index
            
            Col_image <- EBImage::abind(Tiles[Tile_index], along = 2)
            return(Col_image)
          })
          
          #Bind the columns
          Final_image <- EBImage::abind(Cols, along = 1)
          
          #Remove objects manually
          rm(Cols)
          rm(Tiles)
          gc()
          
          #Return the final image
          Final_image <- magick::image_read(Final_image)
        }
      }
      
      #SUBSETTING REQUIRED
      if(any(!is.null(X_crop_coordinates), !is.null(Y_crop_coordinates))){
        
        #If subset coordinates are out of bounds stop the computation
        if(any(X_crop_coordinates[2] > Image_Xmax, Y_crop_coordinates[2] > Image_Ymax)) stop(paste0("Max X_crop_coordinates and y_crop_coordinates allowed are ", Image_Xmax, " and ", Image_Ymax, ", respectively"))
        
        #We will check if the pixel size of the subset is within the reduction factor
        Total_pixel_count <- Image_Xmax * Image_Ymax
        Subset_pixel_count <- (X_crop_coordinates[2] - X_crop_coordinates[1]) * (Y_crop_coordinates[2] - Y_crop_coordinates[1])
        Subset_below_reduction <- Subset_pixel_count / Total_pixel_count <= Reduction_factor
        
        #COMPUTE THE APPROX CROP SIZE (if larger than 1Gb proceed to tiling)
        Subset_approx_size <- Image_size * Subset_pixel_count/Total_pixel_count
        
        Uncomplete_processing <- TRUE #Set the status of the image
        
        
        #If the image subset is smaller than the reduction factor then attempt direct extraction (if crop is <1Gb)
        if(Subset_below_reduction & Subset_approx_size < 1){
          cat("\nAttempting direct image subsetting...")
          try({
            Final_image <- RBioFormats::read.image(Image_directory,
                                                    proprietary.metadata = FALSE,
                                                    read.metadata = TRUE,
                                                    normalize = TRUE,
                                                    subset = list(X = c(X_crop_coordinates[1]:X_crop_coordinates[2]),
                                                                  Y = c(Y_crop_coordinates[1]:Y_crop_coordinates[2]))
            )
            
            #Change color mode as required
            if(!Is_RGB) EBImage::colorMode(Final_image) <- "Grayscale"
            if(Is_RGB) EBImage::colorMode(Final_image) <- "Color"
            
            #Keep desired channels only
            if(any(!is.null(Ordered_Channels), !is.null(Channels_to_keep))) Final_image <- Final_image[,,Channels_to_keep_index]
            
            Final_image <- magick::image_read(Final_image)
            
            Uncomplete_processing <- FALSE #Turn to FALSE
          })
        }
        
        #If image is larger than the pixel requirement or blunt subsetting was unsuccessful then proceed with tiling
        if(any(!Subset_below_reduction | Uncomplete_processing)){
          
          #Compute the final reduction factor if required (Crop size related to final pixel count required)
          Final_resize_factor <- 1
          if(!Subset_below_reduction){
            Final_resize_factor <- Reduction_factor / (Subset_pixel_count / Total_pixel_count)
          }
          
          #First try blunt subsetting and pixel resize (if required for !subsjets_below reduction)
          if(!Subset_below_reduction & Subset_approx_size < 1){
            cat("\nAttempting direct image subsetting...")
            try({
              Final_image <- RBioFormats::read.image(Image_directory,
                                                      proprietary.metadata = FALSE,
                                                      read.metadata = TRUE,
                                                      normalize = TRUE,
                                                      subset = list(X = c(X_crop_coordinates[1]:X_crop_coordinates[2]),
                                                                    Y = c(Y_crop_coordinates[1]:Y_crop_coordinates[2]))
              )
              #Change color mode as required
              if(!Is_RGB) EBImage::colorMode(Final_image) <- "Grayscale"
              if(Is_RGB) EBImage::colorMode(Final_image) <- "Color"
              
              #Keep desired channels only
              if(any(!is.null(Ordered_Channels), !is.null(Channels_to_keep))) Final_image <- Final_image[,,Channels_to_keep_index]
              
              #Downsize according to requirements
              Final_image <- EBImage::resize(Final_image, w = ceiling(dim(Final_image)[1]*sqrt(Final_resize_factor)), h = ceiling(dim(Final_image)[2]*sqrt(Final_resize_factor))) #Sqroot to account for double square reduction
              Final_image <- magick::image_read(Final_image)
              
              Uncomplete_processing <- FALSE #Turn to FALSE
            })
          }
          
          #If blunt subsetting is not successfull either for Subsetting below reduction or above reduction where direct attempt has failed
          if(Uncomplete_processing){
            cat("\nAttempting subsetting by tiling...")
            
            #Compute tile size
            N_rows_columns <- ceiling(sqrt(Subset_approx_size/0.5)) #Each tile is aimed to be 0.5Gb of size (OPTIMAL BASED ON BENCHMARK)
            Subset_width <- X_crop_coordinates[2] - X_crop_coordinates[1]
            Subset_height <- Y_crop_coordinates[2] - Y_crop_coordinates[1]
            
            #Compute the tile coordinates
            #Define the X and Y positions to visit
            X_Positions <- seq(from = X_crop_coordinates[1], to = X_crop_coordinates[2], by = ceiling(Subset_width/N_rows_columns))
            Y_Positions <- seq(from = Y_crop_coordinates[1], to = Y_crop_coordinates[2], by = ceiling(Subset_height/N_rows_columns))
            
            #If the Max X and Y is equal to the Image MAX then remove the final X or Y position
            if(X_Positions[length(X_Positions)] == X_crop_coordinates[2]) X_Positions <- X_Positions[-length(X_Positions)]
            if(Y_Positions[length(Y_Positions)] == Y_crop_coordinates[2]) Y_Positions <- Y_Positions[-length(Y_Positions)]
            
            Tile_Position_tibble <- tidyr::expand_grid(X_Positions, Y_Positions)
            names(Tile_Position_tibble) <- c("X_Position_Min", "Y_Position_Min")
            Tile_Position_tibble$X_Position_Max <- Tile_Position_tibble$X_Position_Min + ceiling(Subset_width/N_rows_columns)-1
            Tile_Position_tibble$Y_Position_Max <- Tile_Position_tibble$Y_Position_Min + ceiling(Subset_height/N_rows_columns)-1
            Tile_Position_tibble$Tile_ID <- stringr::str_c("Tile", 1:nrow(Tile_Position_tibble))
            #If the final X or Y position is larger than the Xmax or Ymax of the image make the ammendments
            Tile_Position_tibble$X_Position_Max[Tile_Position_tibble$X_Position_Max > X_crop_coordinates[2]] <- X_crop_coordinates[2]
            Tile_Position_tibble$Y_Position_Max[Tile_Position_tibble$Y_Position_Max > Y_crop_coordinates[2]] <- Y_crop_coordinates[2]
            
            #Generate the column and the row index
            Col_index_tibble <- tibble::tibble(X_Position_Min = unique(Tile_Position_tibble$X_Position_Min),
                                               Column_index = stringr::str_c("Col_", 1:length(unique(Tile_Position_tibble$X_Position_Min))))
            
            Row_index_tibble <- tibble::tibble(Y_Position_Min = unique(Tile_Position_tibble$Y_Position_Min),
                                               Row_index = stringr::str_c("Row_", 1:length(unique(Tile_Position_tibble$Y_Position_Min))))
            #Bind to the tile_Position_tibble
            Tile_Position_tibble <- dplyr::left_join(Tile_Position_tibble, Col_index_tibble, by = "X_Position_Min") 
            Tile_Position_tibble <- dplyr::left_join(Tile_Position_tibble, Row_index_tibble, by = "Y_Position_Min") 
            
            
            #To increase JAVA RAM use
            rJava::.jinit(parameters="-Xmx200g") 
            #Compute all the tiles
            cat("\nTiling the image")
            
            #Plan the multicore session
            future::plan("future::multisession", workers = N_cores)
            options(future.globals.maxSize = Inf, future.rng.onMisuse = "ignore")
            furrr::furrr_options(scheduling = Inf)
            
            Tiles <- 
              furrr::future_map(1:nrow(Tile_Position_tibble), function(Image_index){
                rJava::.jinit(parameters="-Xmx200g") #To increase JAVA RAM use
                Tile_x_min <- Tile_Position_tibble$X_Position_Min[[Image_index]]
                Tile_x_max <- Tile_Position_tibble$X_Position_Max[[Image_index]]
                
                Tile_y_min <- Tile_Position_tibble$Y_Position_Min[[Image_index]]
                Tile_y_max <- Tile_Position_tibble$Y_Position_Max[[Image_index]]
                
                library(RBioFormats)
                
                Tile_image <- RBioFormats::read.image(Image_directory,
                                                      proprietary.metadata = FALSE,
                                                      read.metadata = TRUE,
                                                      normalize = TRUE,
                                                      subset = list(X = c(Tile_x_min:Tile_x_max),
                                                                    Y = c(Tile_y_min:Tile_y_max))
                )
                
                #Change color mode as required
                if(!Is_RGB) EBImage::colorMode(Tile_image) <- "Grayscale"
                if(Is_RGB) EBImage::colorMode(Tile_image) <- "Color"
                
                #Keep desired channels only
                if(any(!is.null(Ordered_Channels), !is.null(Channels_to_keep))) Tile_image <- Tile_image[,,Channels_to_keep_index]
                
                #Reduce the amount required
                if(Final_resize_factor != 1) Tile_image <- EBImage::resize(Tile_image, w = ceiling(dim(Tile_image)[1]*sqrt(Final_resize_factor)), h = ceiling(dim(Tile_image)[2]*sqrt(Final_resize_factor))) #Sqroot to account for double square reduction
                
                return(Tile_image)
                
              },
              .progress = TRUE)
            future::plan("future::sequential")
            gc()
            
            #Build the columns
            cat("\nMerging tiles")
            Cols <- purrr::map(unique(Tile_Position_tibble$Column_index), function(Current_column_index){
              Tile_index <- Tile_Position_tibble$Column_index == Current_column_index
              
              Col_image <- EBImage::abind(Tiles[Tile_index], along = 2)
              return(Col_image)
            })
            
            #Bind the columns
            Final_image <- EBImage::abind(Cols, along = 1)
            
            #Remove objects manually
            rm(Cols)
            rm(Tiles)
            gc()
            
            Final_image <- magick::image_read(Final_image)
          }
          
        }
      }
    }
    
    ####RETURN THE FINAL PRODUCT####
    if(Save_processed_images){
      cat("\nWriting image in directory...")
      New_image_name_elements <- list(sub("(.*)\\.[^.]*$", "\\1", Image_name),
                                      stringr::str_c("10^", Log10_pixel_output,"Pixels", collapse = ","))
      if(any(!is.null(X_crop_coordinates), !is.null(Y_crop_coordinates))) New_image_name_elements$X_crop_limits <- stringr::str_c("SubsetX", stringr::str_c(X_crop_coordinates, collapse = "_"))
      if(any(!is.null(X_crop_coordinates), !is.null(Y_crop_coordinates))) New_image_name_elements$Y_crop_limits <- stringr::str_c("Y", stringr::str_c(Y_crop_coordinates, collapse = "_"))
      New_image_name_elements$Extension <- ".tif"
      
      
      New_image_name <- paste(unlist(New_image_name_elements), collapse = "_")
      
      try({
        magick::image_write(Final_image, paste0(Output_Directory, "/", New_image_name), format = "tif", depth = 16, compression = "LZW")
        cat(paste0("\n", New_image_name,  " has been successfully saved at ", Output_Directory))
        })
    }
    
    cat("\n DONE!")
    
    #Generate the final list
    Final_list <- list(Image = Final_image,
                       Original_Dims = c(Image_Xmax, Image_Ymax),
                       Current_Dims = c(magick::image_info(Final_image[1])$width, magick::image_info(Final_image[1])$height),
                       Reduction_factor = round(Reduction_factor, digits = 2)
    )
    
    if(any(!is.null(X_crop_coordinates), !is.null(Y_crop_coordinates))) Final_list$X_crop_coords <- X_crop_coordinates
    if(any(!is.null(X_crop_coordinates), !is.null(Y_crop_coordinates))) Final_list$Y_crop_coords <- Y_crop_coordinates
    #If cropping has been performed the reduction factor is the amount of pixels of the final image divided by the amount of pixels of the original subset query
    if(any(!is.null(X_crop_coordinates), !is.null(Y_crop_coordinates))) Final_list$Reduction_factor <- round(
        (Final_list$Current_Dims[1] * Final_list$Current_Dims[2]) / 
        ((X_crop_coordinates[2] - X_crop_coordinates[1]) * (Y_crop_coordinates[2] - Y_crop_coordinates[1])),
        digits = 2
        )
    
    
    return(Final_list)
  }
