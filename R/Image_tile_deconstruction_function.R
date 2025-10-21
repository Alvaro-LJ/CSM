#' Divide images into tile chuncks
#'
#' `Image_tile_deconstruction_function()` divides into smaller chunks all the images in a directory.
#' This allows faster computations when using CSM apps, performing image thresholding and cell segmentation.
#' The names of the resulting images should not be modified. If modified this can
#' result in unexpected results when using companion functions [Image_from_tile_rebuilder()], [Tile_to_image_cell_arrange_function()].
#' Tile overlap may be useful to avoid segmentation errors in cells close to the tile edge.
#'
#' Allowed image formats are tiff and ome.tiff. The ome.tiff files have been less extensively tested than tiff files.
#'
#' Warning!! The function uses the RBioFormats Bioconductor package that parses the Bioformats Java library to R. Sometimes even if the user
#' has installed the RBioFormats package from bioconductor the function will still return an error. If this happens, the issue can be solved
#' checking the installed Java version in the computer and manually installing the CRAN rJava package.
#'
#' @param Directory Character string specifying the path to the folder where images to be chopped into tiles are present.
#' @param Output_directory Character string specifying the path to the folder where  output images are written. It must be an empty folder.
#' @param RGB_Color_images Logical. Is the image a RGB color image?
#' @param Ordered_Channels Character vector specifying image channels. If RGB_Color_images is TRUE it will be set automatically to c("R", "G", "B").
#' @param Channels_to_keep Character vector specifying image channels to be kept after image chopping.
#' @param Tile_pixel_size Integer specifying the size of the tiles in pixels.
#' @param Tile_Overlap Integer specifying the overlap amount between tiles in pixels.
#' @param N_cores Integer. Number of cores to parallelize your computation.
#' @returns The function writes the tiles in the ouput directory.
#'
#' @examples
#' \dontrun{
#' Image_tile_deconstruction_function(
#' Directory = "Input_image_directory",
#' Output_directory = "Output_directory",
#' Ordered_Channels = c("Channel_1", "Channel_2", "Channel_3", "Channel_4", "Channel_5"),
#' Channels_to_keep = c("Channel_1", "Channel_2", "Channel_3", "Channel_4", "Channel_5"),
#' RGB_Color_images = FALSE,
#' Tile_pixel_size = 500,
#' Tile_Overlap = 0,
#' N_cores = 2
#')
#' }
#' @export

Image_tile_deconstruction_function <-
  function(Directory = NULL,
           Output_directory = NULL,
           RGB_Color_images = FALSE,
           Ordered_Channels = NULL,
           Channels_to_keep = NULL,
           Tile_pixel_size = NULL,
           Tile_Overlap = 0,
           N_cores = 1){

    #Check installation of suggested packages
    {
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


    on.exit({
      future::plan("future::sequential")
      gc()
    }
    )

    #Check arguments
    if(!length(dir(Directory)) >= 1) stop("Directory does not contain any files") #Check directory files
    if(berryFunctions::is.error(dir(Output_directory))) stop("Invalid output directory") #Check output directory
    if(length(dir(Output_directory)) != 0) stop("Output directory should be an empty folder") #Check that output directory is empty
    if(!is.logical(RGB_Color_images)) stop("RGB_Color_images must be a logical value")
    #If dealing with RGB images assing channels to RGB
    if(RGB_Color_images){
      Ordered_Channels <- c("R", "G", "B")
      Channels_to_keep <- c("R", "G", "B")
    }
    if(!all(Channels_to_keep %in% Ordered_Channels)){
      Problematic_Channels <- Channels_to_keep[!Channels_to_keep %in% Ordered_Channels]
      stop(paste0("The following Channels_to_keep are not present in Ordered_Channels: ", stringr::str_c(Problematic_Channels, collapse = ", ")))
    }
    if(!all(N_cores >= 1 & N_cores%%1 == 0)) stop("N_cores must be an integer value > 0")
    if(!all(Tile_pixel_size >= 1, Tile_pixel_size%%1 == 0)) stop("Tile_pixel_size must be an integer value > 0")
    if(!all(Tile_Overlap >= 0, Tile_Overlap%%1 == 0,  Tile_Overlap < Tile_pixel_size)) stop("Tile_Overlap must be an integer value >= 0 and smaller than Tile_pixel_size")


    print("Obtaining Image metadata")
    #Build the lookup table with the info
    Image_names <- dir(Directory, full.names = TRUE) #Get the full directory
    Image_names_short <- dir(Directory, full.names = FALSE) #Get image names
    Channels_to_keep_index <- match(Channels_to_keep, Ordered_Channels) #Define the channels to keep
    Image_info <-purrr::map_dfr(Image_names, function(Image){
      Metadata <- RBioFormats::read.metadata(Image,
                                             filter.metadata = TRUE,
                                             proprietary.metadata = FALSE)
      #if the names of the metadata is null then we are dealing with image series (like the OME.tiff)
      if(is.null(names(Metadata))){
        Series_N <- length(Metadata@.Data)
        #ALWAYS select the first series
        Metadata <- Metadata@.Data[[1]]
      }
      else Series_N <- 1

      tibble(N_series = Series_N,
             X_size = Metadata$coreMetadata$sizeX,
             Y_size = Metadata$coreMetadata$sizeY,
             N_Channels = Metadata$coreMetadata$sizeC,
             Z_Planes = Metadata$coreMetadata$sizeZ,
             T_Planes = Metadata$coreMetadata$sizeT,
             Pixel_depth = Metadata$coreMetadata$bitsPerPixel,
             Pixel_type = Metadata$coreMetadata$pixelType)
    })
    Image_info <- dplyr::bind_cols(tibble(Image_names = Image_names_short), Image_info)
    Image_info$Tile_pixel_size <- Tile_pixel_size
    if(Tile_Overlap == 0) Image_info$Approx_n_tiles <- ceiling(Image_info$X_size / Image_info$Tile_pixel_size) * ceiling(Image_info$Y_size / Image_info$Tile_pixel_size)
    if(Tile_Overlap > 0) Image_info$Approx_n_tiles <- ceiling(Image_info$X_size / (Image_info$Tile_pixel_size-Tile_Overlap)) * ceiling(Image_info$Y_size / (Image_info$Tile_pixel_size-Tile_Overlap))
    Image_info$Channels_to_keep <- length(Channels_to_keep_index)
    Image_info$Whole_path <- Image_names

    #If there are single tile Images, then print a message and proceed
    if(any(Image_info$Approx_n_tiles == 1)){
      Single_tile_images <- unlist(Image_info %>% dplyr::filter(Approx_n_tiles == 1) %>% dplyr::select(Image_names))
      message(paste0("The following images are smaller or equal to the tile size, hence they will be left unchanged: ", stringr::str_c(Single_tile_images, collapse = ", ")))
    }

    print("Running size approximation test")
    #Run a small test to make the user understand the final size of the individual tiles
    {
      #Select the random sample and obtain their information
      Random_sample <- sample(unlist(Image_info %>% dplyr::filter(Approx_n_tiles != 1) %>% dplyr::select(Whole_path)),
                              size = 1)
      N_series_Random <- Image_info[Image_info$Whole_path == Random_sample, 2]
      N_Channels_Random <- Image_info[Image_info$Whole_path == Random_sample, 5]
      Z_Planes_Random <- Image_info[Image_info$Whole_path == Random_sample, 6]
      T_Planes_Random <- Image_info[Image_info$Whole_path == Random_sample, 7]

      #Increase JAVA RAM use to 200Gb
      rJava::.jinit(parameters="-Xmx200g") #To increase JAVA RAM use

      #IF single series then run as usual
      if(N_series_Random == 1){
        #If channels codes for the actual channels then run the following
        if(N_Channels_Random == length(Ordered_Channels)){
          Image_subset <- RBioFormats::read.image(Random_sample,
                                                  proprietary.metadata = FALSE,
                                                  read.metadata = TRUE,
                                                  normalize = TRUE,
                                                  subset = list(X = c(1:Tile_pixel_size),
                                                                Y = c(1:Tile_pixel_size))
          )
        }
        #If not try the algorithm with the Z planes
        else if(Z_Planes_Random == length(Ordered_Channels)){
          Image_subset <- RBioFormats::read.image(Random_sample,
                                                  proprietary.metadata = FALSE,
                                                  read.metadata = TRUE,
                                                  normalize = TRUE,
                                                  subset = list(X = c(1:Tile_pixel_size),
                                                                Y = c(1:Tile_pixel_size))
          )
        }
        else if(T_Planes_Random == length(Ordered_Channels)){
          Image_subset <- RBioFormats::read.image(Random_sample,
                                                  proprietary.metadata = FALSE,
                                                  read.metadata = TRUE,
                                                  normalize = TRUE,
                                                  subset = list(X = c(1:Tile_pixel_size),
                                                                Y = c(1:Tile_pixel_size))
          )
        }
        #If not this means that the channels are not well specified
        else stop("Length of Ordered_Channels does not match the number of channels in the image")
      }

      #If multiseries then select the first one
      if(N_series_Random > 1){
        #If channels codes for the actual channels then run the following
        if(N_Channels_Random == length(Ordered_Channels)){
          Image_subset <- RBioFormats::read.image(Random_sample,
                                                  proprietary.metadata = FALSE,
                                                  read.metadata = TRUE,
                                                  normalize = TRUE,
                                                  subset = list(X = c(1:Tile_pixel_size),
                                                                Y = c(1:Tile_pixel_size))
          )[[1]]
        }
        #If not try the algorithm with the Z planes
        else if(Z_Planes_Random == length(Ordered_Channels)){
          Image_subset <- RBioFormats::read.image(Random_sample,
                                                  proprietary.metadata = FALSE,
                                                  read.metadata = TRUE,
                                                  normalize = TRUE,
                                                  subset = list(X = c(1:Tile_pixel_size),
                                                                Y = c(1:Tile_pixel_size))
          )[[1]]
        }
        else if(T_Planes_Random == length(Ordered_Channels)){
          Image_subset <- RBioFormats::read.image(Random_sample,
                                                  proprietary.metadata = FALSE,
                                                  read.metadata = TRUE,
                                                  normalize = TRUE,
                                                  subset = list(X = c(1:Tile_pixel_size),
                                                                Y = c(1:Tile_pixel_size))
          )[[1]]
        }
        else stop("Length of Ordered_Channels does not match the number of channels in the image")
      }
      #Lets generate a temporary directory to store the image
      # Create a temporary directory
      temp_dir <- tempdir()
      # Define a file path within the temp directory
      file_path <- file.path(temp_dir, "Test_image.tiff")

      #Define the color mode
      if(!RGB_Color_images) EBImage::colorMode(Image_subset) <- "Grayscale"
      if(RGB_Color_images) EBImage::colorMode(Image_subset) <- "Color"
      #Subset the channels and write the image
      Image_subset <- Image_subset[,,Channels_to_keep_index]
      EBImage::writeImage(Image_subset, file_path)
      #Retrieve information
      File_info <- file.info(file_path)
      #Remove the temporary directory
      unlink(temp_dir, recursive = FALSE)

      #Print the message of the approximate image size and the number of tiles generated
      print(paste0("This is the approximate tile image size: ", round(File_info$size/1000000, 2), " Mb"))
      print(Image_info %>% dplyr::select(Image_names, Approx_n_tiles) %>% dplyr::arrange(desc(Approx_n_tiles)))
      #If tile overlap is positive print a message
      if(Tile_Overlap > 0) message("Tile overlap is implemented to deal with edge bias during cell segmentation. However, current CSM version does not support image from tile rebuild using overlapping tiles.")
    }

    answer <- menu(c("Proceed", "Abort"), title = "Should the tiling process proceed")
    #If user decides to stop then abort function and return stop message
    if(answer == 2) stop("The function has been stopped.")

    #Proceed with the function iterating by Image
    print("Tiling images")
    purrr::map(1:nrow(Image_info), function(Tibble_row){
      #Generate the tile Info tibble
      Image_name <- Image_info$Image_names[Tibble_row]
      Image_Series <- Image_info$N_series[Tibble_row]
      Image_N_Channels <- Image_info$N_Channels[Tibble_row]
      Image_Z_Planes <- Image_info$Z_Planes[Tibble_row]
      Image_T_Planes <- Image_info$T_Planes[Tibble_row]
      Image_Xmax <- Image_info$X_size[Tibble_row]
      Image_Ymax <- Image_info$Y_size[Tibble_row]
      Image_tile_size <- Image_info$Tile_pixel_size[Tibble_row]
      Image_path <- Image_info$Whole_path[Tibble_row]
      Pixel_depth <- Image_info$Pixel_depth[Tibble_row]


      #If Number of tiles
      if(Image_info$Approx_n_tiles[Tibble_row] > 1){

        #Define the X and Y positions to visit
        X_Positions <- seq(from = 1, to = Image_Xmax, by = Tile_pixel_size-Tile_Overlap)
        Y_Positions <- seq(from = 1, to = Image_Ymax, by = Tile_pixel_size-Tile_Overlap)

        #If the Max X and Y is equal to the Image MAX then remove the final X or Y position
        if(X_Positions[length(X_Positions)] == Image_Xmax) X_Positions <- X_Positions[-length(X_Positions)]
        if(Y_Positions[length(Y_Positions)] == Image_Ymax) Y_Positions <- Y_Positions[-length(Y_Positions)]

        Tile_Position_tibble <- tidyr::expand_grid(X_Positions, Y_Positions)
        names(Tile_Position_tibble) <- c("X_Position_Min", "Y_Position_Min")
        Tile_Position_tibble$X_Position_Max <- Tile_Position_tibble$X_Position_Min + Image_tile_size - 1
        Tile_Position_tibble$Y_Position_Max <- Tile_Position_tibble$Y_Position_Min + Image_tile_size - 1
        Tile_Position_tibble$Tile_ID <- stringr::str_c("Tile", 1:nrow(Tile_Position_tibble))
        #If the final X or Y position is larger than the Xmax or Ymax of the image make the ammendments
        Tile_Position_tibble$X_Position_Max[Tile_Position_tibble$X_Position_Max > Image_Xmax] <- Image_Xmax
        Tile_Position_tibble$Y_Position_Max[Tile_Position_tibble$Y_Position_Max > Image_Ymax] <- Image_Ymax

        #Now we use this Tibble to run the tiling process
        future::plan("future::multisession", workers = N_cores)
        options(future.globals.maxSize = Inf, future.rng.onMisuse = "ignore")
        furrr::furrr_options(scheduling = Inf)
        furrr::future_map(1:nrow(Tile_Position_tibble), function(Tile_index){
          library(RBioFormats)
          library(rJava)
          #Generate a name for every tile
          Tile_Image_name <- stringr::str_c(sub("(.*)\\.[^.]+$", "\\1", Image_name),
                                            "(",
                                            Tile_Position_tibble$Tile_ID[Tile_index], "-",
                                            Tile_Position_tibble$X_Position_Min[Tile_index], "-",
                                            Tile_Position_tibble$X_Position_Max[Tile_index], "-",
                                            Tile_Position_tibble$Y_Position_Min[Tile_index], "-",
                                            Tile_Position_tibble$Y_Position_Max[Tile_index],
                                            ")",
                                            ".tiff")

          rJava::.jinit(parameters="-Xmx200g") #To increase JAVA RAM use

          #Try it
          tryCatch(
            {
              if(Image_Series == 1){
              Tile_image <- RBioFormats::read.image(Image_path,
                                                    proprietary.metadata = FALSE,
                                                    read.metadata = TRUE,
                                                    normalize = TRUE,
                                                    subset = list(X = c(Tile_Position_tibble$X_Position_Min[Tile_index]:Tile_Position_tibble$X_Position_Max[Tile_index]),
                                                                  Y = c(Tile_Position_tibble$Y_Position_Min[Tile_index]:Tile_Position_tibble$Y_Position_Max[Tile_index]))
              )
              }

              #If the image is composed of various series then select the first one
              if(Image_Series > 1){
                #If the Number of channels match the ordered channels
                Tile_image <- RBioFormats::read.image(Image_path,
                                                      proprietary.metadata = FALSE,
                                                      read.metadata = TRUE,
                                                      normalize = TRUE,
                                                      subset = list(X = c(Tile_Position_tibble$X_Position_Min[Tile_index]:Tile_Position_tibble$X_Position_Max[Tile_index]),
                                                                    Y = c(Tile_Position_tibble$Y_Position_Min[Tile_index]:Tile_Position_tibble$Y_Position_Max[Tile_index]))
                )[[1]]
              }
              #write it in the output directory
              if(!RGB_Color_images) EBImage::colorMode(Tile_image) <- "Grayscale"
              if(RGB_Color_images) EBImage::colorMode(Tile_image) <- "Color"
              Tile_image <- Tile_image[,,Channels_to_keep_index]
              EBImage::writeImage(Tile_image, paste0(Output_directory, "/", Tile_Image_name))


            },
            error = function(e) {
              message(paste0("Error in ", sub("(.*)\\.[^.]+$", "\\1", Image_name), "-", Tile_Position_tibble$Tile_ID[Tile_index], ":",  e$message))
            }
          )

        },
        .progress = TRUE)
        future::plan("future::sequential")
        gc()
      }

      #if a single slide is required then run a simple code
      if(Image_info$Approx_n_tiles[Tibble_row] == 1){
        library(RBioFormats)

        Tile_Image_name <- stringr::str_c(sub("(.*)\\.[^.]+$", "\\1", Image_name),
                                          "(",
                                          "Tile1", "-",
                                          "1", "-",
                                          Image_Xmax, "-",
                                          "1", "-",
                                          Image_Ymax,
                                          ")",
                                          ".tiff")

        rJava::.jinit(parameters="-Xmx200g") #To increase JAVA RAM use

        tryCatch(
          {
            #If the image is of a single series proceed as usual
            if(Image_Series == 1){
              Tile_image <- RBioFormats::read.image(Image_path,
                                                    proprietary.metadata = FALSE,
                                                    read.metadata = TRUE,
                                                    normalize = TRUE
              )
            }

            #If the image is composed of various series then select the first one
            if(Image_Series > 1){
              Tile_image <- RBioFormats::read.image(Image_path,
                                                    proprietary.metadata = FALSE,
                                                    read.metadata = TRUE,
                                                    normalize = TRUE
              )[[1]]
            }

            #write it in the output directory
            if(!RGB_Color_images) EBImage::colorMode(Tile_image) <- "Grayscale"
            if(RGB_Color_images) EBImage::colorMode(Tile_image) <- "Color"
            Tile_image <- Tile_image[,,Channels_to_keep_index]
            EBImage::writeImage(Tile_image, paste0(Output_directory, "/", Tile_Image_name))


          },
          #If gives an error print the error type
          error = function(e) {
            message(paste0("Error in ", sub("(.*)\\.[^.]+$", "\\1", Image_name), "-Tile1: ",  e$message))
          }
        )
      }

      #Print the name of the image that has been tiled
      sub("(.*)\\.[^.]+$", "\\1", Image_name)
    }, .progress = list(clear = F,
                        name = "Generating tiles",
                        show_after = 1,
                        type = "iterator"))
  }
