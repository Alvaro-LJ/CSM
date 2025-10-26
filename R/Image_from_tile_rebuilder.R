#' Build images from tile image chunks
#'
#' `Image_from_tile_rebuilder()` builds a complete image from it's tiles.
#' Warning: Image reconstruction is limited to ≈ 2-3Gb.
#' Image reconstruction is designed to rebuild complete images from tiles once these have undergone pixel thresholding using [Pixel_Threshold_calculator()].
#' @param Directory Character string specifying the path to the folder where image tiles are present.
#' @param Output_directory Character string specifying the path to the folder where  output images are written. It must be an empty folder.
#' @param RGB_Color_images Logical. Is the image a RGB color image?
#' @param N_cores Integer. Number of cores to parallelize your computation.
#' @returns The function writes the images in the ouput directory.
#'
#' @examples
#' \dontrun{
# Create temporary input and output directories------------------------------
#' Input_Dir <- tempfile(pattern = "tempdir1_Input")
#' Output_Dir_tiles <- tempfile(pattern = "tempdir2_Output")
#' Output_Dir_reconstruction <- tempfile(pattern = "tempdir3_Output")
#' dir.create(Input_Dir, recursive = TRUE)
#' dir.create(Output_Dir_tiles, recursive = TRUE)
#' dir.create(Output_Dir_reconstruction, recursive = TRUE)
#'
#' #Save images in Input directory
#' purrr::map(1:2,
#' function(Image){
#'    EBImage::writeImage(CSM_MiniMultiTiff_test[[Image]], file.path(Input_Dir, names(CSM_MiniMultiTiff_test)[Image]))
#' })
#'
#' #Divide images into tiles------------------------------------------------------
#'Image_tile_deconstruction_function(
#'    Directory = Input_Dir,
#'    Output_directory = Output_Dir_tiles,
#'    Ordered_Channels = c("DAPI", "PDL1", "GZMB", "PD1", "CK-EPCAM", "CD8a", "FOXP3"),
#'    Channels_to_keep = c("DAPI", "GZMB", "CK-EPCAM", "CD8a"),
#'    RGB_Color_images = FALSE,
#'    Tile_pixel_size = 250,
#'    Tile_Overlap = 0,
#'    N_cores = 1
#')
#'
#'#Rebuild images from tiles-----------------------
#'Image_from_tile_rebuilder(
#'    Directory = Output_Dir_tiles,
#'    Output_directory = Output_Dir_reconstruction,
#'    RGB_Color_images = FALSE,
#'    N_cores = 1
#')
#'
#'#Check the files created-------------------------------------------------
#'list.files(Output_Dir_reconstruction)
#'
#'#Remove directories---------------------------------------------------------
#'unlink(c(Input_Dir, Output_Dir_tiles, Output_Dir_reconstruction), recursive = TRUE)
#'
#' }
#' @export

Image_from_tile_rebuilder <-
  function(Directory = NULL,
           Output_directory = NULL,
           RGB_Color_images = FALSE,
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

    #check that the directory contains at least 1 images
    if(!length(dir(Directory)) >= 1) stop("Directory does not contain any files") #Check directory files
    if(berryFunctions::is.error(dir(Output_directory))) stop("Invalid output directory") #Check output directory
    if(length(dir(Output_directory)) != 0) stop("Output directory should be an empty folder") #Check that output directory is empty
    if(!is.logical(RGB_Color_images)) stop("RGB_Color_images must be a logical value")
    if(!all(N_cores >= 1 & N_cores%%1 == 0)) stop("N_cores must be an integer value > 0")

    #Start building the lookup table of the images contained in the director
    Image_names <- dir(Directory, full.names = TRUE)
    Image_names_short <- dir(Directory, full.names = FALSE)

    print("Checking if any image is above 16bit depth")
    #If any image is above 16bits of depth print a message
    Above_16 <-purrr::map_lgl(Image_names, function(Path){
      RBioFormats::read.metadata(Path)$coreMetadata$bitsPerPixel > 16
    })
    if(any(Above_16)) message(paste0(sum(Above_16), "/", length(Above_16),
                                     " images have a bit depth above 16. The final rebuild image will be limited to 16bit depth.")
    )

    #Now build the actual look-up table
    Look_up_table <- tibble(File_names = sub("(.*)\\.[^.]+$", "\\1", Image_names_short))
    Look_up_table$Image_name <- gsub("\\([^)]*\\)", "", Look_up_table$File_names)
    Look_up_table$Tile_info <- regmatches(Look_up_table$File_names, regexpr("(?<=\\()[^)]*(?=\\))", Look_up_table$File_names, perl = TRUE))
    Tile_info_list <- stringr::str_split(Look_up_table$Tile_info, "-")
    Look_up_table$Tile_ID <-purrr::map_chr(Tile_info_list, ~.[[1]])
    Look_up_table$X_Position_Min <-purrr::map_int(Tile_info_list, ~as.integer(.[[2]]))
    Look_up_table$X_Position_Max <-purrr::map_int(Tile_info_list, ~as.integer(.[[3]]))
    Look_up_table$Y_Position_Min <-purrr::map_int(Tile_info_list, ~as.integer(.[[4]]))
    Look_up_table$Y_Position_Max <-purrr::map_int(Tile_info_list, ~as.integer(.[[5]]))
    Look_up_table$Whole_path <- Image_names

    print("Rebuilding images from tiles")
    #Rebuild the images iterating by image
    future::plan("future::multisession", workers = N_cores)
    options(future.globals.maxSize = Inf, future.rng.onMisuse = "ignore")
    furrr::furrr_options(scheduling = Inf)

    furrr::future_map(unique(Look_up_table$Image_name), function(Image){
      #Obtain individual images and arrange by X tile position and Y tile position
      Image_tile_data <- Look_up_table %>% dplyr::filter(Image_name == Image)
      Image_tile_data <- Image_tile_data %>% dplyr::arrange(X_Position_Min, Y_Position_Min)


      #The final image will be binding the generated columns
      Final_image <-
        EBImage::abind(
          purrr::map(unique(Image_tile_data$X_Position_Min), function(Column_position){
            #Generate the row data per column
            Column_information <- Image_tile_data %>% dplyr::filter(X_Position_Min == Column_position)

            rJava::.jinit(parameters="-Xmx200g") #To increase JAVA RAM use

            #Build the actual column
            Column_image <- EBImage::abind(map(Column_information$Whole_path, ~RBioFormats::read.image(.,
                                                                                                       proprietary.metadata = FALSE,
                                                                                                       read.metadata = TRUE,
                                                                                                       normalize = TRUE)),
                                           along = 2) #Generate the columns
            return(Column_image)
          }),
          along = 1 #Bind the columns
        )

      #Modify image type to  Grayscale
      if(!RGB_Color_images) EBImage::colorMode(Final_image) <- "Grayscale"
      if(RGB_Color_images) EBImage::colorMode(Final_image) <- "Color"

      EBImage::writeImage(Final_image, paste0(Output_directory, "/", Image, ".tiff"))
    }, .progress = TRUE)

    future::plan("future::sequential")
    gc()

    print("Done!")
  }
