#' Generates cell segmentation mask using the watershed algorithm
#'
#' Performs watershed cell segmentation using the simpleSeg Bioconductor R package.
#' Parameters can then be obtained using [Watershed_tester_app()].
#'
#' @param Directory Character specifying the path to the folder where images to be segmented are present.
#' @param Parameter_list List obtained using the [Watershed_tester_app()] containing segmentation parameters.
#' @param Ordered_Channels Character vector specifying image channels in their exact order.
#' @param N_cores Integer. Number of cores to parallelize your computation.
#'
#' @param Output_directory A character value indicating the path to the empty folder where output images should be written.
#'
#' @param Nuclear_marker (Used if Parameters_list is NULL) Character vector of markers corresponding to the nuclei.
#' @param Cell_body_method (Used if Parameters_list is NULL) Method of cytoplasm identification. Can be 'none', 'dilate' or 'discModel'.
#' @param Min_pixel (Used if Parameters_list is NULL) Integer value specifying the minimum pixels for an object to be recognized as a cell and not noise.
#' @param Smooth_amount (Used if Parameters_list is NULL) Numeric value specifying the amount of Gaussian smoothing to be applied to the image.
#' @param Normalization (Used if Parameters_list is NULL) Single value or vector specifying the transformations from "sqrt", "asinh", "norm99", "maxThresh" and "tissueMask.
#' @param Watershed_type (Used if Parameters_list is NULL) Method used to perform watersheding. Accepted values: "intensity", "distance" or "combine".
#' @param Tolerance_value (Used if Parameters_list is NULL) Numeric value specifying the minimum height of the object in the units of image intensity between its highest point (seed) and the point where it contacts another object (MAY BE NULL).
#' @param Neighborhood_distance (Used if Parameters_list is NULL) Radius of the neighborhood in pixels for the detection of neighboring objects. Higher value smooths out small objects.
#' @param Disc_size (Used if Parameters_list is NULL) The size of dilation around nuclei to create cell disc or capture cytoplasm.
#' @param Tissue_mask_markers (Used if Parameters_list is NULL) A vector specifying the channels to be used to create the tissue mask if specified in transformation.
#' @param Perform_PCA (Used if Parameters_list is NULL) A logical value specifying wether to run PCA on aggregated nucleus markers in order to detect the cellular nucclei.
#'
#' @param Perform_nuclear_channel_processing (Used if Parameters_list is NULL) A logical value specifying if nuclear channel pre-processing should be performed.
#' @param Black_level (Used if Parameters_list is NULL)(Used if Perform_nuclear_channel_processing is TRUE) Numeric indicating the % below the max intensity level to be set to absolute black.
#' @param White_level (Used if Parameters_list is NULL)(Used if Perform_nuclear_channel_processing is TRUE) Numeric indicating the % above the max intensity level to be set to absolute white.
#' @param Gamma_level (Used if Parameters_list is NULL)(Used if Perform_nuclear_channel_processing is TRUE) Numeric value between -3 and +3 to indicate channel gamma.
#' @param Equalize (Used if Parameters_list is NULL)(Used if Perform_nuclear_channel_processing is TRUE) A logical value specifying if channel should be equalized.
#' @param Opening_kernel_size (Used if Parameters_list is NULL)(Used if Perform_nuclear_channel_processing is TRUE) Opening kernel size (set to 1 if no opening is required).
#' @param Closing_kernel_size (Used if Parameters_list is NULL)(Used if Perform_nuclear_channel_processing is TRUE) Closing kernel size (set to 1 if no closing is required).
#'
#' @returns Writes the cell mask images in the output directory
#'
#' @seealso [Watershed_tester_app()]
#'
#' @examples
#' \dontrun{
#'#Create temporary input directory----------------------------------------
#'Input_Dir <- tempfile(pattern = "tempdir1_Input")
#'dir.create(Input_Dir, recursive = TRUE)
#'
#'Output_Dir <- tempfile(pattern = "tempdir2_Output")
#'dir.create(Output_Dir, recursive = TRUE)
#'
#'#Save images in Input directory
#'purrr::map(1:2,
#'           function(Image){
#'             EBImage::writeImage(CSM_MiniMultiTiff_test[[Image]],
#'                                 file.path(Input_Dir, names(CSM_MiniMultiTiff_test)[Image]))
#'          })
#'
#'#Check a segmentation parameters list obtained using the dedicated function----------
#'print(CSM_SegmentParams_test)
#'
#'
#'#Generate the whole cell masks using watershed---------------------------------------
#'Watershed_cell_mask_generator(
#'  Directory = Input_Dir,
#'  Output_directory = Output_Dir,
#'  Parameter_list = CSM_SegmentParams_test
#')
#'
#'#Prepare parameters to obtain subcellular nuclear and cytoplasmic compartments---
#'Subcellular_compartment_parameters <-
#'  list(Nuclear = list(Channel_name = "DAPI",
#'                      Type = "Positive",
#'                      Threshold_type = "Arbitrary",
#'                      Threshold_value = 0.05,
#'                      Blurr = TRUE,
#'                      Sigma = 0.4),
#'       Cytoplasmic = list(Channel_name = "DAPI",
#'                          Type = "Negative",
#'                          Threshold_type = "Arbitrary",
#'                          Threshold_value = 0.05,
#'                          Blurr = TRUE,
#'                          Sigma = 0.4)
#'  )
#'
#'#Obtain the cell feature matrix--------------------------------------------
#'Cell_feature_extractor(Image_directory = Input_Dir,
#'                       Mask_directory = Output_Dir,
#'                       Ordered_Channels = CSM_SegmentParams_test[["Channels_to_keep"]],
#'                       Subcellular_compartment_parameter_list = Subcellular_compartment_parameters)
#'
#'#Remove directories---------------------------------------------------------
#'unlink(Input_Dir, recursive = TRUE)
#'unlink(Output_Dir, recursive = TRUE)
#' }
#'
#'
#' @export

Watershed_cell_mask_generator <-
  function(

    Directory,
    Parameter_list = NULL,
    Ordered_Channels = NULL,
    N_cores = 1,

    Output_directory = NULL,

    Nuclear_marker = NULL, #marker or list of markers corresponding to the nuclei
    Cell_body_method = NULL, #Method of cytoplasm identification. Can be 'none', dilate', 'discModel' or the name of a dedicated cytoplasm marker
    Min_pixel = NULL, # Minimum pixels for an object to be recognized as a cell and not noise
    Smooth_amount = NULL, #The amount of Gaussian smoothing to be applied to the image
    Normalization = NULL, #Single value or list specifying the transformations from "sqrt", "asinh", "norm99", "maxThresh" and "tissueMask
    Watershed_type = NULL, #Method used to perform watersheding. Accepted values: "intensity", "distance" or "combine"
    Tolerance_value = NULL, #minimum height of the object in the units of image intensity between its highest point (seed) and the point where it contacts another object (MAY BE NULL)
    Neighborhood_distance = NULL, #Radius of the neighborhood in pixels for the detection of neighboring objects. Higher value smooths out small objects.
    Disc_size = NULL, #The size of dilation around nuclei to create cell disc or capture cytoplasm
    Tissue_mask_markers = NULL, #A vector specifying the channels to be used to create the tissue mask if specified in transforms
    Perform_PCA = FALSE, #Whether to run PCA on aggregated nucleus markers in order to detect the cellular nucclei

    Perform_nuclear_channel_processing = NULL,
    Black_level = NULL,
    White_level = NULL,
    Gamma_level = NULL,
    Equalize = NULL,
    Opening_kernel_size = NULL,
    Closing_kernel_size = NULL

    ){
    #####Check suggested packages####
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
      if(!requireNamespace("magick", quietly = FALSE)) stop(
        paste0("magick CRAN package is required to execute the function. Please install using the following code: ",
               expression(install.packages("benchmarkme")))
      )
    }

    #####Specify what to do on exit####
    on.exit(gc())

    #Check regular arguments
    if(length(dir(Directory)) < 1) stop("No files found at the directory provided. Please check out the path.")
    if(!all(N_cores >= 1, N_cores%%1 == 0)) stop("N_cores must be an integer value > 0")
    if(!dir.exists(Output_directory)) stop("Issue with output directory. Please check the path provided.")
    if(length(list.files(Output_directory)) != 0) stop("Output_directory must be an empty directory")

    #If parameter list is provided obtain the parameters
    if(!is.null(Parameter_list)){
      cat("\n Parameter_list provided. Any other segmentation related parameters will be ignored") #Print a message

      Ordered_Channels <- Parameter_list[["Ordered_Channels"]]
      Nuclear_marker <- Parameter_list[["Nuclear_marker"]]
      Cell_body_method <- Parameter_list[["Cell_body_method"]]
      Min_pixel <- Parameter_list[["Min_pixel"]]
      Smooth_amount <- Parameter_list[["Smooth_amount"]]
      Normalization <- Parameter_list[["Normalization"]]
      Watershed_type <- Parameter_list[["Watershed_type"]]
      Tolerance_value <- Parameter_list[["Tolerance_value"]]
      Neighborhood_distance <- Parameter_list[["Neighborhood_distance"]]
      Disc_size <- Parameter_list[["Disc_size"]]
      Tissue_mask_markers <- Parameter_list[["Tissue_mask_markers"]]
      Perform_PCA <- Parameter_list[["Perform_PCA"]]

      Perform_nuclear_channel_processing <- Parameter_list[["Perform_nuclear_channel_processing"]]
      Black_level <- Parameter_list[["Black_level"]]
      White_level <- Parameter_list[["White_level"]]
      Gamma_level <- Parameter_list[["Gamma_level"]]
      Equalize <- Parameter_list[["Equalize"]]
      Opening_kernel_size <- Parameter_list[["Opening_kernel_size"]]
      Closing_kernel_size <- Parameter_list[["Closing_kernel_size"]]
    }

    #If null import from arguments
    else if(is.null(Parameter_list)){
      Ordered_Channels <- Ordered_Channels
      Nuclear_marker <- Nuclear_marker
      Cell_body_method <- Cell_body_method
      Min_pixel <- Min_pixel
      Smooth_amount <- Smooth_amount
      Normalization <- Normalization
      Watershed_type <- Watershed_type
      Tolerance_value <- Tolerance_value
      Neighborhood_distance <- Neighborhood_distance
      Disc_size <- Disc_size
      Tissue_mask_markers <- Tissue_mask_markers
      Perform_PCA <- Perform_PCA
      Perform_nuclear_channel_processing <- Perform_nuclear_channel_processing
      Black_level <- Black_level
      White_level <- White_level
      Gamma_level <- Gamma_level
      Equalize <- Equalize
      Opening_kernel_size <- Opening_kernel_size
      Closing_kernel_size <- Closing_kernel_size
    }

    #Check arguments by generating a argument check vector and message vector
    Argument_checker <- c(Empty_directory = length(dir(Directory)) >= 1,
                          N_cores_OK = (N_cores >= 1 & N_cores%%1 == 0),
                          Nuclear_OK = unlist(Nuclear_marker) %in% Ordered_Channels,
                          Cell_body_OK = Cell_body_method %in% c("none", "dilate", "discModel"),
                          Min_pixel_OK = all(is.numeric(Min_pixel), Min_pixel%%1 == 0, Min_pixel >= 0),
                          Smooth_amount_OK = all(is.numeric(Smooth_amount), Smooth_amount >= 0),
                          Normalization_OK = Normalization %in% c("sqrt", "asinh", "norm99", "maxThresh", "tissueMask"),
                          Watershed_type_OK = Watershed_type %in% c("intensity", "distance", "combine"),
                          Tolerance_value_OK = if(!is.null(Tolerance_value)) {
                            all(is.numeric(Tolerance_value), Tolerance_value >= 0)
                          } else(TRUE),
                          Neighborhood_distance_OK = all(Neighborhood_distance%%1 == 0, Neighborhood_distance >= 0),
                          Disc_size_OK = all(Disc_size%%1 == 0, Disc_size >= 0),
                          Tissue_mask_markers_OK =  if(!is.null(Tissue_mask_markers)) {
                            all(Tissue_mask_markers %in% Ordered_Channels)
                          } else(TRUE),
                          Perform_PCA_OK = is.logical(Perform_PCA),
                          Perform_nuclear_channel_processing_OK = is.logical(Perform_nuclear_channel_processing)
    )

    Stop_messages <- c(Empty_directory = "No files found at the directory provided. Please check out the path.",
                       N_cores_OK = "N_cores must be a positive integer value",
                       Nuclear_OK = stringr::str_c("Nuclear marker specified not found in ", stringr::str_c(Ordered_Channels, collapse = ", ")),
                       Cell_body_OK = "Cell body method must be one of the following: none, dilate, discModel",
                       Min_pixel_OK = "Min pixel must be a positive integer value",
                       Smooth_amount_OK = "Smooth amount must be a positive numeric value",
                       Normalization_OK = "Normalization method must be one of the following: sqrt, asinh, norm99, maxThresh, tissueMask",
                       Watershed_type_OK = "Watershed type must be NULL or one of the following: intensity, distance, combine",
                       Tolerance_value_OK = "Tolerance value must be NULL or a positive numeric value",
                       Neighborhood_distance_OK = "Neighborhood distance must be a positive integer value",
                       Disc_size_OK = "Disc size must be a positive integer value",
                       Tissue_mask_markers_OK =  stringr::str_c("Tissue mask must be NULL or one of the following: ", stringr::str_c(Ordered_Channels, collapse = ", ")),
                       Perform_PCA_OK = "Perform PCA must be a logical value",
                       Perform_nuclear_channel_processing_OK = "Perform_nuclear_channel_processing must be a logical value")


    #Check arguments and stop if necessary
    if(!all(Argument_checker)){
      stop(cat(Stop_messages[!Argument_checker],
               fill = sum(!Argument_checker)))
    }

    #Check specifically nuclear processing arguments if Perform_nuclear_channel_processing is true
    if(Perform_nuclear_channel_processing){
      if(!all(is.numeric(Black_level), Black_level >= 0, Black_level <= 100, Black_level < White_level)) stop("Black_level must be a numeric value between 0 - 100 and smaller than White_level")
      if(!all(is.numeric(White_level), White_level >= 0, White_level <= 100)) stop("White_level must be a numeric value between 0 - 100")
      if(!all(is.numeric(Gamma_level), Gamma_level >= -3, Gamma_level <= 3)) stop ("Gamma_level must be a numeric value between -3 and +3")
      if(!is.logical(Equalize)) stop("Equalize must be a logical value")
      if(!all(Opening_kernel_size%%1 == 0, Opening_kernel_size >= 1)) stop("Opening_kernel_size must be a positive integer value >= 1")
      if(!all(Closing_kernel_size%%1 == 0, Closing_kernel_size >= 1)) stop("Closing_kernel_size must be a positive integer value >= 1")
    }

    #Obtain the channels and the names of the directory
    Channels <- Ordered_Channels
    Simple_Image_Names <- dir(Directory, full.names = FALSE)
    Complete_Image_Names <- dir(Directory, full.names = TRUE)
    Channels_in_use <- unique(c(Nuclear_marker, Tissue_mask_markers))
    Channels_in_use_index <- match(Channels_in_use, Ordered_Channels)

    #save exit function if parallelization fails
    on.exit({
      future::plan("future::sequential")
      gc()
    })

    #We make the clusters
    future::plan("future::multisession", workers = N_cores)
    options(future.globals.maxSize = Inf, future.rng.onMisuse = "ignore")
    furrr::furrr_options(scheduling = Inf)

    #Execute the function
    furrr::future_walk(seq_along(1:length(Complete_Image_Names)), function(Index){
      RESULT_TIBBLE <- try({
        #Pre-process nuclear channels if required
        if(Perform_nuclear_channel_processing){
          Image <- magick::image_read(Complete_Image_Names[Index]) #Import image

          Nuclear_channels_number <- match(Nuclear_marker, Ordered_Channels)

          #Apply changes to every nuclear channel
          for(index in Nuclear_channels_number){
            Image_Modified <- Image[index]

            if(Equalize) Image_Modified <- Image_Modified %>% magick::image_equalize() #Equalize if necesary
            Image_Modified <- Image_Modified %>% magick::image_level(black_point = Black_level,
                                                                     white_point = White_level,
                                                                     mid_point = 10^Gamma_level) #Chane withe, black and gamma
            Image_Modified <- Image_Modified %>% magick::as_EBImage() #turn to EBImage object
            Image_Modified  <-
              Image_Modified %>% EBImage::opening(EBImage::makeBrush(size = Opening_kernel_size, shape = "disc")) %>%
              EBImage::closing(EBImage::makeBrush(size = Closing_kernel_size, shape = "disc")) #opening and closing

            Image_Modified <- magick::image_read(Image_Modified) #again as Magick image

            Image[index] <- Image_Modified
          }
          #Returna a EBImage object
          Image <- Image %>% magick::as_EBImage()
        }

        #Import image with the EBImage importer if no processing is required
        else{
          Image <- EBImage::readImage(Complete_Image_Names[Index])
        }

        #Retain only the required channels
        Image <- Image[,,Channels_in_use_index]

        #Transform it to cytoImage object
        Image <- cytomapper::CytoImageList(Image)
        cytomapper::channelNames(Image) <- Channels_in_use #define channel names
        S4Vectors::mcols(Image)$imageID <- as.character(Simple_Image_Names[Index])#Modify name

        #Perform cell segmentation
        Seg_results <- try(
          simpleSeg::simpleSeg(Image,
                               nucleus = Nuclear_marker,
                               cellBody = Cell_body_method,
                               sizeSelection = Min_pixel,
                               smooth = Smooth_amount,
                               transform = Normalization,
                               watershed = Watershed_type,
                               tolerance = Tolerance_value,
                               ext = Neighborhood_distance,
                               discSize = Disc_size,
                               tissue = Tissue_mask_markers,
                               pca = Perform_PCA,
                               cores = 1)
          )



        #Generate the name
        Mask_name <- stringr::str_c(Simple_Image_Names[Index], "_CellMask.tif")
        File_name <- stringr::str_c(Output_directory, "/", Mask_name)

        #Modify the values of the image
        Seg_results <- try(Seg_results[[1]]/max(Seg_results[[1]]))

        #Print a message if error and write empty segmentation mask
        if(berryFunctions::is.error(Seg_results) | max(Seg_results) == 0){
          Colored_print(paste0("\nIssue with ", Simple_Image_Names[Index], ". Unable to generate cell mask or no cells identified. \n Generating blank cell mask"), color = "red")

          #Get dimensions
          Image_width <- dim(Image[[1]])[2]
          Image_height <- dim(Image[[1]])[1]

          #Generate matrix
          Image_matrix <- matrix(0, nrow = Image_height, ncol = Image_width)
          Seg_results <- EBImage::as.Image(Image_matrix)

          #Write the matrix
          EBImage::writeImage(Seg_results, File_name, bits.per.sample = 16, compression = "LZW")

        } else{
          #Save modify the file
          EBImage::writeImage(Seg_results, File_name, bits.per.sample = 16, compression = "LZW")
        }

        #Remove Image to save RAM space
        rm(Image)
        rm(Seg_results)
        gc()
      }
      )
    },
    .progress = TRUE)
    future::plan("future::sequential")
    gc()

    #Print done
    print("DONE!")

  }
