#' Calculates MFI for a set of images
#'
#' Calculates the Mean Fluorescence Intensity of a target channel for every image in a directory. MFI can be calculated inside compartments using the Target_mask argument.
#' @param N_cores Integer. Number of cores to parallelize your computation.
#' @param Directory Character string specifying the path to the folder where images are present.
#' @param Ordered_Channels Character vector specifying image channels in their exact order.
#' @param Channels_to_keep Character vector indicating channels used in tissue mask generation.
#' @param Target_channel Character value indicating the target channel to be thresholded.
#'
#' @param Target_masks (OPTIONAL) A named list of the parameters used to calculate the compartment where MFI will be calculated (see details).
#'
#' @param Threshold_type_tissueMask Type of threshold to performe tissue mask. Either 'Otsu', 'Arbitrary' or 'Absolute'.
#' @param Threshold_value_tissueMask Numeric value used if Arbitrary is the threshold type of choice.
#' @param Blurr_tissueMask Logical value indicating if image blurring be performed before tissue mask generation.
#' @param Sigma_tissueMask Numeric value indicating the sigma value to perform Gaussian blurring.
#'
#'
#' @returns Returns a tibble with the MFI per image.
#'
#' @details
#'Target_masks should be a named lists. The names of the list should be any of the Channels_to_keep. Each element should be a list with the following named elements:
#'Mask_name, Threshold_type, Threshold_value, Blurr, Sigma. These parameters should follow the same rules as the ones used to calculate tissueMask (see examples).
#'
#' @seealso [Image_thresholding_app_launcher()], [Binary_threshold_image_combinator()], [MFI_Experimet_Calculator()]
#'
#' @examples
#' \dontrun{
#'#Create temporary input and output directories------------------------------
#' Input_Dir <- tempfile(pattern = "tempdir1_Input")
#' dir.create(Input_Dir, recursive = TRUE)
#'
#' #Save images in Input directory--------------------------------------------
#' purrr::map(1:2,
#' function(Image){
#'    EBImage::writeImage(CSM_MiniMultiTiff_test[[Image]],
#'    file.path(Input_Dir, names(CSM_MiniMultiTiff_test)[Image]))
#' })
#'
#' #Optionally generate marker masks to obtain the target MFI----------------------------
#' #Calculate the target in a combined CK, CD8 mask that will be overlayed to the tissue mask
#'Target_mask_list <- list(CK = list(Mask_name = "CK-EPCAM",
#'                                   Threshold_type = "Arbitrary",
#'                                   Threshold_value = 0.1,
#'                                   Blurr = TRUE,
#'                                   Sigma = 1),
#'                       CD8 = list(Mask_name = "CD8a",
#'                                  Threshold_type = "Arbitrary",
#'                                  Threshold_value = 0.001,
#'                                  Blurr = TRUE,
#'                                  Sigma = 1)
#')
#'
#'
#' #Then calculate the MFI---------------------------
#'MFI_Experimet_Calculator(
#'    N_cores = 1,
#'    Directory = Input_Dir,
#'    Ordered_Channels = c("DAPI", "PDL1", "GZMB", "PD1", "CK-EPCAM", "CD8a", "FOXP3"),
#'    Channels_to_keep = c("DAPI", "PDL1", "GZMB", "PD1", "CK-EPCAM", "CD8a", "FOXP3"),
#'    Target_channel = "GZMB",
#'
#'    Target_masks = Target_mask_list,
#'
#'    Threshold_type_tissueMask = "Arbitrary",
#'    Threshold_value_tissueMask = 0.01,
#'    Blurr_tissueMask = TRUE,
#'    Sigma_tissueMask = 0.5
#')
#'
#'#Remove directories---------------------------------------------------------
#'unlink(Input_Dir, recursive = TRUE)
#' }
#'
#' @export

MFI_Experimet_Calculator <-
  function(N_cores = 1,
           Directory,
           Ordered_Channels,
           Channels_to_keep,
           Target_channel,

           Target_masks = NULL,

           Threshold_type_tissueMask,
           Threshold_value_tissueMask = NULL,
           Blurr_tissueMask = FALSE,
           Sigma_tissueMask = NULL
  ){

    #Check suggested packages
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
      if(!requireNamespace("magick", quietly = FALSE)) stop(
        paste0("magick CRAN package is required to execute the function. Please install using the following code: ",
               expression(install.packages("magick")))
      )
    }


    on.exit({
      future::plan("future::sequential")
      gc()
    })

    #Argument check general arguments
    Argument_checker <- c(N_cores_OK = (N_cores >= 1 & N_cores%%1 == 0),
                          Empty_directory = length(dir(Directory)) >= 1,
                          Channels_OK = all(Channels_to_keep %in% Ordered_Channels),
                          Target_channel_OK = Target_channel %in% Channels_to_keep,

                          Threshold_type_tissueMask_OK = Threshold_type_tissueMask %in% c("Otsu", "Arbitrary", "Absolute"),
                          Threshold_value_tissueMask_OK = if(Threshold_type_tissueMask == "Arbitrary"){
                            all(is.numeric(Threshold_value_tissueMask), Threshold_value_tissueMask >=0, Threshold_value_tissueMask <= 1)
                          } else(TRUE),
                          Blurr_tissueMask_OK = is.logical(Blurr_tissueMask),
                          Sigma_tissueMask_OK = if(Blurr_tissueMask){
                            all(is.numeric(Sigma_tissueMask), Sigma_tissueMask > 0)
                          } else(TRUE)
    )

    Stop_messages <- c(N_cores_OK = "N_cores must be an integer value > 0",
                       Empty_directory = "No files found at the directory provided. Please check out the path.",
                       Channels_OK = stringr::str_c(
                         "The following channels are not present the channel names provided: ",
                         stringr::str_c(Channels_to_keep[!(Channels_to_keep %in% Ordered_Channels)], collapse = ", "),
                         sep = ""),
                       Target_channel_OK = stringr::str_c(Target_channel, " not present in ", stringr::str_c(Channels_to_keep, collapse = ", "), collapse = ""),


                       Threshold_type_tissueMask_OK = "Threshold_type_TissueMask must be one of the following: Otsu, Arbitrary, Absolute",
                       Threshold_value_tissueMask_OK = "Threshold_value_tissueMask must be a single numeric value between 0 and 1",
                       Blurr_tissueMask_OK = "Blurr_tissueMask must be a logical value",
                       Sigma_tissueMask_OK = "Sigma_tissueMask must be a positive numeric value > 0"
    )
    #Check arguments and stop if necessary
    if(!all(Argument_checker)){
      stop(cat(Stop_messages[!Argument_checker],
               fill = sum(!Argument_checker)))
    }

    #Check specifically the Target masks
    #If no target mask required proceed appropriately
    if(is.null(Target_masks)) print("Calculating MFI using the tissue mask alone")
    #If complex target masks then check arguments
    else{
      #Check that target masks are a list
      if(!is.list(Target_masks)) stop("Target_masks must be a list containing mask parameters")

      #Check that target masks items are a list
      if(!all(
        purrr::map_lgl(Target_masks, function(Individual_mask){
          is.list(Individual_mask)
        }))) stop("Individual Items in the Target_mask must be a list")

      #Check names of masks (Mask_name, Threshold_type, Threshold_value, Blurr, Sigma)
      Adequate_names <-
        purrr::map_lgl(Target_masks, function(Individual_mask){
          identical(names(Individual_mask), c("Mask_name", "Threshold_type", "Threshold_value", "Blurr", "Sigma"))
        })
      if(!all(Adequate_names)){
        stop(paste0("Names of Target_masks must be the folloing: Mask_name, Threshold_type, Threshold_value, Blurr, Sigma",
                    "The following masks have inadequate names: ",
                    stringr::str_c(which(!Adequate_names), collapse = ", ")))
      }
      #Check the actual arguments within the list of lists
      purrr::walk(Target_masks, function(Individual_mask){
        #Check mask name present in channels to keep
        if(!Individual_mask[["Mask_name"]] %in% Channels_to_keep) stop(paste0(Individual_mask[["Mask_name"]], ": Invalid mask name.",
                                                                              "It must be present in channels_to_keep: ",
                                                                              stringr::str_c(Channels_to_keep, collapse = ", ")))

        #Check Threshold type
        if(!Individual_mask[["Threshold_type"]] %in% c("Arbitrary", "Otsu")) stop(paste0(Individual_mask[["Mask_name"]], ": Invalid Threshold type.",
                                                                                         "It must be one of the following: Arbitrary, Otsu"))

        #Check Threshold value if Arbitrary
        if(Individual_mask[["Threshold_type"]] == "Arbitrary"){
          if(!all(is.numeric(Individual_mask[["Threshold_value"]]), Individual_mask[["Threshold_value"]] >=0, Individual_mask[["Threshold_value"]] <= 1)){
            stop(paste0(Individual_mask[["Mask_name"]], ": Invalid Threshold_value.",
                        "It must be a numeric value between 0 and 1"))
          }
        }

        #Check blur
        if(!is.logical(Individual_mask[["Blurr"]])) stop(paste0(Individual_mask[["Mask_name"]], ": Invalid Blurr argument.",
                                                                "It must be a logical value"))
        #Check Sigma value
        if(Individual_mask[["Blurr"]]){
          if(!all(is.numeric(Individual_mask[["Sigma"]]), Individual_mask[["Sigma"]] > 0)) stop(paste0(Individual_mask[["Mask_name"]], ": Invalid Sigma argument.",
                                                                                                       "Sigma must be a numeric value > 0"))
        }
      })
    }

    #Get the full image directory and the short image names
    Full_directory <- dir(Directory, full.names = TRUE)
    Image_names <- dir(Directory, full.names = FALSE)

    #Will iterate for every image to obtain tissue mask and MFI
    future::plan("future::multisession", workers = N_cores)
    options(future.globals.maxSize = Inf, future.rng.onMisuse = "ignore")
    furrr::furrr_options(scheduling = Inf)
    RESULTS <-
      suppressMessages(
        furrr::future_map(seq_along(1:length(dir(Directory))), function(Image){
          #Import the image and turn it into EBI image format
          Image <- magick::image_read(Full_directory[Image])
          Image <- Image[which(Channels_to_keep %in% Ordered_Channels)]
          Image <- Image %>% magick::as_EBImage()


          #First generate a tissue mask
          Tissue_mask <- Tissue_mask_generator(Image = Image,
                                               Threshold_type = Threshold_type_tissueMask,
                                               Threshold_value = Threshold_value_tissueMask,
                                               Blurr = Blurr_tissueMask,
                                               Sigma = Sigma_tissueMask)

          #Second generate potential target tissue masks and generate the final combined mask (Tissue + target)
          if(!is.null(Target_masks)){
            #Generate a list containing all the masks images
            Target_mask_list <-
              purrr::map(Target_masks, function(Mask_parameters){
                return(
                  Pixel_thresholder(Target = EBImage::getFrame(Image, which(Channels_to_keep == Mask_parameters[["Mask_name"]])),
                                    Tissue_mask = Tissue_mask,
                                    Threshold_type = Mask_parameters[["Threshold_type"]],
                                    Threshold_value = Mask_parameters[["Threshold_value"]],
                                    Blurr = Mask_parameters[["Blurr"]],
                                    Sigma = Mask_parameters[["Sigma"]])[["Image"]]
                )
              })
            #Add tissue mask
            Target_mask_list$Tissue_mask <- Tissue_mask

            #Generate the final Tissue_mask including all masks generated
            Tissue_mask <- Multi_mask_generator(Target_mask_list)

            #Remove target mask list (large object) and run gc()
            rm(Target_mask_list)
            gc()
          }

          #Get the target channel and remove Image (large object) and run gc()
          Target_Image <- EBImage::getFrame(Image, which(Channels_to_keep == Target_channel))
          rm(Image)
          gc()

          #Finally calculate MFI in mask
          Result <- MFI_calculator(Target = Target_Image,
                                   Tissue_mask = Tissue_mask)
          return(Result)
        }, .progress = TRUE)
      )
    #Return to single core
    future::plan("future::sequential")
    gc()

    #Generate the name of the mask
    if(is.null(Target_masks)) Mask_name <- "Tissue"
    if(!is.null(Target_masks)){
      Mask_name <- c(purrr::map_chr(Target_masks, ~.[["Mask_name"]]), "Tissue")
      Mask_name <- stringr::str_c(Mask_name, collapse = "_")
    }

    #Generate the names of the MFI_score and mask name variables
    MFI_name <- stringr::str_c(c(Target_channel, "_MFI_in_", Mask_name), collapse = "")
    Mask_name <- stringr::str_c(Mask_name, "_Pixels", collapse = "")

    #Generate the final tibble
    RESULTS_tibble <- tibble(Subject_Names = Image_names)
    RESULTS_tibble$Value <-purrr::map_dbl(RESULTS, ~.[["MFI"]])
    RESULTS_tibble$Area <-purrr::map_dbl(RESULTS, ~.[["Total_foreground_pixels"]])
    #Change the names
    names(RESULTS_tibble)[c(2,3)] <- c(MFI_name, Mask_name)

    return(RESULTS_tibble)

  }
