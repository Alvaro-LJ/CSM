#' Performs image pixel thresholding according to user preferences
#'
#' Generates binary or multi-level images based on pixel thresholding.
#' @param N_cores Integer. Number of cores to parallelize your computation.
#' @param Directory Character string specifying the path to the folder where images are present.
#' @param Ordered_Channels Character vector specifying image channels in their exact order.
#' @param Channels_to_keep Character vector indicating channels used in tissue mask generation.
#' @param Target_channel Character value indicating the target channel to be thresholded.
#' @param Save_processed_images Logical value indicating if thresholded images should be written (see details).
#' @param Output_directory Character specifying the path to the folder where  output images are written. It must be an empty folder.
#'
#' @param Local_thresholding A logical value indicating if local (per image) threshold should be calculated (see details).
#' @param Threshold_type Type of threshold should be one of the following: 'Arbitrary', 'Otsu' or 'Multilevel' (see details).
#' @param Threshold_value Numeric value to be used as threshold cut-off point. Must be provided for Arbitrary and CAN be provided for Multilevel (as a numeric vector).
#' @param Levels A integer indicating the number of desired levels to be calculated for multilevel thresholding. Only used if Threshold_value is NULL.
#'
#' @param Threshold_type_tissueMask Type of threshold to performe tissue mask. Either 'Otsu', 'Arbitrary' or 'Absolute'.
#' @param Threshold_value_tissueMask Numeric value used if Arbitrary is the threshold type of choice
#' @param Blurr_tissueMask Logical value indicating if image blurring be performed before tissue mask generation
#' @param Sigma_tissueMask Numeric value indicating the sigma value to perform Gaussian blurring
#'
#' @param Blurr_target Logical value indicating if image blurring be performed before target thresholding
#' @param Sigma_target Numeric value indicating the sigma value to perform Gaussian blurring
#'
#' @returns Returns a tibble with the total foreground pixels and pixels above threshold per image
#'
#' @details
#' If processed images are saved, these can be further combined with cell position data using [Cell_to_pixel_distance_calculator()].
#' Local thresholding calculates a threshold value for every image in the experiment. If set to FALSE, global thresholding will be applied, and a unique threshold value will be calculated for all images.
#' Otsu threshold is calculated using the EBImage::otsu function. Automatic multilevel threshold is calculated using the imagerExtra::ThresholdML function.
#'
#'
#' @seealso [Image_thresholding_app_launcher()], [Binary_threshold_image_combinator()], [MFI_Experimet_Calculator()]
#'
#' @examples
#' \dontrun{
#' Pixel_Threshold_calculator(
#' N_cores = 2,
#' Directory = "Image_directory",
#' Ordered_Channels = c("Channel_1", "Channel_2", "Channel_3", "Channel_4", "Channel_5"),
#' Channels_to_keep = c("Channel_1", "Channel_2", "Channel_3", "Channel_4", "Channel_5"),
#' Target_channel = "Channel_1",
#'
#' Save_processed_images = TRUE,
#' Output_Directory = "Output_directory",
#'
#' Local_thresholding = FALSE,
#' Threshold_type = "Arbitrary",
#' Threshold_value = 0.02,
#' Levels = 3,
#'
#' Threshold_type_tissueMask = "Absolute",
#' Threshold_value_tissueMask = 0.001,
#' Blurr_tissueMask = TRUE,
#' Sigma_tissueMask = 0.5,
#'
#' Blurr_target = FALSE,
#' Sigma_target = NULL
#' )
#' }
#'
#' @export

Pixel_Threshold_calculator <-
  function(
    N_cores = 1,
    Directory = NULL,
    Ordered_Channels = NULL,
    Channels_to_keep = NULL,
    Target_channel = NULL,
    Save_processed_images = FALSE,
    Output_Directory = NULL,

    Local_thresholding = NULL,
    Threshold_type = NULL,
    Threshold_value = NULL,
    Levels = NULL,

    Threshold_type_tissueMask = NULL,
    Threshold_value_tissueMask = NULL,
    Blurr_tissueMask = NULL,
    Sigma_tissueMask = NULL,

    Blurr_target = NULL,
    Sigma_target = NULL
  ){
    #Check suggested packages
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
      if(all(Threshold_type == "Multilevel", is.null(Threshold_value))){
        if(!requireNamespace("imagerExtra", quietly = FALSE)) stop(
          paste0("imagerExtra CRAN package is required to execute the function. Please install using the following code: ",
                 expression(install.packages("imagerExtra")))
        )
      }
    }


    on.exit({
      future::plan("future::sequential")
      gc()
    })

    #Argument check
    Argument_checker <- c(N_cores_OK = (N_cores >= 1 & N_cores%%1 == 0),
                          Empty_directory = length(dir(Directory)) >= 1,
                          Channels_OK = all(Channels_to_keep %in% Ordered_Channels),
                          Target_channel_OK = Target_channel %in% Channels_to_keep,
                          Save_processed_images_OK = is.logical(Save_processed_images),
                          Output_Directory_OK = if(Save_processed_images){
                            length(dir(Output_Directory)) == 0
                          }else(TRUE),

                          Local_thresholding_OK = is.logical(Local_thresholding),
                          Threshold_type_OK = Threshold_type %in% c("Otsu", "Arbitrary", "Multilevel"),
                          Threshold_value_OK = if(Threshold_type == "Arbitrary"){
                            all(is.numeric(Threshold_value), Threshold_value >=0, Threshold_value <= 1)
                          } else if(Threshold_type == "Multilevel"){
                            any(is.null(Threshold_value),
                                all(length(Threshold_value) > 1, is.numeric(Threshold_value), dplyr::if_else(Threshold_value >= 0 & Threshold_value <= 1, TRUE, FALSE))
                            )
                          } else(TRUE),
                          Levels_OK = if(Threshold_type == "Multilevel"){
                            all(is.numeric(Levels), Levels%%1 == 0, Levels >= 2)
                          }else(TRUE),
                          Threshold_type_tissueMask_OK = Threshold_type_tissueMask %in% c("Otsu", "Arbitrary", "Absolute"),
                          Threshold_value_tissueMask_OK = if(Threshold_type_tissueMask == "Arbitrary"){
                            all(is.numeric(Threshold_value_tissueMask), Threshold_value_tissueMask >=0, Threshold_value_tissueMask <= 1)
                          }else(TRUE),
                          Blurr_tissueMask_OK = is.logical(Blurr_tissueMask),
                          Sigma_tissueMask_OK = if(Blurr_tissueMask){
                            all(is.numeric(Sigma_tissueMask), Sigma_tissueMask > 0)
                          }else(TRUE),
                          Blurr_target_OK = is.logical(Blurr_target),
                          Sigma_target_OK = if(Blurr_target){
                            all(is.numeric(Sigma_target), Sigma_target > 0)
                          }else(TRUE)
    )

    Stop_messages <- c(N_cores_OK = "N_cores must be an integer value > 0",
                       Empty_directory = "No files found at the directory provided. Please check out the path.",
                       Channels_OK = stringr::str_c(
                         "The following channels are not present the channel names provided: ",
                         stringr::str_c(Channels_to_keep[!(Channels_to_keep %in% Ordered_Channels)], collapse = ", "),
                         sep = ""),
                       Target_channel_OK = stringr::str_c(Target_channel, " not present in ", stringr::str_c(Channels_to_keep, collapse = ", "), collapse = ""),
                       Save_processed_images_OK = "Save_processed_images must be a logical value",
                       Output_Directory_OK = "Output_Directory must be an empty folder",
                       Local_thresholding_OK = "Local_thresholding must be a logical value",
                       Threshold_type_OK = "Threshold_type must be one of the following: Otsu, Arbitrary, Multilevel",
                       Threshold_value_OK = "Threshold_value must be a single numeric value between 0 and 1 if Arbitrary or multiple numeric values between 0 and 1 for Multilevel",
                       Levels_OK = "Levels must be an integer value > 1",
                       Threshold_type_tissueMask_OK = "Threshold_type_TissueMask must be one of the following: Otsu, Arbitrary, Absolute",
                       Threshold_value_tissueMask_OK = "Threshold_value_tissueMask must be a single numeric value between 0 and 1",
                       Blurr_tissueMask_OK = "Blurr_tissueMask must be a logical value",
                       Sigma_tissueMask_OK = "Sigma_tissueMask must be a positive numeric value > 0",
                       Blurr_target_OK = "Blurr_target must be a logical value",
                       Sigma_target_OK = "Sigma_target must be a positive numeric value > 0"
    )
    #Check arguments and stop if necessary
    if(!all(Argument_checker)){
      stop(cat(Stop_messages[!Argument_checker],
               fill = sum(!Argument_checker)))
    }

    #Get the image full directory
    Image_names <- dir(Directory, full.names = TRUE)
    Image_names_short <- dir(Directory, full.names = FALSE)
    Channels_to_keep_index <- which(Channels_to_keep %in% Ordered_Channels)

    #If user has decided to use Multilevel and has supplied LEVELS and Thresholds value, Threshold values will prevail
    if(all(Threshold_type == "Multilevel", !is.null(Threshold_value), !is.null(Levels))){
      message("Both Threshold_value and Levels provided. Threshold_value will be used in image thresholding")
    }

    #The function is branched according to Local vs Global thresholding
    #First local thresholding
    if(Local_thresholding){
      print("Performing Local thresholding")

      #Will iterate for every image
      future::plan("future::multisession", workers = N_cores)
      options(future.globals.maxSize = Inf, future.rng.onMisuse = "ignore")
      furrr::furrr_options(scheduling = Inf)

      RESULTS <- suppressMessages(
        furrr::future_map(seq_along(1:length(Image_names)), function(Image_index){
          #Import the image and keep only required channels and the target channel
          Image <- magick::image_read(Image_names[Image_index])[Channels_to_keep_index]
          Target_Image <- Image[which(Target_channel == Channels_to_keep)]

          Image <- magick::as_EBImage(Image)
          Target_Image <- magick::as_EBImage(Target_Image)

          #First generate the tissue mask and then remove the original image (no longer required)
          Tissue_mask <- Tissue_mask_generator(Image = Image,
                                               Threshold_type = Threshold_type_tissueMask,
                                               Threshold_value = Threshold_value_tissueMask,
                                               Blurr = Blurr_tissueMask,
                                               Sigma = Sigma_tissueMask)
          rm(Image)
          gc()

          #Now we proceed with threshold calculation
          #First binary thresholding
          if(Threshold_type != "Multilevel"){
            Target_Image <- Pixel_thresholder(Target = Target_Image,
                                              Tissue_mask = Tissue_mask,
                                              Threshold_type = Threshold_type,
                                              Threshold_value = Threshold_value,
                                              Blurr = Blurr_target,
                                              Sigma = Sigma_target)
            #save images and tissue masks if required by user
            if(Save_processed_images){
              #Mask
              EBImage::writeImage(Tissue_mask, paste0(Output_Directory, "/", "Processed_", Image_names_short[Image_index], "_Tissue_mask", ".tiff"))
              #Result
              EBImage::writeImage(Target_Image$Image, paste0(Output_Directory, "/", "Processed_",
                                                             stringr::str_replace_all(Image_names_short[Image_index], pattern = "_", replacement = "."),
                                                             "_",
                                                             stringr::str_replace_all(Target_channel, pattern = "_", replacement = "."),  "_BinaryThresholded", ".tiff"))
            }
            #Return the actual image
            return(Target_Image)
          }

          #Multilevel
          if(Threshold_type == "Multilevel"){
            #If threshold values have been supplied use them
            if(!is.null(Threshold_value)){
              Threshold_levels <- base::sort(Threshold_value)
            }
            #Else compute thresholds with imagerExtra
            else{
              Target_Image[!Tissue_mask] <- 0
              Threshold_levels <- imagerExtra::ThresholdML(imager::cimg(array(as.vector(Target_Image), dim = c(1, length(as.vector(Target_Image)), 1, 1))),
                                                           k = (Levels-1),
                                                           returnvalue = TRUE)
            }
            Target_Image <- Pixel_Multilevel_thresholder(Target = Target_Image,
                                                         Tissue_mask = Tissue_mask,
                                                         Threshold_values = Threshold_levels,
                                                         Blurr = Blurr_target,
                                                         Sigma = Sigma_target)
            #save images and tissue masks if required by user
            if(Save_processed_images){
              #Mask
              EBImage::writeImage(Tissue_mask, paste0(Output_Directory, "/","Processed_", Image_names_short[Image_index], "_Tissue_mask", ".tiff"))
              #Result (divided by the number of breaks to get a graylevel image)
              EBImage::writeImage(Target_Image$Image/length(Threshold_levels), paste0(Output_Directory, "/", "Processed_",
                                                                                      stringr::str_replace_all(Image_names_short[Image_index], pattern = "_", replacement = "."),
                                                                                      "_",
                                                                                      stringr::str_replace_all(Target_channel, pattern = "_", replacement = "."),
                                                                                      "_MultiThresholded", ".tiff"))
            }

            return(Target_Image)
          }
        }, .progress = TRUE)
      )

      #Return to single core
      future::plan("future::sequential")
      gc()
    }

    if(!Local_thresholding){
      print("Performing Global thresholding")

      #We will need to calculate and store temporarily the tissue masks, these will be used to calculate thresholds
      future::plan("future::multisession", workers = N_cores)
      options(future.globals.maxSize = Inf, future.rng.onMisuse = "ignore")
      furrr::furrr_options(scheduling = Inf)
      print("Generating tissue masks")

      Tissue_mask_list <- suppressMessages(
        furrr::future_map(seq_along(1:length(Image_names)), function(Image_index){
          Image <- magick::image_read(Image_names[Image_index])[Channels_to_keep_index]
          Image <- magick::as_EBImage(Image)

          #First generate the tissue mask and then remove the original image (no longer required)
          Tissue_mask <- Tissue_mask_generator(Image = Image,
                                               Threshold_type = Threshold_type_tissueMask,
                                               Threshold_value = Threshold_value_tissueMask,
                                               Blurr = Blurr_tissueMask,
                                               Sigma = Sigma_tissueMask)
          rm(Image)
          gc()
          return(Tissue_mask)
        }, .progress = TRUE)
      )
      gc()


      #If threshold type is otsu or arbitrary, then proceed accordingly
      if(Threshold_type != "Multilevel"){

        #Proceed with otsu
        if(Threshold_type == "Otsu"){
          print("Calculating global Otsu threshold")
          #Generate a list that contains the image and the vectorized version of the image where tissue masks have been applied
          Composite_Image_list <- suppressMessages(
            furrr::future_map2(.x = 1:length(Image_names), .y = Tissue_mask_list, function(.x, .y){
              #Import target image
              Image <- magick::image_read(Image_names[.x])[Channels_to_keep_index]
              Target_Image <- Image[which(Target_channel == Channels_to_keep)]
              Target_Image <- magick::as_EBImage(Target_Image)


              #Turn target image values outside tissue mask to 0
              Target_Image[!.y] <- 0

              #return as vector
              return(list(Image = Target_Image,
                          Vector = as.vector(Target_Image)))
            }, .progress = TRUE)
          )
          gc()

          #Generare a unified vector using all images
          Common_vector <- unlist(purrr::map(Composite_Image_list, function(Image) Image[["Vector"]]))
          #Apply otsu algorithm
          Threshold_global_otsu <- EBImage::otsu(array(Common_vector, dim = c(1, length(Common_vector))), range = c(min(Common_vector), max(Common_vector)), levels = length(unique(Common_vector)))

          #Obtain results (arbitrary with global Otsu threshold)
          print("Thresholding images")
          RESULTS <- suppressMessages(
            furrr::future_map2(.x = Composite_Image_list, .y = Tissue_mask_list, function(.x, .y){
              Pixel_thresholder(Target = .x[["Image"]],
                                Tissue_mask = .y,
                                Threshold_type = "Arbitrary",
                                Threshold_value = Threshold_global_otsu,
                                Blurr = Blurr_target,
                                Sigma = Sigma_target)
            }, .progress = TRUE)
          )
          gc()
          #If images need to be stored
          if(Save_processed_images){
            print("Writing tissue mask images")
            #Tissue mask
            suppressMessages(
              furrr::future_map(seq_along(1:length(Tissue_mask_list)), function(Tissue_mask_index){
                EBImage::writeImage(Tissue_mask_list[[Tissue_mask_index]],
                                    paste0(Output_Directory, "/", "Processed_", Image_names_short[Tissue_mask_index], "_Tissue_mask", ".tiff"))
              }, .progress = TRUE)
            )
            gc()
            #Target image
            print("Writing target images")
            suppressMessages(
              furrr::future_map(seq_along(1:length(RESULTS)), function(Result_Image_index){
                EBImage::writeImage(RESULTS[[Result_Image_index]][["Image"]],
                                    paste0(Output_Directory, "/", "Processed_",
                                           stringr::str_replace_all(Image_names_short[Result_Image_index], pattern = "_", replacement = "."),
                                           "_",
                                           stringr::str_replace_all(Target_channel, pattern = "_", replacement = "."), "_BinaryThresholded", ".tiff"))
              }, .progress = TRUE)
            )
            gc()
          }
        }

        #Proceed with user defined threshold
        if(Threshold_type == "Arbitrary"){
          #Get user defined threshold
          Threshold_arbitrary <- Threshold_value
          #Calculate the results by iterating along every image
          print("Thresholding images")
          RESULTS <- suppressMessages(
            furrr::future_map2(.x = seq_along(1:length(Image_names)), .y = Tissue_mask_list, function(.x, .y){
              Target_Image <- magick::image_read(Image_names[.x])[Channels_to_keep_index]
              Target_Image <- Target_Image[which(Target_channel == Channels_to_keep)]
              Target_Image <- magick::as_EBImage(Target_Image)
              Pixel_thresholder(Target = Target_Image,
                                Tissue_mask = .y,
                                Threshold_type = "Arbitrary",
                                Threshold_value = Threshold_arbitrary,
                                Blurr = Blurr_target,
                                Sigma = Sigma_target)
            }, .progress = TRUE)
          )
          gc()
          #If save images is required then write them in the output directory
          if(Save_processed_images){
            print("Writing tissue mask images")
            #Tissue mask
            suppressMessages(
              furrr::future_map(seq_along(1:length(Tissue_mask_list)), function(Tissue_mask_index){
                EBImage::writeImage(Tissue_mask_list[[Tissue_mask_index]],
                                    paste0(Output_Directory, "/", "Processed_", Image_names_short[Tissue_mask_index], "_Tissue_mask", ".tiff"))
              }, .progress = TRUE)
            )
            gc()
            print("Writing target images")
            #Target image
            suppressMessages(
              furrr::future_map(seq_along(1:length(RESULTS)), function(Result_Image_index){
                EBImage::writeImage(RESULTS[[Result_Image_index]][["Image"]],
                                    paste0(Output_Directory, "/", "Processed_",
                                           stringr::str_replace_all(Image_names_short[Result_Image_index], pattern = "_", replacement = "."),
                                           "_",
                                           stringr::str_replace_all(Target_channel, pattern = "_", replacement = "."), "_BinaryThresholded", ".tiff"))
              }, .progress = TRUE)
            )
            gc()
          }
        }
      }
      #If multilevel is required
      if(Threshold_type == "Multilevel"){
        #Generate a list that contains the image and the vectorized version of the image where tissue masks have been applied
        Composite_Image_list <- suppressMessages(
          furrr::future_map2(.x = 1:length(Image_names), .y = Tissue_mask_list, function(.x, .y){
            #Import target image
            Image <- magick::image_read(Image_names[.x])[Channels_to_keep_index]
            Target_Image <- Image[which(Target_channel == Channels_to_keep)]
            Target_Image <- magick::as_EBImage(Target_Image)

            #Turn target image values outside tissue mask to 0
            Target_Image[!.y] <- 0

            #return as vector
            return(list(Image = Target_Image,
                        Vector = as.vector(Target_Image)))
          }, .progress = TRUE)
        )
        gc()
        #Generare a unified vector using all images
        Common_vector <- unlist(purrr::map(Composite_Image_list, function(Image) Image[["Vector"]]))
        #If user has provided thresholds use them, if not compute using Imagerextra
        if(!is.null(Threshold_value)){
          Threshold_global_multilevel <- Threshold_value
        } else{
          print("Calculating Multi level thresholds")
          Threshold_global_multilevel <- imagerExtra::ThresholdML(imager::cimg(array(Common_vector, dim = c(1, length(Common_vector), 1, 1))),
                                                                  k = (Levels-1),
                                                                  returnvalue = TRUE)
        }
        #Obtain results (arbitrary with global Otsu threshold)
        print("Thresholding images")
        RESULTS <- suppressMessages(
          furrr::future_map2(.x = Composite_Image_list, .y = Tissue_mask_list, function(.x, .y){
            Pixel_Multilevel_thresholder(Target = .x[["Image"]],
                                         Tissue_mask = .y,
                                         Threshold_values = Threshold_global_multilevel,
                                         Blurr = Blurr_target,
                                         Sigma = Sigma_target)
          }, .progress = TRUE)
        )
        gc()
        #If images need to be stored
        if(Save_processed_images){
          print("Writing tissue mask images")
          #Tissue mask
          suppressMessages(
            furrr::future_map(seq_along(1:length(Tissue_mask_list)), function(Tissue_mask_index){
              EBImage::writeImage(Tissue_mask_list[[Tissue_mask_index]],
                                  paste0(Output_Directory, "/", "Processed_", Image_names_short[Tissue_mask_index], "_Tissue_mask", ".tiff"))
            }, .progress = TRUE)
          )
          gc()
          print("Writing target images")
          #Target image divided by the number of breaks to get a graylevel image
          suppressMessages(
            furrr::future_map(seq_along(1:length(RESULTS)), function(Result_Image_index){
              EBImage::writeImage(RESULTS[[Result_Image_index]][["Image"]]/length(Threshold_global_multilevel),
                                  paste0(Output_Directory, "/", "Processed_",
                                         stringr::str_replace_all(Image_names_short[Result_Image_index], pattern = "_", replacement = "."),
                                         "_",
                                         stringr::str_replace_all(Target_channel, pattern = "_", replacement = "."),
                                         "_MultiThresholded", ".tiff"))
            }, .progress = TRUE)
          )
          gc()
        }
      }

      #Return to single core
      future::plan("future::sequential")
      gc()
    }

    print("Generating results summary")
    #Generate the final tibble
    Final_tibble <- tibble(Subject_Names = Image_names_short,
                           Total_foreground_pixels =purrr::map_dbl(RESULTS, ~.[["Total_foreground_pixels"]]))

    #Generate the the result thresholds
    if(Threshold_type == "Multilevel") Target_threshold_vector <-purrr::map_chr(RESULTS, ~.[["Threshold_value"]])
    if(Threshold_type != "Multilevel") Target_threshold_vector <-purrr::map_dbl(RESULTS, ~.[["Threshold_value"]])

    #Generate summary
    if(Threshold_type == "Multilevel"){
      #Get the results
      RESULTS <-purrr::map_dfr(RESULTS, function(Image) Image$Pixel_count)
      names(RESULTS) <- stringr::str_c("Value_", names(RESULTS), sep = "")
      #If na turn to 0
      RESULTS[is.na(RESULTS)] <- 0

      #Calculate the multiScore n_pixels*value / total foreground pixels
      Multi_score <-purrr::map2_dfc(.x = RESULTS[-1], .y = 1:ncol(RESULTS[-1]), function(.x, .y) .x*.y)
      Multi_score <- apply(Multi_score, MARGIN = 1, function(Row) sum(Row))/Final_tibble$Total_foreground_pixels

      #Generate the proportion
      RESULTS_PROP <-purrr::map_dfc(RESULTS, function(column) column/Final_tibble$Total_foreground_pixels)
      names(RESULTS_PROP) <- stringr::str_c("PER_", names(RESULTS_PROP), sep = "")
      RESULTS <-dplyr::bind_cols(RESULTS, RESULTS_PROP)

      #Add the multiScore
      RESULTS$Multi_score <- Multi_score
    }
    if(Threshold_type != "Multilevel"){
      RESULTS <- tibble(Positive_pixels =purrr::map_dbl(RESULTS, ~.[["Pixel_count"]]))
      RESULTS$Prop_positive <- RESULTS$Positive_pixels/Final_tibble$Total_foreground_pixels
    }

    #Add the thresholds
    RESULTS$Target_threshold_value <- Target_threshold_vector

    #Return the final value
    return(dplyr::bind_cols(Final_tibble, RESULTS))
  }
