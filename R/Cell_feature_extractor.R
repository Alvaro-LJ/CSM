#' Generates a cell feature matrix based on images and cell mask
#'
#' The function computes a cell feature matrix containing spatial coordinates, morphologic features and expression values of cells given an image and its corresponding cell segmentation mask.
#'
#' @param N_cores Integer. Number of cores to parallelize your computation.
#'
#' @param Image_directory Character specifying the path to the folder where images are present.
#' @param Mask_directory Character specifying the path to the folder where cell segmentation masks are present.
#' @param Ordered_Channels Character vector specifying image channels in their exact order.
#'
#' @param Compute_basic_features A logical value indicating if expression values should be computed.
#' @param quantiles_to_calculate A numeric vector containing values between 0 and 1 indicating the expression quantiles to be computed (0.25, 0.5 and 0.75 by default).
#' @param Compute_texture_features A logical value indicating if texture features should be computed for all the markers measured.
#' @param Texture_pixel_distance If Compute_texture_features is TRUE, de pixel distance taken into account to compute the GLCM and texture features.
#'
#' @param Subcellular_compartment_parameter_list A list containing parameters to be employed for subcellular feature expression computation (see details).
#'
#'
#' @details
#' The function first performs name matching between Images and their cell masks. If non-unique matches are found, the function will
#' return an error. By default, the function will compute basic cell shape and location. In addition, the user may request the computation of feature expression data for
#' every cell, the texture features of marker expression and subcellular marker expression. In order to perform subcellular expression computation, the user needs to provide
#' the Subcellular_compartment_parameter_list.
#'
#' The Subcellular_compartment_parameter_list should be a named list. The names of the list should be the name of the subcellular compartment to be analyzed.
#' Each element should be a list with the following named elements: Channel_name, Type, Threshold_type, Threshold_value, Blurr, Sigma.
#' These parameters should follow the same rules as the ones used to calculate tissueMask in the [Pixel_Threshold_calculator()] or the [MFI_Experimet_Calculator()] (see examples).
#'
#' \itemize{
#' \item{Channel_name: A single character value or a character vector indicating the channels present In Ordered_Channels that should be used for subcellular mask generation.}
#' \item{Type: One of the following: Positive - the subcellular mask will be computed and applied directly. Negative - the inverse (negative) of the subcellular mask will be used.}
#' \item{Threshold_type: One of the following: Arbitrary, Otsu.}
#' \item{Threshold_value: If Threshold_type is arbitrary, the threshold value to be applied.}
#' \item{Blurr: A logical value indicating if image blurring should be performed before thresholding.}
#' \item{Sigma: If Blurr is set to TRUE, the Sigma value indicates the sigma for Gaussian blurr filtering.}
#'}
#'
#' @seealso [Watershed_cell_mask_generator()], [Watershed_tester_app()]
#'
#' @returns A tibble containing cell data.
#'
#' @examples
#' \dontrun{
#' #Create temporary input directory----------------------------------------
#' Input_Dir <- tempfile(pattern = "tempdir1_Input")
#' dir.create(Input_Dir, recursive = TRUE)
#'
#' Output_Dir <- tempfile(pattern = "tempdir2_Output")
#' dir.create(Output_Dir, recursive = TRUE)
#'
#' #Save images in Input directory
#' purrr::map(1:2,
#'           function(Image){
#'             EBImage::writeImage(CSM_MiniMultiTiff_test[[Image]],
#'                                 file.path(Input_Dir, names(CSM_MiniMultiTiff_test)[Image]))
#'          })
#'
#' #Check a segmentation parameters list obtained using the dedicated function----------
#' print(CSM_SegmentParams_test)
#'
#'
#' #Generate the whole cell masks using watershed---------------------------------------
#' Watershed_cell_mask_generator(
#'  Directory = Input_Dir,
#'  Output_directory = Output_Dir,
#'  Parameter_list = CSM_SegmentParams_test
#' )
#'
#' #Prepare parameters to obtain subcellular nuclear and cytoplasmic compartments---
#' Subcellular_compartment_parameters <-
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
#' #Obtain the cell feature matrix--------------------------------------------
#' Cell_feature_extractor(Image_directory = Input_Dir,
#'                       Mask_directory = Output_Dir,
#'                       Ordered_Channels = CSM_SegmentParams_test[["Channels_to_keep"]],
#'                       Subcellular_compartment_parameter_list = Subcellular_compartment_parameters)
#'
#' #Remove directories---------------------------------------------------------
#' unlink(Input_Dir, recursive = TRUE)
#' unlink(Output_Dir, recursive = TRUE)
#' }
#'
#' @import patchwork
#' @export


Cell_feature_extractor <-
  function(
    N_cores = 1,

    Image_directory = NULL,
    Mask_directory = NULL,

    Ordered_Channels = NULL,

    Compute_basic_features = TRUE,
    quantiles_to_calculate = c(0.25, 0.5, 0.75),
    Compute_texture_features = FALSE,
    Texture_pixel_distance = 1,

    Subcellular_compartment_parameter_list = NULL
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
    }

    #####Check arguments####
    cat("\n Checking provided arguments")

    if(!all(N_cores >= 1, N_cores%%1 == 0)) stop("N_cores must be an integer value > 0")

    if(!all(dir.exists(Image_directory), length(dir(Image_directory)) >= 1)) stop("Issue with Image_directory provided. Please review.")
    if(!all(dir.exists(Mask_directory), length(dir(Mask_directory)) >= 1)) stop("Issue with Mask_directory provided. Please review.")

    if(!all(length(Compute_basic_features) == 1, is.logical(Compute_basic_features))) stop("Compute_basic_features must be a logical value")
    if(!all(is.numeric(quantiles_to_calculate), quantiles_to_calculate > 0, quantiles_to_calculate < 1)) stop("quantiles_to_calculate must be numeric values between 0 and 1")
    if(!all(length(Compute_texture_features) == 1, is.logical(Compute_texture_features))) stop("Compute_texture_features must be a logical value")
    if(Compute_texture_features){
      if(!all(length(Texture_pixel_distance) == 1, is.numeric(Texture_pixel_distance), Texture_pixel_distance%%1 == 0, Texture_pixel_distance >= 1)) stop("Texture_pixel_distance must be a integer value >= 1")
    }

    #Check specifically the Subcellular_compartment_parameter_list
    if(!is.null(Subcellular_compartment_parameter_list)){
      #Check that the Subcellular_compartment_parameter_list is a list
      if(!is.list(Subcellular_compartment_parameter_list)) stop("Subcellular_compartment_parameter_list must be a named list")

      #Check that the names of the list are unique
      if(any(duplicated(names(Subcellular_compartment_parameter_list)))) stop("Names of Subcellular_compartment_parameter_list must be unique")

      #Check that every item in the list is a list
      if(!all(
        purrr::map_lgl(Subcellular_compartment_parameter_list, function(Individual_mask){
          is.list(Individual_mask)
        }))) stop("Individual Items in the Subcellular_compartment_parameter_list must be a list")

      #Check names of masks (Mask_name, Threshold_type, Threshold_value, Blurr, Sigma)
      Adequate_names <-
        purrr::map_lgl(Subcellular_compartment_parameter_list, function(Individual_mask){
          identical(names(Individual_mask), c("Channel_name", "Type", "Threshold_type", "Threshold_value", "Blurr", "Sigma"))
        })
      if(!all(Adequate_names)){
        stop(paste0("Each item in Subcellular_compartment_parameter_list must contain the folloing: Channel_name, Type, Threshold_type, Threshold_value, Blurr, Sigma",
                    "\nThe following masks have inadequate names: ",
                    stringr::str_c(names(Subcellular_compartment_parameter_list)[!Adequate_names], collapse = ", ")))
      }


      #Check the actual arguments within the list of lists
      purrr::walk(Subcellular_compartment_parameter_list, function(Individual_mask){
        #Check mask name present in channels to keep
        if(!Individual_mask[["Channel_name"]] %in% Ordered_Channels) stop(paste0(Individual_mask[["Channel_name"]], ": Invalid mask name.",
                                                                                 "It must be present in Ordered_Channels: ",
                                                                                 stringr::str_c(Ordered_Channels, collapse = ", ")))

        #Check Type of mask
        if(!Individual_mask[["Type"]] %in% c("Positive", "Negative")) stop(paste0(Individual_mask[["Mask_name"]], ": Invalid Type.",
                                                                                  "It must be one of the following: Positive, Negative"))

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


    #####Generate the look up table####
    cat("\n Matching images and masks")

    Names_Images <- dir(Image_directory, full.names = FALSE)
    FullNames_Images <- dir(Image_directory, full.names = TRUE)
    Names_Mask <- dir(Mask_directory, full.names = FALSE)
    FullNames_Mask <- dir(Mask_directory, full.names = TRUE)

    #Compute the closest Mask name for every Image name
    Closest_name_vector <- purrr::map_dbl(Names_Images, function(Image_name){
      Distance_vector <- adist(Image_name, Names_Mask, fixed = TRUE, ignore.case = TRUE)
      which.min(Distance_vector)
    })

    Look_up_table <- tibble::tibble(Image_names = Names_Images,
                                    Image_paths = FullNames_Images,

                                    Mask_names = Names_Mask[Closest_name_vector],
                                    Mask_paths = FullNames_Mask[Closest_name_vector])

    #If any mask name is duplicated, then print an error
    if(any(duplicated(Look_up_table$Mask_names))){
      Conflictive_masks <- Look_up_table$Mask_names[duplicated(Look_up_table$Mask_names)]

      Colored_print("ERROR!!! Multiple matches between image and mask names. Please review the following conflictive Image/mask multiple matches", color = "red")
      print(Look_up_table[Look_up_table$Mask_names %in% Conflictive_masks, c(1,3)], nrow = 20)

      return(Look_up_table[Look_up_table$Mask_names %in% Conflictive_masks, c(1,3)])
    }

    #####Test the Subcellular_compartment_parameter_list if provided####
    if(!is.null(Subcellular_compartment_parameter_list)){
      cat("\n Checking the Subcellular_compartment_parameter_list provided")

      purrr::walk(seq_along(1:length(Subcellular_compartment_parameter_list)),
                  function(Compartment_test_index){

                    #Print a message
                    cat(paste0("\n Testing ", names(Subcellular_compartment_parameter_list)[Compartment_test_index], " subcellular cell masks"))

                    #Obtain Mask parameters
                    Channels_invoveld_in_mask <- Subcellular_compartment_parameter_list[[Compartment_test_index]][["Channel_name"]]
                    Type_of_mask <- Subcellular_compartment_parameter_list[[Compartment_test_index]][["Type"]]
                    Threshold_type_mask <- Subcellular_compartment_parameter_list[[Compartment_test_index]][["Threshold_type"]]
                    Threshold_value_mask <- Subcellular_compartment_parameter_list[[Compartment_test_index]][["Threshold_value"]]
                    Blurr_mask <- Subcellular_compartment_parameter_list[[Compartment_test_index]][["Blurr"]]
                    Sigma_mask <- Subcellular_compartment_parameter_list[[Compartment_test_index]][["Sigma"]]
                    Channels_to_be_thresholded <- match(x = Channels_invoveld_in_mask,
                                                        table = Ordered_Channels)

                    #Prepare for the while loop
                    All_potential_image_indexes <- 1:nrow(Look_up_table)
                    Prepare_to_exit <- "NO"

                    #Enter the loop
                    while(Prepare_to_exit == "NO" & length(All_potential_image_indexes) >= 1){
                      #Pick a random number
                      Random_index_image <- sample(All_potential_image_indexes, size = 1)
                      #Remove the number from the vector
                      All_potential_image_indexes <- All_potential_image_indexes[-which(All_potential_image_indexes == Random_index_image)]
                      #If no random images left then re start the list
                      if(length(All_potential_image_indexes) < 1){
                        Colored_print("All potential images have been used. The testing image list will be restored")
                        All_potential_image_indexes <- 1:nrow(Look_up_table)
                      }

                      #Select the Image and mask path
                      Image <- EBImage::readImage(Look_up_table$Image_paths[Random_index_image])
                      Mask <- EBImage::readImage(Look_up_table$Mask_paths[Random_index_image])

                      #Check that both have the same dimensions
                      if(!identical(dim(Image)[c(1,2)], #Only the first two dimensions of the image
                                    dim(Mask))){
                        Problematic_image <- Look_up_table$Image_names[Random_index_image]
                        Problematic_mask <- Look_up_table$Mask_names[Random_index_image]

                        stop(paste0(Problematic_image, " dimensions are not equal to ", Problematic_mask, " dimensions. Please review."))
                      }

                      #Modify the mask to harbor only integer values
                      Object_id_tibble <- tibble::tibble(value = unique(sort(as.vector(Mask))))
                      Object_id_tibble <- Object_id_tibble[Object_id_tibble$value != 0,]
                      Object_id_tibble$Object_id <- 1:nrow(Object_id_tibble)
                      Match <- match(Mask, Object_id_tibble$value)
                      Mask[!is.na(Match)] <- Object_id_tibble$Object_id[Match[!is.na(Match)]]


                      #Generate the binary thresholded_mask
                      Subcellular_mask <-
                        Tissue_mask_generator(
                          Image = Image[,,Channels_to_be_thresholded],
                          Threshold_type = Threshold_type_mask,
                          Threshold_value = Threshold_value_mask,
                          Blurr = Blurr_mask,
                          Sigma = Sigma_mask
                        )
                      #If the mask is negative then invert
                      if(Type_of_mask == "Negative") Subcellular_mask <- !Subcellular_mask

                      #Obtain the final mask
                      Final_mask <- Mask * Subcellular_mask

                      #print the result
                      Plot_original <-
                        ggplot() +
                        annotation_raster(EBImage::colorLabels(Mask), xmin = 1, xmax = dim(Mask)[1], ymin = 1, ymax = dim(Mask)[2]) +
                        coord_cartesian(xlim = c(1, dim(Mask)[1]), ylim = c(1, dim(Mask)[2]), expand = FALSE) +
                        ggtitle(paste0(Look_up_table$Image_names[Random_index_image], " - Original cell mask")) +
                        theme(plot.title = element_text(hjust = 0.5))

                      Plot_final <-
                        ggplot() +
                        annotation_raster(EBImage::colorLabels(Final_mask), xmin = 1, xmax = dim(Mask)[1], ymin = 1, ymax = dim(Mask)[2]) +
                        coord_cartesian(xlim = c(1, dim(Mask)[1]), ylim = c(1, dim(Mask)[2]), expand = FALSE) +
                        ggtitle(paste0(Look_up_table$Image_names[Random_index_image], " - ", names(Subcellular_compartment_parameter_list)[Compartment_test_index])) +
                        theme(plot.title = element_text(hjust = 0.5))

                      #Plot the final result
                      plot(Plot_original + Plot_final + patchwork::plot_layout(ncol = 2,nrow = 1))

                      #Ask the user if the algorihtm should proceed
                      answer <- menu(c("Proceed", "Abort", "Test another"), title = "\nShould the analysis proceed?")

                      #If user decides to stop then abort function and return stop message
                      if(answer == 2) stop("The function has been stopped. Please tune parameters for a better result")
                      if(answer == 1) Prepare_to_exit <- "YES!!!"
                      if(answer == 3) cat("\n  Testing another random sample")
                    }
                  })
    }

    ####COMPUTE CELL FEATURE TIBBLE FOR EVERY IMAGE####
    print("Computing cell features...")

    #Set the on exit statement
    on.exit({
      future::plan("future::sequential")
      gc()
    })

    #Set the cores
    future::plan("future::multisession", workers = N_cores)
    options(future.globals.maxSize = Inf, future.rng.onMisuse = "ignore")
    furrr::furrr_options(scheduling = Inf)

    Final_result <-
      furrr::future_map_dfr(1:nrow(Look_up_table),
                            function(Image_index){

                              #Select the Image and mask path
                              Image <- EBImage::readImage(Look_up_table$Image_paths[Image_index])
                              Mask <- EBImage::readImage(Look_up_table$Mask_paths[Image_index])

                              #Check that both have the same dimensions
                              if(!identical(dim(Image)[c(1,2)], #Only the first two dimensions of the image
                                            dim(Mask))){
                                Problematic_image <- Look_up_table$Image_names[Image_index]
                                Problematic_mask <- Look_up_table$Mask_names[Image_index]

                                Colored_print(paste0("Issue with ", Problematic_image, ". Dimensions between Image and mask do not match"), color = "red")
                                return(tibble::tibble()) #Return an empty tibble
                              }

                              #Check if image mask is empty, if its empty return a message and an empty tibble
                              if(max(Mask) == 0){
                                Colored_print(paste0("Issue with ", Look_up_table$Mask_names[Image_index], ". no cell mask identified"), color = "red")
                                return(tibble::tibble()) #Return an empty tibble
                              }

                              #Modify the mask to harbor only integer values
                              Object_id_tibble <- tibble::tibble(value = unique(sort(as.vector(Mask))))
                              Object_id_tibble <- Object_id_tibble[Object_id_tibble$value != 0,]
                              Object_id_tibble$Object_id <- 1:nrow(Object_id_tibble)
                              Match <- match(Mask, Object_id_tibble$value)
                              Mask[!is.na(Match)] <- Object_id_tibble$Object_id[Match[!is.na(Match)]]

                              #If subcellular parameter list provided then generate the list with the new masks
                              if(is.null(Subcellular_compartment_parameter_list)) SubCellular_compartment_image_list <- NULL

                              if(!is.null(Subcellular_compartment_parameter_list)){
                                #Generate a list containing the new masks
                                SubCellular_compartment_image_list <-
                                  purrr::map(seq_along(1:length(Subcellular_compartment_parameter_list)),
                                             function(Compartment_test_index){
                                               #Obtain Mask parameters
                                               Channels_invoveld_in_mask <- Subcellular_compartment_parameter_list[[Compartment_test_index]][["Channel_name"]]
                                               Type_of_mask <- Subcellular_compartment_parameter_list[[Compartment_test_index]][["Type"]]
                                               Threshold_type_mask <- Subcellular_compartment_parameter_list[[Compartment_test_index]][["Threshold_type"]]
                                               Threshold_value_mask <- Subcellular_compartment_parameter_list[[Compartment_test_index]][["Threshold_value"]]
                                               Blurr_mask <- Subcellular_compartment_parameter_list[[Compartment_test_index]][["Blurr"]]
                                               Sigma_mask <- Subcellular_compartment_parameter_list[[Compartment_test_index]][["Sigma"]]
                                               Channels_to_be_thresholded <- match(x = Channels_invoveld_in_mask,
                                                                                   table = Ordered_Channels)

                                               #Generate the binary thresholded_mask
                                               Subcellular_mask <-
                                                 Tissue_mask_generator(
                                                   Image = Image[,,Channels_to_be_thresholded],
                                                   Threshold_type = Threshold_type_mask,
                                                   Threshold_value = Threshold_value_mask,
                                                   Blurr = Blurr_mask,
                                                   Sigma = Sigma_mask
                                                 )
                                               #If the mask is negative then invert
                                               if(Type_of_mask == "Negative") Subcellular_mask <- !Subcellular_mask

                                               #Obtain the final mask and return it
                                               Final_mask <- Mask * Subcellular_mask
                                               return(Final_mask)
                                             })
                                #change the names
                                names(SubCellular_compartment_image_list) <- names(Subcellular_compartment_parameter_list)

                                #If any of the masks is empty, remove it from the analysis
                                Any_Subcellular_empty <- purrr::map_lgl(SubCellular_compartment_image_list, ~max(.) == 0)

                                if(any(Any_Subcellular_empty)){
                                  Problematic_image <- Look_up_table$Image_names[Image_index]
                                  Problematic_subcellular <- names(SubCellular_compartment_image_list)[Any_Subcellular_empty]

                                  Colored_print(paste0(Problematic_image, " - Empty mask for ", Problematic_subcellular, ". No subcompartment features will be computed"), color = "red")

                                  #Keep the adequate compartments
                                  SubCellular_compartment_image_list <- SubCellular_compartment_image_list[!Any_Subcellular_empty]
                                }
                              }


                              #Obtain the tibble with Cell_image_mask_profiler
                              Final_tibble <-
                                Cell_image_mask_profiler(
                                  Image = Image,
                                  Mask = Mask,

                                  Ordered_Channels = Ordered_Channels,

                                  Compute_basic_features = Compute_basic_features,
                                  quantiles_to_calculate = quantiles_to_calculate,
                                  Compute_texture_features = Compute_texture_features,
                                  Texture_pixel_distance = Texture_pixel_distance,

                                  SubCellular_compartment_list = SubCellular_compartment_image_list
                                )

                              #Remove elements just in case
                              rm(Image)
                              rm(Mask)
                              rm(SubCellular_compartment_image_list)
                              gc()


                              #Add to the final tibble the image id column
                              Final_tibble <- dplyr::bind_cols(tibble::tibble(image_id = Look_up_table$Image_names[Image_index]),
                                                               Final_tibble)

                              #Return the result
                              return(Final_tibble)

                            }, .progress = TRUE)

    future::plan("future::sequential")
    gc()

    #Return the result
    return(Final_result)
  }




