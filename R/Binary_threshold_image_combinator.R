#' Performs combination of multiple binary thresholded images using boolean operators
#'
#' This function can be used to generate binary thresholded images that result from combining several thresholded images using the operators AND, OR, NOT
#'
#' @param ... Strings indicating the path to folders containing binary thresholded images. Images from different directories will be matched by name.
#' @param Save_processed_images Logical value indicating if thresholded images should be written (see details).
#' @param Output_directory Character specifying the path to the folder where  output images are written. It must be an empty folder.
#' @param Operators Character vector indicating how the provided directories should be combined. Accepted values are 'AND', 'OR', 'NOT'
#' @param N_cores Integer. Number of cores to parallelize your computation.
#'
#' @returns Returns a tibble with the total positive pixels per image after combining binary thresholded images
#'
#' @details
#' If processed images are saved, these can be further combined with cell position data using [Cell_to_pixel_distance_calculator()].
#'
#' @seealso [Image_thresholding_app_launcher()], [Binary_threshold_image_combinator()], [MFI_Experimet_Calculator()]
#' @export

Binary_threshold_image_combinator <-
  function(...,
           Save_processed_images = NULL, #Should thresholded images saved? These thresholded images can be used in further analyses
           Output_Directory = NULL, #Path where output images should be stored
           Operators = NULL,
           N_cores = NULL){

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
    }


    on.exit({
      future::plan("future::sequential")
      gc()
    })

    #Check the cores
    if(!all(N_cores >= 1 & N_cores%%1 == 0)) stop("N_cores must be an integer value > 0")

    #Generate the directory list
    Directory_list <- list(...)

    #Check Directory list
    #At least two directories provided
    if(length(Directory_list) < 2) stop("At least 2 directories containing images must be present")

    #Check that all directories contain BinaryThresholded images
    if(!all(purrr::map_lgl(Directory_list, function(Directory){
      length(dir(Directory, full.names = FALSE, pattern = "BinaryThresholded")) > 0
    }))){
      Problematic_directories <- Directory_list[map_lgl(Directory_list, function(Directory){
        length(dir(Directory, full.names = FALSE, pattern = "BinaryThresholded")) == 0
      })]
      stop(paste0("The following directories do not contain BinaryThresholded images: ", stringr::str_c(Problematic_directories, collapse = ", ")))
    }

    #Get the image names from every directory
    File_names <-purrr::map(Directory_list, function(Directory){
      dir(Directory, full.names = FALSE, pattern = "BinaryThresholded")
    })
    Full_File_names <-purrr::map(Directory_list, function(Directory){
      dir(Directory, full.names = TRUE, pattern = "BinaryThresholded")
    })

    #Obtain Image_names
    Image_names <-purrr::map(File_names, function(Names){
      Interim_names <- gsub("BinaryThresholded.*", "", Names) #Remove anything after BinaryThresholded
      gsub("(_[^_]*)_?$", "", Interim_names) #Remove the marker information contained between "_"
    })

    #Obtain marker names
    Marker_names <-purrr::map(File_names, function(Names){
      Interim_names <- gsub("BinaryThresholded.*", "", Names) #Remove anything after BinaryThresholded
      Interim_names <- sub(".*_([^_]+)_*$", "\\1", Interim_names) #Obtain the marker information contained between "_"
      unique(Interim_names)
    })
    #Check that each directory contains images from a single marker
    if(!all(purrr::map_lgl(Marker_names, ~length(unique(.)) == 1))){
      Problematic_directories <- Directory_list[map_lgl(Marker_names, ~length(unique(.)) != 1)]
      stop(paste0("The following directories contain thresholded images from several markers: ", stringr::str_c(Problematic_directories, collapse = ", ")))
    }

    #Calculate the intersect of Image_names
    Image_names_intersect <- purrr::reduce(Image_names, function(x, y) intersect(x, y))
    #Check that all intersect names are in all directories provided if not generate a message
    if(!all(purrr::map_lgl(Image_names, ~all(. %in% Image_names_intersect)))){
      Problematic_directories <- Directory_list[!map_lgl(Image_names, ~all(. %in% Image_names_intersect))]
      message("Absence of matching images accross all directories have been found. Only images present in all directories will be used")
    }

    #Check the Operators
    if(!all(Operators %in% c("AND", "OR", "NOT"))) stop("Operators must be any of the following: AND, OR, NOT")
    if(!any(length(Operators) == 1, length(Operators) == length(Directory_list)-1)) stop(paste0("Operators must be a unique value or a vector of length ", length(Directory_list)-1))

    #Check processing image arguments
    if(!is.logical(Save_processed_images)) stop("Save_processed_images must be a logical value")
    if(Save_processed_images){
      if(length(dir(Output_Directory)) >= 1) stop("Output_Directory must be an empty folder")
    }

    #Generate an intercalated vector of markers and operators
    #First increase Operators length to the adequate value if single value provided
    if(length(Operators) == 1) {
      Operators <- rep(Operators, length(Marker_names) - 1)
    }
    Intercalated_vector <- vector("list", length(Marker_names) + length(Operators))
    Intercalated_vector[seq(1, length(Intercalated_vector), by = 2)] <- unlist(Marker_names)
    Intercalated_vector[seq(2, length(Intercalated_vector) - 1, by = 2)] <- Operators
    Intercalated_vector <- unlist(Intercalated_vector)

    #return a message and provide a menu
    message(paste0("The following combination strategy will be applied sequentially: ",
                   stringr::str_c(Intercalated_vector, collapse = "  ")))
    answer <- menu(c("Proceed", "Abort"), title = "Should the thresholded image combination process proceed")
    if(answer == 2) stop("The function has been stopped.")

    #Generate the look up table
    Look_up_table <- suppressMessages(purrr::map_dfc(1:length(File_names), function(List_index){
      Full_File_names[[List_index]][match(Image_names_intersect, Image_names[[List_index]])]
    }))
    names(Look_up_table) <- unlist(Marker_names)

    #Proceed with combination
    #Iterate over look_up table rows
    future::plan("future::multisession", workers = N_cores)
    options(future.globals.maxSize = Inf, future.rng.onMisuse = "ignore")
    furrr::furrr_options(scheduling = Inf)
    RESULTS <- furrr::future_map_dfr(1:length(Image_names_intersect), function(Row){
      Directory_vector <- unlist(Look_up_table[Row,])

      #Obtain the image list
      Images <-purrr::map(Directory_vector, function(Directory){
        EBImage::readImage(Directory) > 0
      })

      #If single operator then we can use reduce
      if(length(unique(Operators)) == 1){
        #Reduce the image accroding to Operator
        Images <- purrr::reduce(Images, function(Image_1, Image_2){
          if(unique(Operators) == "AND") Final_image <- Image_1 & Image_2
          if(unique(Operators) == "OR") Final_image <- Image_1 | Image_2
          if(unique(Operators) == "NOT") Final_image <- Image_1 & !Image_2
          return(Final_image)
        })
      }
      #If multiple operators then proceed according to operators
      else{
        Images <- purrr::reduce2(.x = Images, .y = Operators, function(out, .x, .y){
          if(.y == "AND") Final_image <- out & .x
          if(.y == "OR") Final_image <- out | .x
          if(.y == "NOT") Final_image <- out & !.x
          return(Final_image)
        })
      }

      #If images need to be written, then write
      if(Save_processed_images){
        EBImage::writeImage(Images, paste0(Output_Directory, "/", Image_names_intersect[Row], "_",
                                           "CombinedThreshold", stringr::str_c(gsub("-", ".", Intercalated_vector), collapse = "-"),
                                           ".tiff"))
      }

      #Generate the summary tibble
      return(tibble(Subject_Names = Image_names_intersect[Row],
                    Total_positive_pixels = sum(Images),
                    Combination = stringr::str_c(Intercalated_vector, collapse = "  ")))
    }, .progress = TRUE)
    #Return to single core
    future::plan("future::sequential")
    gc()

    #Return the summary of results
    return(RESULTS)
  }
