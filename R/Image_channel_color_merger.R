#' Merges image channels using different colors
#'
#' The function generates and RGB image resulting from merging various image channels.
#'
#' @param Image A magick image object generated using magick::image_read().
#' @param Parameter_list A list containing merging parameters (see details and examples below).
#'
#' @returns A magick image object with the channels merged.
#'
#' @details
#' The parameter list must be a list of lists. Every element in parameter list must be named. Individual list elements within the Parameter_list must contain
#' the parameters used to process that specific channel before merging. This include indicating the channel being processed (a numeric value indicating the channel position within the image),
#' the color to be applied to the channel (indicated by name or HEX code), the gamma correction (a numeric value between 0.001 and 1000), the min and max pixel value (a numeric value between 0 and 100 for both of them) and
#' equalize (a logical value).
#'
#' @examples
#' \dontrun{
#' #Read the image to be processed----------------------------------------------
#' Image <- magick::image_read(CSM_MiniMultiTiff_test[[2]])
#' 
#' #Generate the merging parameters (in this case 4 channels will be merged)----
#' Merge_parameter_list <- 
#'   list(
#'       DAPI = list(channel_index = 1,
#'                   color = "blue",
#'                   gamma = 1,
#'                   min = 0,
#'                   max = 80,
#'                   equalize = FALSE),
#'
#'       CK = list(channel_index = 5,
#'                 color = "white",
#'                 gamma = 1.1,
#'                 min = 5,
#'                 max = 20,
#'                 equalize = FALSE),
#'  
#'       CD8 = list(channel_index = 6,
#'                  color = "#f0703a",
#'                  gamma = 1,
#'                  min = 10,
#'                  max = 30,
#'                  equalize = FALSE),
#'
#'       GZMB = list(channel_index = 3,
#'                   color = "cyan",
#'                   gamma = 1,
#'                   min = 5,
#'                   max = 20,
#'                   equalize = FALSE)
#'
#' )
#'
#' #Perform the merging
#' Image_channel_color_merger(Image = Image,
#'                           Parameter_list = Merge_parameter_list)
#' 
#'}
#' @export

Image_channel_color_merger <- 
  function(Image,
           Parameter_list){
    
    ########Check suggested packages########
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
      
    }
    
    ########Check arguments########
    if(class(Image) != "magick-image") stop("Image must be imported using magick::image_read()")
    
    #Check the Parameter_list
    if(!is.list(Parameter_list)) stop("Parameter_list must be a list")
    if(length(Parameter_list) > 7) warning("Merging more than 7 channels is strongly discouraged")
    if(any(names(Parameter_list) == "")) stop("Lists in Parameter_list must be named")
    if(length(unique(names(Parameter_list))) != length(Parameter_list)) stop("Names in Parameter_list must not be duplicated")
    #Check that all elements in parameter list contain the adequate elements
    Conflictive_elements <- purrr::map_lgl(Parameter_list, ~identical(names(.), c("channel_index", "color", "gamma", "min", "max", "equalize")))
    if(any(!Conflictive_elements)){
      Conflictive_items <- names(Parameter_list)[!Conflictive_elements]
      stop(paste0("All list elements in Parameter_list must contain the following elements channel_index, color, gamma, min, max, equalize", "\n",
                  "The following elements do not contain all the required elements: ", "\n",
                  str_c(Conflictive_items, collapse = "\n")))
    }
    
    #Check if channel index is possible
    Conflictive_channel_index <- purrr::map_lgl(Parameter_list, ~.[["channel_index"]] <= length(Image))
    if(any(!Conflictive_channel_index)){
      Conflictive_items <- names(Parameter_list)[!Conflictive_channel_index]
      stop(paste0("channel_index must be a number between 1 and ", length(Image),  "\n",
                  "The following elements do not contain adequate channel_index: ", "\n",
                  str_c(Conflictive_items, collapse = "\n")))
    }
    
    #Check that colors are adequate
    Colors_chosen <- purrr::map_chr(Parameter_list, ~.[["color"]])
    Error_in_color <- purrr::map_lgl(Colors_chosen, function(color) berryFunctions::is.error(grDevices::col2rgb(color)))
    if(any(Error_in_color)) {
      Conflictive_items <- names(Parameter_list)[Error_in_color]
      stop(paste0("color must have a valid name or HEX value", "\n",
                  "The following elements do not contain adequate color: ", "\n",
                  str_c(Conflictive_items, collapse = "\n")))
    }
    #Check that colors are not duplicated
    if(any(duplicated(Colors_chosen))){
      Conflictive_items <- names(Parameter_list)[duplicated(Colors_chosen)]
      stop(paste0("color must be unique", "\n",
                  "The following elements do not contain unique colors: ", "\n",
                  str_c(Conflictive_items, collapse = "\n")))
    }
    
    #Check gamma value
    Gamma_chosen_numeric <- purrr::map_lgl(Parameter_list, ~is.numeric(.[["gamma"]]))
    if(any(!Gamma_chosen_numeric)){
      Conflictive_items <- names(Parameter_list)[!Gamma_chosen_numeric]
      stop(paste0("gamma must be a numeric value between 0.001 and 1000 ",  "\n",
                  "The following elements do not contain adequate gamma value: ", "\n",
                  str_c(Conflictive_items, collapse = "\n")))
    }
    Gamma_chosen <- purrr::map_dbl(Parameter_list, ~.[["gamma"]])
    if(any(c(Gamma_chosen < 0.001, Gamma_chosen > 1000))){
      Conflictive_items <- names(Parameter_list)[Gamma_chosen < 0.001 | Gamma_chosen > 1000]
      stop(paste0("gamma must be a numeric value between 0.001 and 1000 ",  "\n",
                  "The following elements do not contain adequate gamma value: ", "\n",
                  str_c(Conflictive_items, collapse = "\n")))
    }
    
    #Check min and max
    Min_chosen_numeric <- purrr::map_lgl(Parameter_list, ~is.numeric(.[["min"]]))
    if(any(!Min_chosen_numeric)){
      Conflictive_items <- names(Parameter_list)[!Min_chosen_numeric]
      stop(paste0("min must be a numeric value between 0 and 100 and smaller than max",  "\n",
                  "The following elements do not contain adequate min value: ", "\n",
                  str_c(Conflictive_items, collapse = "\n")))
    }
    Max_chosen_numeric <- purrr::map_lgl(Parameter_list, ~is.numeric(.[["max"]]))
    if(any(!Max_chosen_numeric)){
      Conflictive_items <- names(Parameter_list)[!Max_chosen_numeric]
      stop(paste0("max must be a numeric value between 0 and 100",  "\n",
                  "The following elements do not contain adequate max value: ", "\n",
                  str_c(Conflictive_items, collapse = "\n")))
    }
    Min_Max_tibble <- tibble::tibble(Channel = names(Parameter_list),
                                    Min = purrr::map_dbl(Parameter_list, ~.[["min"]]),
                                    Max = purrr::map_dbl(Parameter_list, ~.[["max"]])
    )
    if(any(Min_Max_tibble$Min > 100 | Min_Max_tibble$Min < 0)){
      Conflictive_items <- names(Parameter_list)[Min_Max_tibble$Min > 100 | Min_Max_tibble$Min < 0]
      stop(paste0("min must be a numeric value between 0 and 100 and smaller than max",  "\n",
                  "The following elements do not contain adequate min value: ", "\n",
                  str_c(Conflictive_items, collapse = "\n")))
    }
    if(any(Min_Max_tibble$Max > 100 | Min_Max_tibble$Max < 0)){
      Conflictive_items <- names(Parameter_list)[Min_Max_tibble$Max > 100 | Min_Max_tibble$Max < 0]
      stop(paste0("max must be a numeric value between 0 and 100", "\n",
                  "The following elements do not contain adequate max value: ", "\n",
                  str_c(Conflictive_items, collapse = "\n")))
    }
    if(any(Min_Max_tibble$Max <= Min_Max_tibble$Min)){
      Conflictive_items <- names(Parameter_list)[Min_Max_tibble$Max <= Min_Max_tibble$Min]
      stop(paste0("min must be a numeric value between 0 and 100 and smaller than max", "\n",
                  "The following elements do not contain adequate min or max value: ", "\n",
                  str_c(Conflictive_items, collapse = "\n")))
    }
    
    #Check equalize
    Equalize_logical <- purrr::map_lgl(Parameter_list, ~is.logical(.[["equalize"]]))
    if(any(!Equalize_logical)){
      Conflictive_items <- names(Parameter_list)[!Equalize_logical]
      stop(paste0("equalize must be a logical value",  "\n",
                  "The following elements do not contain adequate equalize value: ", "\n",
                  str_c(Conflictive_items, collapse = "\n")))
    }
    
    ########Reduce the targets########
    Image <- purrr::map(Parameter_list, function(Channel_parameters){
      Channel <- Image[Channel_parameters[["channel_index"]]]
      
      #Perform image equalization as requested by user
      if(Channel_parameters[["equalize"]]) Channel <- Channel %>% magick::image_equalize()
      #Perform image white adjustment
      Channel <- Channel %>%
        magick::image_level(black_point = Channel_parameters[["min"]],
                            white_point = Channel_parameters[["max"]],
                            mid_point = Channel_parameters[["gamma"]]
        )
      #Change to EBImage object
      Channel <- Channel %>% magick::as_EBImage()
      
      #Transform to grayscale just in case
      EBImage::colorMode(Channel) <- "Grayscale"
      
      #Get the color RGB parameters
      RGB_values <- grDevices::col2rgb(Channel_parameters[["color"]])
      
      Red_coeff <- RGB_values[1]/255
      Green_coeff <- RGB_values[2]/255
      Blue_coeff <- RGB_values[3]/255
      
      
      #Modify the channel color
      Channel <- EBImage::channel(Channel, "gray")
      Channel <- EBImage::rgbImage(
        red  = Red_coeff * Channel,
        green = Green_coeff * Channel,
        blue  = Blue_coeff * Channel
      )
      
      return(Channel)
    })
    
    #Add the channels into a single image
    Image <- purrr::reduce(Image, function(Channel_1, Channel_2) Channel_1 + Channel_2)
    
    #Return as a magick image
    return(magick::image_read(Image))
  }






