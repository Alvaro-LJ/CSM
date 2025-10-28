#' Thresholds RGB channels of a color image
#'
#' Intended for internal use only

#' @param Image A EBImage object.
#' @param Target A numeric vector of length = 3. Represents the target RGB channels.
#' @param Color_Tolerance A numeric vector of length = 3. Represents the individual tolerances of target RGB channels.
#' @param Final_Tolerance A numeric value specifying the final tolerance.
#'
#' @details
#' Used in [Channel_deconvolution_function()], [Color_deconvolution_App_launcher()]
#'
#'
#' @returns The function returns a greyscale EBImage object
#' @keywords Internal

EBImage_color_thresholder <-
  function(Image = NULL,
           Target = c(1, 1, 1),
           Color_Tolerance = NULL,
           Final_Tolerance = NULL){
    #Check arguments
    if(!EBImage::is.Image(Image)) stop("Image must be of class EBImage")
    if(!all(is.integer(as.integer(Target)), length(Target) == 3)) stop("Target must contain the target color information in RGB format")
    if(!all(purrr::map_lgl(Target, function(value) value >= 0 & value <= 255))) stop("Target must be of the 8-bit RGB format")
    if(!all(is.numeric(Color_Tolerance), length(Color_Tolerance) == 3, all(Color_Tolerance >= 0), all(Color_Tolerance <= 1))) stop("Color_Tolerance must be a numeric vector of length = 3 with numeric values between 0 and 1")
    if(!all(is.numeric(Final_Tolerance), Final_Tolerance >= 0, Final_Tolerance <= 1))stop("Final_Tolerance must be a single numeric value between 0 and 1")

    #Target colors are transformed to 0-1 scale
    RED <- Target[[1]]/255
    GREEN <- Target[[2]]/255
    BLUE <- Target[[3]]/255

    #Extract each channel
    RED_channel <- EBImage::channel(Image, "red")
    GREEN_channel <- EBImage::channel(Image, "green")
    BLUE_channel <- EBImage::channel(Image, "blue")

    #Each channel will be the difference between the value and the target according to tolerance
    RED_channel <- 1 - abs(RED_channel - RED)
    RED_channel[RED_channel < 1-Color_Tolerance[[1]]] <- 0

    GREEN_channel <- 1 - abs(GREEN_channel - GREEN)
    GREEN_channel[GREEN_channel < 1-Color_Tolerance[[2]]] <- 0

    BLUE_channel <- 1 - abs(BLUE_channel - BLUE)
    BLUE_channel[BLUE_channel < 1-Color_Tolerance[[3]]] <- 0


    #Merge channels to a grey colorscale
    Final_image <-
      EBImage::channel(EBImage::rgbImage(red = RED_channel,
                                         green = GREEN_channel,
                                         blue = BLUE_channel),
                       mode = "grey")
    Final_image[Final_image < 1-Final_Tolerance] <- 0

    #return the final EBImage object
    return(Final_image)
  }
