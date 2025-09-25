#' Thresholds RGB channels of a color image and performs pre and post-processing steps
#'
#' Intended for internal use only

#' @param Image A EBImage object.
#' @param Parameters A list of image pre and post-processing
#' @returns The function returns a greyscale EBImage object
#' @keywords internal
#' @export

Channel_deconvolution_function <-
  function(Image = NULL,
           Parameters = NULL){
    Photo <- Image

    #Apply required processing
    Photo <- Photo %>% magick::image_modulate(brightness = Parameters[["Brightness"]], saturation = Parameters[["Saturation"]], hue = Parameters[["Hue"]])
    if(Parameters[["Equalize"]]) Photo <- Photo %>% magick::image_equalize()
    if(Parameters[["Normalize"]]) Photo <- Photo %>% magick::image_normalize()
    if(Parameters[["Contrast"]]) Photo <- Photo %>% magick::image_contrast(sharpen = Parameters[["Sharpen"]])
    if(Parameters[["Reduce_color"]]) Photo <- Photo %>% magick::image_quantize(max = Parameters[["Max_color"]], dither = TRUE)

    #Exctract color
    #Transform to EBI object
    Photo <- Photo %>% magick::image_convert(type = NULL,
                                             colorspace = "RGB",
                                             depth = 8,
                                             antialias = NULL,
                                             matte = FALSE) %>%
      magick::as_EBImage() %>%
      EBImage::toRGB()
    #Extract Image
    Photo <- EBImage_color_thresholder(Image = Photo,
                                       Target = c(Parameters[["RED_value"]], Parameters[["GREEN_value"]], Parameters[["BLUE_value"]]),
                                       Color_Tolerance = Parameters[["Color_Tolerance"]],
                                       Final_Tolerance = Parameters[["Final_Tolerance"]]
    )
    #Apply post processing if required
    #Apply transformations before erosion and dilation
    if(Parameters[["Post_normalize"]]) Photo <- magick::image_read(Photo) %>% magick::image_normalize() %>% magick::as_EBImage()
    if(Parameters[["Post_equalize"]]) Photo <- magick::image_read(Photo) %>% magick::image_equalize() %>% magick::as_EBImage()

    #if erosion is required then execute the rounds of erosion the user requires
    if(Parameters[["Erode"]]){
      Photo <- purrr::reduce(
        .x = seq_along(1:Parameters[["Erode_rounds"]]),
        .f = function(img, ...) EBImage::erode(img, kern = EBImage::makeBrush(Parameters[["Erode_kern"]], "disc")),
        .init = Photo
      )
    }


    #If dilation is required then execute the dilation rounds
    if(Parameters[["Dilate"]]){
      Photo <- purrr::reduce(
        .x = seq_along(1:Parameters[["Dilate_rounds"]]),
        .f = function(img, ...) EBImage::dilate(img, kern = EBImage::makeBrush(Parameters[["Dilate_kern"]], "disc")),
        .init = Photo
      )
    }

    return(Photo)
  }
