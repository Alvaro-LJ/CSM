#' Thresholds the foreground pixels of an image according to user parameters into several levels
#'
#' Internal use only
#' @param Target An EBImage type object with a single channel to be thresholded
#' @param Tissue_mask An EBImage type object of the binary tissue mask
#' @param Threshold_type Type of threshold method to be performed: 'Arbitrary' or 'Otsu'
#' @param Threshold_value If Threshold_type = 'Arbitrary', the value applied to threshold the image
#' @param Blurr A logical value indicating if blurring should be applied before thresholding
#' @param Sigma Sigma applied to the Gaussian kernel for blurring
#' @returns An EBImage object containing the binary thresholded image
#' @keywords Internal
#' Thresholds the foreground pixels of an image according to user parameters into several levels
#'
#' Internal use only
#' @param Target An EBImage type object with a single channel to be thresholded
#' @param Tissue_mask An EBImage type object of the binary tissue mask
#' @param Threshold_values The value applied to threshold the image
#' @param Blurr A logical value indicating if blurring should be applied before thresholding
#' @param Sigma Sigma applied to the Gaussian kernel for blurring
#' @returns An EBImage object containing the multi-level thresholded image
#'
#' @details
#' Used in [Image_thresholding_app_launcher()], [Pixel_Threshold_calculator()]
#'
#'
#' @keywords Internal

Pixel_Multilevel_thresholder <-
  function(Target = NULL,
           Tissue_mask = NULL,
           Threshold_values = NULL,
           Blurr = NULL,
           Sigma = NULL){

    #Import Image
    Image <- Target

    #Apply blurring if required
    if(Blurr){
      Image <- EBImage::gblur(Image, sigma = Sigma, boundary = "replicate")
    }

    #Turn to 0 pixels outside tissue mask
    Target[!Tissue_mask] <- 0

    #For every threshold value, calculate a binary image (logical) then sum up the results
    Image <- purrr::reduce(map(Threshold_values, function(Threshold) Image > Threshold),
                           function(Image1, Image2) Image1 + Image2
    )

    Pixel_counts <-purrr::map_dbl(unique(as.vector(Image)), function(Value) sum(Image == Value))
    names(Pixel_counts) <- as.character(unique(as.vector(Image)))

    #If pixel count number is above total pixels then remove pixels from the lowest group (the 0 one)
    if(sum(Pixel_counts) > sum(Tissue_mask)){
      Diff <- sum(Pixel_counts) - sum(Tissue_mask)
      Pixel_counts[1] <- Pixel_counts[1] - Diff
    }

    return(list(
      Image = Image,
      Pixel_count = Pixel_counts,
      Total_foreground_pixels = sum(Tissue_mask),
      Threshold_value = stringr::str_c(Threshold_values, collapse = "_")
    ))
  }
