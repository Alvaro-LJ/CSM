#' Thresholds the foreground pixels of an image according to user parameters
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

Pixel_thresholder <-
  function(Target = NULL,
           Tissue_mask = NULL,
           Threshold_type = NULL,
           Threshold_value = NULL,
           Blurr = NULL,
           Sigma = NULL){

    #Import Image
    Image <- Target

    #Apply blurring if required
    #Perform image blurring if required
    if(Blurr){
      Image <- EBImage::gblur(Image, sigma = Sigma, boundary = "replicate")
    }

    #Turn to 0 pixels outside tissue mask
    Target[!Tissue_mask] <- 0

    #Perform Otsu thresholding if required
    if(Threshold_type == "Otsu"){
      Otsu_thres <- EBImage::otsu(Image)
      Image <- Image >= Otsu_thres
      return(list(Image = Image,
                  Pixel_count = sum(Image),
                  Total_foreground_pixels = sum(Tissue_mask),
                  Threshold_value = Otsu_thres)
      )
    }

    #Perform thresholding according to arbitrary value
    if(Threshold_type == "Arbitrary"){
      Image <- Image >= Threshold_value
      return(list(Image = Image,
                  Pixel_count = sum(Image),
                  Total_foreground_pixels = sum(Tissue_mask),
                  Threshold_value = Threshold_value)
      )
    }
  }
