#' Generates a tissue mask based on all the channels of an image and thresholding parameters
#'
#' Internal use only
#' @param Image An EBImage type object
#' @param Threshold_type Type of threshold method to be performed: 'Absolute', 'Otsu' or 'Arbitrary'
#' @param Threshold_value If Threshold_type = 'Arbitrary', the value applied to threshold the image
#' @param Blurr A logical value indicating if blurring should be applied before thresholding
#' @param Sigma Sigma applied to the Gaussian kernel for blurring
#' @returns An EBImage object containing the binary tissue mask of foreground pixels
#' @keywords Internal

Tissue_mask_generator <-
  function(Image = NULL,
           Threshold_type = NULL,
           Threshold_value = NULL,
           Blurr = NULL,
           Sigma = NULL){

    #First if Image has more than 1 frame sum all frames into 1 and normalize result
    if(EBImage::numberOfFrames(Image) > 1) {
      Image <- purrr::reduce(EBImage::getFrames(Image),
                             function(Image1, Image2) Image1 + Image2)
      Image <- Image / max(Image)
    }
    #If not normalized to max Image value
    else{
      Image <- Image / max(Image)
    }

    #Perform image blurring if required
    if(Blurr){
      Image <- EBImage::gblur(Image, sigma = Sigma, boundary = "replicate")
    }


    #Perform Otsu thresholding if required
    if(Threshold_type == "Otsu"){
      Otsu_thres <- EBImage::otsu(Image)
      Image <- Image >= Otsu_thres
      return(Image)
    }

    #Perform thresholding according to arbitrary value
    if(Threshold_type == "Arbitrary"){
      Image <- Image >= Threshold_value
      return(Image)
    }

    #Perform thresholding if any pixel is above 0
    if(Threshold_type == "Absolute"){
      Image <- Image > 0
      return(Image)
    }
  }
