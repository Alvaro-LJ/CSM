#' Calculates the mean fluorescence intensity of the foreground pixels of an image
#'
#' Internal use only
#' @param Target An EBImage type object with a single channel to be thresholded
#' @param Tissue_mask An EBImage type object of the binary tissue mask
#' @returns A list containing the MFI and the total foreground pixels
#' @keywords Internal

MFI_calculator <-
  function(Target = NULL,
           Tissue_mask = NULL){
    #Replace non tissue pixels for 0
    Target[!Tissue_mask] <- 0
    #Calculate MFI (sum of intensity divided by total foreground pixels)
    return(list(MFI = sum(Target) / sum(Tissue_mask),
                Total_foreground_pixels = sum(Tissue_mask))
    )
  }
