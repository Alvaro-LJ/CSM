#' Calculates a combined tissue mask combining different binary thresholded images
#'
#' Internal use only
#' @param ... EBImage type objects containing binary thresholded tissue masks
#'
#' @details
#' Used in [MFI_Experimet_Calculator()]
#'
#' @returns The result of the intersection of all images provided
#' @keywords Internal

Multi_mask_generator <-
  function(...){
    purrr::reduce(...,
                  function(Image1, Image2) Image1 & Image2)

  }
