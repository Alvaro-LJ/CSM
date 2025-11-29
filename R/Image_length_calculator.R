#' Generate a summary of the width, length and surface of images in a dataset
#'
#' The function can be used  to plan the image tiling strategy.
#'
#' @param DATA A dataframe or tibble containing cell feature data.
#'
#' @returns A tibble containing image dimension information.
#'
#' @seealso [Suggested_Tile_Size_Calculator()], [Image_tiling_processing_function()]
#'
#' @examples
#' \dontrun{
#' Image_length_calculator(DATA = CSM_Phenotypecell_test)
#' }
#'
#' @export

Image_length_calculator <-
  function(DATA) {
    DF <- dplyr::bind_cols(unique(DATA$Subject_Names),
                          purrr::map_dfr(unique(DATA$Subject_Names), function(Image) {
                            Interim <- DATA %>% dplyr::filter(Subject_Names == Image)
                            c(Width = max(Interim$X)- min(Interim$X),
                              Height = max(Interim$Y) - min(Interim$Y),
                              Surface = (max(Interim$X)- min(Interim$X))*(max(Interim$Y) - min(Interim$Y))
                            )
                          })
    )
    plot(DF %>% tidyr::pivot_longer(-1) %>% ggplot(aes(x = value)) + facet_wrap(~name, "free", ncol = 1, nrow = 3) + geom_histogram(bins = 50) +
           cowplot::theme_cowplot()+
           scale_x_continuous("Size"))
    names(DF)[1] <- "Subject_Names"
    return(DF)
  }
