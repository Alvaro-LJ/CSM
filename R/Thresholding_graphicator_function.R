#' Generates thresholding results summary plots
#'
#' The function generates plots comparing feature expression values and thresholding results for individual images.
#'
#' @param DATA A dataframe or tibble containing cell feature data.
#' @param DATA_thresholded A dataframe or tibble containing cell feature data that has been thresholded.
#' @param Marker_names A character vector indicating features that need to be plotted.
#' @param Image_name A character value indicating the image to be plotted.
#'
#' @returns Plots individual cells feature expression values and cells that are above threshold
#'
#' @examples
#' \dontrun{
#' #Threshold data
#'DATA_thresholded <- Thresholding_function(
#'   DATA = CSM_Arrangedcellfeaturedata_test,
#'   Strategy = "EBI_Otsu",
#'   Local_thresholding = FALSE,
#'   Method_autothreshold = "Otsu",
#'   number_iterations_TriClass = 20,
#'   Percentile = 0.5,
#'   Defined_threshold = 0.1,
#'   Levels = 3
#' )
#'
#' #Plot results
#' Thresholding_graphicator_function(
#'   DATA = CSM_Arrangedcellfeaturedata_test,
#'   DATA_thresholded = DATA_thresholded,
#'   Marker_names = names(DATA_thresholded)[-c(1:4)],
#'   Image_name = "ABCM22001_B09_MiniCrop.tif"
#')
#' }
#'
#'
#' @export

Thresholding_graphicator_function <-
  function(DATA = NULL,
           DATA_thresholded = NULL,
           Marker_names = NULL,
           Image_name = NULL) {
    #Test the supplied function arguments
    if(!all(c(Image_name %in% DATA$Subject_Names, Image_name %in% DATA_thresholded$Subject_Names))){
      stop("Image_name provided not present in DATA")
    }

    if(!all(c(Marker_names %in% names(DATA), Marker_names %in% names(DATA_thresholded)))){
      stop("Marker_names not present in DATA")
    }

    #First we need to select our image and our markers
    DATA <- DATA %>% dplyr::filter(Subject_Names == Image_name) %>% dplyr::select(1:4, all_of(Marker_names))
    DATA_thresholded <- DATA_thresholded %>% dplyr::filter(Subject_Names == Image_name) %>% dplyr::select(1:4, all_of(Marker_names))

    #Generate the plot of the values by marker
    print("Generating plots")
    PLOT_1 <-
      DATA %>% tidyr::pivot_longer(-c(1:4)) %>%
      ggplot(aes(x = X, y = Y, color = value)) +
      scale_colour_gradient("", low = alpha("white", 0), high = alpha("black", 0.95))+
      facet_wrap(~name, ncol = 1, nrow = ncol(DATA)-4, "free") +
      geom_point(size = 1.5) +
      cowplot::theme_cowplot()+
      theme(axis.line = element_blank(),
            axis.ticks = element_blank(),
            axis.text = element_blank(),
            axis.title = element_blank(),
            legend.text = element_text(size = 10))

    #Generate the plot of positive cells
    PLOT_2 <-
      DATA_thresholded %>% tidyr::pivot_longer(-c(1:4)) %>%
      ggplot(aes(x = X, y = Y, color = as.character(value))) +
      scale_colour_manual(values = c(alpha("white", 0), alpha("black", 0.95)))+
      facet_wrap(~name, ncol = 1, nrow = ncol(DATA)-4, "free") +
      geom_point(size = 1.5) +
      cowplot::theme_cowplot()+
      guides(color = "none")+
      theme(axis.line = element_blank(),
            axis.ticks = element_blank(),
            axis.text = element_blank(),
            axis.title = element_blank())

    #Combine both plots with patchwork
    patchwork::wrap_plots(PLOT_1, PLOT_2, ncol = 2)
  }
