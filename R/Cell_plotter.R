#' Plot the cells of a single image
#'
#' The function generates a plot of the cells of a single image. It is similar to [Cell_image_plot_generator()] but the original image is not drawn in the background.
#'
#' @param DATA A dataframe or tibble containing cell feature data.
#' @param Image_name A character value indicating the image to be plotted.
#' @param Phenotypes_included A character vector indicating the phenotype labels that will be included in the plot.
#'
#' @returns Generates a plot of the image cells.
#'
#' @seealso [Cell_image_plot_generator()]
#'
#' @export

Cell_plotter <-
  function(DATA = NULL,
           Image_name = NULL,
           Phenotypes_included = NULL) {
    #Check arguments
    if(!Image_name %in% DATA$Subject_Names){
      stop("Image_name not present in Subject_Names")
    }
    if(!all(Phenotypes_included %in% unique(DATA$Phenotype))){
      Missing_phenotypes <- Phenotypes_included[!Phenotypes_included %in% unique(DATA$Phenotype)]
      stop(paste0(stringr::str_c(Missing_phenotypes, collapse = ", "), " not present in the phenotypes of DATA"))
    }

    DATA %>% dplyr::select(Subject_Names, X, Y, Phenotype) %>% dplyr::filter(Subject_Names == Image_name, Phenotype %in% Phenotypes_included) %>%
      ggplot(aes(x = X, y = Y, color = Phenotype)) +
      geom_point(size = 2, alpha = 0.95)+
      cowplot::theme_cowplot()+
      scale_x_continuous("", labels = NULL) +
      scale_y_continuous("", labels = NULL) +
      scale_color_manual("Cell Type", values = unname(pals::polychrome(length(Phenotypes_included))))+

      guides(color = guide_legend(override.aes = list(size = 12)))+
      theme(axis.line = element_blank(),
            axis.ticks = element_blank(),
            legend.text = element_text(size = 20),
            legend.title = element_text(size = 20))
  }
