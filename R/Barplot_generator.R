#' Generate a summary plot of phenotype composition across images
#'
#' After quantifying cell phenotypes in every image using [Phenotype_quantifier()] function, the user can plot the results using this function.
#'
#' @param DATA A dataframe or tibble containing a summary of phenotypes by sample like the one obtained after running [Phenotype_quantifier()].
#' @param Phenotypes_included A character vector indicating the phenotype labels that will be included in the plot.
#' @param Ordering_phenotype A character value indicating the phenotype that will guide barplot ordering.
#'
#' @returns Generates a barplot of the phenotypes composition of every image in DATA.
#'
#' @examples
#' \dontrun{
#' #Generate cell count summary----------------------------------
#' Phenotypes_by_Sample <-
#' Phenotype_quantifier(
#'     DATA = CSM_Phenotypecell_test,
#'     Calculate_Density = FALSE
#')
#'
#' #Generate barplot with desired cell phenotype labels-----------
#' Barplot_generator(
#'     DATA = Phenotypes_by_Sample,
#'     Phenotypes_included = c("TUMOR", "CD8_GZMBpos", "CD8_GZMBneg"),
#'     Ordering_phenotype = "CD8_GZMBpos"
#')
#' }
#'
#' @export

Barplot_generator <-
  function(DATA,
           Phenotypes_included,
           Ordering_phenotype) {

    if(!all(Phenotypes_included %in% names(DATA))){
      stop(paste0("Phenotypes not correctly stated. Choose from ", stringr::str_c(names(DATA)[-1], collapse = ", ")))
    }

    else if(!all(Ordering_phenotype %in% Phenotypes_included)) {
      stop(paste0("Ordering phenotypes not correctly stated. Choose from ", stringr::str_c(names(DATA)[-1], collapse = ", ")))
    }

    else {

      DATA <- DATA  %>% dplyr::select(Subject_Names, all_of(Phenotypes_included))
      DATA <-dplyr::bind_cols(DATA[1],
                              (DATA[-1]/rowSums(DATA[-1]))*100)
      DATA <- DATA[order(DATA[[which(names(DATA) == Ordering_phenotype)]],decreasing=TRUE),]
      DATA <- DATA %>%dplyr::mutate(Subject_Names = factor(Subject_Names, levels = DATA$Subject_Names))

      PLOT <- DATA %>% tidyr::pivot_longer(-1) %>%
        ggplot(aes(x = Subject_Names, y = value, fill = name)) + geom_col(position = "fill", color = "black", width = 0.5) +
        cowplot::theme_cowplot() +
        scale_x_discrete("Image") +
        scale_y_continuous("%")+
        scale_fill_manual("Cell Type", values = unname(pals::polychrome(ncol(DATA)-1)))+
        theme(axis.text.x = element_text(angle = -90, vjust = 0.5, size = 8, color = "black"),
              legend.text = element_text(size = 20),
              legend.title = element_text(size = 20))

      plot(PLOT)
      return(invisible(PLOT))
    }
  }
