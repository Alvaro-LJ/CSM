#' Generate a summary plot of global heterogeneity across images.
#'
#' After quantifying heterogeneity for every image using [Global_heterogeneity_calculator()], the user can plot the results using this function.
#'
#' @param DATA A dataframe or tibble containing a summary of global heterogeneity by sample like the one obtained after running Global_heterogeneity_calculator()].
#' @param Metric A character vector indicating the metric that will be displayed in the plot.
#'
#' @returns Generates a barplot of the heterogeneity of every image in DATA.
#'
#' @examples
#' #Calculate the global heterogeneity by sample----------------------------
#'Global_Heterogeneity_by_sample <-
#'  Global_heterogeneity_calculator(
#'    DATA = CSM_Phenotypecell_test,
#'    Phenotypes_included = unique(CSM_Phenotypecell_test$Phenotype)
#' )
#'
#' #Generate the barplot----------------------------
#'Barplot_Heterogeneity_generator(
#'    DATA = Global_Heterogeneity_by_sample,
#'    Metric = "Gini"
#')
#'
#' @export

Barplot_Heterogeneity_generator <-
  function(DATA = NULL,
           Metric = NULL) {
    if(!(Metric %in% names(DATA))){
      stop(paste0("Metrics not correctly stated. Choose from: ", stringr::str_c(names(DATA)[-1], collapse = ", ")))
    }

    else {

      DATA <- DATA  %>% dplyr::select(Subject_Names, all_of(Metric))
      names(DATA)[-1] <- "metric"

      PLOT <-
        DATA %>%
        ggplot(aes(x = forcats::fct_reorder(Subject_Names, metric), y = metric, fill = metric)) + geom_col(color = "black", width = 0.5) +
        cowplot::theme_cowplot() +
        scale_x_discrete("Image") +
        scale_y_continuous("Global heterogeneity")+
        scale_fill_viridis_c(as.character(Metric)) +
        theme(axis.text.x = element_text(angle = -90, vjust = 0.5, size = 8, color = "black"),
              legend.text = element_text(size = 20),
              legend.title = element_text(size = 20))

      plot(PLOT)
      return(PLOT)
    }

  }
