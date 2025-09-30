#' Thresholds cell features applying a unique approach for every marker
#'
#'  In contrast to `Thresholding_function()` this function allows selecting a feature specific approach for thresholding. The user can test this approaches using the [Thresholding_tester_app()].
#'  Thresholded features can be further processed using the [Marker_combinator_generator()] and [Phenotype_assigner_function()] to obtain cell phenotypes.
#'
#' @param DATA A dataframe or tibble containing cell feature data.
#' @param Variables_tibble A dataframe or tibble summarizing thresholding approach for every marker. It must contain the following columns: variable, Strategy, Local_thresholding, Method_autothreshold,
#' number_iterations_TriClass, Percentile, Defined_threshold, Levels. Values for every of this columns follow the same rules as for [Thresholding_function()].
#'
#' @details
#' Local thresholding is generally discouraged. It may capture image-specific artifacts or patterns and not true feature positivity.
#'
#' EBI_Otsu calculates Otsu threshold using EBImage::otsu function. Kmeans calculates threshold using the mmand::threshold function. Kmeans_Otsu calculates the threhold using imager::threshold function.
#' Autothreshold calculates thresholds using the autothresholdr::auto_thresh_mask. TriClass_Otsu calculates thresholds using imagerExtra::ThresholdTriclass function.
#' Multi_level thresholds are calculated using imagerExtra::ThresholdML function.
#'
#' @seealso [Thresholding_tester_app()], [Thresholding_function()], [Thresholding_exploration_function()], [Thresholding_summary_function()], [Thresholding_graphicator_function()].
#'
#' @returns Returns a tibble with thresholded cell features. Every feature is thresholded using specific parameters. For binary thresholding features are converted to logical vectors, for multi-level are converted to numeric vectors.
#'
#' @export

Thresholding_function_tailored <-
  function(DATA = NULL,
           Variables_tibble = NULL) {
    if(!identical(c("variable", "Strategy", "Local_thresholding", "Method_autothreshold",
                    "number_iterations_TriClass", "Percentile", "Defined_threshold", "Levels"),
                  names(Variables_tibble))) {
      stop("Names and order of the Variables_tibble must be variable, Strategy, Local_thresholding, Method_autothreshold, number_iterations_TriClass, Percentile, Defined_threshold, Levels")
    }
    if(!all(unique(Variables_tibble$variable) %in% names(DATA))) {
      stop(paste0("Names provided by the Variable column in Variables_tibble should be one of: ", stringr::str_c(names(DATA), collapse = ", ")))
    }

    #Check suggested packages
    {
      if("EBI_Otsu" %in% Variables_tibble$Strategy){
        if(!requireNamespace("EBImage", quietly = TRUE)) stop(
          paste0("EBImage Bioconductor package is required to execute the function. Please install using the following code: ",
                 expression({
                   if (!require("BiocManager", quietly = TRUE))
                     install.packages("BiocManager")

                   BiocManager::install("EBImage")
                 })
          )
        )
      }
      if("Kmeans" %in% Variables_tibble$Strategy){
        if(!requireNamespace("mmand", quietly = FALSE)) stop(
          paste0("mmand CRAN package is required to execute the function. Please install using the following code: ",
                 expression(install.packages("mmand")))
        )
      }
      if("Kmeans_Otsu" %in% Variables_tibble$Strategy){
        if(!requireNamespace("imager", quietly = FALSE)) stop(
          paste0("imager CRAN package is required to execute the function. Please install using the following code: ",
                 expression(install.packages("imager")))
        )
      }
      if("Autothreshold" %in% Variables_tibble$Strategy){
        if(!requireNamespace("autothresholdr", quietly = FALSE)) stop(
          paste0("autothresholdr CRAN package is required to execute the function. Please install using the following code: ",
                 expression(install.packages("autothresholdr")))
        )
      }
      if(any("TriClass_Otsu" %in% Variables_tibble$Strategy, "Multi_level" %in% Variables_tibble$Strategy)){
        if(!requireNamespace("imagerExtra", quietly = FALSE)) stop(
          paste0("imagerExtra CRAN package is required to execute the function. Please install using the following code: ",
                 expression(install.packages("imagerExtra")))
        )
      }
    }

    DATA <- DATA
    dplyr::bind_cols(DATA[1:4],
                     purrr::pmap_dfc(Variables_tibble,
                                     function(variable, Strategy, Local_thresholding, Method_autothreshold,
                                              number_iterations_TriClass, Percentile, Defined_threshold, Levels) {
                                       DATA <- DATA %>% dplyr::select(1:4, all_of(variable))

                                       Thresholding_function(DATA = DATA,
                                                             Strategy = Strategy,
                                                             Local_thresholding = Local_thresholding,
                                                             Method_autothreshold = Method_autothreshold,
                                                             number_iterations_TriClass = number_iterations_TriClass,
                                                             Percentile = Percentile,
                                                             Defined_threshold = Defined_threshold,
                                                             Levels = Levels)[5]
                                     },
                                     .progress = list(clear = F,
                                                      name = "Calculating thresholds",
                                                      show_after = 2,
                                                      type = "iterator"))
    )
  }
