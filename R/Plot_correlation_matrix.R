#' Draws a correlation matrix
#'
#' The function calculates correlation between variables and draws a correlation matrix. Highly correlated or anti-correlated variables can cause trouble in subsequent analyses.
#' Therefore, checking for correlation is a good practice before proceeding with analyses.
#'
#' @param DATA A dataframe or tibble containing cell feature data.
#' @param Variables_included A character vector indicating the names of the features to be used.
#' @param Correlation_method Either "pearson" or "spearman".
#' @returns Generates a correlation matrix plot
#'
#' @examples
#' \dontrun{
#' Plot_correlation_matrix(
#'   DATA = CSM_Arrangedcellfeaturedata_test,
#'   Variables_included = names(CSM_Arrangedcellfeaturedata_test)[-c(1:4)],
#'   Correlation_method = "pearson"
#' )
#' }
#'
#' @export

Plot_correlation_matrix <-
  function(DATA,
           Variables_included,
           Correlation_method) {
    #Check suggested packages
    if(!requireNamespace("corrplot", quietly = FALSE)) stop(
      paste0("corrplot CRAN package is required to execute the function. Please install using the following code: ",
             expression(install.packages("corrplot")))
    )

    #Check Variables included in analysis
    if(!all(Variables_included %in% names(DATA))) {
      Missing_arguments <- Variables_included[!Variables_included %in% names(DATA)]
      stop(paste0(stringr::str_c(Missing_arguments, collapse = ", "), " not found in DATA"))
    }

    if(!all(c(length(Correlation_method) == 1,
              Correlation_method %in% c("pearson", "spearman"))))
    {
      stop("Correlation_method must be one of pearson or spearman")
    }

    #Import data and select only desired variabels
    DATA <- DATA %>% dplyr::select(all_of(Variables_included))
    #Compute correlation matrix
    cor_DATA <- stats::cor(DATA, method = Correlation_method)
    #Compute correlation plot
    corrplot::corrplot(cor_DATA, method = "shade", type = "lower", order = "hclust", addCoef.col = "black", number.cex = 0.8, tl.cex = 0.8,
                       tl.pos = "lt", tl.col = "black")
  }
