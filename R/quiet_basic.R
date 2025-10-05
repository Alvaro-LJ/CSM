#' Interal function to RUN SPIAT safely and without bugs
#'
#' Intended for internal use only
#'
#' @details
#' Used in [Own_grid_metrics()], [SPIAT_Heterogeneity_Analyzer()]
#'
#'
#' @keywords Internal

quiet_basic <- function(x) {
  sink(tempfile())
  on.exit(sink())
  invisible(force(x))
}
