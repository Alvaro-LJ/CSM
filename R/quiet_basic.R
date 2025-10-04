#' Interal function to RUN SPIAT safely and without bugs
#'
#' Intended for internal use only
#'
#' @keywords Internal

quiet_basic <- function(x) {
  sink(tempfile())
  on.exit(sink())
  invisible(force(x))
}
