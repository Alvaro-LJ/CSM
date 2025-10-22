#' Splits cell feature matrix
#'
#' The function is used to split cell feature matrix into a tibble containing features that are desired to be used and features that are required to be set aside.
#' Feature that are set aside are stored in a tibble containing cell_no and spatial coordinates.
#' Features set aside can be used in fine tuning of cell phenotypes using [ReClustering_function()] or functional analysis using [Cell_functional_assessment()]
#'
#' @param DATA A dataframe or tibble containing cell feature data.
#' @param Markers_to_set Character vector indicating which cell features will be placed aside.
#' @returns A list with two elements: Aside (containing markers removed from main cell feature dataset) and DATA containing remaining features
#'
#' @examples
#' Data_set_aside(
#'   DATA = CSM_Arrangedcellfeaturedata_test,
#'   Markers_to_set = "GZMB_AVERAGE"
#' )
#'
#' @export

Data_set_aside <-
  function(DATA = NULL,
           Markers_to_set = NULL) {

    if(!identical(c("Cell_no", "X", "Y", "Subject_Names"), names(DATA)[c(1:4)])) {
      stop("Your data does not contain adequate format (Cell_no, X, Y, Subject_Names). Please format using the Data_arrange_function.")
    }

    if(!all(Markers_to_set %in% names(DATA))) {
      Missing_arguments <- c(Markers_to_set)[!Markers_to_set %in% names(DATA)]
      stop(paste0(stringr::str_c(Missing_arguments, collapse = ", "), " not found in DATA"))
    }

    else{list(Aside = DATA %>% dplyr::select(1:4, all_of(Markers_to_set)),
              DATA = DATA %>% dplyr::select(-any_of(Markers_to_set))
    )}
  }
