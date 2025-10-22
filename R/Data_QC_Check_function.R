#' Checks for potential flaws in cell feature data
#'
#' `Data_QC_Check_function()` checks that there are no missing values in the cell feature dataset and that there are no constant value variables.
#' @param DATA A dataframe or tibble containing cell feature data.
#' @returns A message indicating potential issues in cell feature dataset.
#'
#' @examples
#' #No QC flags-------------------------
#' Data_QC_Check_function(
#' DATA = CSM_Arrangedcellfeaturedata_test
#' )
#'
#' #Constant value flag----------------
#' Data_QC_Check_function(
#' DATA = CSM_Arrangedcellfeaturedata_test %>% dplyr::mutate(Constant_variable = 1)
#' )
#'
#' #Missing value flag----------------
#' Data_QC_Check_function(
#' DATA = CSM_Arrangedcellfeaturedata_test %>% dplyr::mutate(Missing_variable = NA)
#' )
#'
#' @export

Data_QC_Check_function <-
 function(DATA) {
  if( (sum(is.na(DATA)) == 0) && (sum(purrr::map_lgl(DATA[-c(1:4)], function(x) length(unique(x)) >1)) == ncol(DATA)-4) ) {
    print("No missing values in data. No constant value markers. Proceed with analysis")
  }
  else if(sum(is.na(DATA)) != 0 && (sum(purrr::map_lgl(DATA[-c(1:4)], function(x) length(unique(x)) >1)) == ncol(DATA)-4)){
    paste0("Missing values in col: ", names(DATA)[purrr::map_dbl(DATA, function(x) sum(is.na(x)))>0])
  }
  else if(sum(is.na(DATA)) == 0 && sum(purrr::map_lgl(DATA[-c(1:4)], function(x) length(unique(x)) >1)) < ncol(DATA)-4) {
    paste0("constant values in col: ", names(DATA)[-c(1:4)][purrr::map_lgl(DATA[-c(1:4)], function(x) length(unique(x)) == 1)])
  }
  else {
    c(
      paste0("Missing values in col: ", names(DATA)[purrr::map_dbl(DATA, function(x) sum(is.na(x)))>0]),
      paste0("constant values in col: ", names(DATA)[-c(1:4)][purrr::map_lgl(DATA[-c(1:4)], function(x) length(unique(x)) == 1)])
    )
  }
}
