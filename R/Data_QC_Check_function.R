#' Checks for potential flaws in cell feature data
#'
#' `Data_QC_Check_function()` checks that there are no missing values in the cell feature dataset and that there are no constant value variables.
#' @param DATA A dataframe or tibble containing cell feature data.
#' @returns A message indicating potential issues in cell feature dataset.
#'
#' @examples
#' \dontrun{
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
#' Data_QC_Check_function(
#'DATA = CSM_Arrangedcellfeaturedata_test %>% dplyr::mutate(X = -X)
#')
#' }
#'
#' @export

Data_QC_Check_function <-
  function(DATA) {
    if(all(sum(is.na(DATA)) == 0,
            sum(purrr::map_lgl(DATA[-c(1:4)], function(x) length(unique(x)) >1)) == ncol(DATA)-4,
            DATA[["X"]] < 0,
            DATA[["Y"]] < 0)) {
      Colored_print("No missing values in data. No constant value markers. X and Y coordinates are positive. Proceed with analysis", color  = "green")
    }
    if(sum(is.na(DATA)) != 0){
      Colored_print(paste0("WARNING!!!!! Missing values in col: ", paste(names(DATA)[purrr::map_dbl(DATA, function(x) sum(is.na(x)))>0], sep = ", ")), 
                    color = "red")
      
    }
    if(sum(purrr::map_lgl(DATA[-c(1:4)], function(x) length(unique(x)) >1)) < ncol(DATA)-4) {
      Colored_print(paste0("WARNING!!!!! constant values in col: ", paste(names(DATA)[-c(1:4)][purrr::map_lgl(DATA[-c(1:4)], function(x) length(unique(x)) == 1)], sep = ", ")),
                    color = "red")
      
    }
    if(any(DATA[["X"]] < 0)){
      Conflictive_samples <- unique(DATA$Subject_Names[DATA[["X"]] < 0])
      Colored_print(paste0("WARNING!!!!! Negative X coord values present in the following samples: ", stringr::str_c(Conflictive_samples, collapse = ", ")),
                    color = "red")
    }
    if(any(DATA[["Y"]] < 0)){
      Conflictive_samples <- unique(DATA$Subject_Names[DATA[["Y"]] < 0])
      Colored_print(paste0("WARNING!!!!! Negative Y coord values present in the following samples: ", stringr::str_c(Conflictive_samples, collapse = ", ")),
                    color = "red")
    }
  }
