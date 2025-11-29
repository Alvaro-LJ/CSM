#' Formats image metadata
#'
#' `Clinical_Data_arrange_function()` must be run in order to format image metadata. Image metadata variables can be either numeric or strings. Time to event data must be coded in two separate features (one for event and other for time).
#' This function is similar to [Data_arrange_function()], without adequate formatting, image metadata will no be able to work with [Clinical_Data_analyzer()].
#'
#' @param DATA A dataframe or tibble containing image metadata.
#' @param Subject_Names A character indicating the column name containing image name.
#' @param Outcomes_to_keep A character vector indicating the names of the columns to be kept in the analysis.
#' @returns Returns a tibble with image features.
#'
#' @examples
#' \dontrun{
#' Clinical_Data_arrange_function(
#'      DATA = CSM_ClinicalTMA_test,
#'      Subject_Names = "Sample",
#'      Outcomes_to_keep = c("AGE", "MMRP_status", "DEATH", "OS_m")
#' )
#' }
#'
#' @export

Clinical_Data_arrange_function <-
  function(DATA,
           Subject_Names,
           Outcomes_to_keep) {
    #Check arguments
    if(!all(c(Subject_Names, Outcomes_to_keep) %in% names(DATA))) {
      Missing_arguments <- c(Subject_Names, Outcomes_to_keep)[!c(Subject_Names, Outcomes_to_keep) %in% names(DATA)]
      stop(paste0(stringr::str_c(Missing_arguments, collapse = ", "), " not found in DATA"))
    }

    else{
      #Import X Y ID data
      DATA_interim <- DATA[as.character(Subject_Names)]
      names(DATA_interim) <- "Subject_Names"

      #Bind coordinates and ID with markers
      DATA_interim <- dplyr::bind_cols(DATA_interim, DATA[unique(Outcomes_to_keep)])
      names(DATA_interim) <- stringr::str_replace_all(names(DATA_interim), "-", "_")
      names(DATA_interim) <- stringr::str_replace_all(names(DATA_interim), " ", "_")
      return(DATA_interim)
    }
  }
