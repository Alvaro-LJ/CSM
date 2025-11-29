#' Formats cell feature data
#'
#' `Data_arrange_function()` must be run in order to format a cell feature matrix. Many CSM functions won't work if the data has not been adequately formatted.
#'
#' @param DATA A dataframe or tibble containing cell feature data.
#' @param X A character indicating the column name containing X coordinates.
#' @param Y A character indicating the column name containing Y coordinates.
#' @param Subject_Names A character indicating the column name containing image names.
#' @param Markers_to_keep A character vector indicating the names of the columns to be kept in the analysis.
#' @returns Returns a tibble with cell features. The function adds a unique Cell_no to every cell in the dataset.
#'
#' @examples
#' \dontrun{
#' Data_arrange_function(
#'   DATA = CSM_RAWcellfeaturedata_test,
#'   X = 'm.cx',
#'   Y = "m.cy",
#'   Subject_Names = "imageID",
#'   Markers_to_keep = c("CK-EPCAM_AVERAGE", "CD8a_AVERAGE", "GZMB_AVERAGE")
#' )
#' }
#'
#'
#' @export

Data_arrange_function <-
  function(DATA,
           X,
           Y,
           Subject_Names,
           Markers_to_keep) {

  #Check arguments by generating a argument check vector and message vector
  if(!all(c(X, Y, Subject_Names, Markers_to_keep) %in% names(DATA))) {
    Missing_arguments <- c(X, Y, Subject_Names, Markers_to_keep)[!c(X, Y, Subject_Names, Markers_to_keep) %in% names(DATA)]
    stop(paste0(str_c(Missing_arguments, collapse = ", "), " not found in DATA"))
  }

  else{
    #Import X Y ID data
    DATA_interim <- DATA[as.character(c(X, Y, Subject_Names))]
    names(DATA_interim) <- c("X", "Y", "Subject_Names")

    #Bind coordinates and ID with markers
    DATA_interim <- dplyr::bind_cols(DATA_interim, DATA[unique(Markers_to_keep)])
    names(DATA_interim) <- stringr::str_replace_all(names(DATA_interim), "-", "_")

    #Arrange data according to subject Names (required in order to work in the cell ID labeling)
    DATA_interim <- DATA_interim %>% dplyr::arrange(Subject_Names)

    #Assign cell Specific ID
    DATA_interim <- DATA_interim %>% dplyr::mutate(Cell_no = stringr::str_c("CELL", as.character(unlist(purrr::map(
      purrr::map_dbl(unique(DATA_interim$Subject_Names), function(x) {
        nrow(DATA_interim %>% dplyr::filter(Subject_Names == x))
      }), function(x) {
        1:x
      }, .progress = list(clear = F,
                          name = "Adding Cell ID",
                          show_after = 1,
                          type = "iterator")))), sep ="_"))
    #Include the Subject_Name in the Cell_no
    DATA_interim$Cell_no <- stringr::str_c(DATA_interim$Cell_no, DATA_interim$Subject_Names, sep = "__")

    #Final DATA
    DATA_interim <- DATA_interim[c(ncol(DATA_interim), 1:(ncol(DATA_interim)-1))]
    names(DATA_interim) <- stringr::str_replace_all(names(DATA_interim), "-", "_")
    names(DATA_interim) <- stringr::str_replace_all(names(DATA_interim), " ", "_")
    return(DATA_interim)
  }
}
