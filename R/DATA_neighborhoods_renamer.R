#' Modifies neighborhood labels
#'
#' The function allows the user to modify neighborhood labels. The function works with any dataframe or tibble containing a column named 'Neighborhood_assignment'.
#'
#' @param DATA A dataframe or tibble containing neighborhood data, including a column named 'Neighborhood_assignment'.
#' @param New_names A character vector indicating the new names. The length must be equal to the number of unique neighborhood labels.
#' @param Old_names (OPTIONAL) A character vector indicating the neighborhood labels that need to be modified. The length of this vector must be equal to the length of New_names.
#'
#' @returns Returns a tibble with cell features and the Neighborhood_assignment column with the modified labels.
#'
#' @examples
#' \dontrun{
#' print(unique(CSM_Neighborhoods_test$Neighborhood_assignment))
#'
#' New_names_DATA <-
#'  DATA_neighborhoods_renamer(
#'     DATA = CSM_Neighborhoods_test,
#'     New_names = c("Tumor_low_density", "Tumor_high_density"),
#'     Old_names = c("Tumor_Other", "Tumor_rich")
#' )
#'
#' print(unique(New_names_DATA$Neighborhood_assignment))
#' }
#'
#' @export

DATA_neighborhoods_renamer <-
  function(DATA = NULL,
           New_names = NULL,
           Old_names = NULL) {
    #Check that DATA has a neighborhood_assignment variable
    if(!"Neighborhood_assignment" %in% names(DATA)) stop("DATA must be generated using Neighborhood_discovery_function or UTAG_Neighborhood_identifier")


    #If new names need to be assigned then proceed as usual
    if(is.null(Old_names)){
      #Check if provided names are equal to number of hoods
      if(length(New_names) != length(unique(DATA$Neighborhood_assignment))) {
        stop(paste0("Provided New_names should match the number of Neighborhoods in the analysis. Number of neighborhoods: ", length(unique(DATA$Neighborhood_assignment)),
                    ". Names provided: ", length(New_names)))
      }
      #Proceed with renaming
      DATA_neighborhoods <- DATA
      names_tibble <- tibble(Neighborhood_assignment = factor(1:length(unique(DATA_neighborhoods$Neighborhood_assignment))),
                             New_names = New_names)
      Result <- dplyr::left_join(DATA_neighborhoods, names_tibble, by = "Neighborhood_assignment")
      return(Result %>%dplyr::mutate(Neighborhood_assignment = New_names) %>% dplyr::select(-New_names))
    }

    #If new names and old names are provided then check first that the lengths are equal and check that old names are present in data. Then proceed
    if(!is.null(Old_names)){
      if(!all(Old_names %in% unique(DATA$Neighborhood_assignment))){
        Conflictive_names <- Old_names[!Old_names %in% unique(DATA$Neighborhood_assignment)]
        stop(paste0("The following Old_names are not present in Neighborhood_assignment of DATA: ", stringr::str_c(Conflictive_names, collapse = ", ")))
      }
      if(length(New_names) != length(Old_names)) stop("Provided New_names and Old_names should be of equal length")

      DATA_neighborhoods <- DATA

      for(Index in 1:length(New_names)){
        Replace_old <- Old_names[Index]
        Replacement <- New_names[Index]

        DATA_neighborhoods$Neighborhood_assignment[DATA_neighborhoods$Neighborhood_assignment == Replace_old] <- Replacement
      }
      return(DATA_neighborhoods)
    }
  }
