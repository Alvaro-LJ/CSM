#' Modifies cell phenotype labels
#'
#' The function allows the user to modify cell phenotype labels. The function works with any dataframe or tibble containing a column named 'Phenotype'.
#'
#' @param DATA A dataframe or tibble containing cell feature data, including a column named 'Phenotype'.
#' @param New_names A character vector indicating the new names. The length must be equal to the number of unique cell phenotype labels.
#' @param Old_names (OPTIONAL) A character vector indicating the phenotype labels that need to be modified. The length of this vector must be equal to the length of New_names.
#'
#' @returns Returns a tibble with cell features and the phenotype column with the modified labels.
#'
#' @export

DATA_Phenotype_renamer <-
  function(DATA = NULL,
           New_names = NULL,
           Old_names = NULL) {

    #If all names need to be replaced then apply as usual
    if(is.null(Old_names)){
      #Check if provided names are equal to number of hoods
      if(length(New_names) != length(unique(DATA$Phenotype))) {
        stop(paste0("Provided New_names should match the number of Phenotypes in the analysis. Number of Phenotypes: ", length(unique(DATA$Phenotype)),
                    ". Names provided: ", length(New_names)))
      }

      #Proceed with renaming
      DATA_Phenotypes <- DATA
      names_tibble <- tibble(Phenotype = factor(1:length(unique(DATA_Phenotypes$Phenotype))),
                             New_names = New_names)
      Result <-dplyr::left_join(DATA_Phenotypes, names_tibble, by = "Phenotype")
      Result %>%dplyr::mutate(Phenotype = New_names) %>% dplyr::select(-New_names)
    }

    #If new names and old names are provided then check first that the lengths are equal and check that old names are present in data. Then proceed
    if(!is.null(Old_names)){
      if(!all(Old_names %in% unique(DATA$Phenotype))){
        Conflictive_names <- Old_names[!Old_names %in% unique(DATA$Phenotype)]
        stop(paste0("The following Old_names are not present in Phenotypes of DATA: ", stringr::str_c(Conflictive_names, collapse = ", ")))
      }
      if(length(New_names) != length(Old_names)) stop("Provided New_names and Old_names should be of equal length")

      DATA_Phenotypes <- DATA

      for(Index in 1:length(New_names)){
        Replace_old <- Old_names[Index]
        Replacement <- New_names[Index]

        DATA_Phenotypes$Phenotype[DATA_Phenotypes$Phenotype == Replace_old] <- Replacement
      }
      return(DATA_Phenotypes)
    }
  }
