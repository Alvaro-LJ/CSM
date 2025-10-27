#' Segments any cell feature according to user provided cut-off values
#'
#' Given a numeric cell feature and user specified cut-off values, a character vector is generated to summarize the feature.
#'
#' @param DATA A dataframe or tibble containing cell feature data.
#' @param DATA_variable A character indicating the feature name to be segmented.
#' @param DATA_cutoff A numeric vector specifying the cut-off points.
#' @param New_labels Character vector specifying the new labels (must be N of cut off - 1).
#' @param Merge A logical value indicating if the new variable should be merged with another variable.
#' @param Var_to_Merge (OPTIONAL) if Merge is TRUE, a character indicating the feature name to be combined with the segmented feature.
#' @returns Returns a tibble with the new segmented feature
#'
#' @examples
#' \dontrun{
#' Marker_segmentator(
#'     DATA = CSM_DistanceToPixelcell_test,
#'     DATA_variable = "CK_DIST",
#'     DATA_cutoff = c(0, 1, 5, 500),
#'     New_labels = c("Inside", "Border", "Outside"),
#'     Merge = TRUE,
#'     Var_to_Merge = "Phenotype"
#' )
#' }
#'
#' @export

Marker_segmentator <-
  function(DATA = NULL,
           DATA_variable = NULL,
           DATA_cutoff = NULL,
           New_labels = NULL,
           Merge = NULL,
           Var_to_Merge = NULL){

    #Check arguments
    if(!DATA_variable %in% names(DATA)) {
      stop(paste0(DATA_variable, " not present in the DATA provided."))
    }

    if(!is.numeric(DATA_cutoff)) stop("DATA_cutoff must be a numeric vector")
    if(!length(New_labels) == length(DATA_cutoff)-1) stop(paste0(length(DATA_cutoff)-1, " New_labels must be provided"))
    if(!is.logical(Merge)) stop("Merge must be a logical value")
    if(Merge){
      if(!Var_to_Merge %in% names(DATA)) stop(paste0(Var_to_Merge, " not present in DATA provided"))
    }

    #Generate the segmented version of the variable
    Segmented_distance <- as.character(cut(DATA[[DATA_variable]], breaks = DATA_cutoff, labels = New_labels))

    #Add it to data or to the desired cell column
    if(!Merge) DATA[[DATA_variable]] <- Segmented_distance
    if(Merge) DATA[[Var_to_Merge]] <- stringr::str_c(DATA[[Var_to_Merge]], Segmented_distance, sep = "_")

    return(DATA)
  }
