#' Generates all thresholded feature combinations present in a thresholded cell feature dataset
#'
#' Given a thresholded cell feature dataset, the function will find all unique combinations of positive features present. These combinations are actually cell phenotypes.
#' The result can then be used to feed the [Phenotype_assigner_function()] to obtain cell phenotypes.
#'
#' @param DATA A dataframe or tibble containing thresholded cell feature data.
#' @param Markers A character vector indicating which features need to be used.
#'
#' @returns Returns a tibble where each row are unique positive feature combinations found in DATA.
#'
#' @seealso [Thresholding_function()], [Thresholding_function_tailored()] [Phenotype_assigner_function()]
#'
#' @examples
#' \dontrun{
#' #Threshold data-------------------------------------
#'DATA_thresholded <- Thresholding_function(
#'   DATA = CSM_Arrangedcellfeaturedata_test,
#'   Strategy = "EBI_Otsu",
#'   Local_thresholding = FALSE,
#'   Method_autothreshold = "Otsu",
#'   number_iterations_TriClass = 20,
#'   Percentile = 0.5,
#'   Defined_threshold = 0.1,
#'   Levels = 3
#' )
#' #Find unique feature positivity combinations----------
#'Phenotype_possibilities <- Marker_combinator_generator(
#'   DATA = DATA_thresholded,
#'   Markers = names(DATA_thresholded)[-c(1:4)]
#')
#'#Assign a phenotype to each combination------------
#'Phenotype_possibilities$Phenotype <- c('TUMOR', 'OTHER', 'CD8_GZMBneg', 'CD8_GZMBneg', 'OTHER', 'CD8_GZMBpos', 'CD8_GZMBpos')
#'
#'#Perform phenotyping-------------------------------
#'Phenotype_assigner_function(
#'   DATA = DATA_thresholded,
#'   Phenotype_possibilities = Phenotype_possibilities
#')
#' }

#'
#' @export

Marker_combinator_generator <-
  function(DATA = NULL, Markers = NULL) {

    #Check arguments
    if(!identical(names(DATA)[1:4], c("Cell_no", "X", "Y", "Subject_Names"))){
      stop("DATA not correctly specified, please format appropiatetly (see Step 0)")
    }
    if(!all(Markers %in% names(DATA))){
      Missing_value <- Markers[!Markers %in% names(DATA)]
      stop(paste0("The following marker is not present in the DATA: ", stringr::str_c(Missing_value, collapse = ", ")))
    }
    if(!all(purrr::map_lgl(DATA %>% dplyr::select(all_of(Markers)), function(Var) {
      is.logical(Var) | length(unique(Var)) <= 10
    }))) stop("Markers supplied should be logical (Generated after thresholding) or Multi-level with less than 10 Levels")

    print("Generating marker combination possibilities")
    Pheno_patterns <- DATA %>% dplyr::select(all_of(Markers)) %>% distinct()
    Pheno_patterns <-  Pheno_patterns %>% dplyr::mutate(Pattern = apply(Pheno_patterns, MARGIN = 1, function(x) stringr::str_c(x, collapse = "_")))

    print("Counting the different possibilities")
    Pattern_count <- DATA %>%dplyr::mutate(Pattern = apply((DATA %>% dplyr::select(all_of(Markers))), MARGIN = 1, function(x) stringr::str_c(x, collapse = "_"))) %>%
      group_by(Pattern) %>% dplyr::count()

    dplyr::left_join(Pheno_patterns, Pattern_count, by = "Pattern") %>% dplyr::arrange(desc(n)) %>% dplyr::rename(Number_of_cells = n) %>%
      dplyr::mutate(Phenotype_no =stringr::str_c("Cell_type_", as.character(1:nrow(Pheno_patterns))))
  }
