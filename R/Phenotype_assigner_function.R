#' Assigns a phenotype to every cell according to user specified labels and thresholding results.
#'
#' The function takes a thresholded cell feature matrix and a tibble obtained using [Marker_combinator_generator()]. The resulting tibble obtained using [Marker_combinator_generator()]
#' must have a new column named 'Phenotype' where the user specifies the name of the phenotype according to the feature positivity pattern. The function uses this information to assign cell phenotypes to all cells in DATA.
#'
#' @param DATA A dataframe or tibble containing thresholded cell feature data.
#' @param Phenotype_possibilities A tibble obtained using [Marker_combinator_generator()] containing a column named 'Phenotype'.
#'
#' @returns Returns a cell feature tibble containing a new column named 'Phenotype' where cell labels are present.
#'
#' @seealso [Thresholding_function()], [Thresholding_function_tailored()], [Marker_combinator_generator()]
#'
#' @examples
#' \dontrun{
#' #Threshold data
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
#' #Find unique feature positivity combinations
#'Phenotype_possibilities <- Marker_combinator_generator(
#'   DATA = DATA_thresholded,
#'   Markers = names(DATA_thresholded)[-c(1:4)]
#')
#'#Assign a phenotype to each combination
#'Phenotype_possibilities$Phenotype <- c('TUMOR', 'OTHER', 'CD8_GZMBneg', 'CD8_GZMBneg', 'OTHER', 'CD8_GZMBpos', 'CD8_GZMBpos')
#'
#'#Perform phenotyping
#'Phenotype_assigner_function(
#'   DATA = DATA_thresholded,
#'   Phenotype_possibilities = Phenotype_possibilities
#')
#' }
#'
#' @export

Phenotype_assigner_function <-
  function(DATA = NULL, Phenotype_possibilities = NULL) {
    if(sum(names(Phenotype_possibilities)[1:(ncol(Phenotype_possibilities)-4)] %in% names(DATA)) != length(names(Phenotype_possibilities)[1:(ncol(Phenotype_possibilities)-4)])) {
      stop("Phenotype markers not matched in thresholded data")
    }

    else if(!identical(names(Phenotype_possibilities)[(ncol(Phenotype_possibilities)-3):ncol(Phenotype_possibilities)],
                       c("Pattern",  "Number_of_cells", "Phenotype_no", "Phenotype"))) {
      stop("Phenotype_possibilities not correctly specified. Please use marker_combinator_generator and specify phenotypes in a column named: 'Phenotype'")
    }

    else {
      print("Performing phenotype assignment")
      Interim <- DATA %>%dplyr::mutate(Pattern = apply((DATA %>% dplyr::select(all_of(names(Phenotype_possibilities)[1:(ncol(Phenotype_possibilities)-4)]))),
                                                       MARGIN = 1,
                                                       function(x) stringr::str_c(x, collapse = "_")))
      Phenotypes <- Phenotype_possibilities %>% dplyr::select(Pattern, Phenotype)
      dplyr::left_join(Interim, Phenotypes, by = "Pattern") %>% dplyr::select(-Pattern)
    }
  }
