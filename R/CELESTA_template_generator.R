#' Generates a CELESTA cell phenotype template
#'
#' In order to run CELESTA, a template indicating expected cell phenotypes present in the tissue must be prepared. This function generates the template that can be modified manually by the user.
#'
#' @param DATA A dataframe or tibble containing cell feature data.
#' @param Markers_to_keep A character vector indicating the names of the columns to be kept in the template.
#' @param Template_name The name of the file that will be written (optionally including the path where the file will be writen).
#'
#' @returns Creates a CSV file containing the template in the working directory (or the path specified).
#'
#' @seealso [CELESTA_phenotyper()]. For more information on how to fullfill the template please visit https://github.com/plevritis-lab/CELESTA.
#'
#' @examples
#' \dontrun{
#' CELESTA_template_generator(
#'    DATA = CSM_Arrangedcellfeaturedata_test,
#'    Markers_to_keep = c("CK_EPCAM_AVERAGE", "CD8a_AVERAGE", "GZMB_AVERAGE"),
#'    Template_name = "Mypath/name_of_the_template"
#'    )
#' }
#'
#'
#' @export

CELESTA_template_generator <-
  function(DATA = NULL,
           Markers_to_keep = NULL,
           Template_name = NULL){
    #Check required packages
    if(!requireNamespace("readr", quietly = FALSE)) stop(
      paste0("readr CRAN package is required to execute the function. Please install using the following code: ",
             expression(install.packages("readr")))
    )


    DATA <- DATA
    #Check that data and Markers_to_keep have been correctly specified
    if(!all(c("X", "Y", "Subject_Names", Markers_to_keep) %in% names(DATA))) {
      Missing_arguments <- c("X", "Y", "Subject_Names", Markers_to_keep)[!c("X", "Y", "Subject_Names", Markers_to_keep) %in% names(DATA)]
      stop(paste0(stringr::str_c(Missing_arguments, collapse = ", "), " not found in DATA"))
    }
    if(!is.character(as.character(Template_name))) stop("Template_name must be a character value")

    #Generate the basic tibble
    Basic <- tibble(Phenotype = NA, Phenotype_Number = NA, Assignment_round = NA, Phenotype_dependency = NA,
                    high_expression_threshold_anchor = NA, high_expression_threshold_index = NA,
                    low_expression_threshold_anchor = 1, low_expression_threshold_index = 1)

    #Generate the marker tibble
    Markers <- as_tibble(matrix(NA, nrow = 1, ncol = length(Markers_to_keep)), .name_repair = "minimal")
    Markers <- setNames(Markers, Markers_to_keep)

    Final_tibble <- dplyr::bind_cols(Basic, Markers)
    Final_tibble

    print(paste0("The CELESTA prior information template will be saved at: ", getwd()))

    readr::write_excel_csv(Final_tibble, paste0(Template_name, ".csv"), na = "NA")

    print("For more information on how to fill in the template please visit: https://github.com/plevritis-lab/CELESTA")
    print("Complete the template and save it as a csv file")
  }
