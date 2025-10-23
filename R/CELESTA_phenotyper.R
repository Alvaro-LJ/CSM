#' Assigns phenotype labels using CELESTA algorithm
#'
#' The function runs the CELESTA algorithm using cell feature data and a template indicating expected cell phenotypes.
#'
#' @param DATA A dataframe or tibble containing cell feature data.
#' @param Template_path A character indicating the path to the CSV file containing the fulfilled template.
#' @param Alternative_CSV_locale A logical value indicating if the CSV has an alternative locale format (like a European vs USA formatting).
#' @param N_cores Integer. Number of cores to parallelize your computation.
#'
#' @param Apply_filters A logical value indicating if cells should be pre-filtered before running the algorithm.
#' @param high_marker_threshold A numeric value indicating the max quantile allowed. Cells with expression above this quantile for all markers will not be phenotyped.
#' @param low_marker_threshold A numeric value indicating the min quantile allowed. Cells with expression below this quantile for all markers will not be phenotyped.
#'
#' @param max_iteration An integer. Number of iterations to run the GMM algorithm used in CELESTA
#' @param cell_change_threshold A numeric value that indicates when the iterative cell-type assignment should stop. The default value is 0.01,
#' which means that if the percentage of additional assigned cells is smaller than 1% of the unassigned cells, then cell-type assignment will stop.
#' The recommended range is 0.01 - 0.05. Note that the higher the cell change threshold, the more cells are left unassigned.
#'
#' @returns Returns a tibble with cell features and a column named 'Phenotype' containing cell labels.
#'
#' @seealso [CELESTA_template_generator()]. For more information on how to fulfill the template please visit https://github.com/plevritis-lab/CELESTA.
#'
#'
#' @examples
#' \dontrun{
#' CELESTA_phenotyper(
#'   DATA = CSM_Arrangedcellfeaturedata_test,
#'   Template_path = "Completed_template_path",
#'   Alternative_CSV_locale = TRUE,
#'   N_cores = 2,
#'   Apply_filters = TRUE,
#'   high_marker_threshold = 0.90,
#'   low_marker_threshold = 0.10,
#'   max_iteration = 10,
#'   cell_change_threshold = 0.01
#')
#' }
#'
#'
#'
#' @export

CELESTA_phenotyper <-
  function(DATA = NULL,
           Template_path = NULL,
           Alternative_CSV_locale = FALSE,
           N_cores = 1,

           Apply_filters = FALSE,
           high_marker_threshold = NULL,
           low_marker_threshold = NULL,

           max_iteration = 10,
           cell_change_threshold = 0.01
  ){

    #Check required packages
    {
      if(!requireNamespace("readr", quietly = FALSE)) stop(
        paste0("readr CRAN package is required to execute the function. Please install using the following code: ",
               expression(install.packages("readr")))
      )
      if(!requireNamespace("CELESTA", quietly = FALSE)) stop(
        paste0("CELESTA GitHub package is required to execute the function. Please install using the following code: ",
               expression(devtools::install_github("plevritis/CELESTA")))
      )
    }


    DATA <- DATA
    #Check arguments
    if(!identical(names(DATA)[1:4], c("Cell_no", "X", "Y", "Subject_Names"))) {
      stop("DATA not correctly specified, please format appropiatetly (see Step 0)")
    }

    if(!is.logical(Alternative_CSV_locale)) stop("Alternative_CSV_locale must be a logical value")
    if(!is.logical(Apply_filters)) stop("Apply_filters must be a logical value")
    if(!all(N_cores >= 1, N_cores%%1 == 0)) stop("N_cores must be a positive integer value")
    if(Apply_filters){
      if(!all(is.numeric(high_marker_threshold), high_marker_threshold >= 0, high_marker_threshold <= 1)) stop("high_marker_threshold should be a numeric value between 0 and 1")
      if(!all(is.numeric(low_marker_threshold), low_marker_threshold >= 0, low_marker_threshold <= 1)) stop("low_marker_threshold should be a numeric value between 0 and 1")
    }
    if(!all(is.numeric(max_iteration), max_iteration > 0, max_iteration%%1 == 0)) stop("max_iteration should be a integer value > 0")
    if(!all(is.numeric(cell_change_threshold), cell_change_threshold >= 0, cell_change_threshold <= 1)) stop("cell_change_threshold should be a numeric value between 0 and 1")

    #Import Prior_marker_info
    if(!Alternative_CSV_locale) Prior_info <- readr::read_csv(Template_path, na = c("NA", ""))
    if(Alternative_CSV_locale) Prior_info <- readr::read_csv2(Template_path, na = c("NA", ""))

    #Check the imported Prior_info
    if(!identical(c("Phenotype",
                    "Phenotype_Number",
                    "Assignment_round",
                    "Phenotype_dependency",
                    "high_expression_threshold_anchor",
                    "high_expression_threshold_index",
                    "low_expression_threshold_anchor",
                    "low_expression_threshold_index"), names(Prior_info)[1:8])) stop("Template not well specified. Please review the file, the pathway provided and the CSV locale")

    #Check that there are no bugs in the Prior_info
    if(sum(is.na(Prior_info$Phenotype)) >= 1){
      stop("The following rows contain unintended NA values present in Template, please review the file provided: ",
           stringr::str_c(which(is.na(Prior_info$Phenotype)), collapse = ", "))
    }
    #Check if any variable has all NA in round 1 phenotypes
    if(any(purrr::map_lgl(Prior_info %>% dplyr::filter(Assignment_round == 1), function(var){
      sum(is.na(var)) == nrow(Prior_info %>% dplyr::filter(Assignment_round == 1))
    }))){
      Conflictive_vars <- names(Prior_info)[map_lgl(Prior_info %>% dplyr::filter(Assignment_round == 1), function(var){
        sum(is.na(var)) == nrow(Prior_info %>% dplyr::filter(Assignment_round == 1))
      })]
      stop(paste0("The following variables have NA values for all phenotypes identified during round 1: ",
                  stringr::str_c(Conflictive_vars, collapse = ", "),
                  ". This can cause errors during CELESTA execution. Please review the file provided."))
    }


    #Build the final Prior_info dataframe
    #Collapse the lineage level, bind it to marker info and add row names
    Lineage_level <- stringr::str_c(Prior_info$Assignment_round, Prior_info$Phenotype_dependency, Prior_info$Phenotype_Number, sep = "_")
    Final_template <- cbind(Prior_info[["Phenotype"]], Lineage_level, Prior_info[-c(1:8)])
    names(Final_template)[1] <- ""

    print(head(Final_template))
    print(paste0("The following markers will be used in the cell phenotyping process", stringr::str_c(names(Final_template)[-1], collapse = ", ")))
    print(paste0(nrow(DATA), " cells will be phenotyped in the process"))
    print(paste0(nrow(Final_template), " phenotypes will be identified"))

    Proceed_CELESTA <- menu(choices = c("Proceed", "Abort"), title = "Should the CELESTA algorithm be executed?")
    if(Proceed_CELESTA == 2) stop("Phenotyping has been aborted")

    #save exit function if parallelization fails
    on.exit({
      future::plan("future::sequential")
      gc()
    })
    #We make the clusters and load required packages
    future::plan("future::multisession", workers = N_cores)
    options(future.globals.maxSize = Inf, future.rng.onMisuse = "ignore")
    furrr::furrr_options(scheduling = Inf)

    Phenotyping_result <- suppressMessages({
      furrr::future_map(unique(DATA$Subject_Names), function(Index){
        library(zeallot)
        library(Rmixmod)
        library(spdep)

        #Select individual images
        DATA <- DATA %>% dplyr::filter(Subject_Names == Index)


        Celesta_object <- CELESTA::CreateCelestaObject(project_title = "Celesta_object",
                                                       prior_marker_info = Final_template,
                                                       imaging_data_file = DATA)
        if(Apply_filters){
          Celesta_object <- CELESTA::FilterCells(Celesta_object,
                                                 high_marker_threshold = high_marker_threshold,
                                                 low_marker_threshold = low_marker_threshold)
        }
        Celesta_object <- CELESTA::AssignCells(Celesta_object,
                                               max_iteration = max_iteration,
                                               cell_change_threshold = cell_change_threshold,
                                               high_expression_threshold_anchor = Prior_info[["high_expression_threshold_anchor"]],
                                               low_expression_threshold_anchor = Prior_info[["low_expression_threshold_anchor"]],
                                               high_expression_threshold_index = Prior_info[["high_expression_threshold_index"]],
                                               low_expression_threshold_index = Prior_info[["low_expression_threshold_index"]],
                                               save_result = FALSE)

        #Obtain final result, link it to original data and exit
        Phenotype_tibble <- as_tibble(Celesta_object@final_cell_type_assignment)$'Final cell type'
        Phenotype_tibble <- tibble(Phenotype = Phenotype_tibble)
        FINAL <-dplyr::bind_cols(DATA, Phenotype_tibble)
        return(FINAL)
      }, .progress = TRUE)
    })

    future::plan("future::sequential")
    gc()

    #return the complete dataset
    return(purrr::map_dfr(Phenotyping_result,dplyr::bind_rows))
  }
