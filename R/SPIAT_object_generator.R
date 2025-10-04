#' Generates SPIAT objects for every image in a dataset
#'
#' The function transforms image information into a SPIAT object that can be used in subsequent analysis that depend on SPIAT package, like
#' [SPIAT_Heterogeneity_Analyzer()], [SPIAT_entropy_gradient_generator()], [SPIAT_Tissue_structuring_function()] and [SPIAT_neighborhood_identifier()].
#'
#' @param DATA_Intensities A dataframe or tibble containing cell feature data with feature expression by cell.
#' @param DATA_Phenotypes A dataframe or tibble containing cell feature data with a column named 'Phenotype' containing cell phenotype labels.
#'
#' @returns A list with SPIAT objects.
#'
#' @seealso [SPIAT_Heterogeneity_Analyzer()], [SPIAT_entropy_gradient_generator()], [SPIAT_Tissue_structuring_function()],[SPIAT_neighborhood_identifier()].
#' @export

SPIAT_object_generator <-
  function(DATA_Intensities = NULL,
           DATA_Phenotypes = NULL) {

    #Check suggested packages
    if(!requireNamespace("SPIAT", quietly = TRUE)) stop(
      paste0("SPIAT Bioconductor package is required to execute the function. Please install using the following code: ",
             expression({
               if (!require("BiocManager", quietly = TRUE))
                 install.packages("BiocManager")

               BiocManager::install("SPIAT")
             })
      )
    )

    #Check arguments
    if(!identical(names(DATA_Intensities)[c(1:4)], c("Cell_no", "X", "Y", "Subject_Names"))) {
      stop("DATA_Intensities should be adequately formated using Data_arrange_function from Step 0")
    }
    if(!identical(names(DATA_Phenotypes)[c(1:4)], c("Cell_no", "X", "Y", "Subject_Names"))) {
      stop("Column Names in Data_Phenotypes are not correctly stated (should contain Cell_no, X, Y and Subject_Names")
    }
    if(!"Phenotype" %in% names(DATA_Phenotypes)) stop("DATA_Phenotypes must contain a Phenotype column")


    #Select data
    DATA_Phenotypes <- DATA_Phenotypes
    DATA_Intensities <- DATA_Intensities[names(DATA_Phenotypes)[-ncol(DATA_Phenotypes)]]

    Formatted_image <-
      purrr::map(unique(DATA_Phenotypes$Subject_Names), function(Image) {

        #Select individual images
        RESULTS_INTERIM <- DATA_Phenotypes %>% dplyr::filter(Subject_Names == Image)
        Results_matrix <- DATA_Intensities %>% dplyr::filter(Subject_Names == Image)

        #Generate the intensity matrix
        Results_matrix <- t(as.matrix(Results_matrix[-c(1:4)]))
        colnames(Results_matrix) <- RESULTS_INTERIM[[1]]

        #Generate the SPIAT IMAGE
        Formatted_image <- SPIAT::format_image_to_spe(intensity_matrix = Results_matrix,
                                                      phenotypes = RESULTS_INTERIM[["Phenotype"]],
                                                      coord_x = RESULTS_INTERIM[[2]],
                                                      coord_y = RESULTS_INTERIM[[3]])
        Formatted_image

      }, .progress = list(clear = F,
                          name = "Generating SPIAT type objects",
                          show_after = 1,
                          type = "iterator"))

    names(Formatted_image) <- unique(DATA_Phenotypes$Subject_Names)
    return(Formatted_image)

  }
