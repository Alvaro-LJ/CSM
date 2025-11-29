#' Analyzes the spatial interaction between two cell types using Multi-level modelling
#'
#' The function calculates a multilevel model to analyze the spatial interaction between two cell types. The spatial cumulative interaction between two cell types,
#' one being the Cell of Origin (COO) and the other being the target cell, is dependent on two main variables: the target cell density and distance from the cell of origin.
#' If the tissue target cell density is higher, or the distance from the COO is larger, the spatial interaction is going to increase.
#' The function tries to model the influence of image metadata on this spatial interaction. It does so by considering cells in the dataset as nested observations,
#' where cells within the same image share specific characteristics (like exposure to the same target cell density or image metadata condition).
#' The model formula will look something like this: Number of Target cells ~ Image_Metadata*Distance_COO + Density_Target + (1|Cell_Of_Origin_ID)
#' The multilevel model is calculated using lme4::lmer function. The contribution of each variables to the total variance explanation is calculated using the partR2::partR2 function.
#'
#' @param DATA_cumulative A list containing cumulative interaction data. It can be computed using [Cumulative_Interaction_generator()] function.
#' @param DATA_Clinical A dataframe or tibble containing image metadata.
#' @param Clinical_var A character value indicating the column name used in the clinical data analysis.
#' @param DATA_Densities A dataframe or tibble containind the cell phenotype density by image. It can be computed using [Phenotype_quantifier()].
#' @param Cell_Of_Origin A character value indicating the cell phenotype label of the Cell of Origin.
#' @param Target_cell A character value indicating the cell phenotype label of the Target cell.
#' @param Calculate_R2 A logical value indicating if R2 should be calculated after model fitting.
#' @param N_bootstrap A integer value indicating the number of permutations to calculate R2 values.
#'
#' @returns Summary of the fitted multi-level model.
#'
#' @examples
#'
#' \dontrun{
#' #Obtain the Cumulative interaction between two cells types accross the desired distance--------
#'DATA_Distances <-
#' Distance_matrix_generator(
#'     N_cores = 1,
#'     DATA = CSM_PhenotypeTMA_test,
#'     Cell_Of_Origin = "CD8",
#'     Target_Cell = "MACROPHAGE",
#'     Allow_Cero_Distance = FALSE,
#'     Perform_edge_correction = FALSE
#')
#'DATA_Cumulative <-
#' Cumulative_Interaction_generator(
#'     N_cores = 1,
#'     DATA = DATA_Distances,
#'     Start_from = 10,
#'     Stop_at = 210,
#'     Sampling_frequency = 50
#')
#'
#'#Obtain cell density by sample to account for density of target cells in the model-------------
#'DATA_AREA <-
#'  Image_size_calculator(
#'     DATA = CSM_PhenotypeTMA_test,
#'     Strategy = "Concave_hull",
#'     Hull_ratio = 0.4
#')
#'Cells_by_sample <-
#'  Phenotype_quantifier(
#'     CSM_PhenotypeTMA_test,
#'     Calculate_Density = TRUE,
#'     DATA_Area = DATA_AREA
#' )
#'
#' #Arrange clinical data------------------------------------------------------------------------
#'DATA_CLINICAL <-
#' Clinical_Data_arrange_function(
#'      DATA = CSM_ClinicalTMA_test,
#'      Subject_Names = "Sample",
#'      Outcomes_to_keep = c("AGE", "MMRP_status", "DEATH", "OS_m")
#' )
#'
#' #Calculate the model--------------------------------------------------------------------------
#' Multi_level_modelling_function(
#'     DATA_cumulative = DATA_Cumulative,
#'     DATA_Clinical = DATA_CLINICAL,
#'     Clinical_var = "MMRP_status",
#'     DATA_Densities = Cells_by_sample,
#'     Cell_Of_Origin = "CD8",
#'     Target_cell = "MACROPHAGE",
#'     Calculate_R2 = FALSE
#' )
#' }
#'
#' @export

Multi_level_modelling_function <-
  function(DATA_cumulative,
           DATA_Clinical,
           Clinical_var,
           DATA_Densities,
           Cell_Of_Origin,
           Target_cell,
           Calculate_R2 = FALSE,
           N_bootstrap = NULL) {

    #Check suggested packages
    {
      if(!requireNamespace("lme4", quietly = FALSE)) stop(
        paste0("lme4 CRAN package is required to execute the function. Please install using the following code: ",
               expression(install.packages("lme4")))
      )
      if(!requireNamespace("partR2", quietly = FALSE)) stop(
        paste0("partR2 CRAN package is required to execute the function. Please install using the following code: ",
               expression(install.packages("partR2")))
      )

    }

    if(length(Clinical_var) > 1) {
      stop("Select a single clinical variable")
    }
    if(!(Clinical_var %in% names(DATA_Clinical))) {
      stop(paste0("Clinical var should be one of: ", stringr::str_c(names(DATA_Clinical)[-1], collapse = ", ")))
    }
    if(is.double(DATA_Clinical[[Clinical_var]])) {
      stop("Clinical variable should be a character vector, not a double vector")
    }
    if(!any(stringr::str_detect(names(DATA_Densities), Target_cell))){
      stop("Target cell density should be present in DATA_Densities")
    }
    if(!is.logical(Calculate_R2)) stop("Calculate_R2 must be a logical value")
    if(Calculate_R2){
      if(!all(N_bootstrap >= 1 & N_bootstrap%%1 == 0)) stop("N_bootstrap must be an integer value > 0")
    }


    #Import the data frames to the function
    DATA_cumulative <- DATA_cumulative
    DATA_Clinical <- DATA_Clinical %>% dplyr::select(Subject_Names, all_of(Clinical_var))
    DATA_Densities <- DATA_Densities %>% dplyr::select(Subject_Names, contains("Density_")) %>% dplyr::select(Subject_Names,
                                                                                                              contains(stringr::str_c("Density_", Target_cell, sep = "")))

    #First we need to create a dataframe with in an adequate format
    print("Formatting data")
    Interim <-purrr::map_dfr(1:length(DATA_cumulative), function(Image) {
      DATA_cumulative[[Image]][[2]] %>% tidyr::pivot_longer(-1) %>%
        dplyr::mutate(name = as.double(name),
                      Cell_Of_Origin_ID = stringr::str_c(names(DATA_cumulative)[Image], Cell_Of_Origin_no, sep = "_"),
                      Subject_Names = names(DATA_cumulative)[Image]) %>%
        dplyr::select(Subject_Names, Cell_Of_Origin_ID, name, value)
    })


    #Add the target density info
    Interim <-dplyr::left_join(Interim, DATA_Densities, by = "Subject_Names")
    #Add the clinical information
    Interim <-dplyr::left_join(Interim, DATA_Clinical, by = "Subject_Names")
    #Modify col names
    names(Interim)[-c(1:2)] <- c("DIST", "N_Target", "Density_Target", "Clin_Group")

    Interim <- Interim %>%dplyr::mutate(DIST = as.double(scale(DIST)),
                                        Density_Target = as.double(scale(Density_Target)),
                                        N_Target = as.double(scale(N_Target)))

    #Fit the model
    print("Fitting the multi-level model")
    Model <- lme4::lmer(N_Target ~ Clin_Group*DIST + Density_Target + (1|Cell_Of_Origin_ID), data = Interim)
    print(summary(Model))

    #Make predictions
    tryCatch({
      Prediction <- tidyr::expand_grid(unique(Interim$DIST),
                                       unique(Interim$Clin_Group))
      names(Prediction) <- c("DIST", "Clin_Group")
      Prediction$Density_Target <- 0
      Prediction$Cell_Of_Origin_ID <- unique(Interim$Cell_Of_Origin_ID)[1]
      Prediction$Prediction <- predict(Model, newdata = Prediction)
      #plot results
      plot(Prediction %>%
             ggplot(aes(x = DIST, y  = Prediction, group = Clin_Group, color = Clin_Group)) +
             geom_line(linewidth = 1.5) +
             cowplot::theme_cowplot() +
             scale_x_continuous(stringr::str_c("Scaled Distance from ", Cell_Of_Origin, sep = "")) +
             scale_y_continuous(stringr::str_c("Scaled number of ", Target_cell, " cells", sep = "")) +
             scale_color_viridis_d(Clinical_var))
    },
    error = function(cond) print("Unable to calculate Multi-level model predictions. Consider subsetting samples or modifying cumulative interaction parameters.")
    )


    #Calculate R2 values if required by user
    if(Calculate_R2){
      print("Calculating model R2 - This can take some time...")
      try(
        Model_R2 <-
          partR2::partR2(Model, partvars = c("DIST", "Clin_Group:DIST", "Density_Target"), R2_type = "conditional", nboot = N_bootstrap, data = Interim)
      )
      if(berryFunctions::is.error(Model_R2)) warning("Unable to calculate predictos R^2 for the current model. This may occur with very large datasets")
      else print(summary(Model_R2))
    }


  }
