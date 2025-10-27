#' Analyzes the spatial interaction using Ripleys K and L functions
#'
#' The function calculates the spatial aggregation of a single cell type using Ripleys K and L functions. Function results are integrated to obtain an AUC,
#' that can serve as a metric of spatial interaction. The function is calculated using the spatstat R metapackage.
#' (specifically the spatstat.explore::envelope). If integration fails, the user may try different Max_distance thresholds to make it work.
#'
#' @param DATA A dataframe or tibble containing a column named 'Phenotype' containing cell phenotype labels.
#' @param Cell_type A character value indicating the cell phenotype to be analyzed.
#' @param Max_distance A numeric value indicating the maximum distance to be analyzed.
#' @param Strategy A character value indicating the function type: either "Ripleys_K" or "Ripleys_L"
#' @param N_simulations A integer value indicating the number of simulations to calculate the random background.
#'
#' @returns A tibble containing results by image.
#'
#' @examples
#' \dontrun{
#'Ripley_function_calculator(
#'    DATA = CSM_Phenotypecell_test,
#'    Cell_type = "TUMOR",
#'    Max_distance = 15,
#'    Strategy = "Ripleys_K",
#'    N_simulations = 10
#')
#' }
#'
#' @export

Ripley_function_calculator <-
  function(DATA = NULL,
           Cell_type = NULL,
           Max_distance = NULL,
           Strategy = NULL,
           N_simulations = NULL){

    #Check suggested packages
    {
      if(!requireNamespace("spatstat", quietly = FALSE)) stop(
        paste0("spatstat CRAN package is required to execute the function. Please install using the following code: ",
               expression(install.packages("spatstat")))
      )
    }

    #Argument validation (Here all the functino has been optimized by ChatpGPT)
    if(is.null(Cell_type) || !is.character(Cell_type)) {
      stop("Cell_type must be a non-null character value.")
    }
    if(!Cell_type %in% unique((DATA %>% dplyr::select(Phenotype))[[1]])) stop(paste0("Cell_type must be one of the following: ",
                                                                                     stringr::str_c(unique((DATA %>% dplyr::select(Phenotype))[[1]]), collapse = ", ")
    )
    )
    if(is.null(Max_distance) || !is.numeric(Max_distance) || Max_distance <= 0) {
      stop("Max_distance must be a positive numeric value.")
    }
    if(is.null(Strategy) || !Strategy %in% c("Ripleys_K", "Ripleys_L")) {
      stop("Strategy must be one of 'Ripleys_K' or 'Ripleys_L'.")
    }
    if(is.null(N_simulations) || !is.numeric(N_simulations) || N_simulations <= 0) {
      stop("N_simulations must be a positive numeric value.")
    }

    # Filter the dataset based on cell type
    DATA_Phenotypes <- DATA %>% dplyr::filter(Phenotype == Cell_type)

    # Map across unique subject names
    purrr::map_dfr(unique(DATA_Phenotypes$Subject_Names), function(Image) {
      Interim <- DATA_Phenotypes %>% dplyr::filter(Subject_Names == Image)  # Select cell type

      # Prepare point pattern
      point_pattern <- spatstat.geom::ppp(
        x = Interim$X, y = Interim$Y,
        n = nrow(Interim),
        window = spatstat.geom::owin(xrange = c(min(Interim$X), max(Interim$X)),
                                     yrange = c(min(Interim$Y), max(Interim$Y)))
      )

      # Define function parameters and run envelope calculations
      fun_name <- if(Strategy == "Ripleys_K") 'Kest' else 'Lest'

      Simulation <- spatstat.explore::envelope(
        point_pattern,
        fun = fun_name,
        funargs = list(rmax = Max_distance, correction = c("Ripley"),
                       nlarge = 3000, var.approx = T, ratio = FALSE),
        nsim = N_simulations,
        verbose = FALSE
      )

      Results <- as_tibble(Simulation)

      if(any(is.na(Results))) {
        warning(paste(Strategy, "function could not be calculated. Try switching Max_distance argument"))
        RESULT_tibble <- tibble(
          Subject_Names = Image,
          Observed_AUC = NA,
          CSR_AUC = NA,
          Observed_limit_AUC = NA,
          Ratio = NA,
          Significant = NA
        )
      } else {
        # Calculate AUC for observed and CSR values
        Observed_AUC <- integrate(approxfun(Results$r, Results$obs), lower = 0, upper = max(Results$r), subdivisions = 2000)$value
        CSR_AUC <- integrate(approxfun(Results$r, Results$theo), lower = 0, upper = max(Results$r), subdivisions = 2000)$value
        Ratio <- Observed_AUC / CSR_AUC

        # Determine significance
        if(all(Ratio > 1 & !any(is.nan(Results$lo)))) {
          Significant <- CSR_AUC < integrate(approxfun(Results$r, Results$lo), lower = 0, upper = max(Results$r), subdivisions = 2000)$value
          Observed_limit_AUC <- integrate(approxfun(Results$r, Results$lo), lower = 0, upper = max(Results$r), subdivisions = 2000)$value
        } else if(all(Ratio < 1 & !any(is.nan(Results$hi)))) {
          Significant <- Observed_AUC > integrate(approxfun(Results$r, Results$hi), lower = 0, upper = max(Results$r), subdivisions = 2000)$value
          Observed_limit_AUC <- integrate(approxfun(Results$r, Results$hi), lower = 0, upper = max(Results$r), subdivisions = 2000)$value
        } else {
          Significant <- NA
          Observed_limit_AUC <- NA
        }

        RESULT_tibble <- tibble(
          Subject_Names = Image,
          Observed_AUC = Observed_AUC,
          CSR_AUC = CSR_AUC,
          Observed_limit_AUC = Observed_limit_AUC,
          Ratio = Ratio,
          Significant = Significant
        )
      }

      return(RESULT_tibble)

    }, .progress = list(clear = FALSE,
                        name = "Calculating spatial aggregation function",
                        show_after = 1,
                        type = "iterator"))
  }
