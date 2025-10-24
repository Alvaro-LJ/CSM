#' Calculates the average target cells within a distance from Cell of Origin for every image
#'
#' The function calculates the average number of target cells within a user defined distance for every image in a dataset.
#'
#' @param DATA A list containing cumulative interaction data. It can be computed using [Cumulative_Interaction_generator()] function.
#' @param DATA_RANDOM (OPTIONAL) A list containing random cumulative interaction data. It can be computed using [Cumulative_Interaction_generator()] function.
#' @param Radius A numeric value indicating the size of the radius to perform analysis. It must have been computed during the generation of the cumulative interaction matrix.
#' @param Include_Random A logical value indicating if the Random background should be used.
#' @param By_Sample_Random If a random background needs to be included, a logical value indicating if the random background should be calculated by sample or by experiment (see details).
#'
#' @seealso [Distance_matrix_generator()], [Random_Distance_matrix_generator()], [Cumulative_Interaction_generator()]
#'
#' @details
#' If By_Sample_Random is TRUE, for every sample, a specific random background will be calculated to compute expected spatial interaction metrics. If By_Sample_Random is FALSE, the random background will be calculated using all image information. The same random background will be used for all images.
#'
#'
#' @returns A tibble containing a summary by sample of spatial interactions.
#'
#' @examples
#' #Generate distance matrix and random distance matrix------------
#'DATA_Distances <-
#' Distance_matrix_generator(
#'     N_cores = 1,
#'     DATA = CSM_Phenotypecell_test,
#'     Cell_Of_Origin = "CD8_GZMBneg",
#'     Target_Cell = "TUMOR",
#'     Allow_Cero_Distance = FALSE,
#'     Perform_edge_correction = FALSE
#')
#'
#'RANDOM_Distances <-
#' Random_Distance_matrix_generator(
#'    N_cores = 1,
#'    DATA = CSM_Phenotypecell_test,
#'    Cell_Of_Origin = "CD8_GZMBneg",
#'    Target_Cell = "TUMOR",
#'    Random_cells_per_sample = 10,
#'    Allow_Cero_Distance = FALSE,
#'    Perform_edge_correction = FALSE
#')
#' #Generate cumulative interactions (must contain the actual radius distance)--
#'DATA_Cumulative <-
#'Cumulative_Interaction_generator(
#'    N_cores = 1,
#'    DATA = DATA_Distances,
#'    Start_from = 25,
#'    Stop_at = 100,
#'    Sampling_frequency = 25
#')
#'
#'RANDOM_Cumulative <-
#'Cumulative_Interaction_generator(
#'    N_cores = 1,
#'    DATA = RANDOM_Distances,
#'    Start_from = 25,
#'    Stop_at = 100,
#'    Sampling_frequency = 25
#')
#'#Calculate the cells in radius-----------------------------------------------
#'Cells_in_Radius_analyzer(
#'    DATA = DATA_Cumulative ,
#'    DATA_RANDOM = RANDOM_Cumulative,
#'    Radius = 50,
#'    Include_Random = T,
#'    By_Sample_Random = T
#')
#'
#'
#'
#' @export

Cells_in_Radius_analyzer <-
  function(DATA = NULL,
           DATA_RANDOM = NULL,
           Radius = NULL,
           Include_Random = NULL,
           By_Sample_Random = NULL) {
    #Check arguments
    if(!is.logical(Include_Random)) stop("Include_Random should be a logical value")
    if(!(as.character(Radius) %in% names(DATA[[1]][[2]]))){
      stop(paste0("Radius should be one of: ", stringr::str_c(names(DATA[[1]][[2]])[-1], collapse = ", ")))
    }

    #If everything is OK proceed with analysis
    if(Include_Random) {
      if(!is.logical(By_Sample_Random)) stop("By_Sample_Random should be a logical value")
      if(!length(intersect(names(DATA), names(DATA_RANDOM))) == length(unique(c(names(DATA), names(DATA_RANDOM))))){
        outersect <- function(x, y) {
          sort(c(x[!x%in%y],
                 y[!y%in%x]))
        }
        Removed_cases <- outersect(names(DATA), names(DATA_RANDOM))

        message(paste0("Only samples present in DATA and DATA_RANDOM will be used. The following samples will be removed: ",
                       stringr::str_c(Removed_cases, collapse = ", ")))
        DATA <- DATA[intersect(names(DATA), names(DATA_RANDOM))]
        DATA_RANDOM <- DATA_RANDOM[intersect(names(DATA), names(DATA_RANDOM))]
      }

      print("Calulating Results")
      #Generate Results
      RESULTS <-purrr::map_dbl(DATA, function(x) {
        Interim <- x[[2]]
        Interim <- Interim %>% dplyr::select(as.character(Radius))  #select the radius length (required to pick up choice from length list)
        mean(Interim[[1]])
      })

      RESULTS_SE <-purrr::map_dbl(DATA, function(x) {
        Interim <- x[[2]]
        Interim <- Interim %>% dplyr::select(as.character(Radius))  #select the radius length (required to pick up choice from length list)
        sd(Interim[[1]]) / sqrt(length(Interim[[1]]))
      })

      #Calculate sample-wise random distribution
      print("Calculating Random Background")
      if(By_Sample_Random){
        #Generate Random results
        RESULTS_RANDOM <-purrr::map_dbl(DATA_RANDOM, function(x) {
          Interim <- x[[2]]
          Interim <- Interim %>% dplyr::select(as.character(Radius))  #select the radius length (required to pick up choice from length list)
          mean(Interim[[1]])
        })

        RESULTS_RANDOM_SE <-purrr::map_dbl(DATA_RANDOM, function(x) {
          Interim <- x[[2]]
          Interim <- Interim %>% dplyr::select(as.character(Radius))
          sd(Interim[[1]]) / sqrt(length(Interim[[1]]))
        })
      }

      #Calculate experiment-wise random distribution
      if(!By_Sample_Random){
        RESULTS_RANDOM <- mean(unlist(
          purrr::map(DATA_RANDOM, function(x){
            Interim <- x[[2]]
            Interim <- Interim %>% dplyr::select(as.character(Radius))
          })
        ))

        RESULTS_RANDOM_SE <- sd(unlist(
          purrr::map(DATA_RANDOM, function(x){
            Interim <- x[[2]]
            Interim <- Interim %>% dplyr::select(as.character(Radius))
          })
        )) / sqrt(length(unlist(
          purrr::map(DATA_RANDOM, function(x){
            Interim <- x[[2]]
            Interim <- Interim %>% dplyr::select(as.character(Radius))
          })
        )))
      }

      #Generate the number of cells
      Cell_A <-purrr::map_dbl(DATA, function(x) x[[1]][[1,2]])
      Cell_B <-purrr::map_dbl(DATA, function(x) x[[1]][[2,2]])
      Name_A <- unique(purrr::map_chr(DATA, function(x) x[[1]][[1,1]]))
      Name_B <- unique(purrr::map_chr(DATA, function(x) x[[1]][[2,1]]))

      #Generate results tibble
      RESULTS_tibble <- tibble(Subject_Names = names(RESULTS), Average_cells = RESULTS)

      RESULTS_tibble$SE <- RESULTS_SE
      RESULTS_tibble <- RESULTS_tibble %>%dplyr::mutate(Result_05CI = Average_cells - 1.96*SE,
                                                        Result_95CI = Average_cells + 1.96*SE)
      RESULTS_tibble$Radius <- Radius

      RESULTS_tibble$Cell_A <- Cell_A
      names(RESULTS_tibble)[[7]] <- Name_A
      RESULTS_tibble$Cell_B <- Cell_B
      names(RESULTS_tibble)[[8]] <- Name_B


      RESULTS_tibble$Random <- RESULTS_RANDOM
      RESULTS_tibble$Random_SE <- RESULTS_RANDOM_SE
      RESULTS_tibble <- RESULTS_tibble %>%dplyr::mutate(Random_05CI = Random - 1.96*Random_SE,
                                                        Random_95CI = Random + 1.96*Random_SE) %>%
        dplyr::mutate(Outside_CI = case_when(Average_cells > Random_95CI | Average_cells < Random_05CI ~ "Significant",
                                             TRUE ~ "Not Significant"),
                      Sign = case_when(Average_cells > Random ~ "Above expected",
                                       TRUE ~ "Below expected"))
      RESULTS_tibble$N_Random_COO <-purrr::map_dbl(DATA_RANDOM, function(DATA_RANDOM) nrow(DATA_RANDOM[[2]]))
      #Establish a lower limit of 0 for CI
      RESULTS_tibble$Result_05CI[RESULTS_tibble$Result_05CI < 0] <- 0
      RESULTS_tibble$Random_05CI[RESULTS_tibble$Random_05CI < 0] <- 0

      print("Generating plot")
      plot(RESULTS_tibble %>% ggplot(aes(x = forcats::fct_reorder(Subject_Names, Average_cells))) +
             geom_col(aes(y = Average_cells), width = 0.5, color = "black", fill = "white", linewidth = 0.7)+
             geom_errorbar(aes(ymin = Result_05CI, ymax = Result_95CI), color = "black", linewidth = 0.9, width = 0.3)+
             geom_errorbar(aes(ymin = Random_05CI, ymax = Random_95CI), color = "red", linewidth = 0.9, width = 0.3)+
             cowplot::theme_cowplot() +
             scale_x_discrete("") +
             scale_y_continuous(stringr::str_c("Average cells within a radius of ", as.character(Radius))) +
             theme(axis.text.x = element_text(angle = -90, vjust = 0.5))
      )

      return(RESULTS_tibble)
    }

    else if(!Include_Random){
      print("Calulating Results")
      #Generate Results
      RESULTS <-purrr::map_dbl(DATA, function(x) {
        Interim <- x[[2]]
        Interim <- Interim %>% dplyr::select(as.character(Radius))  #select the radius length (required to pick up choice from length list)
        mean(Interim[[1]])
      })

      RESULTS_SE <-purrr::map_dbl(DATA, function(x) {
        Interim <- x[[2]]
        Interim <- Interim %>% dplyr::select(as.character(Radius))  #select the radius length (required to pick up choice from length list)
        sd(Interim[[1]]) / sqrt(length(Interim[[1]]))
      })

      #Generate the number of cells
      Cell_A <-purrr::map_dbl(DATA, function(x) x[[1]][[1,2]])
      Cell_B <-purrr::map_dbl(DATA, function(x) x[[1]][[2,2]])
      Name_A <- unique(purrr::map_chr(DATA, function(x) x[[1]][[1,1]]))
      Name_B <- unique(purrr::map_chr(DATA, function(x) x[[1]][[2,1]]))

      #Generate results tibble
      RESULTS_tibble <- tibble(Subject_Names = names(RESULTS), Average_cells = RESULTS)

      RESULTS_tibble$SE <- RESULTS_SE
      RESULTS_tibble <- RESULTS_tibble %>%dplyr::mutate(Result_05CI = Average_cells - 1.96*SE,
                                                        Result_95CI = Average_cells + 1.96*SE)
      RESULTS_tibble$Radius <- Radius

      RESULTS_tibble$Cell_A <- Cell_A
      names(RESULTS_tibble)[[7]] <- Name_A
      RESULTS_tibble$Cell_B <- Cell_B
      names(RESULTS_tibble)[[8]] <- Name_B
      #Establish a lower limit of 0 for CI
      RESULTS_tibble$Result_05CI[RESULTS_tibble$Result_05CI < 0] <- 0

      print("Generating plot")
      plot(RESULTS_tibble %>% ggplot(aes(x = forcats::fct_reorder(Subject_Names, Average_cells))) +
             geom_col(aes(y = Average_cells), width = 0.5, color = "black", fill = "white", linewidth = 0.7)+
             geom_errorbar(aes(ymin = Result_05CI, ymax = Result_95CI), color = "black", linewidth = 0.9, width = 0.3)+
             cowplot::theme_cowplot() +
             scale_x_discrete("") +
             scale_y_continuous(stringr::str_c("Average cells within a radius of ", as.character(Radius))) +
             theme(axis.text.x = element_text(angle = -90, vjust = 0.5))
      )

      return(RESULTS_tibble)
    }
  }
