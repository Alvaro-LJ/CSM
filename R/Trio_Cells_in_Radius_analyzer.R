#' Calculates the average TRIO score for every image
#'
#' The function calculates the average TRIO score for every image given a radius size. For more information on the TRIO score please see [Trio_Cumulative_Interaction_generator()].
#'
#' @param DATA A list containing cumulative interaction data and TRIO scores. It can be computed using [Trio_Cumulative_Interaction_generator()] function.
#' @param DATA_RANDOM (OPTIONAL) A list containing random cumulative interaction data and TRIO scores. It can be computed using [Trio_Cumulative_Interaction_generator()] function.
#' @param Radius A numeric value indicating the size of the radius to perform analysis. It must have been computed during the generation of the cumulative interaction matrix.
#' @param Include_Random A logical value indicating if the Random background should be used.
#' @param By_Sample_Random If a random background needs to be included, a logical value indicating if the random background should be calculated by sample or by experiment (see details).
#'
#' @seealso [Trio_Random_Distance_matrix_generator()], [Trio_Cumulative_Interaction_generator()], [Trio_Min_Distance_analyzer()], [Trio_graph_maker()].
#'
#' @details
#' If By_Sample_Random is TRUE, for every sample, a specific random background will be calculated to compute expected spatial interaction metrics. If By_Sample_Random is FALSE, the random background will be calculated using all image information. The same random background will be used for all images.
#'
#' @returns A tibble containing a summary by sample of spatial interactions.
#'
#' @export

Trio_Cells_in_Radius_analyzer <-
  function(DATA = NULL,
           DATA_RANDOM = NULL,
           Radius = NULL,
           Include_Random = NULL,
           By_Sample_Random = NULL
  ) {

    #Check arguments
    if(!(as.character(Radius) %in% names(DATA[[1]][[4]]))){
      stop(paste0("Radius should be one of: ", stringr::str_c(names(DATA[[1]][[4]])[-1], collapse = ", ")))
    }
    if(!is.logical(Include_Random)) stop("Include_Random should be a logical value")

    print("Proceeding with analysis")
    #If everything is OK proceed with analysis
    if(Include_Random) {
      #Check arguments
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

      print("Proceeding with analysis")
      #Generate Results
      RESULTS <-purrr::map_dbl(DATA, function(x) {
        Interim <- x[[4]]
        Interim <- Interim %>% dplyr::select(as.character(Radius))  #select the radius length (required to pick up choice from length list)
        mean(Interim[[1]])
      })

      RESULTS_SE <-purrr::map_dbl(DATA, function(x) {
        Interim <- x[[4]]
        Interim <- Interim %>% dplyr::select(as.character(Radius))  #select the radius length (required to pick up choice from length list)
        sd(Interim[[1]]) / sqrt(length(Interim[[1]]))
      })

      #Generate Sample-wise random distributions
      if(By_Sample_Random){
        #Generate Random results (Average TRIO score)
        RESULTS_RANDOM <-purrr::map_dbl(DATA_RANDOM, function(x) {
          Interim <- x[[4]]
          Interim <- Interim %>% dplyr::select(as.character(Radius))  #select the radius length (required to pick up choice from length list)
          mean(Interim[[1]])
        })
        #Generate the standard error
        RESULTS_RANDOM_SE <-purrr::map_dbl(DATA_RANDOM, function(x) {
          Interim <- x[[4]]
          Interim <- Interim %>% dplyr::select(as.character(Radius))
          sd(Interim[[1]]) / sqrt(length(Interim[[1]]))
        })
      }

      #Generate Experiment-wise random distribution
      if(!By_Sample_Random){
        RESULTS_RANDOM <-  mean(unlist(
          purrr::map(DATA_RANDOM, function(x) {
            Interim <- x[[4]]
            Interim <- Interim %>% dplyr::select(as.character(Radius))  #select the radius length (required to pick up choice from length list)
            Interim
          })
        ))
        RESULTS_RANDOM_SE <- sd(unlist(
          purrr::map(DATA_RANDOM, function(x) {
            Interim <- x[[4]]
            Interim <- Interim %>% dplyr::select(as.character(Radius))  #select the radius length (required to pick up choice from length list)
            Interim
          })
        )) / sqrt(length(unlist(
          purrr::map(DATA_RANDOM, function(x) {
            Interim <- x[[4]]
            Interim <- Interim %>% dplyr::select(as.character(Radius))  #select the radius length (required to pick up choice from length list)
            Interim
          })
        )))
      }

      #Generate the number of cells
      Cell_A <-purrr::map_dbl(DATA, function(x) x[[1]][[1,2]])
      Cell_B <-purrr::map_dbl(DATA, function(x) x[[1]][[2,2]])
      Cell_C <-purrr::map_dbl(DATA, function(x) x[[1]][[3,2]])
      Name_A <- unique(purrr::map_chr(DATA, function(x) x[[1]][[1,1]]))
      Name_B <- unique(purrr::map_chr(DATA, function(x) x[[1]][[2,1]]))
      Name_C <- unique(purrr::map_chr(DATA, function(x) x[[1]][[3,1]]))

      print("Generating output")
      #Generate results tibble
      RESULTS_tibble <- tibble(Subject_Names = names(RESULTS), Average_TRIO_score = RESULTS)
      RESULTS_tibble$SE <- RESULTS_SE
      RESULTS_tibble <- RESULTS_tibble %>%dplyr::mutate(Result_05CI = Average_TRIO_score - 1.96*SE,
                                                        Result_95CI = Average_TRIO_score + 1.96*SE)
      RESULTS_tibble$Radius <- Radius
      RESULTS_tibble$Cell_A <- Cell_A
      names(RESULTS_tibble)[[7]] <- Name_A
      RESULTS_tibble$Cell_B <- Cell_B
      names(RESULTS_tibble)[[8]] <- Name_B
      RESULTS_tibble$Cell_C <- Cell_C
      names(RESULTS_tibble)[[9]] <- Name_C

      RESULTS_tibble$Random <- RESULTS_RANDOM
      RESULTS_tibble$Random_SE <- RESULTS_RANDOM_SE
      RESULTS_tibble <- RESULTS_tibble %>%dplyr::mutate(Random_05CI = Random - 1.96*Random_SE,
                                                        Random_95CI = Random + 1.96*Random_SE) %>%
        dplyr::mutate(Outside_CI = case_when(Average_TRIO_score > Random_95CI | Average_TRIO_score < Random_05CI ~ "Significant",
                                             TRUE ~ "Not Significant"),
                      Sign = case_when(Average_TRIO_score > Random ~ "Above expected",
                                       TRUE ~ "Below expected"))
      RESULTS_tibble$N_Random_COO <-purrr::map_dbl(DATA_RANDOM, function(DATA_RANDOM) nrow(DATA_RANDOM[[2]]))
      #Establish a lower limit of 0 for CI
      RESULTS_tibble$Result_05CI[RESULTS_tibble$Result_05CI < 0] <- 0
      RESULTS_tibble$Random_05CI[RESULTS_tibble$Random_05CI < 0] <- 0

      print("Generating plot")
      #Plot the results
      plot(RESULTS_tibble %>% ggplot(aes(x = fct_reorder(Subject_Names, Average_TRIO_score))) +
             geom_col(aes(y = Average_TRIO_score), width = 0.5, color = "black", fill = "white", linewidth = 0.7)+
             geom_errorbar(aes(ymin = Result_05CI, ymax = Result_95CI), color = "black", linewidth = 0.9, width = 0.3)+
             geom_errorbar(aes(ymin = Random_05CI, ymax = Random_95CI), color = "red", linewidth = 0.9, width = 0.3)+
             cowplot::theme_cowplot() +
             scale_x_discrete("") +
             scale_y_continuous(str_c("Average TRIO score within a radius of ", as.character(Radius))) +
             theme(axis.text.x = element_text(angle = -90, vjust = 0.5))
      )
      #Return the results
      return(RESULTS_tibble)
    }

    else if(!Include_Random){
      print("Proceeding with analysis")
      #Generate Results
      RESULTS <-purrr::map_dbl(DATA, function(x) {
        Interim <- x[[4]]
        Interim <- Interim %>% dplyr::select(as.character(Radius))  #select the radius length (required to pick up choice from length list)
        mean(Interim[[1]])
      })

      RESULTS_SE <-purrr::map_dbl(DATA, function(x) {
        Interim <- x[[4]]
        Interim <- Interim %>% dplyr::select(as.character(Radius))  #select the radius length (required to pick up choice from length list)
        sd(Interim[[1]]) / sqrt(length(Interim[[1]]))
      })

      #Generate the number of cells
      Cell_A <-purrr::map_dbl(DATA, function(x) x[[1]][[1,2]])
      Cell_B <-purrr::map_dbl(DATA, function(x) x[[1]][[2,2]])
      Cell_C <-purrr::map_dbl(DATA, function(x) x[[1]][[3,2]])
      Name_A <- unique(purrr::map_chr(DATA, function(x) x[[1]][[1,1]]))
      Name_B <- unique(purrr::map_chr(DATA, function(x) x[[1]][[2,1]]))
      Name_C <- unique(purrr::map_chr(DATA, function(x) x[[1]][[3,1]]))

      print("Generating output")
      #Generate results tibble
      RESULTS_tibble <- tibble(Subject_Names = names(RESULTS), Average_TRIO_score = RESULTS)
      RESULTS_tibble$SE <- RESULTS_SE
      RESULTS_tibble <- RESULTS_tibble %>%dplyr::mutate(Result_05CI = Average_TRIO_score - 1.96*SE,
                                                        Result_95CI = Average_TRIO_score + 1.96*SE)
      RESULTS_tibble$Radius <- Radius
      RESULTS_tibble$Cell_A <- Cell_A
      names(RESULTS_tibble)[[7]] <- Name_A
      RESULTS_tibble$Cell_B <- Cell_B
      names(RESULTS_tibble)[[8]] <- Name_B
      RESULTS_tibble$Cell_C <- Cell_C
      names(RESULTS_tibble)[[9]] <- Name_C

      #Establish a lower limit of 0 for CI
      RESULTS_tibble$Result_05CI[RESULTS_tibble$Result_05CI < 0] <- 0

      print("Generating plot")
      #Plot the results
      plot(RESULTS_tibble %>% ggplot(aes(x = fct_reorder(Subject_Names, Average_TRIO_score))) +
             geom_col(aes(y = Average_TRIO_score), width = 0.5, color = "black", fill = "white", linewidth = 0.7)+
             geom_errorbar(aes(ymin = Result_05CI, ymax = Result_95CI), color = "black", linewidth = 0.9, width = 0.3)+
             cowplot::theme_cowplot() +
             scale_x_discrete("") +
             scale_y_continuous(str_c("Average TRIO score within a radius of ", as.character(Radius))) +
             theme(axis.text.x = element_text(angle = -90, vjust = 0.5))
      )
      #Return the results
      return(RESULTS_tibble)
    }
  }
