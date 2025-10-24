#' Calculates the min distance to generate a TRIO between three cell types
#'
#' The function calculates the average min distance to generate a TRIO between three cell types for every image.
#'
#' @param N_cores Integer. Number of cores to parallelize your computation.
#' @param DATA A list containing trio cell distance matrix. It can be computed using [Trio_Distance_matrix_generator()] function.
#' @param DATA_RANDOM (OPTIONAL) A list containing random trio cell distance matrix. It can be computed using [Trio_Random_Distance_matrix_generator()] function.
#' @param Include_Random A logical value indicating if the Random background should be used.
#' @param By_Sample_Random If a random background needs to be included, a logical value indicating if the random background should be calculated by sample or by experiment (see details).
#'
#' @details
#' If By_Sample_Random is TRUE, for every sample, a specific random background will be calculated to compute expected spatial interaction metrics. If By_Sample_Random is FALSE, the random background will be calculated using all image information. The same random background will be used for all images.
#'
#'
#' @returns A tibble containing a summary by sample of spatial interactions.
#'
#' @examples
#' #Generate distance matrix and random distance matrix------------
#'TRIO_Distance <-
#' Trio_Distance_matrix_generator(
#'     N_cores = 1,
#'     DATA = CSM_Phenotypecell_test,
#'     Cell_Of_Origin = "TUMOR",
#'     Target_Cell_1 = "CD8_GZMBneg",
#'     Target_Cell_2 = "CD8_GZMBpos",
#'     Perform_edge_correction = FALSE
#' )
#'
#'TRIO_RANDOM <-
#'  Trio_Random_Distance_matrix_generator(
#'   N_cores = 1,
#'   DATA = CSM_Phenotypecell_test,
#'   Cell_Of_Origin = "TUMOR",
#'   Target_Cell_1 = "CD8_GZMBneg",
#'   Target_Cell_2 = "CD8_GZMBpos",
#'   Random_cells_per_sample = 10,
#'   Perform_edge_correction = FALSE
#')
#'
#'#Analyze the spatial association------------------------------
#'Trio_Min_Distance_analyzer(
#'    N_cores = 1,
#'    DATA = TRIO_Distance,
#'    DATA_RANDOM = TRIO_RANDOM,
#'    Include_Random = TRUE,
#'    By_Sample_Random = TRUE
#')
#'
#' @export

Trio_Min_Distance_analyzer <-
  function(N_cores = NULL,
           DATA = NULL,
           DATA_RANDOM = NULL,
           Include_Random = NULL,
           By_Sample_Random = NULL
  ) {
    #Check arguments
    if(!is.logical(Include_Random)) stop("Include_Random should be a logical value")
    if(!all(N_cores >= 1 & N_cores%%1 == 0)) stop("N_cores must be an integer value > 0")
    #Specify if random samples should be include in the computation
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

      #save exit function if parallelization fails
      on.exit({
        future::plan("future::sequential")
        gc()
      })

      #Generate results (average min distance to TRIO)
      future::plan("future::multisession", workers = N_cores)
      options(future.globals.maxSize = Inf, future.rng.onMisuse = "ignore")
      furrr::furrr_options(scheduling = Inf)
      RESULTS <- furrr::future_map_dbl(DATA, function(DATA) {
        #Obtain the 2 individual distance matrices
        Interim_A_to_B <- DATA[[2]]
        Interim_A_to_C <- DATA[[3]]
        #Calculate the min distance to TRIO by calculating first the min distance from COO to each Target and then calculating the max of both distances
        Max_min_dist <-
          purrr::map_dbl(1:nrow(Interim_A_to_B), function(Cell_A){
            Min_A_to_B <- min(Interim_A_to_B[Cell_A, 2:ncol(Interim_A_to_B)])
            Min_A_to_C <- min(Interim_A_to_C[Cell_A, 2:ncol(Interim_A_to_C)])
            max(c(Min_A_to_B, Min_A_to_C))
          })
        #Calculate the average min distance to TRIO
        mean(Max_min_dist)
      },
      .progress = TRUE)

      #Generate the SE of the mean
      RESULTS_SE <- furrr::future_map_dbl(DATA, function(DATA) {
        Interim_A_to_B <- DATA[[2]]
        Interim_A_to_C <- DATA[[3]]
        #Calculate the min distance to TRIO by calculating first the min distance from COO to each Target and then calculating the max of both distances
        Max_min_dist <-
          purrr::map_dbl(1:nrow(Interim_A_to_B), function(Cell_A){
            Min_A_to_B <- min(Interim_A_to_B[Cell_A, 2:ncol(Interim_A_to_B)])
            Min_A_to_C <- min(Interim_A_to_C[Cell_A, 2:ncol(Interim_A_to_C)])
            max(c(Min_A_to_B, Min_A_to_C))
          })
        #Calculate the SE
        sd(Max_min_dist) / sqrt(length(Max_min_dist))
      },
      .progress = TRUE)

      #Calculate Sample-wise random distribution
      if(By_Sample_Random){
        #Generate Random Results (average min distance to TRIO)
        RESULTS_RANDOM <- furrr::future_map_dbl(DATA_RANDOM, function(DATA) {
          Interim_A_to_B <- DATA[[2]]
          Interim_A_to_C <- DATA[[3]]
          #Calculate the min distance to TRIO by calculating first the min distance from COO to each Target and then calculating the max of both distances
          Max_min_dist <-
            purrr::map_dbl(1:nrow(Interim_A_to_B), function(Cell_A){
              Min_A_to_B <- min(Interim_A_to_B[Cell_A, 2:ncol(Interim_A_to_B)])
              Min_A_to_C <- min(Interim_A_to_C[Cell_A, 2:ncol(Interim_A_to_C)])
              max(c(Min_A_to_B, Min_A_to_C))
            })
          #Calculate the average min distance to TRIO
          mean(Max_min_dist)
        }, .progress = TRUE)

        #Generate the SE of the RANDOM mean
        RESULTS_RANDOM_SE <- furrr::future_map_dbl(DATA_RANDOM, function(DATA) {
          Interim_A_to_B <- DATA[[2]]
          Interim_A_to_C <- DATA[[3]]
          #Calculate the min distance to TRIO by calculating first the min distance from COO to each Target and then calculating the max of both distances
          Max_min_dist <-
            purrr::map_dbl(1:nrow(Interim_A_to_B), function(Cell_A){
              Min_A_to_B <- min(Interim_A_to_B[Cell_A, 2:ncol(Interim_A_to_B)])
              Min_A_to_C <- min(Interim_A_to_C[Cell_A, 2:ncol(Interim_A_to_C)])
              max(c(Min_A_to_B, Min_A_to_C))
            })
          #Calculate the SE
          sd(Max_min_dist) / sqrt(length(Max_min_dist))
        }, .progress = TRUE)
      }

      #Calculate Experiment-wise random distribution
      if(!By_Sample_Random){
        #Generate Random Results (average min distance to TRIO)
        RESULTS_RANDOM <- mean(unlist(
          furrr::future_map(DATA_RANDOM, function(DATA) {
            Interim_A_to_B <- DATA[[2]]
            Interim_A_to_C <- DATA[[3]]
            #Calculate the min distance to TRIO by calculating first the min distance from COO to each Target and then calculating the max of both distances
            Max_min_dist <-
              purrr::map_dbl(1:nrow(Interim_A_to_B), function(Cell_A){
                Min_A_to_B <- min(Interim_A_to_B[Cell_A, 2:ncol(Interim_A_to_B)])
                Min_A_to_C <- min(Interim_A_to_C[Cell_A, 2:ncol(Interim_A_to_C)])
                max(c(Min_A_to_B, Min_A_to_C))
              })
            Max_min_dist
          }, .progress = TRUE)
        ))

        #Generate the SE of the RANDOM mean
        RESULTS_RANDOM_SE <-
          sd(unlist(
            purrr::map(DATA_RANDOM, function(DATA) {
              Interim_A_to_B <- DATA[[2]]
              Interim_A_to_C <- DATA[[3]]
              #Calculate the min distance to TRIO by calculating first the min distance from COO to each Target and then calculating the max of both distances
              Max_min_dist <-
                purrr::map_dbl(1:nrow(Interim_A_to_B), function(Cell_A){
                  Min_A_to_B <- min(Interim_A_to_B[Cell_A, 2:ncol(Interim_A_to_B)])
                  Min_A_to_C <- min(Interim_A_to_C[Cell_A, 2:ncol(Interim_A_to_C)])
                  max(c(Min_A_to_B, Min_A_to_C))
                })
              Max_min_dist
            }, .progress = list(clear = F,
                                name = "Calculate Random CI - Stage 1",
                                show_after = 1,
                                type = "iterator"))
          )) / sqrt(length(unlist(
            purrr::map(DATA_RANDOM, function(DATA) {
              Interim_A_to_B <- DATA[[2]]
              Interim_A_to_C <- DATA[[3]]
              #Calculate the min distance to TRIO by calculating first the min distance from COO to each Target and then calculating the max of both distances
              Max_min_dist <-
                purrr::map_dbl(1:nrow(Interim_A_to_B), function(Cell_A){
                  Min_A_to_B <- min(Interim_A_to_B[Cell_A, 2:ncol(Interim_A_to_B)])
                  Min_A_to_C <- min(Interim_A_to_C[Cell_A, 2:ncol(Interim_A_to_C)])
                  max(c(Min_A_to_B, Min_A_to_C))
                })
              Max_min_dist
            }, .progress = list(clear = F,
                                name = "Calculate Random CI - Stage 2",
                                show_after = 1,
                                type = "iterator"))
          )))
      }

      #Generate results tibble
      RESULTS_tibble <- tibble(Subject_Names = names(RESULTS), value = RESULTS)
      names(RESULTS_tibble)[2] <- "Min_dist_to_TRIO"

      RESULTS_tibble$SE <- RESULTS_SE
      RESULTS_tibble <- RESULTS_tibble %>%dplyr::mutate(Result_05CI = Min_dist_to_TRIO - 1.96*SE,
                                                        Result_95CI = Min_dist_to_TRIO + 1.96*SE)

      RESULTS_tibble$N_COO <-purrr::map_dbl(DATA, function(DATA) sum(!is.na(DATA[[2]]$Cell_Of_Origin_no)))
      RESULTS_tibble$N_Target_1 <-purrr::map_dbl(DATA, function(DATA) sum(!is.na(DATA[[2]][1, -1])))
      RESULTS_tibble$N_Target_2 <-purrr::map_dbl(DATA, function(DATA) sum(!is.na(DATA[[3]][1, -1])))

      RESULTS_tibble$Random <- RESULTS_RANDOM
      RESULTS_tibble$Random_SE <- RESULTS_RANDOM_SE
      RESULTS_tibble <- RESULTS_tibble %>%dplyr::mutate(Random_05CI = Random - 1.96*Random_SE,
                                                        Random_95CI = Random + 1.96*Random_SE) %>%
        dplyr::mutate(Outside_CI = case_when(Min_dist_to_TRIO > Random_95CI | Min_dist_to_TRIO < Random_05CI ~ "Significant",
                                             TRUE ~ "Not Significant"),
                      Sign = case_when(Min_dist_to_TRIO > Random ~ "Above expected",
                                       TRUE ~ "Below expected"))
      RESULTS_tibble$N_Random_COO <-purrr::map_dbl(DATA_RANDOM, function(DATA_RANDOM) nrow(DATA_RANDOM[[2]]))

      #Establish a lower limit of 0 for CI
      RESULTS_tibble$Result_05CI[RESULTS_tibble$Result_05CI < 0] <- 0
      RESULTS_tibble$Random_05CI[RESULTS_tibble$Random_05CI < 0] <- 0

      #Plot the results
      plot(RESULTS_tibble %>% ggplot(aes(x = forcats::fct_reorder(Subject_Names, Min_dist_to_TRIO))) +
             geom_col(aes(y = Min_dist_to_TRIO), width = 0.5, color = "black", fill = "white", linewidth = 0.7)+
             geom_errorbar(aes(ymin = Result_05CI, ymax = Result_95CI), color = "black", linewidth = 0.9, width = 0.3)+
             geom_errorbar(aes(ymin = Random_05CI, ymax = Random_95CI), color = "red", linewidth = 0.9, width = 0.3)+
             cowplot::theme_cowplot() +
             scale_x_discrete("") +
             scale_y_continuous("Min Distance to TRIO") +
             theme(axis.text.x = element_text(angle = -90, vjust = 0.5))
      )

      future::plan("future::sequential")
      gc()
      #Return the final tibble
      return(RESULTS_tibble)
    }

    else if(!Include_Random){
      #save exit function if parallelization fails
      on.exit({
        future::plan("future::sequential")
        gc()
      })

      future::plan("future::multisession", workers = N_cores)
      options(future.globals.maxSize = Inf, future.rng.onMisuse = "ignore")
      furrr::furrr_options(scheduling = Inf)

      #Generate results (average min distance to TRIO)
      RESULTS <- furrr::future_map_dbl(DATA, function(DATA) {
        #Obtain the 2 individual distance matrices
        Interim_A_to_B <- DATA[[2]]
        Interim_A_to_C <- DATA[[3]]
        #Calculate the min distance to TRIO by calculating first the min distance from COO to each Target and then calculating the max of both distances
        Max_min_dist <-
          purrr::map_dbl(1:nrow(Interim_A_to_B), function(Cell_A){
            Min_A_to_B <- min(Interim_A_to_B[Cell_A, 2:ncol(Interim_A_to_B)])
            Min_A_to_C <- min(Interim_A_to_C[Cell_A, 2:ncol(Interim_A_to_C)])
            max(c(Min_A_to_B, Min_A_to_C))
          })
        #Calculate the average min distance to TRIO
        mean(Max_min_dist)
      }, .progress = TRUE)

      #Generate the SE of the mean
      RESULTS_SE <- furrr::future_map_dbl(DATA, function(DATA) {
        Interim_A_to_B <- DATA[[2]]
        Interim_A_to_C <- DATA[[3]]
        #Calculate the min distance to TRIO by calculating first the min distance from COO to each Target and then calculating the max of both distances
        Max_min_dist <-
          purrr::map_dbl(1:nrow(Interim_A_to_B), function(Cell_A){
            Min_A_to_B <- min(Interim_A_to_B[Cell_A, 2:ncol(Interim_A_to_B)])
            Min_A_to_C <- min(Interim_A_to_C[Cell_A, 2:ncol(Interim_A_to_C)])
            max(c(Min_A_to_B, Min_A_to_C))
          })
        #Calculate the SE
        sd(Max_min_dist) / sqrt(length(Max_min_dist))
      }, .progress = TRUE)

      #Generate results tibble
      RESULTS_tibble <- tibble(Subject_Names = names(RESULTS), value = RESULTS)
      names(RESULTS_tibble)[2] <- "Min_dist_to_TRIO"

      RESULTS_tibble$SE <- RESULTS_SE
      RESULTS_tibble <- RESULTS_tibble %>%dplyr::mutate(Result_05CI = Min_dist_to_TRIO - 1.96*SE,
                                                        Result_95CI = Min_dist_to_TRIO + 1.96*SE)

      RESULTS_tibble$N_COO <-purrr::map_dbl(DATA, function(DATA) sum(!is.na(DATA[[2]]$Cell_Of_Origin_no)))
      RESULTS_tibble$N_Target_1 <-purrr::map_dbl(DATA, function(DATA) sum(!is.na(DATA[[2]][1, -1])))
      RESULTS_tibble$N_Target_2 <-purrr::map_dbl(DATA, function(DATA) sum(!is.na(DATA[[3]][1, -1])))

      #Establish a lower limit of 0 for CI
      RESULTS_tibble$Result_05CI[RESULTS_tibble$Result_05CI < 0] <- 0

      #Plot the results
      plot(RESULTS_tibble %>% ggplot(aes(x = forcats::fct_reorder(Subject_Names, Min_dist_to_TRIO))) +
             geom_col(aes(y = Min_dist_to_TRIO), width = 0.5, color = "black", fill = "white", linewidth = 0.7)+
             geom_errorbar(aes(ymin = Result_05CI, ymax = Result_95CI), color = "black", linewidth = 0.9, width = 0.3)+
             cowplot::theme_cowplot() +
             scale_x_discrete("") +
             scale_y_continuous("Average Min Distance to TRIO") +
             theme(axis.text.x = element_text(angle = -90, vjust = 0.5))
      )

      future::plan("future::sequential")
      gc()
      #Return the final tibble
      return(RESULTS_tibble)
    }
  }
