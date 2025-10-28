#' Analyzes the distance pattern between two cell types for all images
#'
#' The function summarizes the spatial interaction pattern between two cell types given a distance matrix.
#'
#' @param N_cores Integer. Number of cores to parallelize your computation.
#' @param DATA A list containing cell distance matrix. It can be computed using [Distance_matrix_generator()] function.
#' @param DATA_RANDOM (OPTIONAL) A list containing random cell distance matrix. It can be computed using [Random_Distance_matrix_generator()] function.
#' @param Metric A character value indicating the type of analysis strategy. One of the following: Average_Distance, Min_Distance, Max_Distance (see details).
#' @param Include_Random A logical value indicating if the Random background should be used.
#' @param By_Sample_Random If a random background needs to be included, a logical value indicating if the random background should be calculated by sample or by experiment (see details).
#'
#' @details
#' Average distance computes the average distance by image between Cells of origin (COO) and target cells. Min distance averages the minimum distance from every COO to target cells. Max ditance averages
#' the maximum distance between COO and target cells.
#'
#' If By_Sample_Random is TRUE, for every sample, a specific random background will be calculated to compute expected spatial interaction metrics. If By_Sample_Random is FALSE, the random background will be calculated using all image information. The same random background will be used for all images.
#'
#'
#' @seealso [Distance_matrix_generator()], [Random_Distance_matrix_generator()], [Cell_to_Cell_graph_maker()]
#'
#' @returns A tibble containing a summary by sample of spatial interactions.
#'
#' @examples
#' \dontrun{
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
#'
#' #Analyze the spatial association------------------------------
#'Distance_analyzer(
#'    N_cores = 1,
#'    DATA = DATA_Distances,
#'    DATA_RANDOM = RANDOM_Distances,
#'    Metric = "Min_Distance",
#'    Include_Random = TRUE,
#'    By_Sample_Random = TRUE
#')
#' }
#'
#' @export

Distance_analyzer <-
  function(N_cores = 1,
           DATA,
           DATA_RANDOM = NULL,
           Metric,
           Include_Random = FALSE,
           By_Sample_Random = NULL) {
    #Check arguments
    if(!is.logical(Include_Random)) stop("Include_Random should be a logical value")
    if(!Metric %in% c("Average_Distance", "Min_Distance", "Max_Distance")) stop("Metric should be one of the following: Average_Distance, Min_Distance, Max_Distance")
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
      #Generate the future
      future::plan("future::multisession", workers = N_cores)
      options(future.globals.maxSize = Inf, future.rng.onMisuse = "ignore")
      furrr::furrr_options(scheduling = Inf)

      if(Metric == "Average_Distance") {
        print("Calculating Result")
        #Generate results
        RESULTS <- furrr::future_map_dbl(DATA, function(DATA) {
          Interim <- DATA[[2]]
          mean(as.matrix(Interim[-1]), na.rm = T)
        })

        RESULTS_SE <- furrr::future_map_dbl(DATA, function(DATA) {
          Interim <- DATA[[2]]
          sd(as.matrix(Interim[-1]), na.rm = T) / sqrt(sum(!is.na(Interim[-1])))
        })

        print("Calculating Random Background")
        if(By_Sample_Random){
          #Generate Random results
          RESULTS_RANDOM <- furrr::future_map_dbl(DATA_RANDOM, function(DATA) {
            Interim <- DATA[[2]]
            mean(as.matrix(Interim[-1]), na.rm = T)
          })

          #Calculate Standard Error of Random Mean
          RESULTS_RANDOM_SE <- furrr::future_map_dbl(DATA_RANDOM, function(DATA) {
            Interim <- DATA[[2]]
            sd(as.matrix(Interim[-1]), na.rm = T) / sqrt(sum(!is.na(Interim[-1])))
          })
        }

        if(!By_Sample_Random){
          #Generate Random results for all Random samples at the same time
          RESULTS_RANDOM <- mean(unlist(
            purrr::map(DATA_RANDOM, function(DATA){
              DATA[[2]][,-1]
            })),
            na.rm = T)
          #Calculate Standard Error of Random Mean for all Random samples at the same time
          RESULTS_RANDOM_SE <- sd(unlist(
            purrr::map(DATA_RANDOM, function(DATA){
              DATA[[2]][,-1]
            })),
            na.rm = T) / sqrt(sum(!is.na(unlist(
              purrr::map(DATA_RANDOM, function(DATA){
                DATA[[2]][,-1]
              })))))
        }


        #Generate results tibble
        RESULTS_tibble <- tibble(Subject_Names = names(RESULTS), value = RESULTS)
        names(RESULTS_tibble)[2] <- Metric

        RESULTS_tibble$SE <- RESULTS_SE
        RESULTS_tibble <- RESULTS_tibble %>%dplyr::mutate(Result_05CI = Average_Distance - 1.96*SE,
                                                          Result_95CI = Average_Distance + 1.96*SE)

        RESULTS_tibble$N_COO <-purrr::map_dbl(DATA, function(DATA) sum(!is.na(DATA[[2]]$Cell_Of_Origin_no)))
        RESULTS_tibble$N_Target <-purrr::map_dbl(DATA, function(DATA) sum(!is.na(DATA[[2]][1, -1])))

        RESULTS_tibble$Random <- RESULTS_RANDOM
        RESULTS_tibble$Random_SE <- RESULTS_RANDOM_SE
        RESULTS_tibble <- RESULTS_tibble %>%dplyr::mutate(Random_05CI = Random - 1.96*Random_SE,
                                                          Random_95CI = Random + 1.96*Random_SE) %>%
          dplyr::mutate(Outside_CI = case_when(Average_Distance > Random_95CI | Average_Distance < Random_05CI ~ "Significant",
                                               TRUE ~ "Not Significant"),
                        Sign = case_when(Average_Distance > Random ~ "Above expected",
                                         TRUE ~ "Below expected"))
        RESULTS_tibble$N_Random_COO <-purrr::map_dbl(DATA_RANDOM, function(DATA_RANDOM) nrow(DATA_RANDOM[[2]]))

        #Establish a lower limit of 0 for CI
        RESULTS_tibble$Result_05CI[RESULTS_tibble$Result_05CI < 0] <- 0
        RESULTS_tibble$Random_05CI[RESULTS_tibble$Random_05CI < 0] <- 0

        print("Generating plot")
        #Plot the results
        plot(RESULTS_tibble %>% ggplot(aes(x = forcats::fct_reorder(Subject_Names, Average_Distance))) +
               geom_col(aes(y = Average_Distance), width = 0.5, color = "black", fill = "white", linewidth = 0.7)+
               geom_errorbar(aes(ymin = Result_05CI, ymax = Result_95CI), color = "black", linewidth = 0.9, width = 0.3)+
               geom_errorbar(aes(ymin = Random_05CI, ymax = Random_95CI), color = "red", linewidth = 0.9, width = 0.3)+
               cowplot::theme_cowplot() +
               scale_x_discrete("") +
               scale_y_continuous("Average Distance") +
               theme(axis.text.x = element_text(angle = -90, vjust = 0.5))
        )

        future::plan("future::sequential")
        gc()
        #Return the final results
        return(RESULTS_tibble)
      }

      else if(Metric == "Min_Distance") {
        print("Calculating Result")
        #Generate results
        RESULTS <- furrr::future_map_dbl(DATA, function(DATA) {
          Interim <- DATA[[2]]
          mean(apply(Interim[-1], MARGIN = 1, function(x) min(x, na.rm = T)))
        })

        RESULTS_SE <- furrr::future_map_dbl(DATA, function(DATA) {
          Interim <- DATA[[2]]
          sd(unlist(apply(Interim[-1], MARGIN = 1, function(x) min(x, na.rm = T)))) / sqrt(nrow(Interim))
        })

        print("Calculating Random Background")
        if(By_Sample_Random){
          #Generate Random Results
          RESULTS_RANDOM <- furrr::future_map_dbl(DATA_RANDOM, function(DATA) {
            Interim <- DATA[[2]]
            mean(apply(Interim[-1], MARGIN = 1, function(x) min(x, na.rm = T)))
          })

          RESULTS_RANDOM_SE <- furrr::future_map_dbl(DATA_RANDOM, function(DATA) {
            Interim <- DATA[[2]]
            sd(unlist(apply(Interim[-1], MARGIN = 1, function(x) min(x, na.rm = T)))) / sqrt(nrow(Interim))
          })
        }

        if(!By_Sample_Random){
          #Generate the Random Results overall
          RESULTS_RANDOM <- mean(unlist(
            furrr::future_map(DATA_RANDOM, function(DATA){
              Interim <- DATA[[2]]
              apply(Interim[-1], MARGIN = 1, function(x) min(x, na.rm = T))
            })
          ))

          #Generate the Random standard deviation overall
          RESULTS_RANDOM_SE <- sd(unlist(
            furrr::future_map(DATA_RANDOM, function(DATA){
              Interim <- DATA[[2]]
              apply(Interim[-1], MARGIN = 1, function(x) min(x, na.rm = T))
            })
          )) / sqrt(sum(!is.na(unlist(
            furrr::future_map(DATA_RANDOM, function(DATA){
              Interim <- DATA[[2]]
              apply(Interim[-1], MARGIN = 1, function(x) min(x, na.rm = T))
            })
          ))))
        }

        #Generate results tibble
        RESULTS_tibble <- tibble(Subject_Names = names(RESULTS), value = RESULTS)
        names(RESULTS_tibble)[2] <- Metric

        RESULTS_tibble$SE <- RESULTS_SE
        RESULTS_tibble <- RESULTS_tibble %>%dplyr::mutate(Result_05CI = Min_Distance - 1.96*SE,
                                                          Result_95CI = Min_Distance + 1.96*SE)

        RESULTS_tibble$N_COO <-purrr::map_dbl(DATA, function(DATA) sum(!is.na(DATA[[2]]$Cell_Of_Origin_no)))
        RESULTS_tibble$N_Target <-purrr::map_dbl(DATA, function(DATA) sum(!is.na(DATA[[2]][1, -1])))

        RESULTS_tibble$Random <- RESULTS_RANDOM
        RESULTS_tibble$Random_SE <- RESULTS_RANDOM_SE
        RESULTS_tibble <- RESULTS_tibble %>%dplyr::mutate(Random_05CI = Random - 1.96*Random_SE,
                                                          Random_95CI = Random + 1.96*Random_SE) %>%
          dplyr::mutate(Outside_CI = case_when(Min_Distance > Random_95CI | Min_Distance < Random_05CI ~ "Significant",
                                               TRUE ~ "Not Significant"),
                        Sign = case_when(Min_Distance > Random ~ "Above expected",
                                         TRUE ~ "Below expected"))
        RESULTS_tibble$N_Random_COO <-purrr::map_dbl(DATA_RANDOM, function(DATA_RANDOM) nrow(DATA_RANDOM[[2]]))

        #Establish a lower limit of 0 for CI
        RESULTS_tibble$Result_05CI[RESULTS_tibble$Result_05CI < 0] <- 0
        RESULTS_tibble$Random_05CI[RESULTS_tibble$Random_05CI < 0] <- 0

        print("Generating plot")
        #Plot the results
        plot(RESULTS_tibble %>% ggplot(aes(x = forcats::fct_reorder(Subject_Names, Min_Distance))) +
               geom_col(aes(y = Min_Distance), width = 0.5, color = "black", fill = "white", linewidth = 0.7)+
               geom_errorbar(aes(ymin = Result_05CI, ymax = Result_95CI), color = "black", linewidth = 0.9, width = 0.3)+
               geom_errorbar(aes(ymin = Random_05CI, ymax = Random_95CI), color = "red", linewidth = 0.9, width = 0.3)+
               cowplot::theme_cowplot() +
               scale_x_discrete("") +
               scale_y_continuous("Average Min Distance") +
               theme(axis.text.x = element_text(angle = -90, vjust = 0.5))
        )

        future::plan("future::sequential")
        gc()
        return(RESULTS_tibble)
      }

      else if(Metric == "Max_Distance") {
        print("Calculating Result")
        #Generate results
        RESULTS <- furrr::future_map_dbl(DATA, function(DATA) {
          Interim <- DATA[[2]]
          mean(apply(Interim[-1], MARGIN = 1, function(x) max(x, na.rm = T)))
        })

        RESULTS_SE <- furrr::future_map_dbl(DATA, function(DATA) {
          Interim <- DATA[[2]]
          sd(unlist(apply(Interim[-1], MARGIN = 1, function(x) max(x, na.rm = T)))) / sqrt(nrow(Interim))
        })

        print("Calculating Random Background")
        if(By_Sample_Random){
          #Generate Random Results
          RESULTS_RANDOM <- furrr::future_map_dbl(DATA_RANDOM, function(DATA) {
            Interim <- DATA[[2]]
            mean(apply(Interim[-1], MARGIN = 1, function(x) max(x, na.rm = T)))
          })

          RESULTS_RANDOM_SE <- furrr::future_map_dbl(DATA_RANDOM, function(DATA) {
            Interim <- DATA[[2]]
            sd(unlist(apply(Interim[-1], MARGIN = 1, function(x) max(x, na.rm = T)))) / sqrt(nrow(Interim))
          })
        }

        if(!By_Sample_Random){
          #Generate the Random Results overall
          RESULTS_RANDOM <- mean(unlist(
            purrr::map(DATA_RANDOM, function(DATA){
              Interim <- DATA[[2]]
              apply(Interim[-1], MARGIN = 1, function(x) max(x, na.rm = T))
            })
          ))

          #Generate the Random standard deviation overall
          RESULTS_RANDOM_SE <- sd(unlist(
            purrr::map(DATA_RANDOM, function(DATA){
              Interim <- DATA[[2]]
              apply(Interim[-1], MARGIN = 1, function(x) max(x, na.rm = T))
            })
          )) / sqrt(sum(!is.na(unlist(
            purrr::map(DATA_RANDOM, function(DATA){
              Interim <- DATA[[2]]
              apply(Interim[-1], MARGIN = 1, function(x) max(x, na.rm = T))
            })
          ))))
        }


        #Generate results tibble
        RESULTS_tibble <- tibble(Subject_Names = names(RESULTS), value = RESULTS)
        names(RESULTS_tibble)[2] <- Metric

        RESULTS_tibble$SE <- RESULTS_SE
        RESULTS_tibble <- RESULTS_tibble %>%dplyr::mutate(Result_05CI = Max_Distance - 1.96*SE,
                                                          Result_95CI = Max_Distance + 1.96*SE)

        RESULTS_tibble$N_COO <-purrr::map_dbl(DATA, function(DATA) sum(!is.na(DATA[[2]]$Cell_Of_Origin_no)))
        RESULTS_tibble$N_Target <-purrr::map_dbl(DATA, function(DATA) sum(!is.na(DATA[[2]][1, -1])))

        RESULTS_tibble$Random <- RESULTS_RANDOM
        RESULTS_tibble$Random_SE <- RESULTS_RANDOM_SE
        RESULTS_tibble <- RESULTS_tibble %>%dplyr::mutate(Random_05CI = Random - 1.96*Random_SE,
                                                          Random_95CI = Random + 1.96*Random_SE) %>%
          dplyr::mutate(Outside_CI = case_when(Max_Distance > Random_95CI | Max_Distance < Random_05CI ~ "Significant",
                                               TRUE ~ "Not Significant"),
                        Sign = case_when(Max_Distance > Random ~ "Above expected",
                                         TRUE ~ "Below expected"))
        RESULTS_tibble$N_Random_COO <-purrr::map_dbl(DATA_RANDOM, function(DATA_RANDOM) nrow(DATA_RANDOM[[2]]))

        #Establish a lower limit of 0 for CI
        RESULTS_tibble$Result_05CI[RESULTS_tibble$Result_05CI < 0] <- 0
        RESULTS_tibble$Random_05CI[RESULTS_tibble$Random_05CI < 0] <- 0

        print("Generating plot")
        #Plot the results
        plot(RESULTS_tibble %>% ggplot(aes(x = forcats::fct_reorder(Subject_Names, Max_Distance))) +
               geom_col(aes(y = Max_Distance), width = 0.5, color = "black", fill = "white", linewidth = 0.7)+
               geom_errorbar(aes(ymin = Result_05CI, ymax = Result_95CI), color = "black", linewidth = 0.9, width = 0.3)+
               geom_errorbar(aes(ymin = Random_05CI, ymax = Random_95CI), color = "red", linewidth = 0.9, width = 0.3)+
               cowplot::theme_cowplot() +
               scale_x_discrete("") +
               scale_y_continuous("Average Max Distance") +
               theme(axis.text.x = element_text(angle = -90, vjust = 0.5))
        )

        future::plan("future::sequential")
        gc()
        return(RESULTS_tibble)
      }
    }

    #Specify what to do if random samples are not included in the analysis
    else if(!Include_Random){

      #save exit function if parallelization fails
      on.exit({
        future::plan("future::sequential")
        gc()
      })
      #Generate future cluster
      future::plan("future::multisession", workers = N_cores)
      options(future.globals.maxSize = Inf, future.rng.onMisuse = "ignore")
      furrr::furrr_options(scheduling = Inf)

      if(Metric == "Average_Distance") {
        print("Calculating Result")

        #Generate results
        RESULTS <- furrr::future_map_dbl(DATA, function(DATA) {
          Interim <- DATA[[2]]
          mean(as.matrix(Interim[-1]), na.rm = T)
        })

        RESULTS_SE <- furrr::future_map_dbl(DATA, function(DATA) {
          Interim <- DATA[[2]]
          sd(as.matrix(Interim[-1]), na.rm = T) / sqrt(sum(!is.na(Interim[-1])))
        })

        #Generate results tibble
        RESULTS_tibble <- tibble(Subject_Names = names(RESULTS), value = RESULTS)
        names(RESULTS_tibble)[2] <- Metric

        RESULTS_tibble$SE <- RESULTS_SE
        RESULTS_tibble <- RESULTS_tibble %>%dplyr::mutate(Result_05CI = Average_Distance - 1.96*SE,
                                                          Result_95CI = Average_Distance + 1.96*SE)

        RESULTS_tibble$N_COO <-purrr::map_dbl(DATA, function(DATA) sum(!is.na(DATA[[2]]$Cell_Of_Origin_no)))
        RESULTS_tibble$N_Target <-purrr::map_dbl(DATA, function(DATA) sum(!is.na(DATA[[2]][1, -1])))

        #Establish a lower limit of 0 for CI
        RESULTS_tibble$Result_05CI[RESULTS_tibble$Result_05CI < 0] <- 0

        #plot the results
        plot(RESULTS_tibble %>% ggplot(aes(x = forcats::fct_reorder(Subject_Names, Average_Distance))) +
               geom_col(aes(y = Average_Distance), width = 0.5, color = "black", fill = "white", linewidth = 0.7)+
               geom_errorbar(aes(ymin = Result_05CI, ymax = Result_95CI), color = "black", linewidth = 0.9, width = 0.3)+
               cowplot::theme_cowplot() +
               scale_x_discrete("") +
               scale_y_continuous("Average Distance") +
               theme(axis.text.x = element_text(angle = -90, vjust = 0.5))
        )

        future::plan("future::sequential")
        gc()
        return(RESULTS_tibble)
      }

      else if(Metric == "Min_Distance") {
        print("Calculating Result")
        #Generate results
        RESULTS <- furrr::future_map_dbl(DATA, function(DATA) {
          Interim <- DATA[[2]]
          mean(apply(Interim[-1], MARGIN = 1, function(x) min(x, na.rm = T)))
        })

        RESULTS_SE <- furrr::future_map_dbl(DATA, function(DATA) {
          Interim <- DATA[[2]]
          sd(unlist(apply(Interim[-1], MARGIN = 1, function(x) min(x, na.rm = T)))) / sqrt(nrow(Interim))
        })

        #Generate results tibble
        RESULTS_tibble <- tibble(Subject_Names = names(RESULTS), value = RESULTS)
        names(RESULTS_tibble)[2] <- Metric

        RESULTS_tibble$SE <- RESULTS_SE
        RESULTS_tibble <- RESULTS_tibble %>%dplyr::mutate(Result_05CI = Min_Distance - 1.96*SE,
                                                          Result_95CI = Min_Distance + 1.96*SE)

        RESULTS_tibble$N_COO <-purrr::map_dbl(DATA, function(DATA) sum(!is.na(DATA[[2]]$Cell_Of_Origin_no)))
        RESULTS_tibble$N_Target <-purrr::map_dbl(DATA, function(DATA) sum(!is.na(DATA[[2]][1, -1])))

        #Establish a lower limit of 0 for CI
        RESULTS_tibble$Result_05CI[RESULTS_tibble$Result_05CI < 0] <- 0

        #Plot the results
        plot(RESULTS_tibble %>% ggplot(aes(x = forcats::fct_reorder(Subject_Names, Min_Distance))) +
               geom_col(aes(y = Min_Distance), width = 0.5, color = "black", fill = "white", linewidth = 0.7)+
               geom_errorbar(aes(ymin = Result_05CI, ymax = Result_95CI), color = "black", linewidth = 0.9, width = 0.3)+
               cowplot::theme_cowplot() +
               scale_x_discrete("") +
               scale_y_continuous("Average Min Distance") +
               theme(axis.text.x = element_text(angle = -90, vjust = 0.5))
        )

        future::plan("future::sequential")
        gc()
        return(RESULTS_tibble)
      }

      else if(Metric == "Max_Distance") {
        print("Calculating Result")
        #Generate results
        RESULTS <- furrr::future_map_dbl(DATA, function(DATA) {
          Interim <- DATA[[2]]
          mean(apply(Interim[-1], MARGIN = 1, function(x) max(x, na.rm = T)))
        })

        RESULTS_SE <- furrr::future_map_dbl(DATA, function(DATA) {
          Interim <- DATA[[2]]
          sd(unlist(apply(Interim[-1], MARGIN = 1, function(x) max(x, na.rm = T)))) / sqrt(nrow(Interim))
        })

        #Generate results tibble
        RESULTS_tibble <- tibble(Subject_Names = names(RESULTS), value = RESULTS)
        names(RESULTS_tibble)[2] <- Metric

        RESULTS_tibble$SE <- RESULTS_SE
        RESULTS_tibble <- RESULTS_tibble %>%dplyr::mutate(Result_05CI = Max_Distance - 1.96*SE,
                                                          Result_95CI = Max_Distance + 1.96*SE)

        RESULTS_tibble$N_COO <-purrr::map_dbl(DATA, function(DATA) sum(!is.na(DATA[[2]]$Cell_Of_Origin_no)))
        RESULTS_tibble$N_Target <-purrr::map_dbl(DATA, function(DATA) sum(!is.na(DATA[[2]][1, -1])))

        #Establish a lower limit of 0 for CI
        RESULTS_tibble$Result_05CI[RESULTS_tibble$Result_05CI < 0] <- 0

        #plot the results
        plot(RESULTS_tibble %>% ggplot(aes(x = forcats::fct_reorder(Subject_Names, Max_Distance))) +
               geom_col(aes(y = Max_Distance), width = 0.5, color = "black", fill = "white", linewidth = 0.7)+
               geom_errorbar(aes(ymin = Result_05CI, ymax = Result_95CI), color = "black", linewidth = 0.9, width = 0.3)+
               cowplot::theme_cowplot() +
               scale_x_discrete("") +
               scale_y_continuous("Average Max Distance") +
               theme(axis.text.x = element_text(angle = -90, vjust = 0.5))
        )

        future::plan("future::sequential")
        gc()
        return(RESULTS_tibble)
      }
    }
  }
