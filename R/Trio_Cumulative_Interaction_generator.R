#' Calculates the cumulative interaction matrix between three cells
#'
#' The cumulative interaction matrix contains information regarding the total potential interactions between COO and target cells. For a a given set of increasing distances,
#' the number of target cells present from every COO is calculated. The trio score is calculated as √(N Target Cell1*N Target cell 2).
#'
#' @param N_cores Integer. Number of cores to parallelize your computation.
#' @param DATA A list containing cell distance matrix. It can be computed using [Trio_Distance_matrix_generator()] or [Trio_Random_Distance_matrix_generator()] functions.
#' @param Start_from A numeric value indicating the minimum distance to be analyzed.
#' @param Stop_at A numeric value indicating the maximum distance to be analyzed.
#' @param Sampling_frequency A numeric value indicating the distance sampling frequency.
#'
#' @seealso [Trio_Random_Distance_matrix_generator()], [Trio_Distance_matrix_generator()], [Trio_Cells_in_Radius_analyzer()].
#'
#' @returns A list containing the cumulative interaction matrix for every image as well as the trio interaction score.
#'
#' @examples
#' #Generate a data distance or random distance matrix--------------------------
#' TRIO_Distance <-
#' Trio_Distance_matrix_generator(
#'     N_cores = 1,
#'     DATA = CSM_Phenotypecell_test,
#'     Cell_Of_Origin = "TUMOR",
#'     Target_Cell_1 = "CD8_GZMBneg",
#'     Target_Cell_2 = "CD8_GZMBpos",
#'     Perform_edge_correction = FALSE
#' )
#' #Calculate the cumulative interaction every 25 pixels up to 100 pixels------
#'Trio_Cumulative_Interaction_generator(
#'   N_cores = 1,
#'   DATA = TRIO_Distance,
#'   Start_from = 25,
#'   Stop_at = 100,
#'   Sampling_frequency = 25
#')
#'
#' @export

Trio_Cumulative_Interaction_generator <-
  function(N_cores = 1,
           DATA = NULL,
           Start_from = NULL,
           Stop_at = NULL,
           Sampling_frequency = NULL
  ) {
    #Check arguments
    if(!all(N_cores >= 1 & N_cores%%1 == 0)) stop("N_cores must be an integer value > 0")
    if(!all(is.numeric(Start_from), is.numeric(Stop_at), is.numeric(Sampling_frequency), Start_from < Stop_at, Sampling_frequency < (Stop_at - Start_from))) {
      stop("Start_from, Stop_at and Sampling_frequency must be numeric values. Start_from must be smaller than Stop_at. Sampling_frequency must be smaller than the range Start_from - Stop_at")
    }
    if(!all(Start_from >= 0, Stop_at > 0, Sampling_frequency > 0)) stop("Start_from, Stop_at and Sampling_frequency must have positive values")

    #Get our data from the global environment
    DATA <- DATA


    #save exit function if parallelization fails
    on.exit({
      future::plan("future::sequential")
      gc()
    })

    #Prepare our clustering
    future::plan("future::multisession", workers = N_cores)
    options(future.globals.maxSize = Inf, future.rng.onMisuse = "ignore")
    furrr::furrr_options(scheduling = Inf)

    RESULTS <-
      suppressMessages({
        furrr::future_map(DATA, function(x) {
          #Get our tibble and format it adequately and remove NA values if necessary
          Longer_distance_matrix_A_to_B <-  x$A_to_B %>% tidyr::pivot_longer(2:ncol(x$A_to_B)) %>% na.omit()

          Longer_distance_matrix_A_to_C <-  x$A_to_C %>% tidyr::pivot_longer(2:ncol(x$A_to_C)) %>% na.omit()

          #Calculate the cumulative distance tibble with the desired sampling strategy for both A to B
          Cumulative_distance_tibble_A_to_B <-
            purrr::map_dfc(seq(from = Start_from, to = Stop_at, by = Sampling_frequency), function(filter) {
              Interim <- Longer_distance_matrix_A_to_B %>%dplyr::mutate(counts = value <= filter) %>% group_by(Cell_Of_Origin_no, counts) %>% dplyr::count() %>% dplyr::filter(counts == T) %>%
                dplyr::ungroup() %>% dplyr::select(-counts)
              Final <-dplyr::left_join(tibble(Cell_Of_Origin_no = unique(Longer_distance_matrix_A_to_B$Cell_Of_Origin_no)), Interim, by = "Cell_Of_Origin_no") %>%
                dplyr::mutate(number = case_when(is.na(n) ~ 0,
                                                 TRUE ~ n)) %>% dplyr::select(-n)
              names(Final) <- c("Cell_Of_Origin_no", as.character(filter))
              Final
            }) %>% dplyr::select(-contains("Cell")) %>%dplyr::mutate(Cell_Of_Origin_no = unique(Longer_distance_matrix_A_to_B$Cell_Of_Origin_no))
          Cumulative_distance_tibble_A_to_B <- Cumulative_distance_tibble_A_to_B[c(ncol(Cumulative_distance_tibble_A_to_B), 1:(ncol(Cumulative_distance_tibble_A_to_B)-1))]

          #Calculate the cumulative distance tibble with the desired sampling strategy for for both A to C
          Cumulative_distance_tibble_A_to_C <-
            purrr::map_dfc(seq(from = Start_from, to = Stop_at, by = Sampling_frequency), function(filter) {
              Interim <- Longer_distance_matrix_A_to_C %>%dplyr::mutate(counts = value <= filter) %>% group_by(Cell_Of_Origin_no, counts) %>% dplyr::count() %>% dplyr::filter(counts == T) %>%
                dplyr::ungroup() %>% dplyr::select(-counts)
              Final <-dplyr::left_join(tibble(Cell_Of_Origin_no = unique(Longer_distance_matrix_A_to_C$Cell_Of_Origin_no)), Interim, by = "Cell_Of_Origin_no") %>%
                dplyr::mutate(number = case_when(is.na(n) ~ 0,
                                                 TRUE ~ n)) %>% dplyr::select(-n)
              names(Final) <- c("Cell_Of_Origin_no", as.character(filter))
              Final
            }) %>% dplyr::select(-contains("Cell")) %>%dplyr::mutate(Cell_Of_Origin_no = unique(Longer_distance_matrix_A_to_C$Cell_Of_Origin_no))
          Cumulative_distance_tibble_A_to_C <- Cumulative_distance_tibble_A_to_C[c(ncol(Cumulative_distance_tibble_A_to_C), 1:(ncol(Cumulative_distance_tibble_A_to_C)-1))]

          #Now we calculate the TRIO Socore (square root of Cumulative interaction AB * Cumulative interaction AC) combining the cumulative interaction matrices
          Trio_Score <-dplyr::bind_cols(Cumulative_distance_tibble_A_to_C[1], as_tibble(sqrt(Cumulative_distance_tibble_A_to_B[-1] * Cumulative_distance_tibble_A_to_C[-1])))

          #Return the final result
          list(Cell_counts = x[[1]],
               Cum_dist_A_to_B = Cumulative_distance_tibble_A_to_B,
               Cum_dist_A_to_C = Cumulative_distance_tibble_A_to_C,
               Trio_Score = Trio_Score
          )
        },
        .progress = TRUE)
      })
    future::plan("future::sequential")
    gc()

    names(RESULTS) <- names(DATA)
    #Return the results
    return(RESULTS)
  }
