#' Calculates the cumulative interaction matrix
#'
#' The cumulative interaction matrix contains information regarding the total potential interactions between COO and target cells. For a a given set of increasing distances,
#' the number of target cells present from every COO is calculated.
#'
#' @param N_cores Integer. Number of cores to parallelize your computation.
#' @param DATA A list containing cell distance matrix. It can be computed using [Distance_matrix_generator()] or [Random_Distance_matrix_generator()] functions.
#' @param Start_from A numeric value indicating the minimum distance to be analyzed.
#' @param Stop_at A numeric value indicating the maximum distance to be analyzed.
#' @param Sampling_frequency A numeric value indicating the distance sampling frequency.
#'
#' @seealso [Random_Distance_matrix_generator()], [Distance_matrix_generator()], [Distance_analyzer()], [Cells_in_Radius_analyzer()].
#'
#' @returns A list containing the cumulative interaction matrix for every image.
#'
#'
#' @examples
#' #Generate a data distance or random distance matrix--------------------------
#' DATA_Distances <-
#' Distance_matrix_generator(
#'     N_cores = 1,
#'     DATA = CSM_Phenotypecell_test,
#'     Cell_Of_Origin = "CD8_GZMBneg",
#'     Target_Cell = "TUMOR",
#'     Allow_Cero_Distance = FALSE,
#'     Perform_edge_correction = FALSE
#')
#'
#' #Calculate the cumulative interaction every 25 pixels up to 100 pixels------
#'Cumulative_Interaction_generator(
#'    N_cores = 1,
#'    DATA = DATA_Distances,
#'    Start_from = 25,
#'    Stop_at = 100,
#'    Sampling_frequency = 25
#')
#'
#' @export

Cumulative_Interaction_generator <-
  function(N_cores = 1,
           DATA = NULL,
           Start_from = NULL,
           Stop_at = NULL,
           Sampling_frequency = NULL) {
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

    RESULTS <- suppressMessages(
      furrr::future_map(DATA, function(x) {
        #Get our tibble and format it adequately and remove NA values if necessary
        Longer_distance_matrix <-  x[[2]] %>% tidyr::pivot_longer(2:ncol(x[[2]])) %>% na.omit()

        #Calculate the cumulative distance tibble with the desired sampling strategy
        Cumulative_distance_tibble <-
          purrr::map_dfc(seq(from = Start_from, to = Stop_at, by = Sampling_frequency), function(filter) {
            Interim <- Longer_distance_matrix %>% dplyr::mutate(counts = value <= filter) %>% dplyr::group_by(Cell_Of_Origin_no, counts) %>% dplyr::count() %>% dplyr::filter(counts == T) %>%
              dplyr::ungroup() %>% dplyr::select(-counts)
            Final <- dplyr::left_join(tibble(Cell_Of_Origin_no = unique(Longer_distance_matrix$Cell_Of_Origin_no)), Interim, by = "Cell_Of_Origin_no") %>%
              dplyr::mutate(number = case_when(is.na(n) ~ 0,
                                               TRUE ~ n)) %>% dplyr::select(-n)
            names(Final) <- c("Cell_Of_Origin_no", as.character(filter))
            Final
          }) %>% dplyr::select(-contains("Cell")) %>% dplyr::mutate(Cell_Of_Origin_no = unique(Longer_distance_matrix$Cell_Of_Origin_no))
        Cumulative_distance_tibble <- Cumulative_distance_tibble[c(ncol(Cumulative_distance_tibble), 1:(ncol(Cumulative_distance_tibble)-1))]

        #Arrange the results in a list and return it
        list(Cell_counts = x[[1]],
             Cumulative_distance = Cumulative_distance_tibble)
      },
      .progress = TRUE)
    )
    future::plan("future::sequential")
    gc()

    names(RESULTS) <- names(DATA)

    return(RESULTS)
  }
