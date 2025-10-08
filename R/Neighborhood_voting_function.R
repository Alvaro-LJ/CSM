#' Performs elections within tiled images with neighborhood information
#'
#' Given a tiled image, the function calculates neighborhood prevalence by tile. It then assigns neighborhood identity to each tile based on results.
#'
#' @param Tiled_images A list containing tiled images obtained using [Image_tiling_processing_function()] with neighborhood information.
#' @param Minimum_cell_no_per_tile An integer indicating the minimum number of cells that a tile must contain. Tiles below the limit will not be included in the analysis.
#' @param Phenotypes_included A character vector indicating the phenotype labels that will be included in the analysis.
#'
#' @returns A list containing image information with by-tile neighborhood counts and a tibble with a summary of the number of winner per image.
#' @export

Neighborhood_voting_function <-
  function(N_cores = NULL,
           Tiled_Images = NULL,
           Minimum_cell_no_per_tile = NULL,
           Neighborhoods_included = NULL) {
    #Import data from a Tiled_Images list
    Tiled_Images <- Tiled_Images

    #Check arguments
    if(!all(N_cores >= 1 & N_cores%%1 == 0)) stop("N_cores must be an integer value > 0")
    if(!all(Minimum_cell_no_per_tile >= 1 & Minimum_cell_no_per_tile%%1 == 0)) stop("Minimum_cell_no_per_tile must be an integer value > 0")
    #Check that the Neighborhoods_included are present in the data
    if(!all(Neighborhoods_included %in% unique(unlist(purrr::map(Tiled_Images, function(df) df[[2]]$Neighborhood_assignment))))) {
      stop(paste0("Neighborhoods_included included must be any of: ", stringr::str_c(unique(unlist(purrr::map(Tiled_Images, function(df) df[[2]]$Neighborhood_assignment))), collapse = ", ")))
    }

    #save exit function if parallelization fails
    on.exit({
      future::plan("future::sequential")
      gc()
    })

    future::plan("future::multisession", workers = N_cores)
    options(future.globals.maxSize = Inf, future.rng.onMisuse = "ignore")
    furrr::furrr_options(scheduling = Inf)


    #Generate the tile info with the votes for each tile
    RESULTS <-
      furrr::future_map(1:length(Tiled_Images), function(Image) {
        #Generate the count of votes within each tile
        Interim <- Tiled_Images[[Image]][[2]]
        Interim <- Interim %>% dplyr::filter(Neighborhood_assignment %in% Neighborhoods_included) #Select only desired neighborhoods
        Tiles <- Interim %>% dplyr::count(tile_id) %>% dplyr::filter(n >= Minimum_cell_no_per_tile) %>% dplyr::select(tile_id) #Filter tiles without enough voters
        Votes <- Interim %>% dplyr::filter(tile_id %in% Tiles[[1]]) %>% dplyr::group_by(tile_id, Neighborhood_assignment) %>% dplyr::count() %>%
          dplyr::pivot_wider(names_from = Neighborhood_assignment, values_from = n) #Count the votes per tile
        Votes[is.na(Votes)] <- 0
        Votes_prop <- Votes[-1] / apply(Votes[-1], MARGIN = 1, function(x) sum(x)) #Calculate the proportion of votes supporting each neighborhood
        names(Votes_prop) <-stringr::str_c("PROP_", names(Votes_prop))
        Votes$N_votes <- apply(Votes[-1], MARGIN = 1, function(x) sum(x)) #Count the total voter per tile

        Votes <-dplyr::bind_cols(Votes, Votes_prop) #Bind absolute votes information with proportion votes

        Tiles_with_votes <-dplyr::left_join(Tiled_Images[[Image]][[1]], Votes, by = "tile_id") %>% na.omit() #Bind tile info with election results

        #Count vote percentages and find the winner by tile
        Vote_count <- Tiles_with_votes %>% dplyr::select(contains("PROP"))

        #If valid vote count (Vote counts is at least one row in length) proceed with calculating the winner
        if(nrow(Vote_count) > 0) {
          #Generate a function that calculates the winner, 2nd and third candidates in each tile
          Election_results <-purrr::map_dfr(seq_along(1:nrow(Vote_count)), function(Tile_no) {
            Election_results <- Vote_count[Tile_no, ] %>% dplyr::pivot_longer(1:ncol(Vote_count)) %>% dplyr::arrange(desc(value)) %>%
              dplyr::mutate(name = substr(name, start = 6, stop = nchar(name))) %>%
              dplyr::mutate(value = case_when(value == 0 ~ NA,
                                              TRUE ~ value),
                            name = case_when(is.na(value) ~ NA,
                                             TRUE ~ name))
            #Specify the results if absolute winner (only one result returned), super winner or simple winner (more than one result returned)
            if(nrow(Election_results) == 1) {
              Election_results[2,] <- NA
              Election_results[3,] <- NA
            }
            else if(nrow(Election_results == 2)) {
              Election_results[3,] <- NA
            }
            Election_results_tibble <- tibble(Winner = Election_results[[1,1]],
                                              Second = Election_results[[2,1]],
                                              Third = Election_results[[3,1]]) %>%
              dplyr::mutate(Majority_type = case_when(Election_results[[1,2]] == 1 ~ "ABSOLUTE",
                                                      Election_results[[1,2]] >= 0.5 & Election_results[[1,2]] <1 ~ "SUPER",
                                                      Election_results[[1,2]] < 0.5 ~ "SIMPLE"))
          })
        }

        #if invalid then the result is NA
        else{
          Election_results <- tibble(Winner = NA,
                                     Second = NA,
                                     Third = NA,
                                     Majority_type = NA)
        }

        #Combine the results
        return(bind_cols(Tiles_with_votes, Election_results))
      },
      .progress = TRUE)
    future::plan("future::sequential")
    gc()

    names(RESULTS) <- names(Tiled_Images)

    Summary_tibble <-purrr::map_dfr(seq_along(1:length(RESULTS)), function(Image) {
      Subject_Names <- names(RESULTS)[Image]
      N_tiles <- nrow(RESULTS[[Image]])
      Counts <- RESULTS[[Image]] %>% dplyr::count(Winner) %>% dplyr::pivot_wider(names_from = Winner, values_from = n)
      Prop_counts <- Counts / N_tiles
      names(Prop_counts) <- stringr::str_c("PROP_", names(Prop_counts))
      Final_tibble <- tibble(Subject_Names = Subject_Names)
      Final_tibble2 <- tibble(N_tiles = N_tiles)

      dplyr::bind_cols(Final_tibble, Counts, Final_tibble2, Prop_counts)

    })

    Summary_tibble[is.na(Summary_tibble)] <- 0

    PROPORTIONS <- Summary_tibble %>% dplyr::select(dplyr::contains("PROP_"))

    OTHER <- Summary_tibble %>% dplyr::select(-dplyr::contains("PROP_"), -Subject_Names, -N_tiles)

    Winner_summary <- dplyr::bind_cols(Summary_tibble %>% dplyr::select(Subject_Names, N_tiles), OTHER, PROPORTIONS)

    return(list(Images = RESULTS,
                Winner_summary = Winner_summary)
    )
  }
