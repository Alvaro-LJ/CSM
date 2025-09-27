#' Removes cells in edge of tissue
#'
#' Removes cells in the edge of tissue calculating the concave hull. Cells in the edge of tissue can sometimes be subject of technical artifacts and can be worth removing them before analysis.
#'
#' @param N_cores Integer. Number of cores to parallelize your computation.
#' @param DATA A dataframe or tibble containing cell feature data.
#' @param Hull_ratio A numeric value indicating the hull ratio. Smaller values calculate more precise edge silhouettes at the cost of being more computationally demanding.
#' @param Distance_to_edge A numeric value indicating the maximum distance to tissue edge allowed
#' @param Image_preview A character value indicating which Subject_Names will be used in the quick test. If NULL a random image will be selected.
#' @returns Returns a tibble with cell features without cells closer to the edge than Distance_to_edge.
#'
#' @export

Cell_in_edge_remover <-
  function(N_cores = NULL,
           DATA = NULL,
           Hull_ratio = NULL,
           Distance_to_edge = NULL,
           Image_preview = NULL
  ) {
    #Check arguments
    #Obtain the data
    DATA <- DATA

    if(!identical(names(DATA)[1:4],  c("Cell_no", "X", "Y", "Subject_Names"))) { #Check if Data is correctly formatted
      stop("DATA provided should have an adecuate format")
    }
    if(!all(N_cores >= 1 & N_cores%%1 == 0)) stop("N_cores must be an integer value > 0")
    if(!all(is.numeric(Hull_ratio), Hull_ratio >= 0, Hull_ratio <= 1)) stop("Hull_ratio must be a numeric value between 0 and 1")
    if(!all(is.numeric(Distance_to_edge), Distance_to_edge > 0)) stop("Distance_to_edge must be a numeric value > 0")
    if(!any(is.null(Image_preview), Image_preview %in% unique(DATA$Subject_Names))) stop(paste0(Image_preview, " not found in Subject_Names"))


    #Perform a random test to allow user to stop the computation if edge correction parameter is not desired
    if(is.null(Image_preview)){
      print("Running edge correction example on a random sample")
      Image_preview <- sample(unique(DATA$Subject_Names), size = 1)
    }
    Sample <- DATA %>% dplyr::filter(Subject_Names == Image_preview)
    Cells_sf <- sf::st_as_sf(Sample , coords = c("X", "Y"))
    Edge_line <- sf::st_cast((Cells_sf %>% dplyr::summarise() %>% sf::st_concave_hull(ratio = Hull_ratio) %>% dplyr::summarise), "LINESTRING")
    Cells_in_Border_vector <- unlist(sf::st_is_within_distance(Cells_sf, Edge_line, sparse = F, dist = Distance_to_edge))

    plot(Sample %>%
           dplyr::mutate(Removed = Cells_in_Border_vector) %>%
           ggplot(aes(x = X, y = Y, color = Cells_in_Border_vector)) +
           geom_point() +
           scale_color_manual("", labels = c("Included", "Removed"), values = c("black", "grey")) +
           theme_minimal() +
           scale_x_continuous("") +
           scale_y_continuous("") +
           theme(panel.grid = element_blank(),
                 axis.text = element_blank(),
                 legend.position = "bottom",
                 legend.text = element_text(size = 12)))

    #Ask the user if the algorihtm should proceed
    answer <- menu(c("Proceed", "Abort"), title = "Should the analysis proceed")
    #If user decides to stop then abort function and return stop message
    if(answer == 2) stop("The function has been stopped. Please tune edge correction parameters for a better result")


    #Remove cells in border
    #save exit function if parallelization fails
    on.exit({
      future::plan("future::sequential")
      gc()
    })

    #Now we calculate our distance matrix
    future::plan("future::multisession", workers = N_cores)
    options(future.globals.maxSize = Inf, future.rng.onMisuse = "ignore")
    furrr::furrr_options(scheduling = Inf)

    RESULTS <-
      furrr::future_map_dfr(unique(DATA$Subject_Names), function(x){
        #Prepare our data
        Image_tibble <- DATA %>% dplyr::filter(Subject_Names == x)
        Cells_sf <- sf::st_as_sf(Image_tibble , coords = c("X", "Y"))
        Edge_line <- sf::st_cast((Cells_sf %>% dplyr::summarise() %>% sf::st_concave_hull(ratio = Hull_ratio) %>% dplyr::summarise), "LINESTRING")
        Cells_in_Border_vector <- unlist(sf::st_is_within_distance(Cells_sf, Edge_line, sparse = F, dist = Distance_to_edge))



        #Print message to warn COO removed in analysis
        message(paste0("Sample ", as.character(x), ": ", sum(Cells_in_Border_vector), " / ", nrow(Image_tibble), " cell/s will be removed due to edge proximity."))

        #Return the Tibble with the cells that are not in the border
        return(Image_tibble[!Cells_in_Border_vector,])
      }, .progress = TRUE)
    future::plan("future::sequential")
    gc()

    return(RESULTS)
  }
