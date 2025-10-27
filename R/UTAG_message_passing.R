#' Performs message passing between neighboring cells
#'
#' The function performs message passing. Neighboring cells share their features according to user preferences. The resulting tibble can be used to feed the [UTAG_Neighborhood_identifier()] function.
#'
#' @param DATA_name A tibble or dataframe object containing cell features.
#' @param COO_to_visit (OPTIONAL) A logical vector indicating which cells from DATA should be visited by the algorithm.
#' @param Neighbor_strategy A character value indicating how neighboring cells are defined: "Number", "Distance" or "Both".
#' @param Message_strategy A character value indicating how features are shared between neighboring cells: "Averaging" or "Sum".
#' @param N_neighbors If strategy is Number or Both, an integer indicating the number of closest neighbors to calculate for each cell.
#' @param Max_dist_allowed If strategy is Distance or Both, a numeric value indicating the distance from every cell to search for neighbors
#' @param Weighting_Strategy A character value indicating the weighting strategy for message passing: "None", "Proximity", "Disregarded_minority", "Both" (see details).
#' @param COO_weight (OPTIONAL) A numeric value between 0 and 1 stating the relative weight of the COO in the final result. If NULL, the weight will be the same as for neighbors.
#' @param N_cores Integer. Number of cores to parallelize your computation.
#'
#' @details
#' Weighting strategy modifies the way neighbors share information:
#'
#' Proximity weighting allows closer neighbors to have a higher influence compared to distant neighbors. For a set of neighbors each located at a distance d from the cell of origin, and being the sum of all neighbor distances D, weight is defined as 1/(d/D).
#'
#' Disregarded_minority weighting allows low expressed features to have higher influence on message passing: First the -log10 of average feature expression is calculated for every feature. The weight is defined as -log10 average / min(-log10 average)
#'
#' Proximity weighting is available for averaging and sum message passing strategies. Disregarded_minority can only be used with averaging.
#'
#' The function is based on the original UTAG pipeline described in the following reference: https://doi.org/10.1038/s41592-022-01657-2 and https://github.com/ElementoLab/utag.
#'
#' @seealso [UTAG_Neighborhood_identifier()]
#'
#' @examples
#' \dontrun{
#' UTAG_message_passing(
#'     DATA = CSM_Arrangedcellfeaturedata_test,
#'
#'     Neighbor_strategy = "Distance",
#'     Message_strategy = "Sum",
#'     Max_dist_allowed = 50,
#'     Weighting_Strategy = "Proximity",
#'
#'     N_cores = 1
#')
#' }
#'
#' @returns A tibble containing cell feature data that have undergone neighbor message passing.

UTAG_message_passing <-
  function(DATA_name = NULL,
           COO_to_visit = NULL,

           Neighbor_strategy = NULL,
           Message_strategy = NULL,


           N_neighbors = NULL,
           Max_dist_allowed = NULL,

           Weighting_Strategy = NULL,
           COO_weight = NULL,

           N_cores = NULL) {
    #Check suggested packages
    {
      if(!requireNamespace("rtree", quietly = FALSE)) stop(
        paste0("rtree GitHub package is required to execute the function. Please install using the following code: ",
               expression(remotes::install_github("akoyabio/rtree")))
      )
      if(!requireNamespace("Matrix", quietly = FALSE)) stop(
        paste0("Matrix CRAN package is required to execute the function. Please install using the following code: ",
               expression(install.packages("Matrix")))
      )
    }


    #What to do on exit
    on.exit({
      future::plan("sequential")
      gc()})

    #Import data
    DATA <- DATA_name
    #Check data structure
    if(!identical(names(DATA)[c(1:4)], c("Cell_no", "X", "Y", "Subject_Names"))) stop("DATA must be formatted adequately")
    #Check other arguments
    if(!any(is.null(COO_to_visit), all(is.logical(COO_to_visit), length(COO_to_visit) == nrow(DATA)))) stop(str_c("COO_to_visit must be NULL or a logical vector of length ", nrow(DATA), collapse = ""))
    if(!Neighbor_strategy %in% c("Number", "Distance", "Both")) stop("Neighbor_strategy must be one of the following: Number, Distance, Both")
    if(!Message_strategy %in% c("Averaging", "Sum")) stop("Message_strategy must be one of the following: Averaging, Sum")
    if(Neighbor_strategy %in% c("Number", "Both")){
      if(!all(is.numeric(N_neighbors), N_neighbors%%1 == 0, N_neighbors > 0)) stop("N_neighbors must be a positive integer value > 0")
    }
    if(Neighbor_strategy %in% c("Distance", "Both")){
      if(!all(is.numeric(Max_dist_allowed), Max_dist_allowed > 0)) stop("Max_dist_allowed must be a numeric value > 0")
    }
    if(Message_strategy == "Averaging"){
      if(!(Weighting_Strategy %in% c("None", "Proximity", "Disregarded_minority", "Both"))) {
        stop("Weighting_Strategy must be one of None, Proximity, Disregarded_minority or Both")
      }
      if(!any(is.null(COO_weight), all(is.numeric(COO_weight), COO_weight > 0, COO_weight < 1))) stop("COO_weight must be NULL or a numeric value between 0 and 1")
    }
    if(Message_strategy == "Sum"){
      if(!(Weighting_Strategy %in% c("None", "Proximity"))) {
        stop("Weighting_Strategy must be one of None or Proximity for Sum strategy")
      }
      if(!is.null(COO_weight)) message("COO_weight is not supported for Sum strategy. Argument will be ignored")
    }
    if(!all(!is.null(N_cores), N_cores >= 1 & N_cores%%1 == 0)) stop("N_cores must be an integer value > 0")

    #Generate a vector of the images that are going to be visited (used in the function)
    Images_to_visit <- unique(DATA$Subject_Names)

    #Generate a list of COO and Targets cells to be included in the analysis
    if(!is.null(COO_to_visit)){
      COO_info <- tibble(Subject_Names = DATA$Subject_Names,
                         COO_Included = COO_to_visit)
      COO_info <-purrr::map(unique(DATA$Subject_Names), function(Image){
        Vector <- COO_info %>% dplyr::filter(Subject_Names == Image) %>% dplyr::select(COO_Included)
        Vector[[1]]
      })
      names(COO_info) <- unique(DATA$Subject_Names)

      #If any image does not have any COO to visit, then remove from analysis
      Images_to_visit <- names(COO_info)[purrr::map_lgl(COO_info, ~sum(.) >= 1)]
    }


    #Add 1 to N_neighbors
    N_neighbors <- N_neighbors + 1

    #Perform message passing
    RESULTS <-
      purrr::map_dfr(Images_to_visit, function(Image){
        #First select cells from the image to be analyzed
        Interim <- DATA %>% dplyr::filter(Subject_Names == Image)

        #Calculate some initial weights if Disregarded_minority needs to be used
        if(Weighting_Strategy %in% c("Disregarded_minority", "Both")){
          Image_mean_marker_intensity <- purrr::map_dbl(Interim[-c(1:4)], function(Marker) -log(mean(Marker), base = 10))
          Image_mean_marker_intensity_weight <- Image_mean_marker_intensity / min(Image_mean_marker_intensity)
        }

        #Calculate the targets
        Tibble_Targets <- cbind(Interim[[2]],Interim[[3]])
        Tibble_Targets <- rtree::RTree(as.matrix(Tibble_Targets))
        #Calculate the COO
        Tibble_COO <- cbind(Interim[[2]], Interim[[3]])
        #If user has selected COO_to visit, filter the Tibble_COO object
        if(!is.null(COO_to_visit)){
          COO_to_visit <- COO_info[[Image]]
          Tibble_COO <- Tibble_COO[COO_to_visit, ]
        }
        Tibble_COO <- as.matrix(Tibble_COO)

        #Calculate the closest neighbors according to the instructions given
        if(Neighbor_strategy == "Number" | Neighbor_strategy == "Both"){
          Index <- rtree::knn(Tibble_Targets, Tibble_COO, as.integer(N_neighbors))
        }
        if(Neighbor_strategy == "Distance"){
          Index <- rtree::withinDistance(Tibble_Targets, Tibble_COO, Max_dist_allowed)
        }

        #save exit function if parallelization fails
        on.exit({
          future::plan("future::sequential")
          gc()
        })
        #Make the cluster
        future::plan("future::multisession", workers = N_cores)
        options(future.globals.maxSize = Inf, future.rng.onMisuse = "ignore")
        furrr::furrr_options(scheduling = Inf)
        #Calculate message passed cells
        Message_passed_cells <-
          furrr::future_map_dfr(1:length(Index), function(row) {
            #If no cells have been removed from the message passing then proceed as usual
            if(is.null(COO_to_visit)){
              #Extract Target cells coordinates
              cell_X <- Interim[[row,2]]
              cell_Y <- Interim[[row,3]]
            }

            #If cells have been removed from the COO ist then proceed in a special manner
            if(!is.null(COO_to_visit)){
              Interim_COO <- Interim[COO_to_visit, ]
              cell_X <- Interim_COO[[row,2]]
              cell_Y <- Interim_COO[[row,3]]
            }

            #Obtain all neighbors (including COO)
            Neighbors_in_window <- Interim[Index[[row]],]

            #Calculate distances
            Neighbors_in_window <- Neighbors_in_window %>%dplyr::mutate(DIST = sqrt((X - cell_X)^2 + (Y - cell_Y)^2)) %>% dplyr::arrange(DIST)

            #Remove neighbors accounting for distance if required
            if(Neighbor_strategy == "Both"){
              Neighbors_in_window <- Neighbors_in_window %>% dplyr::filter(DIST <= Max_dist_allowed)
            }

            #Split by Cell Of Origin and neighbors
            COO <- Neighbors_in_window[1,]
            Neighbors <- Neighbors_in_window[-1,]

            #If no neighbors are identified return the individual cell
            if(nrow(Neighbors) == 0){
              Final_cell <- COO %>% dplyr::select(-DIST) %>%dplyr::mutate(mean_DIST = NA, max_DIST = NA, N_neighbors = 0) #Cells without neighbors are not influenced
              return(Final_cell)
            }

            #If neighbors are present then proceed with message passing with the user selected weighting strategy
            #AVERAGING
            if(Message_strategy == "Averaging"){
              if(Weighting_Strategy == "None"){
                #Assign a FIXED weight to all neighbors (INDEPENDENT OF DISTANCE)
                Neighbors <- Neighbors %>%dplyr::mutate(Neighboring_weight = 1)
                Neighbors$Neighboring_weight[is.infinite(Neighbors$Neighboring_weight)] <- max(Neighbors$Neighboring_weight[!is.infinite(Neighbors$Neighboring_weight)])

                #Calculate weighted mean for each marker according to cell proximity
                Averaged_markers <-purrr::map_dbl(Neighbors[-c(1:4, (ncol(Neighbors)-1), ncol(Neighbors))], function(x) stats::weighted.mean(x, w = Neighbors$Neighboring_weight))

                #Perform message passing
                #If COO_weight is null then perform the classic algorithm (all cells are weighted equally)
                if(is.null(COO_weight)){
                  COO[-c(1:4, ncol(COO))] <-purrr::map2_dfc(.x = unlist(COO[-c(1:4, ncol(COO))]), .y = Averaged_markers, function(.x, .y) stats::weighted.mean(c(.x, .y), w = c(1, nrow(Neighbors))))
                }
                #If not perform the alternative with user supplied COO_weight
                else{
                  COO[-c(1:4, ncol(COO))] <-purrr::map2_dfc(.x = unlist(COO[-c(1:4, ncol(COO))]), .y = Averaged_markers, function(.x, .y) stats::weighted.mean(c(.x, .y), w = c(COO_weight, 1-COO_weight)))
                }
                Final_cell <- COO %>% dplyr::select(-DIST) %>%dplyr::mutate(mean_DIST = mean(Neighbors$DIST), max_DIST = max(Neighbors$DIST), N_neighbors = nrow(Neighbors))
                return(Final_cell)
              }
              else if(Weighting_Strategy == "Proximity"){
                #Assign a weight by neighbor proximity
                Neighbors <- Neighbors %>%dplyr::mutate(Neighboring_weight =  1/ (DIST/sum(Neighbors$DIST)))
                Neighbors$Neighboring_weight[is.infinite(Neighbors$Neighboring_weight)] <- max(Neighbors$Neighboring_weight[!is.infinite(Neighbors$Neighboring_weight)])

                #Calculate weighted mean for each marker according to cell proximity
                Averaged_markers <-purrr::map_dbl(Neighbors[-c(1:4, (ncol(Neighbors)-1), ncol(Neighbors))], function(x) stats::weighted.mean(x, w = Neighbors$Neighboring_weight))

                #Perform message passing
                #If COO_weight is null then perform the classic algorithm (all cells are weighted equally)
                if(is.null(COO_weight)){
                  COO[-c(1:4, ncol(COO))] <-purrr::map2_dfc(.x = unlist(COO[-c(1:4, ncol(COO))]), .y = Averaged_markers, function(.x, .y) stats::weighted.mean(c(.x, .y), w = c(1, nrow(Neighbors))))
                }
                #If not perform the alternative with user supplied COO_weight
                else{
                  COO[-c(1:4, ncol(COO))] <-purrr::map2_dfc(.x = unlist(COO[-c(1:4, ncol(COO))]), .y = Averaged_markers, function(.x, .y) stats::weighted.mean(c(.x, .y), w = c(COO_weight, 1-COO_weight)))
                }
                Final_cell <- COO %>% dplyr::select(-DIST) %>%dplyr::mutate(mean_DIST = mean(Neighbors$DIST), max_DIST = max(Neighbors$DIST), N_neighbors = nrow(Neighbors))
                return(Final_cell)
              }
              else if(Weighting_Strategy == "Disregarded_minority"){
                #Assign a FIXED weight to all neighbors (INDEPENDENT OF DISTANCE)
                Neighbors <- Neighbors %>%dplyr::mutate(Neighboring_weight = 1)
                Neighbors$Neighboring_weight[is.infinite(Neighbors$Neighboring_weight)] <- max(Neighbors$Neighboring_weight[!is.infinite(Neighbors$Neighboring_weight)])

                #Calculate weighted mean for each marker according to cell proximity
                Averaged_markers <-purrr::map_dbl(Neighbors[-c(1:4, (ncol(Neighbors)-1), ncol(Neighbors))], function(x) stats::weighted.mean(x, w = Neighbors$Neighboring_weight))

                #Weight results according to Image marker abundance
                Averaged_markers_weighted <- Averaged_markers * Image_mean_marker_intensity_weight

                #Perform message passing
                #If COO_weight is null then perform the classic algorithm (all cells are weighted equally)
                if(is.null(COO_weight)){
                  COO[-c(1:4, ncol(COO))] <-purrr::map2_dfc(.x = unlist(COO[-c(1:4, ncol(COO))]), .y = Averaged_markers, function(.x, .y) stats::weighted.mean(c(.x, .y), w = c(1, nrow(Neighbors))))
                }
                #If not perform the alternative with user supplied COO_weight
                else{
                  COO[-c(1:4, ncol(COO))] <-purrr::map2_dfc(.x = unlist(COO[-c(1:4, ncol(COO))]), .y = Averaged_markers, function(.x, .y) stats::weighted.mean(c(.x, .y), w = c(COO_weight, 1-COO_weight)))
                }
                Final_cell <- COO %>% dplyr::select(-DIST) %>%dplyr::mutate(mean_DIST = mean(Neighbors$DIST), max_DIST = max(Neighbors$DIST), N_neighbors = nrow(Neighbors))
                return(Final_cell)
              }
              else if(Weighting_Strategy == "Both"){
                #Assign a weight by neighbor proximity
                Neighbors <- Neighbors %>%dplyr::mutate(Neighboring_weight =  1/ (DIST/sum(Neighbors$DIST)))
                Neighbors$Neighboring_weight[is.infinite(Neighbors$Neighboring_weight)] <- max(Neighbors$Neighboring_weight[!is.infinite(Neighbors$Neighboring_weight)])

                #Calculate weighted mean for each marker according to cell proximity
                Averaged_markers <-purrr::map_dbl(Neighbors[-c(1:4, (ncol(Neighbors)-1), ncol(Neighbors))], function(x) stats::weighted.mean(x, w = Neighbors$Neighboring_weight))

                #Weight results according to Image marker abundance
                Averaged_markers_weighted <- Averaged_markers * Image_mean_marker_intensity_weight

                #Perform message passing adjusted by weight
                #If COO_weight is null then perform the classic algorithm (all cells are weighted equally)
                if(is.null(COO_weight)){
                  COO[-c(1:4, ncol(COO))] <-purrr::map2_dfc(.x = unlist(COO[-c(1:4, ncol(COO))]), .y = Averaged_markers, function(.x, .y) stats::weighted.mean(c(.x, .y), w = c(1, nrow(Neighbors))))
                }
                #If not perform the alternative with user supplied COO_weight
                else{
                  COO[-c(1:4, ncol(COO))] <-purrr::map2_dfc(.x = unlist(COO[-c(1:4, ncol(COO))]), .y = Averaged_markers, function(.x, .y) stats::weighted.mean(c(.x, .y), w = c(COO_weight, 1-COO_weight)))
                }
                Final_cell <- COO %>% dplyr::select(-DIST) %>%dplyr::mutate(mean_DIST = mean(Neighbors$DIST), max_DIST = max(Neighbors$DIST), N_neighbors = nrow(Neighbors))
                return(Final_cell)
              }
            }
            #SUM
            if(Message_strategy == "Sum"){
              #If no correction then perform bare sum
              if(Weighting_Strategy == "None"){
                #Add neighbor markers
                Added_markers <-purrr::map_dbl(Neighbors[-c(1:4, ncol(Neighbors))], sum)
                #Add neighbor markers to COO
                COO[-c(1:4, ncol(COO))] <-purrr::map2_dfc(.x = unlist(COO[-c(1:4, ncol(COO))]), .y = Added_markers, sum)
                #Generate the final cell
                Final_cell <- COO %>% dplyr::select(-DIST) %>%dplyr::mutate(mean_DIST = mean(Neighbors$DIST), max_DIST = max(Neighbors$DIST), N_neighbors = nrow(Neighbors))
                return(Final_cell)
              }

              if(Weighting_Strategy == "Proximity"){
                #First get our neighbors distance and calculate weights
                Distance_weight <- Neighbors$DIST / sum(Neighbors$DIST)
                Distance_weight <- (1/Distance_weight) / sum(1/Distance_weight)

                #Then we sum all the values in our matrix to obtain the total credit to be delivered
                Expression_matrix <- Neighbors %>% dplyr::select(-c(1:4), -DIST) %>% Matrix::as.matrix()
                Total_credit <- sum(Expression_matrix)
                Credit_by_Neighbor <- Total_credit*Distance_weight

                #Now we calculate the relative credit each cell will assign to each marker
                Relative_value <- sweep(Expression_matrix, MARGIN = 1, STATS = Matrix::rowSums(Expression_matrix), FUN = "/", check.margin = FALSE)
                Final_expression <- sweep(Relative_value, MARGIN = 1, STATS = Credit_by_Neighbor, FUN = "*", check.margin = FALSE)

                #Calculate final neighbors
                Neighbors_Weighted <-dplyr::bind_cols(Neighbors[c(1:4)], Final_expression, Neighbors["DIST"])
                #Add neighbor markers
                Added_markers <-purrr::map_dbl(Neighbors_Weighted[-c(1:4, ncol(Neighbors_Weighted))], sum)
                #Add neighbor markers to COO
                COO[-c(1:4, ncol(COO))] <-purrr::map2_dfc(.x = unlist(COO[-c(1:4, ncol(COO))]), .y = Added_markers, sum)
                #Generate the final cell
                Final_cell <- COO %>% dplyr::select(-DIST) %>%dplyr::mutate(mean_DIST = mean(Neighbors$DIST), max_DIST = max(Neighbors$DIST), N_neighbors = nrow(Neighbors))
                return(Final_cell)
              }
            }
          },
          .progress = TRUE)
        future::plan("future::sequential")
        gc()
        return(Message_passed_cells)
      }, .progress = list(clear = F,
                          name = "Performing Message Passing",
                          show_after = 1,
                          type = "iterator")
      )

    #Turn results into a single tibble and return
    return(RESULTS)
  }
