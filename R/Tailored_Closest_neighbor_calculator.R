#' Compute the closest neighbors matrix.
#'
#' The function calculates the closest neighbors matrix according to user preferences. The algorithm will visit all the cells in the dataset and compute their neighbors. During this process, the visited cell is considered the Cell of Origin (COO) and the neighbors are the target cells.
#' The closest neighbor matrix can be used to identify neighborhoods using [Neighborhood_discovery_function()].
#' In addition to the matrix, the function generates a correlation plot. The correlation between neighbors composition can be interpreted as attraction repulsion patterns between cell types.
#' Neighbors may be defined by distance (a cell is considered a neighbor if it's within a defined distance D), by number (neighbors are the closest N cells) or both (cells are N closest neighbors as long as they are within distance D).
#'
#'
#' @param N_cores Number of cores to parallelize your computation
#' @param DATA A dataframe or tibble containing cell feature data that and a column named 'Phenotype' containing cell labels.
#' @param Strategy A character value indicating the definition of neighborhoods. Either: "Number", "Distance", "Both".
#' @param N_neighbors If Strategy is Number or Both, an integer value indicating the N closest neighbors to be calculated.
#' @param Include_COO_in_neighborhood A logical value indicating if the Cell of Origin itself should be included in the neighborhood.
#' @param Max_dist_allowed If Strategy is Distance or Both, a numeric value indicating the maximum distance allowed to be considere a neighbor.
#' @param Cell_Of_Origin A character vector indicating which cell phenotype labels will be visited to compute their neighbors.
#' @param Target_Cell A character vector indicating which cell phenotype labels will be considered when calculating the neighbors.
#'
#' @seealso [Neighborhood_discovery_function()], [DATA_neighborhoods_renamer()], [Neighborhood_Quantifier()], [Neighborhood_voting_function()], [Tiled_neighborhoods_graphicator()]
#'
#' @returns A list containing two tibbles. Percentage: A tibble containing the neighbor matrix expressed as percentages. Absolute: A tibble containing the neighbor matrix expressed as cell counts.
#'
#' @examples
#' \dontrun{
#' Tailored_Closest_neighbor_calculator(
#'     N_cores = 1,
#'     DATA = CSM_Phenotypecell_test,
#'     Strategy = "Distance",
#'     Include_COO_in_neighborhood = TRUE,
#'     Max_dist_allowed = 50,
#'     Cell_Of_Origin = c("TUMOR", "CD8_GZMBneg", "CD8_GZMBpos"),
#'     Target_Cell = c("TUMOR", "CD8_GZMBneg", "CD8_GZMBpos", "OTHER")
#')
#' }
#'
#' @export


Tailored_Closest_neighbor_calculator <-
  function(N_cores = 1,
           DATA,
           Strategy,
           N_neighbors = NULL,
           Include_COO_in_neighborhood = TRUE,
           Max_dist_allowed = NULL,
           Cell_Of_Origin,
           Target_Cell
  ){

    #Check suggested packages
    {
      if(!requireNamespace("corrplot", quietly = FALSE)) stop(
        paste0("corrplot CRAN package is required to execute the function. Please install using the following code: ",
               expression(install.packages("corrplot")))
      )
      if(!requireNamespace("rtree", quietly = FALSE)) stop(
        paste0("rtree GitHub package is required to execute the function. Please install using the following code: ",
               expression(remotes::install_github("akoyabio/rtree")))
      )
    }

    DATA <- DATA
    #Check arguments
    if(!all(N_cores >= 1 & N_cores%%1 == 0)) stop("N_cores must be an integer value > 0")
    if(!identical(names(DATA)[c(1:4)], c("Cell_no", "X", "Y", "Subject_Names"))) stop("DATA must be formatted adequately")
    if(!"Phenotype" %in% names(DATA)) stop("DATA must contain a column named 'Phenotype")
    if(!all(c(Cell_Of_Origin, Target_Cell) %in% unique(DATA$Phenotype))) stop(
      paste0("Cell_Of_Origin and Target cell must be one of the following: ", stringr::str_c(unique(DATA$Phenotype), collapse = ", "))
    )
    if(!Strategy %in% c("Number", "Distance", "Both")) stop("Strategy must be any of the following: Number, Distance, Both")
    if(Strategy == "Number" || Strategy == "Both"){
      if(!all(N_neighbors%%1 == 0, N_neighbors >= 2)) stop("N_neighbors must be an integer value > 1")
    }
    if(Strategy == "Distance" || Strategy == "Both"){
      if(!all(is.numeric(Max_dist_allowed), Max_dist_allowed > 0)) stop("Max_dist_allowed must be a numeric value > 0")
    }
    if(!is.logical(Include_COO_in_neighborhood)) stop("Include_COO_in_neighborhood must be a logical value")

    #Import data to the function
    Function_DATA_Phenotypes <- DATA %>% dplyr::select(1:4, Phenotype)

    #Filter samples that dont have any cell of origio or any target cell
    Samples_to_keep <- purrr::map_lgl(unique(Function_DATA_Phenotypes$Subject_Names), function(Image){
      Interim <- Function_DATA_Phenotypes %>% dplyr::filter(Subject_Names == Image)

      COO_present <- sum(Interim$Phenotype %in% Cell_Of_Origin) >= 1
      Target_Cell_present <- sum(Interim$Phenotype %in% Target_Cell) >= 1

      COO_present & Target_Cell_present
    })

    #If any invalid sample is present remove it from analysis and print results
    if(sum(!Samples_to_keep) > 0){
      warning(paste0("Samples without COO or without Target cells will be excluded. The following samples will be eliminated: ",
                     stringr::str_c(unique(Function_DATA_Phenotypes$Subject_Names)[!Samples_to_keep], collapse = ", ")))
      Function_DATA_Phenotypes <- Function_DATA_Phenotypes %>% dplyr::filter(Subject_Names %in% unique(Function_DATA_Phenotypes$Subject_Names)[Samples_to_keep])
    }

    #Compute the neighbors for each cell in each image
    RESULTS <-
      purrr::map_dfr(unique(Function_DATA_Phenotypes$Subject_Names), function(Image) {
        #Select data from each unique image
        Interim <- Function_DATA_Phenotypes %>% dplyr::filter(Subject_Names == Image)

        #Generate specific tibbles for cells of origin and Target cells
        Tibble_COO <- Interim %>% dplyr::filter(Phenotype %in% Cell_Of_Origin)
        Tibble_Targets <- Interim %>% dplyr::filter(Phenotype %in% Target_Cell)

        #Modify these tibbles according to the requirements of the function
        COO <- cbind(Tibble_COO[[2]], Tibble_COO[[3]])
        Targets <- rtree::RTree(cbind(Tibble_Targets[[2]], Tibble_Targets[[3]]))

        #Calculate the index of the closest neighbors in the Targets according to strategy
        if(Strategy == "Number" | Strategy == "Both"){
          if(Include_COO_in_neighborhood) Index <- rtree::knn(Targets, COO, as.integer(N_neighbors))
          if(!Include_COO_in_neighborhood) Index <- rtree::knn(Targets, COO, as.integer(N_neighbors+1))
        }
        if(Strategy == "Distance"){
          Index <- rtree::withinDistance(Targets, COO, Max_dist_allowed)
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
        #Generate a DF that includes for each cell the N closest neighbors
        closest_neighbors_by_cell <-
          furrr::future_map_dfr(seq_along(1:length(Index)), function(row){
            cell_X <- Tibble_COO[[row,2]]
            cell_Y <- Tibble_COO[[row,3]]
            Cell_no_ID <- Tibble_COO[[row,1]]

            RESULT <- Tibble_Targets[Index[[row]],]
            #If COO is excluded, eliminate accordingly
            if(!Include_COO_in_neighborhood) RESULT <- RESULT %>% dplyr::filter(!Cell_no %in% Cell_no_ID)
            RESULT <- RESULT %>%dplyr::mutate(DIST = sqrt((X - cell_X)^2 + (Y - cell_Y)^2)) %>% dplyr::arrange(DIST)

            #If strategy includes distance and number filter out the neighbors that are not in a allowed distance range
            if(Strategy == "Both"){
              RESULT <- RESULT %>% dplyr::filter(DIST <= Max_dist_allowed)
            }
            #Count the neighbors (relative to all neighbors) and distance metrics
            Neighbors_in_window_result <- RESULT %>% dplyr::count(Phenotype) %>% tidyr::pivot_wider(names_from = Phenotype, values_from = n)
            FINAL <-dplyr::bind_cols(Tibble_COO[row,], Neighbors_in_window_result/nrow(RESULT))
            FINAL <- FINAL %>%dplyr::mutate(N_neighbors = nrow(RESULT),
                                            min_DIST = min(RESULT$DIST),
                                            max_DIST = max(RESULT$DIST),
                                            avg_DIST = mean(RESULT$DIST),
                                            median_DIST = quantile(RESULT$DIST, 0.5))
            return(FINAL)
          },
          .progress = TRUE)
        future::plan("future::sequential")
        gc()
        #Replace NA values by 0
        closest_neighbors_by_cell[is.na(closest_neighbors_by_cell)] <- 0

        #If any element of the Tibble target is missing in the neighbors then add 0 manually
        if(!all(unique(Tibble_Targets$Phenotype) %in% names(closest_neighbors_by_cell))){
          #Find missing targets
          names_var <- unique(Tibble_Targets$Phenotype)[!unique(Tibble_Targets$Phenotype) %in% names(closest_neighbors_by_cell)]
          #Generate a 0 tibble and change names
          Absent_tibble <- as_tibble(matrix(0, nrow = nrow(closest_neighbors_by_cell), ncol = length(names_var)))
          names(Absent_tibble) <- names_var

          closest_neighbors_by_cell <-dplyr::bind_cols(closest_neighbors_by_cell, Absent_tibble)
        }
        return(closest_neighbors_by_cell)
      },
      .progress = list(clear = FALSE,
                       name = "Calculating closest neighbors",
                       show_after = 1,
                       type = "iterator"))

    #Change NA values to 0
    RESULTS[is.na(RESULTS)] <- 0

    #Modify column order
    RESULTS <- RESULTS[c("Cell_no","X","Y","Subject_Names","Phenotype",
                         names(RESULTS)[-which(names(RESULTS) %in%
                                                 c("Cell_no","X","Y","Subject_Names","Phenotype",
                                                   "N_neighbors", "min_DIST","max_DIST","avg_DIST","median_DIST"))],
                         "N_neighbors", "min_DIST","max_DIST","avg_DIST","median_DIST")]

    #Plot correlation matrix and distances
    For_correlation <- RESULTS[names(RESULTS)[-which(names(RESULTS) %in%
                                                       c("Cell_no","X","Y","Subject_Names","Phenotype",
                                                         "N_neighbors", "min_DIST","max_DIST","avg_DIST","median_DIST"))]]
    #We need to remove columns without variability
    Drop_vars <-purrr::map_lgl(For_correlation, ~length(unique(.)) == 1)
    if(sum(Drop_vars >= 1)) message(paste0("The following phenotypes appear at a constant rate in neighbors and will be dropped from the correlation analysis: ",
                                           stringr::str_c(names(Drop_vars)[Drop_vars], collapse = ", ")))
    #Select only variables that are keen for correlation
    For_correlation <- For_correlation[!Drop_vars]
    cor_DATA <- cor(For_correlation, method = "pearson")
    corrplot::corrplot(cor_DATA, method = "shade", type = "lower", order = "hclust", addCoef.col = "black", number.cex = 0.8, tl.cex = 0.8,
                       tl.pos = "lt", tl.col = "black")

    plot(RESULTS %>% dplyr::select(max_DIST, avg_DIST, median_DIST) %>% tidyr::pivot_longer(1:3) %>%
           ggplot(aes(x = value)) + facet_wrap(~name, ncol = 1, nrow = 3, "free") + geom_histogram(binwidth = 2)+
           cowplot::theme_cowplot()+
           scale_x_continuous("Distance"))

    RESULTS_Proportion <- RESULTS
    RESULTS_Absolute <- RESULTS

    RESULTS_Absolute[c(names(RESULTS)[-which(names(RESULTS) %in%
                                               c("Cell_no","X","Y","Subject_Names","Phenotype",
                                                 "N_neighbors", "min_DIST","max_DIST","avg_DIST","median_DIST"))])] <-
      round(RESULTS_Absolute[c(names(RESULTS)[-which(names(RESULTS) %in%
                                                       c("Cell_no","X","Y","Subject_Names","Phenotype",
                                                         "N_neighbors", "min_DIST","max_DIST","avg_DIST","median_DIST"))])]*RESULTS$N_neighbors,
            digits = 0)

    #Return the final data
    return(list(
      Percentage = RESULTS_Proportion,
      Absolute_count = RESULTS_Absolute
    ))
  }
