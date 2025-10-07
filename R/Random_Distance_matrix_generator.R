#' Calculates random distance backgrounds between two cell types
#'
#' Given a dataframe with cell features and phenotype labels, random spatial distances between two cell populations is calculated. This is done by
#' randomly permuting cell labels in an iterative manner. For every image the same number of random COO to target cell distances are calculated.
#' The COO and target cell can be the same cell type.
#'
#' @param N_cores Integer. Number of cores to parallelize your computation.
#' @param DATA A dataframe or tibble containing a column named 'Phenotype' containing cell phenotype labels.
#' @param Cell_Of_Origin A character value indicating the cell phenotype label of the Cell of Origin.
#' @param Target_Cell A character value indicating the cell phenotype label of the Target cell.
#' @param Random_cells_per_sample A integer value indicating the number of random iterations to calculate the distances.
#' @param Allow_Cero_Distance A logical value indicating if zero distances should be allowed. This may be relevant when COO and Target cell labels are the same.
#' @param Perform_edge_correction A logical value indicating if COO close to the edge of the tissue should be removed from the analysis.
#' @param Hull_ratio If edge correction needs to be performed, a numeric value indicating the hull ratio. Smaller values calculate more precise edge silhouettes at the cost of being more computationally demanding.
#' @param Distance_to_edge If edge correction needs to be performed, the distance the maximum distance to the edge allowed. Cells closer to the edge will be removed.
#'
#' @seealso [Distance_matrix_generator()], [Cumulative_Interaction_generator()], [Distance_analyzer()], [Cells_in_Radius_analyzer()], [Cell_to_Cell_graph_maker()]
#'
#' @returns A list containing the random distance matrix for each image. Rows represent Cells of Origin and columns represent Target cells.
#'
#' @export

Random_Distance_matrix_generator <-
  function(N_cores = NULL,
           DATA = NULL,
           Cell_Of_Origin = NULL,
           Target_Cell = NULL,
           Random_cells_per_sample = NULL,
           Allow_Cero_Distance = NULL,
           Perform_edge_correction = NULL,
           Hull_ratio = NULL,
           Distance_to_edge = NULL
  ) {
    #Check suggested packages
    if(Perform_edge_correction){
      if(!requireNamespace("sf", quietly = FALSE)) stop(
        paste0("sf CRAN package is required to execute the function. Please install using the following code: ",
               expression(install.packages("sf")))
      )
    }

    #Check arguments
    if(!exists("Advanced_Distance_function_A_B_single", envir = .GlobalEnv)) stop("Please execute STEP 5 required functions")
    Advanced_Distance_function_A_B_single <- get("Advanced_Distance_function_A_B_single",  envir = .GlobalEnv)
    if(!exists(DATA, envir = .GlobalEnv)){
      stop("A DATA_Distance object must be created before running the tiling analysis. Check name supplied to the DATA argument")
    }
    DATA_Phenotypes <- get(DATA, envir = .GlobalEnv)
    if(!identical(names(DATA_Phenotypes)[1:4],  c("Cell_no", "X", "Y", "Subject_Names"))) { #Check if Data is correctly formatted
      stop("DATA provided should have an adecuate format")
    }
    if(!("Phenotype" %in% names(DATA_Phenotypes))) {
      stop("DATA should contain a column named Phenotype specifying the cell types")
    }
    #Check if provided cell types are in the phenotype variable
    if(!all(c(Cell_Of_Origin %in% unique(DATA_Phenotypes$Phenotype),
              Target_Cell %in% unique(DATA_Phenotypes$Phenotype))
    )) {
      stop(paste0("Cell of origin provided and Target cells should be one of: ", stringr::str_c(unique(DATA_Phenotypes$Phenotype), collapse = ", ")))
    }
    if(!all(Random_cells_per_sample%%1 == 0, Random_cells_per_sample > 0)) stop("Random_cells_per_sample must be an integer value > 0")
    if(!all(N_cores >= 1 & N_cores%%1 == 0)) stop("N_cores must be an integer value > 0")
    if(!is.logical(Allow_Cero_Distance)) stop("Allow_Cero_Distance should be a logical value")
    if(!is.logical(Perform_edge_correction)) stop("Perform_edge_correction must be a logical value")
    if(Perform_edge_correction){
      if(!all(is.numeric(Hull_ratio), Hull_ratio >= 0, Hull_ratio <= 1)) stop("Hull_ratio must be a numeric value between 0 and 1")
      if(!all(is.numeric(Distance_to_edge), Distance_to_edge > 0)) stop("Distance_to_edge must be a numeric value > 0")
    }

    #save exit function if parallelization fails
    on.exit({
      future::plan("future::sequential")
      gc()
    })

    #Perform Randomization of labels and calculate distances
    future::plan("future::multisession", workers = N_cores)
    options(future.globals.maxSize = Inf, future.rng.onMisuse = "ignore")
    furrr::furrr_options(scheduling = Inf)

    RESULTS <-
      furrr::future_map(unique(DATA_Phenotypes$Subject_Names), function(Images){
        #Prepare our data
        Image_tibble <- DATA_Phenotypes %>% dplyr::filter(Subject_Names == Images)
        Pre_Cell_counts <- Image_tibble %>% dplyr::count(Phenotype) %>% dplyr::filter(Phenotype %in% c(Cell_Of_Origin, Target_Cell))#Select cells of interest
        #Build a DF where the first row is the COO cell count and the second row is the Target cell count
        Cell_counts <-dplyr::bind_rows(Pre_Cell_counts %>% dplyr::filter(Phenotype == Cell_Of_Origin),
                                       Pre_Cell_counts %>% dplyr::filter(Phenotype == Target_Cell))

        #Acount for samples without any COO or Target Cell
        if(nrow(Cell_counts) < 2) {
          return(
            list(Cell_counts = Cell_counts,
                 Distance_matrix = tibble(Cell_Of_Origin_no = NA,
                                          Target_cell_1 = NA))
          )
        }

        #Acount for samples with a unique COO that is also the target cell (infrequent but may occur)
        if((Cell_counts[[1,1]] == Cell_counts[[2,1]]) & Cell_counts[[1,2]] == 1){
          return(
            list(Cell_counts = Cell_counts,
                 Distance_matrix = tibble(Cell_Of_Origin_no = NA,
                                          Target_cell_1 = NA))
          )
        }

        #If border cell correction required generate the border polygon
        if(Perform_edge_correction){
          Cells_sf <- sf::st_as_sf(Image_tibble , coords = c("X", "Y"))
          Edge_line <- sf::st_cast((Cells_sf %>%dplyr::summarise() %>% sf::st_concave_hull(ratio = Hull_ratio) %>% summarise), "LINESTRING")
          Cells_in_Border_vector <- unlist(sf::st_is_within_distance(Cells_sf, Edge_line, sparse = F, dist = Distance_to_edge))
        }

        #If all ok proceed with analysis
        #Number of random cell generated by each sample
        Distance_matrix <-purrr::map_dfr(seq_along(1:Random_cells_per_sample), function(y){
          I_dont_care <- y
          #First generate new random cell type labels
          new_Phenotype <- sample(Image_tibble$Phenotype, size = nrow(Image_tibble), replace = F)
          Random_tibble <- Image_tibble %>%dplyr::mutate(Phenotype = new_Phenotype)

          #If needed remove COO that are close to border
          if(Perform_edge_correction){
            Random_tibble <- Random_tibble[!(Random_tibble$Phenotype == Cell_Of_Origin & Cells_in_Border_vector), ]

            #If no random cells are remaining then resample the labels once again
            while(nrow(Random_tibble %>% dplyr::filter(Phenotype == Cell_Of_Origin)) == 0){
              new_Phenotype <- sample(Image_tibble$Phenotype, size = nrow(Image_tibble), replace = F)
              Random_tibble <- Image_tibble %>%dplyr::mutate(Phenotype = new_Phenotype)
              Random_tibble <- Random_tibble[!(Random_tibble$Phenotype == Cell_Of_Origin & Cells_in_Border_vector), ]
            }
          }

          #Then apply distance function
          Distance_matrix <- Advanced_Distance_function_A_B_single(DATA = Random_tibble, cell_A = Cell_Of_Origin, cell_B = Target_Cell)#Define COO and target cell
          names(Distance_matrix)[-1] <-stringr::str_c("Random_Target_cell_", as.character(1:(ncol(Distance_matrix)-1)))
          Distance_matrix
        })
        Distance_matrix <- Distance_matrix %>%dplyr::mutate(Cell_A_no =stringr::str_c("Random_COO_cell_", as.character(1:nrow(Distance_matrix))))
        names(Distance_matrix)[1] <- "Cell_Of_Origin_no"

        if(!Allow_Cero_Distance){
          #If distance is 0, substitute it for NA
          Distance_matrix[Distance_matrix == 0] <- NA
        }
        #Return the final list
        return(
          list(Cell_counts = Cell_counts,
               Distance_matrix = Distance_matrix)
        )
      },
      .progress = TRUE)
    future::plan("future::sequential")
    gc()

    names(RESULTS) <- unique(DATA_Phenotypes$Subject_Names)

    Samples_to_remove <-purrr::map_lgl(RESULTS, function(x) is.na(x[[2]][[1,1]]))

    #If samples to remove are present print a warning and proceed to remove required samples
    if(sum(Samples_to_remove) > 0) {
      warning(paste0("Samples without COO or target cells will be removed from the analysis. ",
                     "The following samples will be removed: ", stringr::str_c(names(RESULTS)[Samples_to_remove], collapse = ", ")
      )
      )
      return(RESULTS[!Samples_to_remove])
    }

    else{
      return(RESULTS)
    }
  }
