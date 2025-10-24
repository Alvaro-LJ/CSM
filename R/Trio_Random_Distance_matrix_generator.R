#' Calculates random distance backgrounds between a cell of origin and two target cell types
#'
#' Given a dataframe with cell features and phenotype labels, random spatial distances between three cell populations is calculated. This is done by
#' randomly permuting cell labels in an iterative manner. For every image the same number of random COO to target cell distances are calculated.
#'
#'
#' @param N_cores Integer. Number of cores to parallelize your computation.
#' @param DATA A dataframe or tibble containing a column named 'Phenotype' containing cell phenotype labels.
#' @param Cell_Of_Origin A character value indicating the cell phenotype label of the Cell of Origin.
#' @param Target_Cell_1 A character value indicating the cell phenotype label of the first Target cell.
#' @param Target_Cell_2 A character value indicating the cell phenotype label of the second Target cell.
#' @param Random_cells_per_sample A integer value indicating the number of random iterations to calculate the distances.
#' @param Perform_edge_correction A logical value indicating if COO close to the edge of the tissue should be removed from the analysis.
#' @param Hull_ratio If edge correction needs to be performed, a numeric value indicating the hull ratio. Smaller values calculate more precise edge silhouettes at the cost of being more computationally demanding.
#' @param Distance_to_edge If edge correction needs to be performed, the distance the maximum distance to the edge allowed. Cells closer to the edge will be removed.
#'
#' @seealso [Trio_Distance_matrix_generator()], [Trio_Cumulative_Interaction_generator()], [Trio_Min_Distance_analyzer()], [Trio_Cells_in_Radius_analyzer()], [Trio_graph_maker()]
#'
#' @returns A list containing the random distance matrix for each image. Rows represent Cells of Origin and columns represent Target cells.
#'
#' @examples
#' Trio_Random_Distance_matrix_generator(
#'   N_cores = 1,
#'   DATA = CSM_Phenotypecell_test,
#'   Cell_Of_Origin = "TUMOR",
#'   Target_Cell_1 = "CD8_GZMBneg",
#'   Target_Cell_2 = "CD8_GZMBpos",
#'   Random_cells_per_sample = 10,
#'   Perform_edge_correction = FALSE
#')
#'
#'
#' @export

Trio_Random_Distance_matrix_generator <-
  function(N_cores = NULL,
           DATA = NULL,
           Cell_Of_Origin = NULL,
           Target_Cell_1 = NULL,
           Target_Cell_2 = NULL,
           Random_cells_per_sample = NULL,
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
    DATA_Phenotypes <- DATA
    if(!identical(names(DATA_Phenotypes)[1:4],  c("Cell_no", "X", "Y", "Subject_Names"))) { #Check if Data is correctly formatted
      stop("DATA provided should have an adecuate format")
    }
    if(!("Phenotype" %in% names(DATA_Phenotypes))) {
      stop("DATA should contain a column named Phenotype specifying the cell types")
    }
    #Check if provided cell types are in the phenotype variable
    if(!all(c(Cell_Of_Origin %in% unique(DATA_Phenotypes$Phenotype),
              Target_Cell_1 %in% unique(DATA_Phenotypes$Phenotype),
              Target_Cell_2 %in% unique(DATA_Phenotypes$Phenotype)))
    ) stop(paste0("Cell of origin provided and Target cells should be one of: ", stringr::str_c(unique(DATA_Phenotypes$Phenotype), collapse = ", ")))
    if(!all(Random_cells_per_sample%%1 == 0, Random_cells_per_sample > 0)) stop("Random_cells_per_sample must be an integer value > 0")
    if(!all(N_cores >= 1 & N_cores%%1 == 0)) stop("N_cores must be an integer value > 0")
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

    RESULTS <-suppressWarnings({
      furrr::future_map(unique(DATA_Phenotypes$Subject_Names), function(Images){
        #Prepare our data
        Image_tibble <- DATA_Phenotypes %>% dplyr::filter(Subject_Names == Images)
        Pre_Cell_counts <- Image_tibble %>% dplyr::count(Phenotype) %>% dplyr::filter(Phenotype %in% c(Cell_Of_Origin, Target_Cell_1, Target_Cell_2))#Select cells of interest
        #Build a DF where the first row is the COO cell count and the second row is the Target cell count
        Cell_counts <-dplyr::bind_rows(Pre_Cell_counts %>% dplyr::filter(Phenotype == Cell_Of_Origin),
                                       Pre_Cell_counts %>% dplyr::filter(Phenotype == Target_Cell_1),
                                       Pre_Cell_counts %>% dplyr::filter(Phenotype == Target_Cell_2))

        #Acount for samples without any COO or Target Cell 1 or 2
        if(nrow(Cell_counts) < 3) {
          return(
            list(
              Cell_counts = Cell_counts,
              A_to_B = tibble(
                Cell_Of_Origin_no = NA,
                Target_cell_1 = NA
              ),
              A_to_C = tibble(
                Cell_Of_Origin_no = NA,
                Target_cell_1 = NA
              )
            )
          )
        }

        #If border cell correction required generate the border polygon
        if(Perform_edge_correction){
          Cells_sf <- sf::st_as_sf(Image_tibble , coords = c("X", "Y"))
          Edge_line <- sf::st_cast((Cells_sf %>%dplyr::summarise() %>% sf::st_concave_hull(ratio = Hull_ratio) %>% summarise), "LINESTRING")
          Cells_in_Border_vector <- as.vector(unlist(sf::st_is_within_distance(Cells_sf, Edge_line, sparse = F, dist = Distance_to_edge)))
        }

        #If all ok proceed with analysis
        Distance_matrix <-purrr::map(seq_along(1:Random_cells_per_sample), function(y){ #Number of random cell generated by each sample
          I_dont_care <- y
          #First generate new random cell type labels
          new_Phenotype <- sample(Image_tibble$Phenotype, size = nrow(Image_tibble), replace = FALSE)
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
          Distance_matrix <- Advanced_Distance_function_A_B_C_single(DATA = Random_tibble, cell_A = Cell_Of_Origin, cell_B = Target_Cell_1, cell_C = Target_Cell_2) #Select the cell of origin (cell_A) and target cells (B and C)
          return(Distance_matrix)
        })

        #Obtain individual distance matrices
        Cell_A_to_B <-purrr::map_dfr(Distance_matrix, function(x) {
          Interim <- x[[1]]
          names(Interim)[-1] <-stringr::str_c("Random_Target_cell_", as.character(1:(ncol(Interim)-1)))
          return(Interim)
        })
        Cell_A_to_C <-purrr::map_dfr(Distance_matrix, function(x) {
          Interim <- x[[2]]
          names(Interim)[-1] <-stringr::str_c("Random_Target_cell_", as.character(1:(ncol(Interim)-1)))
          return(Interim)
        })

        #Return all three elements
        return(
          list(
            Cell_counts = Cell_counts,
            A_to_B = Cell_A_to_B,
            A_to_C = Cell_A_to_C
          )
        )

      },
      .progress = TRUE)
    })
    future::plan("future::sequential")
    gc()

    #Change the names of the RESULTS list
    names(RESULTS) <- unique(DATA_Phenotypes$Subject_Names)

    #Calculate for which samples the distance matrix could not be calculated
    Samples_to_remove <-purrr::map_lgl(RESULTS, function(x) is.na(x[[2]][[1,1]]))

    #If samples to remove are present print a warning and proceed to remove required samples
    if(sum(Samples_to_remove) > 0) {
      warning(paste0("Samples without COO or target cells will be removed from the analysis. ",
                     "The following samples will be removed: ", stringr::str_c(names(RESULTS)[Samples_to_remove], collapse = ", ")
      )
      )
      return(RESULTS[!Samples_to_remove])
    }

    #If no samples need to be removed return the unmodified RESULTS tibble
    else{
      return(RESULTS)
    }
  }
