#' Calculates distance between a single cell of origin and two target cells
#'
#' The function calculates the Euclidean distance between three cell phenotypes. One is considered the Cell of Origin (COO) and two others are considered the target cells.
#' The COO and one of the target cells can be the same cell type.
#'
#' @param N_cores Integer. Number of cores to parallelize your computation.
#' @param DATA A dataframe or tibble containing a column named 'Phenotype' containing cell phenotype labels.
#' @param Cell_Of_Origin A character value indicating the cell phenotype label of the Cell of Origin.
#' @param Target_Cell_1 A character value indicating the cell phenotype label of the first Target cell.
#' @param Target_Cell_2 A character value indicating the cell phenotype label of the second Target cell.
#' @param Perform_edge_correction A logical value indicating if COO close to the edge of the tissue should be removed from the analysis.
#' @param Hull_ratio If edge correction needs to be performed, a numeric value indicating the hull ratio. Smaller values calculate more precise edge silhouettes at the cost of being more computationally demanding.
#' @param Distance_to_edge If edge correction needs to be performed, the distance the maximum distance to the edge allowed. Cells closer to the edge will be removed.
#'
#' @seealso [Trio_Random_Distance_matrix_generator()], [Trio_Cumulative_Interaction_generator()], [Trio_Min_Distance_analyzer()], [Trio_Cells_in_Radius_analyzer()], [Trio_graph_maker()]
#'
#' @returns A list containing the distance matrix for each image. Rows represent Cells of Origin and columns represent Target cells.
#'
#' @examples
#' \dontrun{
#' Trio_Distance_matrix_generator(
#'     N_cores = 1,
#'     DATA = CSM_Phenotypecell_test,
#'     Cell_Of_Origin = "TUMOR",
#'     Target_Cell_1 = "CD8_GZMBneg",
#'     Target_Cell_2 = "CD8_GZMBpos",
#'     Perform_edge_correction = FALSE
#' )
#' }
#'
#' @import dplyr
#'
#' @export

Trio_Distance_matrix_generator <-
  function(N_cores = 1,
           DATA,
           Cell_Of_Origin,
           Target_Cell_1,
           Target_Cell_2,
           Perform_edge_correction = FALSE,
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

    DATA_Phenotypes <- DATA

    if(!identical(names(DATA_Phenotypes)[1:4],  c("Cell_no", "X", "Y", "Subject_Names"))) { #Check if Data is correctly formatted
      stop("DATA provided should have an adecuate format")
    }
    if(!("Phenotype" %in% names(DATA_Phenotypes))) {
      stop("DATA should contain a column named Phenotype specifying the cell types")
    }
    if(!all(c(Cell_Of_Origin %in% unique(DATA_Phenotypes$Phenotype),
              Target_Cell_1 %in% unique(DATA_Phenotypes$Phenotype),
              Target_Cell_2 %in% unique(DATA_Phenotypes$Phenotype)
    )
    )) { #Check if provided cell types are in the phenotype variable
      stop(paste0("Cell of origin provided and Target cells should be one of: ", stringr::str_c(unique(DATA_Phenotypes$Phenotype), collapse = ", ")))
    }
    if(!all(N_cores >= 1 & N_cores%%1 == 0)) stop("N_cores must be an integer value > 0")
    if(Perform_edge_correction){
      if(!all(is.numeric(Hull_ratio), Hull_ratio >= 0, Hull_ratio <= 1)) stop("Hull_ratio must be a numeric value between 0 and 1")
      if(!all(is.numeric(Distance_to_edge), Distance_to_edge > 0)) stop("Distance_to_edge must be a numeric value > 0")
    }

    #Perform a random test to allow user to stop the computation if edge correction parameter is not desired
    if(Perform_edge_correction){
      print("Running edge correction example on a random sample")
      Sample <- DATA_Phenotypes %>% dplyr::filter(Subject_Names == sample(unique(DATA_Phenotypes$Subject_Names), size = 1))
      Cells_sf <- sf::st_as_sf(Sample , coords = c("X", "Y"))
      Edge_line <- sf::st_cast((Cells_sf %>% summarise() %>% sf::st_concave_hull(ratio = Hull_ratio) %>% summarise), "LINESTRING")
      Cells_in_Border_vector <- unlist(sf::st_is_within_distance(Cells_sf, Edge_line, sparse = FALSE, dist = Distance_to_edge))

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

      #Remove cells
      Keep_vector_list <-purrr::map(unique(DATA_Phenotypes$Subject_Names), function(x){
        #Prepare our data
        Image_tibble <- DATA_Phenotypes %>% dplyr::filter(Subject_Names == x)
        Cells_sf <- sf::st_as_sf(Image_tibble , coords = c("X", "Y"))
        Edge_line <- sf::st_cast((Cells_sf %>% summarise() %>% sf::st_concave_hull(ratio = Hull_ratio) %>% summarise), "LINESTRING")
        Cells_in_Border_vector <- unlist(sf::st_is_within_distance(Cells_sf, Edge_line, sparse = F, dist = Distance_to_edge))
        #Calculate cells in border
        COO_in_Border_vector <- Image_tibble$Phenotype == Cell_Of_Origin & Cells_in_Border_vector

        #Print message to warn COO removed in analysis
        message(paste0("Sample ", as.character(x), ": ", sum(COO_in_Border_vector), " / ", sum(Image_tibble$Phenotype == Cell_Of_Origin), " ", Cell_Of_Origin, " cell/s will be removed due to edge proximity."))

        #Return a vector with the cells to keep (either no COO or COO not in edge)
        return(!COO_in_Border_vector)
      })
      names(Keep_vector_list) <- unique(DATA_Phenotypes$Subject_Names)
    }

    #save exit function if parallelization fails
    on.exit({
      future::plan("future::sequential")
      gc()
    })

    #Proceed with computation
    future::plan("future::multisession", workers = N_cores)
    options(future.globals.maxSize = Inf, future.rng.onMisuse = "ignore")
    furrr::furrr_options(scheduling = Inf)
    RESULTS <-
      furrr::future_map(unique(DATA_Phenotypes$Subject_Names),  function(x) {
        #Prepare our data
        Image_tibble <-
          DATA_Phenotypes %>% dplyr::filter(Subject_Names == x)
        #First we will filter cells of origin in the border as desired by the user
        if(Perform_edge_correction){
          Cells_to_keep <- Keep_vector_list[[x]]
          Image_tibble <- Image_tibble[Cells_to_keep,]
        }

        Pre_Cell_counts <-
          Image_tibble %>% dplyr::count(Phenotype) %>% dplyr::filter(Phenotype %in% c(Cell_Of_Origin, Target_Cell_1, Target_Cell_2))#Select cells of interest
        #Build a DF where the first row is the COO cell count and the second row is the Target cell count
        Cell_counts <-
          dplyr::bind_rows(
            Pre_Cell_counts %>% dplyr::filter(Phenotype == Cell_Of_Origin),
            Pre_Cell_counts %>% dplyr::filter(Phenotype == Target_Cell_1),
            Pre_Cell_counts %>% dplyr::filter(Phenotype == Target_Cell_2)
          )

        #Acount for samples without any COO or Target Cell 1 or 2
        if (nrow(Cell_counts) < 3) {
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

        #If everything is OK proceed to calculate the distance matrices
        Distance_matrices <- Advanced_Distance_function_A_B_C(DATA = Image_tibble, cell_A = Cell_Of_Origin, cell_B = Target_Cell_1, cell_C = Target_Cell_2)
        return(
          list(Cell_counts = Cell_counts,
               A_to_B = Distance_matrices[[1]],
               A_to_C = Distance_matrices[[2]])
        )
      },
      .progress = TRUE)
    future::plan("future::sequential")
    gc()
    #Add names
    names(RESULTS) <- unique(DATA_Phenotypes$Subject_Names)

    #Calculate which samples have NA values and should be removed
    Samples_to_remove <-purrr::map_lgl(RESULTS, function(x) is.na(x[[2]][[1,1]]))

    #If present print a warning and remove the samples
    if(sum(Samples_to_remove) > 0) {
      warning(paste0("Samples without COO or target cells will be removed from the analysis. ",
                     "The following samples will be removed: ", stringr::str_c(names(RESULTS)[Samples_to_remove], collapse = ", ")
      )
      )
      return(RESULTS[!Samples_to_remove])
    }

    #Else return the resulsts straight forward
    else{
      return(RESULTS)
    }
  }
