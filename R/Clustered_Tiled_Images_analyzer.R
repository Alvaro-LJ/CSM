#' Generates a summary of neighborhood analysis based on clustering of tiles.
#'
#' The function returns neighborhood counts and percentages. Optionally the function may compute neighborhood heterogeneity metrics.
#'
#' @param Tiled_images A list containing neighborhood information generated using [Tiled_Image_Clustering_function()] function.
#' @param Perform_heterogeneity_analysis A logical value indicating if heterogeneity analysis should be performed.
#' @param Graph_Modularity_resolution The resolution parameter for graph modularity analysis.
#'
#' @details
#' Shannon, Simpson,Inverse Simpson and Renyi entropy indexes are calculated using the vegan package.
#'
#' Rao index is calculated using the picante::raoD function.
#'
#' Gini index is calculated using the DescTools::Gini function.
#'
#' Graph modularity is computed using the igraph::modularity function.
#'
#' @returns Generates a tibble containing neighborhood summary by sample and optionally heterogeneity metrics.
#'
#' @examples
#'\dontrun{
#' #Tile images with cell phenotype information---------------------------------
#' Tiled_Images <-
#'  Image_tiling_processing_function(
#'    N_cores = 1,
#'    DATA = CSM_Phenotypecell_test,
#'    Tile_width = 125,
#'    Tile_height = 125,
#'    Variables_to_keep = "Phenotype"
#' )
#'
#' #Cluster cell composition by tile to find neighborhoods---------------------
#' Clustered_Tiled_Images <-
#' Tiled_Image_Clustering_function(
#'     Tiled_images = Tiled_Images,
#'     Minimum_cell_no_per_tile = 4,
#'     Minimum_valid_tiles_per_image = 4,
#'     Phenotypes_included = unique(CSM_Phenotypecell_test$Phenotype),
#'
#'     Cluster_Data = "Cell_Density",
#'
#'     Perform_Dimension_reduction = FALSE,
#'     Cluster_on_Reduced = FALSE,
#'
#'    Strategy = "Consensus_Clustering",
#'    Max_N_Clusters = 5,
#'    Consensus_reps = 2,
#'    Consensus_p_Items = 1,
#'    Consensus_Cluster_Alg = "pam",
#'    Consensus_Distance = "euclidean",
#'    Consensus_Name = "Consensus_clustering_test"
#')
#'
#' #Run the analysis-------------------------------------------------------------
#'Clustered_Tiled_Images_analyzer(
#'    Tiled_images = Clustered_Tiled_Images,
#'    Perform_heterogeneity_analysis = TRUE,
#'    Graph_Modularity_resolution = 0.5
#')
#' }
#'
#' @export




Clustered_Tiled_Images_analyzer <-
  function(Tiled_images,
           Perform_heterogeneity_analysis = FALSE,
           Graph_Modularity_resolution = NULL) {

    #Check suggested packages
    {
      if(Perform_heterogeneity_analysis){
        if(!requireNamespace("picante", quietly = FALSE)) stop(
          paste0("picante CRAN package is required to execute the function. Please install using the following code: ",
                 expression(install.packages("picante")))
        )
        if(!requireNamespace("vegan", quietly = FALSE)) stop(
          paste0("vegan CRAN package is required to execute the function. Please install using the following code: ",
                 expression(install.packages("vegan")))
        )
        if(!requireNamespace("DescTools", quietly = FALSE)) stop(
          paste0("DescTools CRAN package is required to execute the function. Please install using the following code: ",
                 expression(install.packages("DescTools")))
        )
        if(!requireNamespace("igraph", quietly = FALSE)) stop(
          paste0("igraph CRAN package is required to execute the function. Please install using the following code: ",
                 expression(install.packages("igraph")))
        )
      }

    }

    #Check arguments
    if(!all(purrr::map_lgl(Tiled_images, function(x) "Cluster_assignment" %in% names(x)))) stop("Tiled_images must be created using Tiled_Image_Clustering_function")
    if(!is.logical(Perform_heterogeneity_analysis)) stop("Perform_heterogeneity_analysis must be a logical value")
    if(Perform_heterogeneity_analysis){
      if(!all(is.numeric(Graph_Modularity_resolution), Graph_Modularity_resolution > 0)) stop("Graph_Modularity_resolution must be a numeric value > 0")
    }

    #Generate a tibble with the counts
    RESULTS <-purrr::map_dfr(Tiled_images, function(Image){
      Count_DF <- Image %>% dplyr::count(Cluster_assignment) %>% tidyr::pivot_wider(names_from = Cluster_assignment, values_from = n)
      Count_DF$n_tiles <- nrow(Image)
      return(Count_DF)
    }, .progress = list(clear = F,
                        name = "Counting cells in each tile",
                        show_after = 1,
                        type = "iterator"))

    #Substitute na for 0
    RESULTS[is.na(RESULTS)] <- 0

    #Reorder the tibble columns
    RESULTS <- RESULTS[c(which(names(RESULTS) == "n_tiles"), which(names(RESULTS) != "n_tiles"))]

    #Add subject names and reorder the tibble
    RESULTS$Subject_Names <- names(Tiled_images)
    RESULTS <- RESULTS[c(ncol(RESULTS), 1:(ncol(RESULTS)-1))]

    #Calculate the percentage of tiles belonging to each cluster
    Per_DF <-purrr::map_dfc(RESULTS[-c(1:2)], function(Column) Column/RESULTS$n_tiles)
    names(Per_DF) <-stringr::str_c("PER_", names(Per_DF))

    #If modularity is not calculated return the result
    if(!Perform_heterogeneity_analysis){
      return(bind_cols(RESULTS, Per_DF))
    }

    #If modularity needs to be calculated generate graph and calculate the modularity
    else if(Perform_heterogeneity_analysis){
      #We will calculate a graph for each image

      Modularity_scores <-
        purrr::map_dbl(Tiled_images, function(Image){

          #We define the number and ID of edges of the graph
          Graph_tibble <- as_tibble(expand.grid.unique(Image$tile_id, Image$tile_id), .name_repair = "universal_quiet")
          names(Graph_tibble) <- c("from", "to")
          Graph_tibble <- Graph_tibble %>%dplyr::mutate(ID = stringr::str_c(from, to, sep = "_"))

          #We determine the distance between nodes that will be the features of the edges
          DISTANCE_MATRIX <- as_tibble(as.matrix(dist(Image %>% dplyr::select(tile_X_centroid, tile_Y_centroid), method = "euclidean")))
          names(DISTANCE_MATRIX) <- Image$tile_id
          DISTANCE_MATRIX <- DISTANCE_MATRIX %>%dplyr::mutate(from = Image$tile_id)
          DISTANCE_MATRIX <- DISTANCE_MATRIX[c(ncol(DISTANCE_MATRIX), 1:(ncol(DISTANCE_MATRIX)-1))]
          DISTANCE_MATRIX <- DISTANCE_MATRIX %>% tidyr::pivot_longer(-1, names_to = "to", values_to = "weight") %>%
            dplyr::mutate(ID = stringr::str_c(from, to, sep = "_")) %>% dplyr::select(-from, -to)

          #We join both tibbles to generate the final graph
          GRAPH_DF <-dplyr::left_join(Graph_tibble, DISTANCE_MATRIX, by = "ID") %>% dplyr::select(-ID)
          GRAPH_DF <- GRAPH_DF %>% dplyr::mutate(weight = 1/weight)

          #We add the tile clustering information
          NODE_ID <- Image %>% dplyr::select(tile_id, Cluster_assignment)

          #We generate the graph
          Tile_pattern_graph <- igraph::graph_from_data_frame(GRAPH_DF, directed = F, vertices = NODE_ID)

          #We calculate the graph modularity
          igraph::modularity(Tile_pattern_graph, membership = as.factor(NODE_ID$Cluster_assignment), directed = F, resolution = Graph_Modularity_resolution)

        }, .progress = list(clear = F,
                            name = "Calculating graph modularity",
                            show_after = 1,
                            type = "iterator"))
      #Generate a modularity tibble
      Modularity_tibble <- tibble(Subject_Names = names(Tiled_images), Graph_Modularity = Modularity_scores)


      #Generate an interim tibble to calculate other heterogeneity metrics
      Interim <- RESULTS %>% dplyr::select(-n_tiles)
      #Calculate heterogeneity metrics
      Heterogeneity_results <-
        dplyr::bind_cols(Interim[1],
                         purrr::map_dfc(c("shannon" = "shannon", "simpson" = "simpson", "invsimpson" = "invsimpson"),
                                        function(Metric) vegan::diversity(Interim[-1], index = Metric)),
                         tibble(renyi = as.double(vegan::renyi(Interim[-1], scales = Inf))),
                         tibble(Rao_Dkk = picante::raoD(Interim[-1])$Dkk),
                         tibble(Gini = unlist(apply(Interim[-1], MARGIN = 1, function(Sample) DescTools::Gini(Sample, conf.level = NA))))
        )
      #change names
      names(Heterogeneity_results)[-1] <- c("Shannon", "Simpson", "Inverse_simpson", "Renyi_Scale_Inf", "Rao_Dkk", "Gini")

      #Bind both tibbles
      Heterogeneity_results <-dplyr::left_join(Modularity_tibble, Heterogeneity_results, by = "Subject_Names")

      return(
        list(Tile_count_tibble = dplyr::bind_cols(RESULTS, Per_DF),
             Heterogeneity_analysis = Heterogeneity_results)
      )
    }
  }
