#' Generates plots of tiled neighborhood identity
#'
#' The function generates plots of the neighborhood labels of each tile after running the [Tiled_Image_Clustering_function()].
#'
#' @param Tiled_images A list containing tiled images obtained using [Tiled_Image_Clustering_function()].
#' @param Image_name A character value indicating the image to be plotted.
#'
#' @returns A plot of the neighborhood identity of each tile.
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
#' #Graph the results---------------------------------------------------------
#' Clustered_Tiled_Images_graphicator(
#'     Tiled_images = Clustered_Tiled_Images,
#'     Image_name = "ABCM22001_B14_MiniCrop.tif"
#')
#' }
#'
#'
#' @export

Clustered_Tiled_Images_graphicator <-
  function(Tiled_images,
           Image_name){
    #Check arguments
    if(!Image_name %in% names(Tiled_images)) stop(paste0(Image_name, " not found in Tiled_images"))

    Image_to_graph <- Tiled_images[[Image_name]]
    PLOT <- Image_to_graph %>%
      ggplot() +
      geom_rect(aes(group = tile_id, xmin = tile_xmin, ymin = tile_ymin, xmax = tile_xmax, ymax = tile_ymax, fill = Cluster_assignment), color = NA) +
      theme_minimal() +
      scale_fill_manual("Cluster", values = unname(pals::polychrome(length(unique(Image_to_graph$Cluster_assignment))))) +
      theme(panel.grid = element_blank())
    plot(PLOT)
    return(invisible(PLOT))
  }
