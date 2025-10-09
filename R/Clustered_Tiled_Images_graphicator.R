#' Generates plots of tiled neighborhood identity
#'
#' The function generates plots of the neighborhood labels of each tile after running the [Tiled_Image_Clustering_function()].
#'
#' @param Tiled_images A list containing tiled images obtained using [Tiled_Image_Clustering_function()].
#' @param Image_name A character value indicating the image to be plotted.
#'
#' @returns A plot of the neighborhood identity of each tile.
#' @export

Clustered_Tiled_Images_graphicator <-
  function(Tiled_images = NULL,
           Image_name = NULL){
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
    return(PLOT)
  }
