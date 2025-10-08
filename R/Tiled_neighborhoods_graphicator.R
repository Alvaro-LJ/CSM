#' Generates a plot illustrating heterogeneity by tile.
#'
#' @param DATA_elections An object created using the [Neighborhood_voting_function()] function containing neighborhood election results.
#' @param Image_name A character value indicating the image to be plotted.
#' @param Graph_only_winner_neighborhood A logical value indicating if only the winner neighborhood by tile should be plotted.
#'
#' @returns Returns a plot with the tiles and the neighborhood elections results.
#'
#' @export

Tiled_neighborhoods_graphicator <-
  function(DATA_elections = NULL,
           Image_name = NULL,
           Graph_only_winner_neighborhood = NULL) {
    #Check arguments
    if(!is.logical(Graph_only_winner_neighborhood)) stop("Graph_only_winner_neighborhood must be a logical value")
    if(!Image_name %in% names(DATA_elections$Images)) stop("Image_name not present in DATA_elections")

    #Select individual image to graph
    Image_to_graph <- DATA_elections$Images[[Image_name]]

    #Graph all available neighborhoods
    if(!Graph_only_winner_neighborhood) {

      PLOT <- Image_to_graph %>%
        dplyr::select(dplyr::contains("tile_"), dplyr::contains("PROP")) %>%
        pivot_longer(dplyr::contains("PROP")) %>%
        ggplot() + geom_rect(aes(group = tile_id, xmin = tile_xmin, ymin = tile_ymin, xmax = tile_xmax, ymax = tile_ymax, fill = value), color = "black") +
        facet_wrap(~name) + theme_minimal() + theme(panel.grid = element_blank()) +
        scale_fill_gradient(low = "grey", high = "red")
      plot(PLOT)
      return(PLOT)
    }

    else if(Graph_only_winner_neighborhood) {
      PLOT <- Image_to_graph %>%
        ggplot() +
        geom_rect(aes(group = tile_id, xmin = tile_xmin, ymin = tile_ymin, xmax = tile_xmax, ymax = tile_ymax, fill = Winner), color = NA) +
        theme_minimal() +
        scale_fill_manual("Winner neighborhood", values = unname(pals::polychrome(length(unique(Image_to_graph$Winner))))) +
        theme(panel.grid = element_blank())

      plot(PLOT)
      return(PLOT)
    }

  }
