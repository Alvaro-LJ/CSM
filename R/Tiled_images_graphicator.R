#' Generates plots of cell counts by tile
#'
#' The function generates plots of the number of target cells present within each tile for a given image.
#'
#' @param Tiled_images A list containing tiled images obtained using [Image_tiling_processing_function()].
#' @param Image_name A character value indicating the image to be plotted.
#' @param Phenotypes_included A character vector indicating the phenotype labels that will be included in the plots.
#'
#' @returns A list containing plots (one element per phenotype).
#'
#' @examples
#' \dontrun{
#'#Divide cells into tiles---------
#' Tiled_Images <-
#' Image_tiling_processing_function(
#'    N_cores = 2,
#'    DATA = CSM_Phenotypecell_test,
#'    Tile_width = 125,
#'    Tile_height = 125,
#'    Variables_to_keep = "Phenotype"
#')
#'
#' #Graph target cells
#' Tiled_images_graphicator(
#'    Tiled_images = Tiled_Images,
#'    Image_name = "ABCM22001_B14_MiniCrop.tif",
#'    Phenotypes_included = c("TUMOR", "CD8_GZMBneg", "CD8_GZMBpos", "OTHER")
#')
#' }
#'
#'
#' @export

Tiled_images_graphicator <-
  function(Tiled_images = NULL,
           Image_name = NULL,
           Phenotypes_included = NULL) {
    if(!Image_name %in% names(Tiled_images)) stop("Image_name not found in Tiled_images")

    #Import the image data
    Tiled_Images <- Tiled_images[[Image_name]]

    #Check that the phenotypes included are present in the data
    if(!all(Phenotypes_included %in% unique(Tiled_Images[[2]]$Phenotype))) {
      stop(paste0("Phenotypes included must be any of: ", stringr::str_c(unique(Tiled_Images[[2]]$Phenotype), collapse = ", ")))
    }

    else{
      #Generate a cell count by tile info (we need to change the counting system to account for tiles were there are no cells and dissapear from count)
      Interim <- Tiled_Images[[2]] %>% dplyr::select(tile_id, Phenotype) %>% dplyr::filter(Phenotype %in% Phenotypes_included) %>%
        group_by(tile_id) %>%
        dplyr::count(Phenotype) %>% tidyr::pivot_wider(names_from = Phenotype, values_from = n)

      #Generate a tibble and account for cero values
      Interim <- dplyr::left_join(tibble(tile_id = unique(Tiled_Images[[2]]$tile_id)),
                                  Interim,
                                  by = "tile_id")
      Interim[is.na(Interim)] <- 0

      #Join the cell count tibble with the grid information
      For_graph <- dplyr::left_join(Tiled_Images[[1]], Interim, by = "tile_id")


      #Generate a plot for each of the phenotypes
      Plot_list <-purrr::map(8:ncol(For_graph), function(Variable){
        #Generate a tibble containing tile info and only a single marker
        Graph_tibble <-dplyr::bind_cols(For_graph[1:7], For_graph[Variable])
        names(Graph_tibble)[8] <- "Variable"

        #Generate the individual cell tible
        Individual_cells <- Tiled_Images[[2]] %>% dplyr::filter(Phenotype == names(For_graph[Variable]))

        Plot <- Graph_tibble %>% ggplot() +
          geom_rect(aes(group = tile_id, xmin = tile_xmin, ymin = tile_ymin, xmax = tile_xmax, ymax = tile_ymax, fill = Variable), color = "black",
                    alpha = 0.5) +
          geom_point(aes(x = X, y = Y), alpha = 0.4, size = 1.5, data = Individual_cells) +
          cowplot::theme_cowplot()+
          scale_x_continuous("", labels = NULL) +
          scale_y_continuous("", labels = NULL) +
          scale_fill_viridis_c(names(For_graph)[Variable], na.value = "white") +
          ggtitle(names(For_graph)[Variable])+
          theme(panel.grid = element_blank(),
                axis.line = element_blank(),
                axis.ticks = element_blank(),
                legend.text = element_text(size = 15),
                legend.title = element_text(size = 20),
                plot.title = element_text(size = 25, hjust = 0.5, vjust = -3))
        return(Plot)
      })
      plot(cowplot::plot_grid(plotlist = Plot_list))
      return(invisible(Plot_list))
    }
  }
