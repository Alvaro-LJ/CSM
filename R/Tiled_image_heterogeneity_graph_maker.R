#' Generates a plot illustrating heterogeneity by tile.
#'
#' @param DATA A dataframe or tibble containing cell feature data.
#' @param Tiled_images A list containing tiled images obtained using [Tiled_image_heterogeneity_calculator()].
#' @param Image_name A character value indicating the image to be plotted.
#' @param Metric A character value indicating the heterogeneity metric to be plotted.
#'
#' @returns Returns a plot with the cells and the overlying tiles shaded according to the heterogeneity index.
#'
#' @seealso [Tiled_image_heterogeneity_graph_maker()], [Tiled_image_heterogeneity_analyzer()]
#'
#' @examples
#' \dontrun{
#' #Divide cells into tiles---------
#' Tiled_Images <-
#' Image_tiling_processing_function(
#'    N_cores = 2,
#'    DATA = CSM_Phenotypecell_test,
#'    Tile_width = 125,
#'    Tile_height = 125,
#'    Variables_to_keep = "Phenotype"
#')
#'
#' #Calculate heterogeneity by tile----
#' Tiled_image_heterogeneity_calculator(
#'     Tiled_images = Tiled_Images,
#'     Minimum_cell_no_per_tile = 3,
#'     Phenotypes_included = c("TUMOR", "CD8_GZMBneg", "CD8_GZMBpos", "OTHER")
#')
#'
#' #Graph results----------------------
#' Tiled_image_heterogeneity_graph_maker(
#'     DATA = CSM_Phenotypecell_test,
#'     Tiled_images = Heterogeneity_by_tile,
#'     Image_name = "ABCM22001_B09_MiniCrop.tif",
#'     Metric = "Shannon"
#')
#' }
#'
#'
#' @export

Tiled_image_heterogeneity_graph_maker <-
  function(DATA,
           Tiled_images,
           Image_name,
           Metric){
    #Check arguments
    if(!all(Image_name %in% DATA$Subject_Names, Image_name %in% names(Tiled_images))) stop("Image_name not found in DATA or Tiled_images")
    if(!Metric %in% names(Tiled_images[[1]])) stop("Metric not found in DATA")

    Cells <- DATA %>% dplyr::filter(Subject_Names == Image_name) %>% dplyr::select(1:4, Phenotype)
    Tiles <- Tiled_images[[Image_name]]
    Cells <- Cells %>% dplyr::filter(Phenotype %in% names(Tiles))

    Tiles <- cbind(Tiles[1:7], Tiles[Metric])
    names(Tiles)[8] <- "value"

    PLOT <- Cells  %>%
      ggplot() +
      geom_rect(aes(group = tile_id, xmin = tile_xmin, ymin = tile_ymin, xmax = tile_xmax, ymax = tile_ymax, fill = value), color = "black",
                alpha = 0.5, data = Tiles) +
      geom_point(aes(x = X, y = Y, color = Phenotype), size = 2, alpha = 0.8)+
      cowplot::theme_cowplot()+
      scale_x_continuous("", labels = NULL) +
      scale_y_continuous("", labels = NULL) +
      scale_color_manual("Cell Type", values = unname(pals::polychrome(length(unique(Cells$Phenotype)))))+
      guides(color = guide_legend(override.aes = list(size = 12)))+
      theme(panel.grid = element_blank(),
            axis.line = element_blank(),
            axis.ticks = element_blank(),
            legend.text = element_text(size = 20),
            legend.title = element_text(size = 20)) +
      scale_fill_gradient2(Metric,
                           limits = c(min(Tiles$value), max(Tiles$value)),
                           low = "blue", high = "red", mid = "white",
                           midpoint = quantile(Tiles$value, 0.5))
    plot(PLOT)
    return(invisible(PLOT))
  }
