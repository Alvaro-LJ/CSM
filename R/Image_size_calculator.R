#' Approximates tissue size
#'
#' The function approximates tissue size either by tiling the tissue and analyzing the amount of tiles that present cells or by calculating the silouette of the tissue by calculating the concave hull.
#' The area can then be used by other functions to calculate densities like [Phenotype_quantifier()] or [Neighborhood_Quantifier()].
#'
#' @param DATA A dataframe or tibble containing cell feature data.
#' @param Strategy Either 'Tiling' or 'Concave_hull'
#' @param Image_to_plot (Optional) A character value indicating the image to be used in the preview. If NULL the smallest tissue will be used.
#' @param Tile_accuracy Lower values calculate the area in a more precise manner, with higher computational times. If the value is too low, it can classify areas of stroma as not being tissue.
#' @param Hull_ratio A numeric value indicating the hull ratio. Smaller values calculate more precise edge silhouettes at the cost of being more computationally demanding.
#' @returns Returns a tibble with tissue size for each image
#'
#' @export

Image_size_calculator <-
  function(DATA = NULL,
           Strategy = NULL,
           Image_to_plot = NULL,

           Tile_accuracy = NULL,

           Hull_ratio = NULL

  ) {

    #Check suggested package
    if(Strategy == "Convave_hull"){
      if(!requireNamespace("sf", quietly = FALSE)) stop(
        paste0("sf CRAN package is required to execute the function. Please install using the following code: ",
               expression(install.packages("sf")))
      )
    }

    DATA <- DATA
    #Check arguments
    if(!all(c("Subject_Names", "X", "Y") %in% names(DATA)[c(1:4)])) {
      stop("Your Data must contain the following variables: Subject_Names, X, Y. Please format using the Data_arrange_function.")
    }
    if(!Strategy %in% c("Tiling", "Concave_hull")) stop("Strategy must be one of the following: Tiling, Concave_hull")
    if(!any(Image_to_plot %in% unique(DATA$Subject_Names), is.null(Image_to_plot))) stop("Image_to_plot is not present in DATA")
    if(Strategy == "Tiling"){
      if(!all(is.numeric(Tile_accuracy), Tile_accuracy > 0)) stop("Tile_accuracy must be a positive numeric value")
    }
    if(Strategy == "Concave_hull"){
      if(!all(is.numeric(Hull_ratio), Hull_ratio >= 0, Hull_ratio <= 1)) stop("Hull_ratio must be a numeric value between 0 and 1")
    }

    #Tiling strategy
    if(Strategy == "Tiling"){
      Result <- tibble(Subject_Names = unique(DATA$Subject_Names),
                       Area =purrr::map_dbl(unique(DATA$Subject_Names), function(Image) {
                         Interim <- DATA %>% dplyr::filter(Subject_Names == Image) %>%
                           ggplot(aes(x = X, y = Y)) + geom_bin2d(binwidth = Tile_accuracy)
                         N_tiles <- nrow(as_tibble(layer_data(Interim))) #obtain the layer data
                         N_tiles * Tile_accuracy * Tile_accuracy #multiply the tile number by the tile area
                       }, .progress = list(clear = F,
                                           name = "Calculating tissue area",
                                           show_after = 2,
                                           type = "iterator"))
      )

      #If image to plot is NULL select the smallest image
      if(is.null(Image_to_plot)){
        Image_to_plot <- (Result %>% dplyr::arrange(Area))[[1,1]]
      }

      DENSITY_plot <- DATA %>%  dplyr::filter(Subject_Names == Image_to_plot) %>% #plot the smallest sample
        ggplot(aes(x = X, y = Y)) + geom_bin2d(binwidth = Tile_accuracy) +
        cowplot::theme_cowplot() + guides(fill = "none")

      CELL_plot <- DATA %>%  dplyr::filter(Subject_Names == Image_to_plot) %>%
        ggplot(aes(x = X, y = Y)) + geom_point(size = 2) +
        cowplot::theme_cowplot()
      plot(
        patchwork::wrap_plots(CELL_plot,  DENSITY_plot, nrow = 1)
      )
      return(Result)
    }

    #Concave_hull
    if(Strategy == "Concave_hull"){
      Result <- tibble(Subject_Names = unique(DATA$Subject_Names),
                       Area =purrr::map_dbl(unique(DATA$Subject_Names), function(Image){
                         Interim <- DATA %>% dplyr::filter(Subject_Names == Image)
                         Cells_sf <- sf::st_as_sf(Interim, coords = c("X", "Y"))
                         Sample_polygon <- sf::st_cast((Cells_sf %>% dplyr::summarise() %>% sf::st_concave_hull(ratio = Hull_ratio) %>% summarise), "POLYGON")
                         Area <- sf::st_area(Sample_polygon)
                         return(Area)
                       }, .progress = list(clear = F,
                                           name = "Calculating tissue area",
                                           show_after = 2,
                                           type = "iterator")))

      #If image to plot is NULL select the smallest image
      if(is.null(Image_to_plot)){
        Image_to_plot <- (Result %>% dplyr::arrange(Area))[[1,1]]
      }

      #Generate the sample to plot and the plots
      Sample_to_plot <- DATA %>%  dplyr::filter(Subject_Names == Image_to_plot)
      Sample_sf <- sf::st_as_sf(Sample_to_plot, coords = c("X", "Y"))
      Sample_polygon <- sf::st_cast((Sample_sf %>% dplyr::summarise() %>% sf::st_concave_hull(ratio = Hull_ratio) %>% summarise), "POLYGON")
      LAYOUT_plot <-  Sample_polygon %>% ggplot() + geom_sf(linewidth = 1.5, fill = "white", color = "black") +
        cowplot::theme_cowplot() + guides(fill = "none") +
        scale_x_continuous("") +
        scale_y_continuous("") +
        theme(panel.grid = element_blank(),
              axis.text = element_blank())
      CELL_plot <- DATA %>%  dplyr::filter(Subject_Names == Image_to_plot) %>%
        ggplot(aes(x = X, y = Y)) + geom_point(size = 2) +
        cowplot::theme_cowplot()
      plot(
        patchwork::wrap_plots(CELL_plot,  LAYOUT_plot, nrow = 1)
      )
      #Return the final result
      return(Result)
    }
  }
