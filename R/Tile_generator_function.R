#' Generates tiles from any given image
#'
#' Intended for internal use only

#' @param Image_name A character value indicating the name of the image being processed.
#' @param DATA A dataframe or tibble containing cell feature data.
#' @param Tile_width A numeric value indicating the width of the tiles.
#' @param Tile_height A numeric value indicating the height of the tiles.
#' @param Variables_to_keep A character vector indicating the column names to be kept when tiling the image.
#'
#' @details
#' Used in [Image_tiling_processing_function()]
#'
#'
#' @returns A list containing two elements. Tile_info contains information of tiles. Final_tibble contains cell information and the tile where they are located.
#'
#' @seealso [Suggested_Tile_Size_Calculator()], [Image_tiling_processing_function()]
#'
#' @keywords Internal

Tile_generator_function <-
  function(Image_name,
           DATA = NULL,
           Tile_width = NULL,
           Tile_height = NULL,
           Variables_to_keep = NULL) {
    #Check arguments
    Interim <- DATA %>% dplyr::filter(Subject_Names == Image_name)#Select cells from image

    #Calculate X and Y centroids
    X_centroids <- seq(from = (min(Interim$X) + (Tile_width/2)), to = (max(Interim$X) + (Tile_width/2)), by =  Tile_width)
    Y_centroids <- seq(from = (min(Interim$Y) + (Tile_height/2)), to = (max(Interim$Y) + (Tile_height/2)), by =  Tile_height)

    #Find all possible combinations of X with Y to find all tile centroids
    centroids_tibble <- as_tibble(expand.grid(X_centroids, Y_centroids))
    names(centroids_tibble) <- c("tile_X_centroid", "tile_Y_centroid")
    #Describe centroid tiles and assign a tile ID
    Tile_info <- centroids_tibble %>%
      dplyr::mutate(tile_id =stringr::str_c("tile_", 1:nrow(centroids_tibble))) %>%
      dplyr::mutate(tile_xmin = floor(tile_X_centroid - (Tile_width/2)), tile_xmax = ceiling(tile_X_centroid + (Tile_width/2)),
                    tile_ymin = floor(tile_Y_centroid - (Tile_height/2)), tile_ymax = ceiling(tile_Y_centroid + (Tile_height/2)))

    #Now assign the cells of the data a unique tile
    Interim_cells <- Interim
    Interim_tiles <- Tile_info
    Interim_cells <- Interim_cells %>%dplyr::mutate(tile_id =
                                                      purrr::map2_chr(.x = Interim_cells$X, .y = Interim_cells$Y, function(.x, .y) {
                                                        x_position <- (.x >= Interim_tiles$tile_xmin) & (.x < Interim_tiles$tile_xmax)
                                                        y_position <- (.y >= Interim_tiles$tile_ymin) & (.y < Interim_tiles$tile_ymax)
                                                        Final_pos <- x_position == T & y_position == T
                                                        Tile <- Interim_tiles[Final_pos, 3]
                                                        Tile <- Tile[[1]]
                                                        Tile[[1]]
                                                      })
    )
    Final_tibble <-dplyr::left_join(Interim_cells, Interim_tiles, by = "tile_id")
    Final_tibble <- Final_tibble %>% dplyr::select(1:4, dplyr::all_of(Variables_to_keep), dplyr::contains("tile"))
    Final_tibble <- Final_tibble %>%
      dplyr::mutate(Final_status = (Final_tibble$X >= Final_tibble$tile_xmin) & (Final_tibble$X < Final_tibble$tile_xmax) &
                      (Final_tibble$Y >= Final_tibble$tile_ymin) & (Final_tibble$Y < Final_tibble$tile_ymax))
    list(Tile_info,
         Final_tibble)
  }
