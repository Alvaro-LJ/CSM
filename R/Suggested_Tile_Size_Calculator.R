#' Calculates appropriate tile size to satisfy user preferences
#'
#' The function can be used  to plan the image tiling strategy. The function calculates the tile dimensions required to divide an image into a user provided number of rows and columns.
#' The calculation can be performed using the smallest or largest image from the dataset.
#'
#' @param DATA A dataframe or tibble containing cell feature data.
#' @param N_rows An integer value indicating the number of rows to divide into tiles.
#' @param N_cols An integer value indicating the number of columns to divide into tiles.
#' @param Based_on_smaller A logical value indicating if the test should be run using the smallest image. If FALSE the largest image in the dataset will be used.
#' @param Draw_preview A logical value indicating if the test should be plotted.
#'
#' @returns A tibble containing recommended tile dimensions.
#'
#' @seealso [Image_length_calculator()], [Image_tiling_processing_function()]
#'
#' @examples
#' \dontrun{
#' Suggested_Tile_Size_Calculator(
#'    CSM_Phenotypecell_test,
#'    N_cols = 4,
#'    N_rows = 4,
#'    Based_on_smaller = TRUE,
#'    Draw_preview = TRUE
#' )
#' }
#'
#' @export

Suggested_Tile_Size_Calculator <-
  function(DATA = NULL,
           N_rows = NULL,
           N_cols = NULL,
           Based_on_smaller = NULL,
           Draw_preview = NULL){
    #Check arguments
    if(!all(c("Subject_Names", "X", "Y") %in% names(DATA))) {
      stop("Data must contain columns named Subject_Names, X and Y")
    }

    else if(!all(c(N_rows%%1 == 0, N_cols%%1 == 0))) {
      stop("N_rows and N_cols arguments must be integer values")
    }

    else if(!is.logical(Based_on_smaller)) {
      stop("Based_on_smaller should be a logical argument, either TRUE or FALSE")
    }

    else if(!is.logical(Draw_preview)) {
      stop("Draw_preview should be a logical argument, either TRUE or FALSE")
    }

    #Define what to do if smaller image is to be used
    if(Based_on_smaller) {
      Tile_width <- round(
        min(
          purrr::map_dbl(unique(DATA$Subject_Names), function(Image){
            Interim <- DATA %>% dplyr::filter(Subject_Names == Image)
            max(Interim$X)-min(Interim$X)
          })
        )/N_cols,
        digits = 0
      )

      Tile_height <-round(
        min(
          purrr::map_dbl(unique(DATA$Subject_Names), function(Image){
            Interim <- DATA %>% dplyr::filter(Subject_Names == Image)
            max(Interim$Y)-min(Interim$Y)
          })
        )/N_rows,
        digits = 0
      )

      Suggestion_tibble <- tibble('Suggested Tile dimension' = c("Width", "Height", "Squared_tiles"),
                                  value = c(Tile_width, Tile_height, mean(c(Tile_width, Tile_height))))


      #Draw a preview if required
      if(Draw_preview) {
        #Calculate the area of all images
        Area_tibble <-purrr::map_dfr(unique(DATA$Subject_Names), function(Image){
          Interim <- DATA %>% dplyr::filter(Subject_Names == Image)
          Area <- (max(Interim$X) - min(Interim$X)) * (max(Interim$Y) - min(Interim$Y))
          tibble(Subject_Names = Image,
                 Area = Area)
        })
        Area_tibble <- Area_tibble %>% dplyr::arrange(Area)

        Interim <- DATA %>% dplyr::filter(Subject_Names == Area_tibble[[1,1]])
        Tile_width <- Suggestion_tibble[[3,2]]

        Tiled_plot <- Interim %>% ggplot(aes(x = X, y = Y)) + geom_bin2d(fill = "white", color = "black", linewidth = 1.1, binwidth = Tile_width) +
          geom_point(size = 2, alpha = 0.5) + cowplot::theme_cowplot()

        Summary_smaller <- dplyr::as_tibble(ggplot2::layer_data(Tiled_plot)) %>% summarize(N_tiles = length(count),
                                                                                           Min_cells = min(count),
                                                                                           p25_cells = quantile(count, 0.25),
                                                                                           Average_cells = mean(count),
                                                                                           p50_cells = quantile(count, 0.5),
                                                                                           p75_cells = quantile(count, 0.75),
                                                                                           Max_cells = max(count))
        #Plot the tiled version of our data
        plot(
          Tiled_plot
        )
        #Generate the final tibble
        list_final <- list(Tile_dimension = Suggestion_tibble,
                           Image = Summary_smaller)
        names(list_final)[2] <- Area_tibble[[1,1]]

        #Return the final tibble
        return(list_final)
      }

      #Return the final results
      return(Suggestion_tibble)
    }

    if(!Based_on_smaller) {
      Tile_width <- round(
        max(
          purrr::map_dbl(unique(DATA$Subject_Names), function(Image){
            Interim <- DATA %>% dplyr::filter(Subject_Names == Image)
            max(Interim$X)-min(Interim$X)
          })
        )/N_cols,
        digits = 0
      )

      Tile_height <-round(
        max(
          purrr::map_dbl(unique(DATA$Subject_Names), function(Image){
            Interim <- DATA %>% dplyr::filter(Subject_Names == Image)
            max(Interim$Y)-min(Interim$Y)
          })
        )/N_rows,
        digits = 0
      )

      Suggestion_tibble <- tibble('Suggested Tile dimension' = c("Width", "Height", "Squared_tiles"),
                                  value = c(Tile_width, Tile_height, mean(c(Tile_width, Tile_height))))

      #Draw a preview if required
      if(Draw_preview) {
        Area_tibble <-purrr::map_dfr(unique(DATA$Subject_Names), function(Image){
          Interim <- DATA %>% dplyr::filter(Subject_Names == Image)
          Area <- (max(Interim$X) - min(Interim$X)) * (max(Interim$Y) - min(Interim$Y))
          tibble(Subject_Names = Image,
                 Area = Area)
        })
        Area_tibble <- Area_tibble %>% dplyr::arrange(dplyr::desc(Area))

        #Get the interim
        Interim <- DATA %>% dplyr::filter(Subject_Names == Area_tibble[[1,1]])
        Tile_width <- Suggestion_tibble[[3,2]]

        Tiled_plot <- Interim %>% ggplot(aes(x = X, y = Y)) + geom_bin2d(fill = "white", color = "black", linewidth = 1.1, binwidth = Tile_width) +
          geom_point(size = 2, alpha = 0.5) + cowplot::theme_cowplot()
        Summary_smaller <- dplyr::as_tibble(ggplot2::layer_data(Tiled_plot)) %>% summarize(N_tiles = length(count),
                                                                                           Min_cells = min(count),
                                                                                           p25_cells = quantile(count, 0.25),
                                                                                           Average_cells = mean(count),
                                                                                           p50_cells = quantile(count, 0.5),
                                                                                           p75_cells = quantile(count, 0.75),
                                                                                           Max_cells = max(count))
        #Plot the tiled version of our data
        plot(
          Tiled_plot
        )
        #Generate the final tibble
        list_final <- list(Tile_dimension = Suggestion_tibble,
                           Image = Summary_smaller)
        names(list_final)[2] <- Area_tibble[[1,1]]

        #Return the final tibble
        return(list_final)
      }

      #Return the final results
      return(Suggestion_tibble)
    }
  }
