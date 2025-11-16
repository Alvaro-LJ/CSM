#' Generates a plot of the spatial interaction for a single image
#'
#' The function generates a summary plot of the spatial interaction between two cell types.
#' It can graph results for average, min and max distance as well as cell in radius analyses.
#'
#' @param Image_name A character value indicating the image to be plotted.
#' @param DATA_Phenotypes A dataframe or tibble containing a column named 'Phenotype' containing cell phenotype labels.
#' @param Strategy A character value indicating the plot that will be generated. One of the following: "Min_Distance", "Average_Distance", "Max_Distance" or "Cells_in_Radius".
#' @param DATA_Distances If Strategy is "Min_Distance", "Average_Distance" or "Max_Distance",  a distance matrix list generated using [Distance_matrix_generator()] function.
#' @param DATA_Cumulative If Strategy is Cells_in_Radius, a cumulative interaction matrix list generated using [Cumulative_Interaction_generator()] function.
#' @param Radius If strategy is Cells_in_Radius, a numeric value indicating the radius to be plotted. It must have been computed during the generation of the cumulative interaction matrix.
#'
#' @seealso [Distance_matrix_generator()], [Cumulative_Interaction_generator()], [Distance_analyzer()], [Cells_in_Radius_analyzer()]
#'
#' @returns A graph summarizing spatial interactions for the indicated image.
#'
#' @examples
#' \dontrun{
#' #Generate distance matrix----------------------------------------------------
#'DATA_Distances <-
#' Distance_matrix_generator(
#'     N_cores = 1,
#'     DATA = CSM_Phenotypecell_test,
#'     Cell_Of_Origin = "CD8_GZMBneg",
#'     Target_Cell = "TUMOR",
#'     Allow_Cero_Distance = FALSE,
#'     Perform_edge_correction = FALSE
#')
#'
#' #Generate cumulative interactions (must contain the actual radius distance)--
#'DATA_Cumulative <-
#'Cumulative_Interaction_generator(
#'    N_cores = 1,
#'    DATA = DATA_Distances,
#'    Start_from = 25,
#'    Stop_at = 100,
#'    Sampling_frequency = 25
#')
#'
#'#Graph min distance to target-------------------------------------------------
#'Cell_to_Cell_graph_maker(
#'    Image_name = "ABCM22001_B14_MiniCrop.tif",
#'    DATA_Phenotypes = CSM_Phenotypecell_test,
#'    Strategy = "Min_Distance",
#'    DATA_Distances = DATA_Distances
#')
#Graph cells within radius------------------------------------------------------
#'Cell_to_Cell_graph_maker(
#'    Image_name = "ABCM22001_B14_MiniCrop.tif",
#'    DATA_Phenotypes = CSM_Phenotypecell_test,
#'    Strategy = "Cells_in_Radius",
#'    DATA_Cumulative = DATA_Cumulative,
#'    Radius = 50
#')
#' }
#'
#' @export

Cell_to_Cell_graph_maker <-
  function(Image_name,
           DATA_Phenotypes,
           Strategy,
           DATA_Distances = NULL,
           DATA_Cumulative = NULL,
           Radius = NULL){

    # Check that Image_name is provided and is character
    if (is.null(Image_name) || !is.character(Image_name)) {
      stop("Image_name must be provided and must be of type character.")
    }
    # Check that Strategy is provided and valid
    valid_strategies <- c("Min_Distance", "Average_Distance", "Max_Distance", "Cells_in_Radius")
    if (is.null(Strategy) || !Strategy %in% valid_strategies) {
      stop(paste0("Strategy must be one of the following: ", paste(valid_strategies, collapse = ", ")))
    }
    # Check that DATA_Phenotypes is provided and is a data frame
    if (is.null(DATA_Phenotypes) || !is.data.frame(DATA_Phenotypes)) {
      stop("DATA_Phenotypes must be provided and must be a data frame.")
    }
    # Check that DATA_Distances is provided and is a list (for non-"Cells_in_Radius" strategies)
    if (Strategy != "Cells_in_Radius" && (is.null(DATA_Distances) || !is.list(DATA_Distances))) {
      stop("DATA_Distances must be provided and must be a list when using Min_Distance, Average_Distance, or Max_Distance.")
    }
    # Check that DATA_Cumulative is provided and is a list (only for "Cells_in_Radius")
    if (Strategy == "Cells_in_Radius" && (is.null(DATA_Cumulative) || !is.list(DATA_Cumulative))) {
      stop("DATA_Cumulative must be provided and must be a list when using Cells_in_Radius strategy.")
    }
    # Check that Radius is provided and is numeric (only for "Cells_in_Radius")
    if (Strategy == "Cells_in_Radius" && (is.null(Radius) || !is.numeric(Radius))) {
      stop("Radius must be provided and must be numeric when using Cells_in_Radius strategy.")
    }
    # Ensure that Image_name exists in the DATA_Phenotypes and DATA_Distances/Cumulative
    if (Strategy %in% c("Min_Distance", "Average_Distance", "Max_Distance")) {
      if (!all(Image_name %in% unique(DATA_Phenotypes$Subject_Names), Image_name %in% names(DATA_Distances))) {
        stop(paste0("Image_name not present in data. It should be one of: ", stringr::str_c(names(DATA_Distances), collapse = ", ")))
      }
    } else if (Strategy == "Cells_in_Radius") {
      if (!all(Image_name %in% unique(DATA_Phenotypes$Subject_Names), Image_name %in% names(DATA_Cumulative))) {
        stop(paste0("Image_name not present in data. It should be one of: ", stringr::str_c(names(DATA_Cumulative), collapse = ", ")))
      }
    }

    if(Strategy == "Min_Distance") {
      #Import general data frames
      DATA_Phenotypes <- DATA_Phenotypes %>% dplyr::filter(Subject_Names == Image_name)
      DATA_Distances <- DATA_Distances[[Image_name]]

      #Define the COO and Target_Cell
      COO <- DATA_Distances[[1]][[1,1]]
      Target_Cell <- DATA_Distances[[1]][[2,1]]

      #First prepare 3 tibbles, a tibble with COO cells, a tibble with target cells and a tibble with other cells
      Other_tibble <- DATA_Phenotypes %>% dplyr::filter(!Phenotype %in% c(COO, Target_Cell))
      Target_tibble <- DATA_Phenotypes %>% dplyr::filter(Phenotype == Target_Cell)
      COO_tibble <- DATA_Phenotypes %>% dplyr::filter(Phenotype == COO)

      #Calculate the min distance
      For_Join <-dplyr::bind_cols(DATA_Distances[[2]][1],
                                  purrr::map_dbl(1:nrow(DATA_Distances[[2]]), function(Row){
                                    min(DATA_Distances[[2]][Row,-1])
                                  })
      )
      names(For_Join) <- c("Cell_no", "Min_Distance")

      #Join the min distance with the COO tibble
      COO_tibble <-dplyr::left_join(COO_tibble, For_Join, by = "Cell_no") %>% dplyr::filter(!is.na(Min_Distance))

      #Build the plot
      PLOT <-
        ggplot() +
        geom_point(aes(x = X, y = Y), color = "lightgrey", data = Other_tibble) +
        ggforce::geom_circle(aes(x0 = X, y0 = Y, r = Min_Distance, fill = Min_Distance), alpha = 0.3, color = NA, data = COO_tibble)+
        geom_point(aes(x = X, y = Y), size = 2, color = "red", data = Target_tibble) +
        geom_point(aes(x = X, y = Y), size = 2, color = "blue", data = COO_tibble) +
        cowplot::theme_cowplot() +
        scale_x_continuous("", labels = NULL)+
        scale_y_continuous("", labels = NULL)+
        scale_fill_viridis_c("Min dist to target")+
        ggtitle(Image_name)+
        theme(panel.grid = element_blank(),
              axis.line = element_blank(),
              axis.ticks = element_blank(),
              legend.text = element_text(size = 15),
              legend.title = element_text(size = 20),
              legend.position = "bottom",
              plot.title = element_text(size = 25, hjust = 0.5, vjust = -3))
      plot(PLOT)
      return(invisible(PLOT))
    }

    else if(Strategy == "Average_Distance") {
      #Import general data frames
      DATA_Phenotypes <- DATA_Phenotypes %>% dplyr::filter(Subject_Names == Image_name)
      DATA_Distances <- DATA_Distances[[Image_name]]

      #Define the COO and Target_Cell
      COO <- DATA_Distances[[1]][[1,1]]
      Target_Cell <- DATA_Distances[[1]][[2,1]]

      #First prepare 3 tibbles, a tibble with COO cells, a tibble with target cells and a tibble with other cells
      Other_tibble <- DATA_Phenotypes %>% dplyr::filter(!Phenotype %in% c(COO, Target_Cell))
      Target_tibble <- DATA_Phenotypes %>% dplyr::filter(Phenotype == Target_Cell)
      COO_tibble <- DATA_Phenotypes %>% dplyr::filter(Phenotype == COO)

      #Calculate the average distance
      For_Join <-dplyr::bind_cols(DATA_Distances[[2]][1],
                                  purrr::map_dbl(1:nrow(DATA_Distances[[2]]), function(Row){
                                    mean(unlist(DATA_Distances[[2]][Row,-1]))
                                  })
      )
      names(For_Join) <- c("Cell_no", "Average_Distance")

      #Join the average distance with the COO tibble
      COO_tibble <-dplyr::left_join(COO_tibble, For_Join, by = "Cell_no") %>% dplyr::filter(!is.na(Average_Distance))

      #Build the plot
      PLOT <-
        ggplot() +
          geom_point(aes(x = X, y = Y), color = "lightgrey", data = Other_tibble) +
          ggforce::geom_circle(aes(x0 = X, y0 = Y, r = Average_Distance, fill = Average_Distance), alpha = 0.3, color = NA, data = COO_tibble)+
          geom_point(aes(x = X, y = Y), size = 2, color = "red", data = Target_tibble) +
          geom_point(aes(x = X, y = Y), size = 2, color = "blue", data = COO_tibble) +
          cowplot::theme_cowplot() +
          scale_x_continuous("", labels = NULL)+
          scale_y_continuous("", labels = NULL)+
          scale_fill_viridis_c("Average dist to target")+
          ggtitle(Image_name)+
          theme(panel.grid = element_blank(),
                axis.line = element_blank(),
                axis.ticks = element_blank(),
                legend.text = element_text(size = 15),
                legend.title = element_text(size = 20),
                legend.position = "bottom",
                plot.title = element_text(size = 25, hjust = 0.5, vjust = -3))
      plot(PLOT)
      return(invisible(PLOT))
    }

    else if(Strategy == "Max_Distance") {
      #Import general data frames
      DATA_Phenotypes <- DATA_Phenotypes %>% dplyr::filter(Subject_Names == Image_name)
      DATA_Distances <- DATA_Distances[[Image_name]]

      #Define the COO and Target_Cell
      COO <- DATA_Distances[[1]][[1,1]]
      Target_Cell <- DATA_Distances[[1]][[2,1]]

      #First prepare 3 tibbles, a tibble with COO cells, a tibble with target cells and a tibble with other cells
      Other_tibble <- DATA_Phenotypes %>% dplyr::filter(!Phenotype %in% c(COO, Target_Cell))
      Target_tibble <- DATA_Phenotypes %>% dplyr::filter(Phenotype == Target_Cell)
      COO_tibble <- DATA_Phenotypes %>% dplyr::filter(Phenotype == COO)

      #Calculate the max distance
      For_Join <-dplyr::bind_cols(DATA_Distances[[2]][1],
                                  purrr::map_dbl(1:nrow(DATA_Distances[[2]]), function(Row){
                                    max(DATA_Distances[[2]][Row,-1])
                                  })
      )
      names(For_Join) <- c("Cell_no", "Max_Distance")

      #Join the max distance with the COO tibble
      COO_tibble <-dplyr::left_join(COO_tibble, For_Join, by = "Cell_no") %>% dplyr::filter(!is.na(Max_Distance))

      #Build the plot
      PLOT <-
        ggplot() +
          geom_point(aes(x = X, y = Y), color = "lightgrey", data = Other_tibble) +
          ggforce::geom_circle(aes(x0 = X, y0 = Y, r = Max_Distance, fill = Max_Distance), alpha = 0.3, color = NA, data = COO_tibble)+
          geom_point(aes(x = X, y = Y), size = 2, color = "red", data = Target_tibble) +
          geom_point(aes(x = X, y = Y), size = 2, color = "blue", data = COO_tibble) +
          cowplot::theme_cowplot() +
          scale_x_continuous("", labels = NULL)+
          scale_y_continuous("", labels = NULL)+
          scale_fill_viridis_c("Max dist to target")+
          ggtitle(Image_name)+
          theme(panel.grid = element_blank(),
                axis.line = element_blank(),
                axis.ticks = element_blank(),
                legend.text = element_text(size = 15),
                legend.title = element_text(size = 20),
                legend.position = "bottom",
                plot.title = element_text(size = 25, hjust = 0.5, vjust = -3))
      plot(PLOT)
      return(invisible(PLOT))
    }

    else if(Strategy == "Cells_in_Radius"){
      #Import general data frames
      DATA_Phenotypes <- DATA_Phenotypes %>% dplyr::filter(Subject_Names == Image_name)
      DATA_Cumulative_Interaction <- DATA_Cumulative[[Image_name]]

      #Check that radius size is present in the data
      if(!(as.character(Radius) %in% names(DATA_Cumulative_Interaction[[2]]))){
        stop(paste0("Radius size should be one of: ", stringr::str_c(names(DATA_Cumulative_Interaction[[2]])[-1], collapse = ", ")))
      }

      #Define the COO and Target_Cell
      COO <- DATA_Cumulative_Interaction[[1]][[1,1]]
      Target_Cell <- DATA_Cumulative_Interaction[[1]][[2,1]]

      #First prepare 3 tibbles, a tibble with COO cells, a tibble with target cells and a tibble with other cells
      Other_tibble <- DATA_Phenotypes %>% dplyr::filter(!Phenotype %in% c(COO, Target_Cell))
      Target_tibble <- DATA_Phenotypes %>% dplyr::filter(Phenotype == Target_Cell)
      COO_tibble <- DATA_Phenotypes %>% dplyr::filter(Phenotype == COO)

      #Calculate the cells within radius
      For_Join <- DATA_Cumulative_Interaction[[2]] %>% dplyr::select(1, as.character(Radius))
      names(For_Join) <-
        names(For_Join) <- c("Cell_no", "Cells_in_Radius")

      #Join the max distance with the COO tibble
      COO_tibble <-dplyr::left_join(COO_tibble, For_Join, by = "Cell_no") %>% dplyr::filter(!is.na(Cells_in_Radius))

      #Build the plot
      PLOT <-
        ggplot() +
          geom_point(aes(x = X, y = Y), color = "lightgrey", data = Other_tibble) +
          ggforce::geom_circle(aes(x0 = X, y0 = Y, r = Radius, fill = Cells_in_Radius), alpha = 0.3, color = NA, data = COO_tibble)+
          geom_point(aes(x = X, y = Y), size = 2, color = "red", data = Target_tibble) +
          geom_point(aes(x = X, y = Y), size = 2, color = "blue", data = COO_tibble) +
          cowplot::theme_cowplot() +
          scale_x_continuous("", labels = NULL)+
          scale_y_continuous("", labels = NULL)+
          scale_fill_viridis_c(stringr::str_c("Cells in ", as.character(Radius), " radius"))+
          ggtitle(Image_name)+
          theme(panel.grid = element_blank(),
                axis.line = element_blank(),
                axis.ticks = element_blank(),
                legend.text = element_text(size = 15),
                legend.title = element_text(size = 20),
                legend.position = "bottom",
                plot.title = element_text(size = 25, hjust = 0.5, vjust = -3))
      plot(PLOT)
      return(invisible(PLOT))
    }
  }
