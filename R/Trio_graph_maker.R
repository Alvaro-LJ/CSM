#' Generates a plot of the TRIO spatial interaction for a single image
#'
#' The function generates a summary plot of TRIO spatial interactions. It can graph results for min dist to TRIO and TRIO score within radius approaches

#' @param Image_name A character value indicating the image to be plotted.
#' @param DATA_Phenotypes A dataframe or tibble containing a column named 'Phenotype' containing cell phenotype labels.
#' @param Strategy A character value indicating the plot that will be generated. One of the following: "Min_Distance" or "TRIO_in_Radius".
#' @param TRIO_Distances If Strategy is Min_Distance, a trio distance matrix list generated using [Trio_Distance_matrix_generator()] function.
#' @param TRIO_Cumulative If Strategy is TRIO_in_Radius, a trio cumulative interaction matrix list generated using [Trio_Cumulative_Interaction_generator()] function.
#' @param Radius If strategy is TRIO_in_Radius, a numeric value indicating the radius to be plotted. It must have been computed during the generation of the cumulative interaction matrix.
#'
#' @seealso [Trio_Distance_matrix_generator()], [Trio_Cumulative_Interaction_generator()], [Trio_Min_Distance_analyzer()], [Trio_Cells_in_Radius_analyzer()]
#'
#' @returns A graph summarizing spatial interactions for the indicated image.
#'
#' @examples
#' \dontrun{
#'#Generate distance matrix-----------------------------------
#'TRIO_Distance <-
#' Trio_Distance_matrix_generator(
#'     N_cores = 1,
#'     DATA = CSM_Phenotypecell_test,
#'     Cell_Of_Origin = "TUMOR",
#'     Target_Cell_1 = "CD8_GZMBneg",
#'     Target_Cell_2 = "CD8_GZMBpos",
#'     Perform_edge_correction = FALSE
#' )
#'
#'#Generate cumulative interactions (must contain the actual radius distance)--
#'TRIO_Cumulative <-
#'Trio_Cumulative_Interaction_generator(
#'   N_cores = 1,
#'   DATA = TRIO_Distance,
#'   Start_from = 25,
#'   Stop_at = 100,
#'   Sampling_frequency = 25
#')
#'
#'#Analyze min distance
#'Trio_graph_maker(
#'    Image_name = "ABCM22001_B14_MiniCrop.tif",
#'    DATA_Phenotypes = CSM_Phenotypecell_test,
#'    Strategy = "Min_Distance",
#'    TRIO_Distances = TRIO_Distance
#')
#'
#'#Analyze TRIO score in radius
#'Trio_graph_maker(
#'    Image_name = "ABCM22001_B14_MiniCrop.tif",
#'    DATA_Phenotypes = CSM_Phenotypecell_test,
#'    Strategy = "TRIO_in_Radius",
#'    TRIO_Cumulative = TRIO_Cumulative,
#'    Radius = 50
#')
#' }
#'
#' @export

Trio_graph_maker <-
  function(Image_name,
           DATA_Phenotypes,
           Strategy,
           TRIO_Distances = NULL,
           TRIO_Cumulative = NULL,
           Radius = NULL
  ){
    #Check arguments (Generated with ChatGPT)
    if(is.null(Image_name) || !is.character(Image_name)) {
      stop("Image_name must be specified as a non-null character string.")
    }
    # Check if Strategy is provided and matches one of the expected values
    if(is.null(Strategy) || !Strategy %in% c("Min_Distance", "TRIO_in_Radius")) {
      stop("Strategy must be either 'Min_Distance' or 'TRIO_in_Radius'.")
    }
    # Check if TRIO_Distances is provided when Strategy is 'Min_Distance'
    if(Strategy == "Min_Distance" && is.null(TRIO_Distances)) {
      stop("TRIO_Distances must be provided when Strategy is 'Min_Distance'.")
    }
    # Check if TRIO_Cumulative is provided when Strategy is 'TRIO_in_Radius'
    if(Strategy == "TRIO_in_Radius" && is.null(TRIO_Cumulative)) {
      stop("TRIO_Cumulative must be provided when Strategy is 'TRIO_in_Radius'.")
    }
    # Check if Radius is provided when Strategy is 'TRIO_in_Radius'
    if(Strategy == "TRIO_in_Radius" && is.null(Radius)) {
      stop("Radius must be provided as a non-null value when Strategy is 'TRIO_in_Radius'.")
    }

    #First Min Distance analysis
    if(Strategy == "Min_Distance") {
      #Check that Image_name is in the data
      if(!all(Image_name %in% unique(DATA_Phenotypes$Subject_Names),
              Image_name %in% names(TRIO_Distances))){
        stop(paste0("Image_name not present in data. It should be one of: ", stringr::str_c(names(TRIO_Distances), collapse = ", ")))
      }
      #Import general data frames
      DATA_Phenotypes <- DATA_Phenotypes %>% dplyr::filter(Subject_Names == Image_name)
      DATA_Distances <- TRIO_Distances[[Image_name]]

      #Define the COO and Target_Cell
      COO <- DATA_Distances[[1]][[1,1]]
      Target_Cell_1 <- DATA_Distances[[1]][[2,1]]
      Target_Cell_2 <- DATA_Distances[[1]][[3,1]]

      #First prepare 4 tibbles, a tibble with COO cells, two tibbles with one target cell each, and a tibble with other cells
      Other_tibble <- DATA_Phenotypes %>% dplyr::filter(!Phenotype %in% c(COO, Target_Cell_1, Target_Cell_2))
      Target_tibble <- DATA_Phenotypes %>% dplyr::filter(Phenotype %in% c(Target_Cell_1, Target_Cell_2))
      COO_tibble <- DATA_Phenotypes %>% dplyr::filter(Phenotype == COO)

      #Build the tibble with the COO and closest Target coordinates
      COO_Target_Coordinates <- tibble::tibble(COO_no = DATA_Distances[[2]]$Cell_Of_Origin_no)
      COO_Target_Coordinates <-
        dplyr::left_join(COO_Target_Coordinates, DATA_Phenotypes %>% dplyr::select(Cell_no, X, Y),
                         by = dplyr::join_by(COO_no == Cell_no))
      names(COO_Target_Coordinates)[2:3] <- c("X_COO", "Y_COO")


      COO_Target_Coordinates$Target_1_no <- names(DATA_Distances[[2]][-1])[max.col(DATA_Distances[[2]][-1]*-1, ties.method = "random")]
      COO_Target_Coordinates <-
        dplyr::left_join(COO_Target_Coordinates, DATA_Phenotypes %>% dplyr::select(Cell_no, X, Y),
                         by = dplyr::join_by(Target_1_no == Cell_no))
      names(COO_Target_Coordinates)[5:6] <- c("X_Target_1", "Y_Target_1")

      COO_Target_Coordinates$Target_2_no <- names(DATA_Distances[[3]][-1])[max.col(DATA_Distances[[3]][-1]*-1, ties.method = "random")]
      COO_Target_Coordinates <-
        dplyr::left_join(COO_Target_Coordinates, DATA_Phenotypes %>% dplyr::select(Cell_no, X, Y),
                         by = dplyr::join_by(Target_2_no == Cell_no))
      names(COO_Target_Coordinates)[8:9] <- c("X_Target_2", "Y_Target_2")

      #Build the plot
      PLOT <-
        ggplot() +
        geom_point(aes(x = X, y = Y), color = "lightgrey", data = Other_tibble) +
        geom_point(aes(x = X, y = Y, color = Phenotype), size = 2.5, data = Target_tibble) +
        geom_point(aes(x = X, y = Y), size = 2.5, color = "black", data = COO_tibble) +
        geom_segment(aes(x = X_COO, y = Y_COO, xend = X_Target_1, yend = Y_Target_1),
                     arrow = grid::arrow(length = unit(3, "pt"), type = "closed"),
                     color = "black",
                     linewidth = 0.75,
                     alpha = 0.5,
                     data = COO_Target_Coordinates) +
        geom_segment(aes(x = X_COO, y = Y_COO, xend = X_Target_2, yend = Y_Target_2),
                     arrow = grid::arrow(length = unit(3, "pt"), type = "closed"),
                     color = "black",
                     linewidth = 0.75,
                     alpha = 0.5,
                     data = COO_Target_Coordinates) +
        cowplot::theme_cowplot() +
        scale_x_continuous("", labels = NULL)+
        scale_y_continuous("", labels = NULL)+
        scale_fill_viridis_c("Min dist to TRIO")+
        scale_color_discrete("")+
        guides(color = guide_legend(override.aes = list(size = 8)),
               fill = guide_colorbar(theme = theme(legend.key.width = unit(10, "cm"))))+
        ggtitle(Image_name)+
        theme(panel.grid = element_blank(),
              axis.line = element_blank(),
              axis.ticks = element_blank(),
              legend.text = element_text(size = 12),
              legend.title = element_text(size = 20),
              legend.position = "bottom",
              plot.title = element_text(size = 25, hjust = 0.5, vjust = -3))
      plot(PLOT)
      return(invisible(PLOT))
    }

    #Then TRIO in radius analysis
    else if(Strategy == "TRIO_in_Radius"){
      #Check that Image_name is in the data
      if(!all(Image_name %in% unique(DATA_Phenotypes$Subject_Names),
              Image_name %in% names(TRIO_Cumulative))){
        stop(paste0("Image_name not present in data. It should be one of: ", stringr::str_c(names(TRIO_Cumulative), collapse = ", ")))
      }
      #Import general data frames
      DATA_Phenotypes <- DATA_Phenotypes %>% dplyr::filter(Subject_Names == Image_name)
      DATA_Cumulative_Interaction <- TRIO_Cumulative[[Image_name]]

      #Check that radius size is present in the data
      if(!(as.character(Radius) %in% names(DATA_Cumulative_Interaction[[2]]))){
        stop(paste0("Radius size should be one of: ", stringr::str_c(names(DATA_Cumulative_Interaction[[2]])[-1], collapse = ", ")))
      }

      #Define the COO and Target_Cell
      COO <- DATA_Cumulative_Interaction[[1]][[1,1]]
      Target_Cell_1 <- DATA_Cumulative_Interaction[[1]][[2,1]]
      Target_Cell_2 <- DATA_Cumulative_Interaction[[1]][[3,1]]

      #First prepare 3 tibbles, a tibble with COO cells, a tibble with target cells and a tibble with other cells
      Other_tibble <- DATA_Phenotypes %>% dplyr::filter(!Phenotype %in% c(COO, Target_Cell_1, Target_Cell_2))
      Target_tibble <- DATA_Phenotypes %>% dplyr::filter(Phenotype %in% c(Target_Cell_1, Target_Cell_2))
      COO_tibble <- DATA_Phenotypes %>% dplyr::filter(Phenotype == COO)

      #Calculate the cells within radius
      For_Join <- DATA_Cumulative_Interaction[["Trio_Score"]] %>% dplyr::select(1, as.character(Radius))
      names(For_Join) <- c("Cell_no", "Trio_Score_in_Radius")


      #Join the max distance with the COO tibble
      COO_tibble <-dplyr::left_join(COO_tibble, For_Join, by = "Cell_no") %>% dplyr::filter(!is.na(Trio_Score_in_Radius))

      #Build the plot
      PLOT <-
        ggplot() +
        geom_point(aes(x = X, y = Y), color = "lightgrey", data = Other_tibble) +
        ggforce::geom_circle(aes(x0 = X, y0 = Y, r = Radius, fill = Trio_Score_in_Radius), alpha = 0.3, color = NA, data = COO_tibble)+
        geom_point(aes(x = X, y = Y, color = Phenotype), size = 2.5, data = Target_tibble) +
        geom_point(aes(x = X, y = Y), size = 2, color = "black", data = COO_tibble) +
        cowplot::theme_cowplot() +
        scale_x_continuous("", labels = NULL)+
        scale_y_continuous("", labels = NULL)+
        scale_color_discrete("")+
        scale_fill_viridis_c(stringr::str_c("TRIO score ", as.character(Radius), " radius"))+
        guides(color = guide_legend(override.aes = list(size = 8)),
               fill = guide_colorbar(theme = theme(legend.key.width = unit(10, "cm"))))+
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
