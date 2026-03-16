#' Generates a plot comparing the tile pattern and the structures identified
#'
#' Given a tiled image cell data (generated using the [Image_tiling_processing_function()] or the [Tiled_Image_Clustering_function()]) and
#' Morphometry_data (generated using the [Structure_morphometry_calculator()]) the function will graph a specific image comparing the tile pattern
#' and the structures identified. Structures can be identified can be colored by different structure metrics.
#'
#'
#' @param Tiled_Images A list containing tiled images obtained using [Image_tiling_processing_function()] or the [Tiled_Image_Clustering_function()].
#' @param Variable_name A character value indicating the column name containing the cell labels.
#' @param DATA_Morphometry A list containing structure morphometry data (generated using the [Structure_morphometry_calculator()]).
#'
#' @param Image_name A single character value specifying the name of the image to be plotted.
#' @param Color_by A single character value indicating the variable defining the color of the structures.
#'
#' @returns Plots the resulting graphs and returns a list with the graphs.
#'
#' @seealso [Image_tiling_processing_function()], [Tiled_Image_Clustering_function()], [Structure_morphometry_calculator()], [Structure_morphometry_analyzer()]
#'
#' @examples
#' \dontrun{
#' #Divide cells into discrete tiles
#'Tiled_Images <-
#' Image_tiling_processing_function(
#'  N_cores = 1,
#'  DATA = CSM_Phenotypecell_test,
#'  Tile_width = 25,
#'  Tile_height = 25,
#'  Variables_to_keep = "Phenotype")
#'
#' #Try to identify morphological structures for OTHER cells according to tile distribution
#'Structures <-
#' Structure_morphometry_calculator(N_cores = 1,
#'
#'                                 Tiled_Images = Tiled_Images,
#'                                 Variable_name = "Phenotype",
#'                                 Feature_to_analyze = "OTHER",
#'
#'                                 Min_tile_number = 1,
#'                                 Smooth_amount = 0.01,
#'                                 Tolerance_value = 1,
#'                                 Neighborhood_distance = 1,
#'                                 Minimum_cell_no_per_tile = 2)
#'
#' #Analyze the results
#' Structure_morphometry_analyzer(DATA = Structures)
#'
#' #Graph the results
#' Structure_morphometry_graph_maker(Tiled_Images = Tiled_Images,
#'                                   Variable_name = "Phenotype",
#'                                   DATA_Morphometry = Structures,
#'                                   Color_by = "Structure_ID",
#'                                   Image_name = "ABCM22001_B09_MiniCrop.tif")
#'}
#'
#'
#' @export


Structure_morphometry_graph_maker <-
  function(
    Tiled_Images,
    Variable_name,
    DATA_Morphometry,
    Image_name,
    Color_by = "s.area"
  ){

    ########Check suggested packages########
    {
      if(!requireNamespace("EBImage", quietly = TRUE)) stop(
        paste0("EBImage Bioconductor package is required to execute the function. Please install using the following code: ",
               expression({
                 if (!require("BiocManager", quietly = TRUE))
                   install.packages("BiocManager")

                 BiocManager::install("EBImage")
               })
        )
      )
    }

    ########Check Arguments########
    #Check the DATA_Morphometry element
    if(!is.list(DATA_Morphometry)) stop("DATA_Morphometry must be generated using the Structure_morphometry_calculator function")

    #Check content of DATA_Morphometry is adecuate
    Conflictive_elements <-
      !purrr::map_lgl(DATA_Morphometry, function(Image){
        identical(names(Image), c("Morphology_tibble", "Structure_mask"))
      })

    if(any(Conflictive_elements)) stop("DATA_Morphometry must be generated using the Structure_morphometry_calculator function")

    #Check that color by is adequate
    if(!all(length(Color_by) == 1,
            Color_by %in% c("Structure_ID", "s.area", "s.perimeter", "s.radius.mean", "s.radius.sd", "s.radius.min", "s.radius.max",
                            "m.majoraxis", "m.eccentricity", "m.theta_AVERAGE")
            )
    ) stop("Color_by must be one of the following: s.area, s.perimeter, s.radius.mean, s.radius.sd, s.radius.min, s.radius.max, m.majoraxis, m.eccentricity, m.theta_AVERAGE")

    #Check Tiled_Images
    if(!is.list(Tiled_Images)) stop("Tiled_Images should be a list created with the Image_tiling_processing_function or Tiled_Image_Clustering_function")

    List_content_tibble <- unique(purrr::map_lgl(names(Tiled_Images), ~tibble::is_tibble(Tiled_Images[[.]])))
    if(length(List_content_tibble) > 1) stop("Tiled_Images should be a list created with the Image_tiling_processing_function or Tiled_Image_Clustering_function")

    #Check Image_name provided
    if(!all(length(Image_name) == 1,
            Image_name %in% names(DATA_Morphometry),
            Image_name %in% names(Tiled_Images))) stop(paste0("Image name should be a single character value included in: \n", stringr::str_c(names(DATA_Morphometry), collapse = "\n")))#DATA morphometry is the rate limiting data source

    ########GENERATE THE GRAPHS########
    #Select the Tiled_Image and DATA_Morphometry
    Tiled_Images <- Tiled_Images[[Image_name]]
    DATA_Morphometry <- DATA_Morphometry[[Image_name]]
    Selected_cell_type <- unique(DATA_Morphometry[[1]][["Cell_analyzed"]])

    #First the tiled images
    #If the Tiled_Images were generated with the Tiled_image clustering function then perform appropriate corrections
    if(List_content_tibble){
      #If variable name not present stop computation
      if(!Variable_name %in% names(Tiled_Images)) stop(paste0(Variable_name, " not present in Tiled_Images column names."))

      #Else transform the Tiled_images (turning tiles with selected cell type to 1 and 0 otherwise)
      Tiled_Images  <- Tiled_Images %>%
        dplyr::select(Subject_Names, tile_X_centroid, tile_Y_centroid, tile_id, tile_xmin, tile_xmax, tile_ymin, tile_ymax, dplyr::all_of(Variable_name))

      names(Tiled_Images)[9] <- "Target_column"

      Tiled_Images <- Tiled_Images %>% dplyr::mutate(Target_column = dplyr::case_when(Target_column == Selected_cell_type ~ 1,
                                                                                      TRUE ~ 0))
      #Compute all the possible tiles
      Full_tile_tibble <- tidyr::expand_grid(tile_xmin = unique(Tiled_Images$tile_xmin), tile_ymin = unique(Tiled_Images$tile_ymin))
      #Compute the tile width and height
      Tile_width <- unique(Tiled_Images$tile_xmax - Tiled_Images$tile_xmin)
      Tile_height <- unique(Tiled_Images$tile_ymax - Tiled_Images$tile_ymin)

      Full_tile_tibble <- Full_tile_tibble %>% dplyr::mutate(tile_xmax = tile_xmin + Tile_width,
                                                             tile_ymax = tile_ymin + Tile_height)
      #Join with existing data
      Full_tile_tibble <-
        dplyr::left_join(Full_tile_tibble, Tiled_Images, dplyr::join_by(tile_xmin, tile_xmax, tile_ymin, tile_ymax))

      Full_tile_tibble$tile_id[is.na(Full_tile_tibble$tile_id)] <- "ZGenerated_tile"
      Full_tile_tibble$tile_X_centroid[is.na(Full_tile_tibble$tile_X_centroid)] <- Full_tile_tibble$tile_xmin[is.na(Full_tile_tibble$tile_X_centroid)] + (Tile_width/2)
      Full_tile_tibble$tile_Y_centroid[is.na(Full_tile_tibble$tile_Y_centroid)] <- Full_tile_tibble$tile_ymin[is.na(Full_tile_tibble$tile_Y_centroid)] + (Tile_height/2)
      Full_tile_tibble$Subject_Names <- Image_name
      Full_tile_tibble$Target_column[is.na(Full_tile_tibble$Target_column)] <- 0


      Full_tile_tibble <- Full_tile_tibble %>% dplyr::select(Subject_Names, tile_X_centroid, tile_Y_centroid, tile_id, tile_xmin, tile_xmax, tile_ymin, tile_ymax, Target_column)

      Tiled_Images <- Full_tile_tibble
    }

    #If it has not been obtained with the Tile_image clustering function proceed as usual
    if(!List_content_tibble){

      #Check that Variable_name is present in all images
      if(!Variable_name %in% names(Tiled_Images[[2]])) stop(paste0(Variable_name, " not present in Tiled_Images column names."))

      #Compute the counts per tile
      Interim <- Tiled_Images[[2]] %>% dplyr::rename("Target" = dplyr::all_of(Variable_name))

      Count_per_tile <- Interim %>% dplyr::filter(Target == Selected_cell_type) %>% dplyr::count(tile_id)
      Tile_info <- left_join(Tiled_Images[[1]], Count_per_tile, by = "tile_id")
      Tile_info <- Tile_info %>% mutate(Target_column = case_when(n >= 1 ~ n,
                                                                  TRUE ~ 0))
      Tile_info$Subject_Names <- Image_name

      Tile_info <- Tile_info %>%
        dplyr::select(Subject_Names, tile_X_centroid, tile_Y_centroid, tile_id, tile_xmin, tile_xmax, tile_ymin, tile_ymax, Target_column)

      Tiled_Images <- Tile_info
    }

    Tiles_plot <-
      Tiled_Images %>%
      ggplot() +
      geom_rect(aes(xmin = tile_xmin, xmax = tile_xmax, ymin = tile_ymin, ymax = tile_ymax, fill = as.numeric(Target_column))) +
      scale_x_continuous("", labels = NULL) +
      scale_y_continuous("", labels = NULL) +
      scale_fill_gradient2(Selected_cell_type, low = "black", mid = "black", high = "red", midpoint = 0.1) +
      ggtitle("Tile pattern", Image_name) +
      theme_minimal() +
      theme(panel.grid = element_blank(),
            aspect.ratio =
              (max(Tiled_Images$tile_xmax) - min(Tiled_Images$tile_xmin)) /
              (max(Tiled_Images$tile_ymax) - min(Tiled_Images$tile_ymin))
      )

    #Now the structure mask
    Tissue_mask <- DATA_Morphometry[[2]] %>% EBImage::rotate(angle = 90) %>% EBImage::flop()

    Tissue_mask_tibble <- tibble(Row_index = rep(dim(Tissue_mask)[[1]]:1, times= dim(Tissue_mask)[[2]]),
                                 Col_index = rep(1:dim(Tissue_mask)[[2]], each = dim(Tissue_mask)[[1]]),
                                 Structure_ID = as.vector(Tissue_mask))
    Row_max <- max(Tissue_mask_tibble$Row_index)
    Col_max <- max(Tissue_mask_tibble$Col_index)
    Tissue_mask_tibble <- Tissue_mask_tibble %>% dplyr::filter(Structure_ID != 0)

    if(Color_by != "Structure_ID"){
      Tissue_mask_tibble <- left_join(Tissue_mask_tibble, DATA_Morphometry[[1]] %>% dplyr::select(Structure_ID, dplyr::all_of(Color_by)), by = "Structure_ID")
      names(Tissue_mask_tibble)[4] <- "Target"

      Mask_plot <-
        Tissue_mask_tibble %>%
        ggplot(aes(x = Col_index, y = Row_index, fill = Target)) +
        geom_tile(na.rm = TRUE) +
        scale_x_continuous("", labels = NULL, limits = c(1, Col_max)) +
        scale_y_continuous("", labels = NULL, limits = c(1, Row_max)) +
        scale_fill_viridis_c(Color_by, option = "C") +
        ggtitle("Structures", Image_name) +
        theme_minimal() +
        theme(panel.grid = element_blank(),
              panel.background = element_rect(fill = "black"),
              aspect.ratio = Col_max / Row_max
        )
    }
    else{
      Mask_plot <-
        Tissue_mask_tibble %>%
        ggplot(aes(x = Col_index, y = Row_index, fill = as.factor(Structure_ID))) +
        geom_tile(na.rm = TRUE) +
        scale_x_continuous("", labels = NULL, limits = c(1, Col_max)) +
        scale_y_continuous("", labels = NULL, limits = c(1, Row_max)) +
        scale_fill_manual(Color_by, values = unname(pals::polychrome(n = length(unique(Tissue_mask_tibble$Structure_ID))))) +
        ggtitle("Structures", Image_name) +
        theme_minimal() +
        theme(panel.grid = element_blank(),
              panel.background = element_rect(fill = "black"),
              aspect.ratio = Col_max / Row_max
        )
    }


    #plot both graphs
    plot(Tiles_plot + Mask_plot + patchwork::plot_layout(nrow = 1, ncol = 2))

    #return a list with both plots
    return(invisible(
      list(Tiles_plot = Tiles_plot,
           Mask_plot = Mask_plot)
    ))
  }

