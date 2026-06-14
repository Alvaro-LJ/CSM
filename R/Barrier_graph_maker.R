#' Generate a plot of the barrier metrics of single samples
#'
#'
#'
#' @param DATA_Phenotypes A dataframe or tibble containing a column named 'Phenotype' containing cell phenotype labels.
#' @param Phenotypes_included A character vector indicating the cell phenotypes that will be included in the plot.
#'
#' @param Barrier_DATA A tibble containing barrier cell data obtained using the [Barrier_effect_calculator()] function.
#' @param Image_name A character value indicating the image to be plotted (it must be present in DATA_Phenotypes and Barrier_DATA).
#' @param Color_by A character value indicating the column name in Barrier_DATA that should be used to color cell points.
#' @param Point_size A numeric value indicating the size of the points.
#'
#' @param Image_directory (OPTIONAL) Character string specifying the path to the image file that needs to be plotted.
#' @param Pixel_distance_ratio (OPTIONAL) A numeric value indicating the ratio between pixel size / cell coordinates. Use this argument when cell coordinates have been transformed to a distance unit like microns.
#' @param Channel_to_display A integer value indicating the Channel index to be displayed (1 for single channel images).
#' @param Image_rotate (OPTIONAL) A integer value indicating the degrees of rotation of the image.
#' @param Image_x_flip A logical value indicating if X image flip should be performed.
#' @param Image_y_flip A logical value indicating if Y image flip should be performed.
#'
#' @details
#' If pixel_distance_ratio is provided, it will modify the X and Y coordinates of the provided DATA_Phenotypes. The function will assume that
#' X and Y coordinates of the Barrier_DATA has already been modified with the pixel_distance_ratio argument of the [Barrier_effect_calculator()].
#'
#'
#' @returns A plot of the barrier effect for a given sample.
#'
#' @seealso [Barrier_effect_calculator()], [Barrier_effect_analyzer()]
#'
#' @examples
#' \dontrun{
#' #Create temporary input and output directories------------------------------
#' Input_Dir <- tempfile(pattern = "tempdir1_Input")
#' dir.create(Input_Dir, recursive = TRUE)
#'
#' #Save images in Input directory
#' purrr::map(1:2,
#'           function(Image){
#'             EBImage::writeImage(CSM_MiniMultiTiff_test[[Image]],
#'                                file.path(Input_Dir, names(CSM_MiniMultiTiff_test)[Image]))
#'           })
#'
#' #Generate the Image_parameter_list------------------------------
#' Image_Parameters <- list(
#'  PDL1 = list(
#'    Directory = Input_Dir,
#'    Channel = 2,
#'    Image_rotate = NULL,
#'    Image_x_flip = FALSE,
#'    Image_y_flip = TRUE),
#'  GZMB = list(
#'    Directory = Input_Dir,
#'    Channel = 3,
#'    Image_rotate = NULL,
#'    Image_x_flip = FALSE,
#'    Image_y_flip = TRUE))
#'
#' #Compute the cells and pixels in the barrier------------------------------
#' Barrier_DATA <- Barrier_effect_calculator(
#'  N_cores = 1,
#'
#'  DATA = CSM_Phenotypecell_test,
#'
#'
#'  Cell_Of_Origin = "CD8_GZMBneg",
#'  Target_Cell = "TUMOR",
#'  Barrier_Cell = c("OTHER", "CD8_GZMBpos"),
#'  Image_parameter_list = Image_Parameters,
#'
#'  Distance_sampled = 30,
#'
#'  Perform_edge_correction = FALSE,
#'
#'  Polygon_angle = 30,
#'  Target_cell_Angle_tolerance = 30
#' )
#'
#' #Summarize the results by sample------------------------------
#' Barrier_effect_analyzer(
#'  DATA = Barrier_DATA,
#'  Analysis_type = "Area"
#'  )
#'
#' #Plot the results------------------------------
#' Barrier_graph_maker(
#'   DATA_Phenotypes = CSM_Phenotypecell_test,
#'   Phenotypes_included = unique(CSM_Phenotypecell_test$Phenotype),
#'
#'   Barrier_DATA = Barrier_DATA,
#'   Image_name = "ABCM22001_B14_MiniCrop.tif",
#'
#'   Color_by = "GZMB_PixelValueSum_Area_Ratio",
#'   Point_size = 2,
#'
#'   Image_directory = dir(Input_Dir, full.names = TRUE)[[2]],
#'   Pixel_distance_ratio = NULL,
#'   Channel_to_display = 3,
#'   Image_rotate = 0,
#'   Image_x_flip = FALSE,
#'   Image_y_flip = TRUE
#' )
#'}
#'
#' @import dplyr
#'
#' @export

Barrier_graph_maker <- function(
    DATA_Phenotypes,
    Phenotypes_included,
    Barrier_DATA,
    Image_name,
    Color_by,
    Point_size = 1,


    Image_directory = NULL,
    Pixel_distance_ratio = NULL,
    Channel_to_display = 1,
    Image_rotate = 0,
    Image_x_flip = FALSE,
    Image_y_flip = FALSE
){

  #Check suggested packages
  {
    if(!requireNamespace("sf", quietly = FALSE)) stop(
      paste0("sf CRAN package is required to execute the function. Please install using the following code: ",
             expression(install.packages("sf")))
    )
    if(!is.null(Image_directory)){
      if(!requireNamespace("magick", quietly = FALSE)) stop(
        paste0("magick CRAN package is required to execute the function. Please install using the following code: ",
               expression(install.packages("magick")))
      )
    }
  }

  #Check arguments
  if(!identical(names(DATA_Phenotypes)[1:4],  c("Cell_no", "X", "Y", "Subject_Names"))) { #Check if Data is correctly formatted
    stop("DATA_phenotypes provided should have an adecuate format")
  }
  if(!identical(names(Barrier_DATA)[1:5],  c("Cell_no", "COO_X", "COO_Y", "Subject_Names", "Phenotype"))) { #Check if Data is correctly formatted
    stop("Barrier_DATA provided should have been generated using the Barrier_effect_calculator function")
  }
  if(!all(c("Original_N_Targets", "Analyzed_N_Targets", "Cumulative_distance", "Analyzed_Area", "Polygon_objects") %in% names(Barrier_DATA))){
    stop("Barrier_DATA provided should have been generated using the Barrier_effect_calculator function")
  }
  if(!Image_name %in% DATA_Phenotypes[["Subject_Names"]]) stop(paste0(Image_name, " not found in DATA_Phenotypes Subject_Names"))
  if(!Image_name %in% Barrier_DATA[["Subject_Names"]]) stop(paste0(Image_name, " not found in Barrier_DATA Subject_Names"))
  if(!"Phenotype" %in% names(DATA_Phenotypes)) stop("DATA_Phenotypes should contain a column named Phenotype specifying the cell types")
  if(!all(length(Color_by) == 1, Color_by %in% names(Barrier_DATA))) stop("Color_by should be a single column name present in Barrier_DATA")
  if(!all(is.numeric(Point_size), Point_size > 0)) stop("Point_size must be a numeric value > 0")
  if(!all(Phenotypes_included %in% unique(DATA_Phenotypes$Phenotype))) stop(paste0("Phenotypes_included must be any of the following: ",
                                                                                   stringr::str_c(unique(DATA_Phenotypes$Phenotype), collapse = ", ")))

  #Check image parameters if required
  if(!is.null(Image_directory)){
    if(!is.null(Pixel_distance_ratio)) {
      if(!all(is.numeric(Pixel_distance_ratio), Pixel_distance_ratio > 0)) stop("Pixel_distance_ratio must be NULL or a numeric value > 0")
    }
    if(!is.null(Channel_to_display)){
      if(!all(is.numeric(Channel_to_display), Channel_to_display%%1 == 0, Channel_to_display >= 1)) stop("Channel_to_display must be either NULL or an integer value larger than 0")
    }
    if(!is.null(Image_rotate)) {
      if(!all(is.numeric(Image_rotate), Image_rotate >= 0, Image_rotate <= 360, Image_rotate%%1 == 0)) stop("Image_rotate must be NULL or a numeric value between 0 and 360")
    }
    if(!is.logical(Image_x_flip)) stop("Image_x_flip must be a logical value")
    if(!is.logical(Image_y_flip)) stop("Image_y_flip must be a logical value")
  }

  #If Image directory and Pixel_distance_ratio provided, the DATA_Phenotypes coordinates must be modified (Barrier_DATA is asumed to be already modified)
  if(all(!is.null(Image_directory), !is.null(Pixel_distance_ratio))) DATA_Phenotypes <- DATA_Phenotypes %>% dplyr::mutate(X = X*Pixel_distance_ratio,
                                                                                                                          Y = Y*Pixel_distance_ratio)

  #Filter the DATA_Phenotypes and the Barrier_data to obtain data from a single Subject_Names
  DATA_Phenotypes <- DATA_Phenotypes %>% dplyr::filter(Subject_Names == Image_name)
  Barrier_DATA <- Barrier_DATA %>% dplyr::filter(Subject_Names == Image_name)

  #Filter the DATA_Phenotypes to obtain only the phenotypes required and not being the COO phenotype
  DATA_Phenotypes <- DATA_Phenotypes %>% dplyr::filter(Phenotype %in% c(Phenotypes_included, unique(Barrier_DATA$Phenotype)))

  #Obtain the polygon geometries and generate a dataframe with the required parameters to guide color filling
  Polygon_list <- purrr::map(1:nrow(Barrier_DATA), function(Row_index) Barrier_DATA[[Row_index, "Polygon_objects"]][[1]][['Unified_polygons']])
  Polygon_geometries <- do.call(c, lapply(Polygon_list, sf::st_geometry))
  Polygon_geometries_df <- sf::st_as_sf(data.frame(id = seq_along(Polygon_geometries), geometry = Polygon_geometries))
  Polygon_geometries_df <- Polygon_geometries_df %>% dplyr::mutate(Polygon_fill = Barrier_DATA[[Color_by]])

  #If Image provided
  if(!is.null(Image_directory)){
    #Read image
    if(is.null(Channel_to_display)) Image <- magick::image_read(Image_directory)
    if(!is.null(Channel_to_display)) Image <- magick::image_read(Image_directory)[Channel_to_display]

    #Execute rotation if required
    if(!is.null(Image_rotate)) Image <- Image %>% magick::image_rotate(degrees = Image_rotate)
    #Execute x or y flip if required
    if(Image_x_flip) Image <- Image %>% magick::image_flop()
    if(Image_y_flip) Image <- Image %>% magick::image_flip() #Here we will flip the image as required by user as annotation_raster will flip it for plotting purposes

    #obtain the max width and height
    X_max <- magick::image_info(Image)$width
    Y_max <- magick::image_info(Image)$height

    #Generate the plot
    PLOT <- ggplot() +
      #1st plot the image
      annotation_raster(Image, xmin = 1, xmax = X_max, ymin = 0, ymax = Y_max) +

      #2nd layer the sf object
      geom_sf(aes(fill = Polygon_fill), color = NA, alpha = 0.5, data = Polygon_geometries_df) +

      #3rd layer are the general points
      geom_point(aes(X, Y, color = Phenotype), size = Point_size, data = DATA_Phenotypes) +

      #Change_scales
      scale_color_manual("", values = unname(pals::polychrome(n = length(unique(DATA_Phenotypes$Phenotype))))) +
      scale_fill_viridis_c(Color_by) +

      #title and theme
      ggtitle(Image_name) +
      theme_minimal() +
      theme(panel.grid = element_blank(),
            plot.title = element_text(hjust = 0.5))

    plot(PLOT)

    return(invisible(PLOT))
  }

  #If image not provided
  if(is.null(Image_directory)){

    #Generate the graph
    PLOT <- ggplot() +
      #1st layer the sf object
      geom_sf(aes(fill = Polygon_fill), color = NA, alpha = 0.5, data = Polygon_geometries_df) +

      #2nd layer are the general points
      geom_point(aes(X, Y, color = Phenotype), size = Point_size, data = DATA_Phenotypes) +

      #Change_scales
      scale_color_manual("", values = unname(pals::polychrome(n = length(unique(DATA_Phenotypes$Phenotype))))) +
      scale_fill_viridis_c(Color_by) +

      #title and theme
      ggtitle(Image_name) +
      theme_minimal() +
      theme(panel.grid = element_blank(),
            plot.title = element_text(hjust = 0.5))

    plot(PLOT)

    return(invisible(PLOT))
  }
}
















