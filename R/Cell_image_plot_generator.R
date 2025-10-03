#' Plot an image with a cell overlay.
#'
#' The function generates a plot with an image in the background and cells overlayed. Cell information must have been obtained from the image in order to match X Y coordinates with image coordinates.
#'
#' @param Image_directory Character string specifying the path to the image file that needs to be plotted.
#' @param Channel_to_display A integer value indicating the Channel index to be displayed (1 for single channel images).
#' @param Image_rotate (OPTIONAL) A integer value indicating the degrees of rotation of the image.
#' @param Image_x_flip A logical value indicating if X image flip should be performed.
#' @param Image_y_flip A logical value indicating if X image flip should be performed.
#' @param Gamma_level A numeric value between -3 and +3 indicating the gamma level to apply to the image.
#' @param Equalize A logical value indicating if the image should be equalized.
#' @param Black_level A integer value between 0 and 100 indicating the minimum pixel intensity.
#' @param White_level A integer value between 0 and 100 indicating the maximum pixel intensity.
#'
#' @param DATA A dataframe or tibble containing cell feature data.
#' @param Image_name A character value indicating the image to be plotted.
#' @param Color_by (OPTIONAL) A character value indicating the column name in DATA that should be used to color cell points.
#' @param Point_size A numeric value indicating the size of the points.
#' @param Pixel_distance_ratio (OPTIONAL) A numeric value indicating the ratio between pixel size / cell coordinates. Use this argument when cell coordinates have been transformed to a distance unit like microns.
#'
#' @returns Generates a plot with the image and the cells overlay.
#'
#' @seealso [Cell_plotter()]
#'
#' @export

Cell_image_plot_generator <-
  function(Image_directory = NULL,
           Channel_to_display = NULL,
           Image_rotate = NULL,
           Image_x_flip = FALSE,
           Image_y_flip = FALSE,
           Gamma_level = 0,
           Equalize = FALSE,
           Black_level = 0,
           White_level = 100,

           DATA = NULL,
           Image_name = NULL,
           Color_by = NULL,
           Point_size = 1,
           Pixel_distance_ratio = NULL
  ){
    #Check suggested packages
    if(!requireNamespace("magick", quietly = FALSE)) stop(
      paste0("magick CRAN package is required to execute the function. Please install using the following code: ",
             expression(install.packages("magick")))
    )

    on.exit(gc())

    DATA <- DATA
    #Check arguments
    if(!is.null(Image_rotate)) {
      if(!all(is.numeric(Image_rotate), Image_rotate >= 0, Image_rotate <= 360, Image_rotate%%1 == 0)) stop("Image_rotate must be NULL or a numeric value between 0 and 360")
    }
    if(!is.logical(Image_x_flip)) stop("Image_x_flip must be a logical value")
    if(!is.logical(Image_y_flip)) stop("Image_y_flip must be a logical value")
    if(!all(is.numeric(Gamma_level), Gamma_level >= -3, Gamma_level <= 3)) stop ("Gamma_level must be a numeric value between -3 and +3")
    if(!is.logical(Equalize)) stop("Equalize must be a logical value")
    if(!all(is.numeric(Black_level), Black_level >= 0, Black_level <= 100, Black_level < White_level)) stop("Black_level must be a numeric value between 0 - 100 and smaller than White_level")
    if(!all(is.numeric(White_level), White_level >= 0, White_level <= 100)) stop("White_level must be a numeric value between 0 - 100")

    if(!identical(names(DATA)[1:4], c("Cell_no", "X", "Y", "Subject_Names"))) stop("DATA must be adequately formatted")
    if(!Image_name %in% DATA[["Subject_Names"]]) stop(paste0(Image_name, " not found in DATA Subject_Names"))
    if(!any(is.null(Color_by), Color_by %in% names(DATA))) stop(paste0(Color_by, " not found in DATA"))
    if(!all(is.numeric(Point_size), Point_size > 0)) stop("Point_size must be a numeric value > 0")
    if(!is.null(Pixel_distance_ratio)) {
      if(!all(is.numeric(Pixel_distance_ratio), Pixel_distance_ratio > 0)) stop("Pixel_distance_ratio must be NULL or a numeric value > 0")
    }

    #Lets get the data and select only the cells and the selected column specifying the color
    DATA <- DATA %>% dplyr::filter(Subject_Names == Image_name)

    #Select the desired variables
    if(!is.null(Color_by)) DATA <- DATA %>% dplyr::select(1:4, all_of(Color_by)) else DATA <- DATA %>% dplyr::select(1:4)
    if(!is.null(Color_by)) names(DATA)[5] <- "Variable"

    #Adjust to pixel ratio if required
    if(!is.null(Pixel_distance_ratio)){
      DATA$X <- DATA$X * Pixel_distance_ratio
      DATA$Y <- DATA$Y * Pixel_distance_ratio
    }

    #Work on the image
    Image <- magick::image_read(Image_directory)[Channel_to_display] #Get original Image
    if(Image_x_flip) Image <- Image %>% magick::image_flop()
    if(Image_y_flip) Image <- Image %>% magick::image_flip()
    if(!is.null(Image_rotate)) Image <- Image %>% magick::image_rotate(degrees = Image_rotate)
    #Get image info
    Image_data <- magick::image_info(Image)
    Image_width <- Image_data$width
    Image_height <- Image_data$height

    #Perform image modifications
    if(Equalize) Image <- Image %>% magick::image_equalize() #Equalize if necesary
    Image <-  Image %>% magick::image_level(black_point = Black_level,
                                            white_point = White_level,
                                            mid_point = 10^Gamma_level) #Chane withe, black and gamma

    #If color needs to be applied
    if(!is.null(Color_by)){
      #If its numeric
      if(is.numeric(DATA[["Variable"]])){
        PLOT <- ggplot() +
          annotation_raster(Image, xmin = 0, xmax = Image_width, ymin = 0, ymax = Image_height, interpolate = TRUE) +
          geom_point(aes(x = X, y = Y, color = Variable), size = Point_size, data = DATA) +
          scale_color_viridis_c() +
          guides(color = guide_legend(title = Color_by)) +
          theme(axis.title = element_blank(),
                axis.text = element_blank(),
                axis.ticks = element_blank(),
                axis.line = element_blank(),
                panel.grid = element_blank())
        plot(PLOT)
        return(PLOT)
      }

      #If is character
      else{
        PLOT <- ggplot() +
          annotation_raster(Image, xmin = 0, xmax = Image_width, ymin = 0, ymax = Image_height, interpolate = TRUE) +
          geom_point(aes(x = X, y = Y, color = Variable), size = Point_size, data = DATA) +
          scale_color_viridis_d() +
          guides(color = guide_legend(title = Color_by)) +
          theme(axis.title = element_blank(),
                axis.text = element_blank(),
                axis.ticks = element_blank(),
                axis.line = element_blank(),
                panel.grid = element_blank())
        plot(PLOT)
        return(PLOT)
      }
    }

    #If color is not required
    if(is.null(Color_by)){
      PLOT <- ggplot() +
        annotation_raster(Image, xmin = 0, xmax = Image_width, ymin = 0, ymax = Image_height, interpolate = TRUE) +
        geom_point(aes(x = X, y = Y), color = "grey", size = Point_size, data = DATA) +
        theme(axis.title = element_blank(),
              axis.text = element_blank(),
              axis.ticks = element_blank(),
              axis.line = element_blank(),
              panel.grid = element_blank())
      plot(PLOT)
      return(PLOT)
    }

  }
