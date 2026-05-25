#' Generates individual core images from TMA whole slide images
#'
#' The function uses downsized representations of whole slide TMA images to find core position and extract them from the original image (see details)
#'
#' @param N_cores Integer. Number of cores to parallelize your computation.
#' @param Label_text_size Numeric value indicating the size of the text of the labels (used to generate adequately sized images).
#'
#' @param Original_Image_Directory A character value containing the path to the original TMA image.
#' @param Downsized_Image A list created using the [Smart_image_importer()] function. It contains the downsized image representation of the TMA and the original pixel size.
#' @param Output_Directory A character value containing the path to the folder where core images will be saved.
#' @param RGB_Color_images A logical value indicating if TMA contains RGB color images (these will be processed in a different manner to grayscale multi-channel tif images).
#' @param Ordered_Channels (OPTIONAL) If image has multiple channel, the user can provide channel names. These names will be used with Channels_to_keep (see below) to filter the resulting image.
#' @param Channels_to_keep (OPTIONAL) If image has multiple channel, the user can provide channel names. These names will be used with Channels_to_keep to filter the resulting image.
#'
#' @param Threshold_type_tissueMask Type of threshold to performe tissue mask. Either 'Otsu', 'Arbitrary' or 'Absolute'.
#' @param Threshold_value_tissueMask Numeric value used if Arbitrary is the threshold type of choice.
#' @param Blurr_tissueMask Logical value indicating if image blurring be performed before tissue mask generation.
#' @param Sigma_tissueMask Numeric value indicating the sigma value to perform Gaussian blurring.
#'
#' @param Erode_kernel_size Erosion kernel size. Used if Erode_rounds is >= 1.
#' @param Erode_rounds Number of erosion rounds to be performed.
#' @param Dilate_kernel_size Dilation kernel size. Used if Dilate_rounds is >= 1.
#' @param Dilate_rounds Number of dilation rounds to be performed.
#' @param Min_size A integer value indicating the minimum size (in pixels) a TMA core must have to be considered a valid core.
#'
#' @param Watershed_tolerance A numeric value > 0 indicating the tolerance of the watershed algorithm. Higher values reduce the number of individual TMA cores identified.
#' @param Watershed_ext A numeric value indicating the extension watershed parameter. Higher values reduce the number of individual TMA cores identified.
#'
#' @param Tolerance_column A numeric value indicating the tolerance from perfect alignment to consider a new column being identified (see details).
#' @param Tolerance_row A numeric value indicating the tolerance from perfect alignment to consider a new row being identified (see details).
#'
#' @param Core_names_matrix (OPTIONAL BUT HIGHLY RECOMMENDED) A dataframe or matrix containing TMA core names. These names will be used to identify output core images.
#'
#'
#'
#' @returns The function exports individual core images identified to the indicated output directory.
#'
#' @details
#' The function works with downsized whole slide images generated using the [Smart_image_importer()]. First a tissue mask is generated to identify the foreground and background
#' pixels. Afterwards this mask is processed (hole filling, erosion and finally dilation). Then the mask is transformed to a distmap of the distance of every pixel to the closest
#' background pixel. Afterwards, this distmap y segmented using a watershed algorithm to identify individual tumor cores. The bounding box of these cores is used to extract
#' the core images from the originial TMA whole slide image.
#'
#' If the image provided is an RGB image, in order to generate the mask the image is first inverted.
#'
#' After identifying the position of the TMA cores, the number of rows and columns present in the TMA is computed by clustering the X and Y centroid coordinates of the cores.
#' If the tolerance values are high, higher rates of missalignment will be tolerated to consider new rows/columns.
#'
#' After identifying the position of TMA cores and the row and column they belong to, cores are labelled according to the user provided matrix.
#'
#'
#' @examples
#' \dontrun{
#'
#' TMA_image_downsized <-
#'     Smart_image_importer(Image_directory = "Path/to/original/TMA/Image",
#'                          Log10_pixel_output = 6
#'   )
#'
#' TMA_dearrayer_function(Original_Image_Directory = "Path/to/original/TMA/Image",
#'                       Downsized_Image = TMA_image_downsized,
#'                       Output_Directory = "Path/to/ouput/directory")
#'
#'}
#' @export

TMA_dearrayer_function <-
  function(
    Original_Image_Directory,
    Downsized_Image,
    Output_Directory,
    RGB_Color_images = FALSE,
    N_cores = 1,
    Label_text_size = 7,

    Ordered_Channels = NULL,
    Channels_to_keep = NULL,


    Threshold_type_tissueMask = "Arbitrary",
    Threshold_value_tissueMask = 0.1,
    Blurr_tissueMask = TRUE,
    Sigma_tissueMask = 5,

    Erode_kernel_size = 3,
    Erode_rounds = 0,
    Dilate_kernel_size = 3,
    Dilate_rounds = 0,
    Min_size = 100,

    Watershed_tolerance = 1,
    Watershed_ext = 1,

    Tolerance_column = 100,
    Tolerance_row = 100,

    Core_names_matrix = NULL

    ){

    #####WHAT TO DO ON EXIT#####
    on.exit({
      future::plan("future::sequential")
      gc()
    }
    )

    #####CHECK SUGGESTED PACKAGES####
    #Check installation of suggested packages
    {
      if(!requireNamespace("magick", quietly = FALSE)) stop(
        paste0("magick CRAN package is required to execute the function. Please install using the following code: ",
               expression(install.packages("magick"))
        )
      )
      if(!requireNamespace("RBioFormats", quietly = TRUE)) stop(
        paste0("RBioFormats Bioconductor package is required to execute the function. Please install using the following code: ",
               expression({
                 if (!require("BiocManager", quietly = TRUE))
                   install.packages("BiocManager")

                 BiocManager::install("RBioFormats")
               })
        )
      )
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

    ####Argument check####
    #Input and output dir and the Downsized image
    if(!file.exists(Original_Image_Directory)) stop(paste0("Issue with Original_Image_Directory directory: ", Original_Image_Directory))
    if(!dir.exists(Output_Directory)) stop(paste0("Issue with output directory: ", Output_Directory))
    if(!identical(names(Downsized_Image)[1:4], c("Image", "Original_Dims", "Current_Dims", "Reduction_factor"))) stop("Downsized_image must have been generated using the Smart_image_importer")
    if(!is.logical(RGB_Color_images)) stop("RGB_Color_images must be a logical value")

    #Check the cores
    if(!all(N_cores >= 1 & N_cores%%1 == 0)) stop("N_cores must be an integer value > 0")

    #Check tissue thresholding parameters
    if(!Threshold_type_tissueMask %in% c("Otsu", "Arbitrary", "Absolute")) stop("Threshold_type_TissueMask must be one of the following: Otsu, Arbitrary, Absolute")
    if(Threshold_type_tissueMask == "Arbitrary"){
      if(!all(is.numeric(Threshold_value_tissueMask), Threshold_value_tissueMask >=0, Threshold_value_tissueMask <= 1)) stop("Threshold_value_tissueMask must be a single numeric value between 0 and 1")
    }
    if(!is.logical(Blurr_tissueMask)) stop("Blurr_tissueMask must be a logical value")
    if(Blurr_tissueMask){
      if(!all(is.numeric(Sigma_tissueMask), Sigma_tissueMask > 0)) stop("Sigma_tissueMask must be a positive numeric value > 0")
    }

    #Check kernels for opening and closing
    if(!all(Erode_kernel_size%%1 == 0, Erode_kernel_size >= 1)) stop("Erode_kernel_size must be a positive integer value >= 1")
    if(!all(Erode_rounds%%1 == 0, Erode_kernel_size >= 0)) stop("Erode_rounds must be a positive integer value >= 0")
    if(!all(Dilate_kernel_size%%1 == 0, Dilate_kernel_size >= 1)) stop("Dilate_kernel_size must be a positive integer value >= 1")
    if(!all(Dilate_rounds%%1 == 0, Dilate_rounds >= 0)) stop("Dilate_rounds must be a positive integer value >= 0")

    #Check watershed arguments and size
    if(!all(is.numeric(Watershed_tolerance), Watershed_tolerance > 0)) stop("Watershed_tolerance must be a numeric value > 0")
    if(!all(is.numeric(Watershed_ext), Watershed_ext > 0)) stop("Watershed_ext must be a numeric value > 0")
    if(!all(is.numeric(Min_size), Min_size%%1 == 0, Min_size > 0)) stop("Min_size must be a positive integer value > 0")

    #Check tolerances
    if(!all(is.numeric(Tolerance_column), Tolerance_column >= 1, length(Tolerance_column) == 1)) stop("Tolerance_column must be a single numeric value larger than 0")
    if(!all(is.numeric(Tolerance_row), Tolerance_row >= 1, length(Tolerance_row) == 1)) stop("Tolerance_row must be a single numeric value larger than 0")

    #Core_names_matrix will be checked latter, just check if provided and if not print a message
    if(is.null(Core_names_matrix)){
      Colored_print("\n Core_names_matrix not provided. Assigning core names by row/col position. \n PROVIDING Core_names_matrix IS STRONGLY ENCOURAGED!!!!!", color = "red")
    }


    ####Define used functions####
    assign_irregular_grid <- function(x, y, tol_x, tol_y) {

      # columnas (X)
      hc_x <- hclust(dist(x), method = "complete")
      col_raw <- cutree(hc_x, h = tol_x)
      col_centers <- tapply(x, col_raw, median)
      col <- match(col_raw, order(col_centers))

      # filas (Y)
      hc_y <- hclust(dist(y), method = "complete")
      row_raw <- cutree(hc_y, h = tol_y)
      row_centers <- tapply(y, row_raw, median)
      row <- match(row_raw, order(row_centers))

      return(tibble::as_tibble(data.frame(row = row, col = col)))
    }

    ggplot_to_EBImage <- function(p) {
      tf <- tempfile(fileext = ".png")

      #scale it to 1
      ggplot2::ggsave(tf, p, scale = 1)

      img <- EBImage::readImage(tf)
      unlink(tf)

      img
    }

    flip_order <- function(x,
                           what = c("rows", "cols", "both")) {
      what <- match.arg(what)

      if (what %in% c("rows", "both")) {
        x <- x[rev(seq_len(nrow(x))), , drop = FALSE]
      }
      if (what %in% c("cols", "both")) {
        x <- x[, rev(seq_len(ncol(x))), drop = FALSE]
      }
      return(x)
    }


    ####Obtain the image####
    Image_downsized <- Downsized_Image

    if(!RGB_Color_images) Image_for_mask <- Image_downsized$Image %>% magick::as_EBImage()
    if(RGB_Color_images){
      Image_for_mask <- Image_downsized$Image %>% magick::image_negate()  %>% magick::as_EBImage()
      EBImage::colorMode(Image_for_mask) <- "Grayscale"
    }

    ####Generate the core segmentation mask and ask user####
    cat("\n Generating TMA core mask")

    #Generate the basic tissue mask
    Tissue_mask <-
      Tissue_mask_generator(
        Image = Image_for_mask,
        Threshold_type = Threshold_type_tissueMask,
        Threshold_value = Threshold_value_tissueMask,
        Blurr = Blurr_tissueMask,
        Sigma = Sigma_tissueMask
        )

    #Fill holes at the start of the process
    Tissue_mask <- EBImage::fillHull(Tissue_mask)

    #Erode the number of rounds required
    if(Erode_rounds > 0){
      Tissue_mask <- purrr::reduce(
        .x = seq_along(1:Erode_rounds),
        .f = function(img, ...) EBImage::erode(img, kern = EBImage::makeBrush(Erode_kernel_size, "disc")),
        .init = Tissue_mask
      )
    }


    #Dilate the number of rounds required
    if(Dilate_rounds > 0){
      Tissue_mask <- purrr::reduce(
        .x = seq_along(1:Dilate_rounds),
        .f = function(img, ...) EBImage::dilate(img, kern = EBImage::makeBrush(Dilate_kernel_size, "disc")),
        .init = Tissue_mask
      )
    }

    #Fill holes at the end of the process
    Tissue_mask <- EBImage::fillHull(Tissue_mask)

    #Turn to distmap to closest background pixel
    Tissue_mask <- EBImage::distmap(Tissue_mask)
    Tissue_mask <- EBImage::watershed(EBImage::normalize(Tissue_mask), tolerance = Watershed_tolerance, ext = Watershed_ext)


    #Plot the result
    Initial_plot <-
      ggplot() +
      annotation_raster(EBImage::colorLabels(EBImage::flip(Tissue_mask)), xmin = 1, xmax = Image_downsized$Current_Dims[1],  ymin = 1, ymax = Image_downsized$Current_Dims[2]) +
      scale_x_continuous("", limits = c(0, Image_downsized$Current_Dims[1]+1)) +
      scale_y_continuous("", limits = c(0, ymax = Image_downsized$Current_Dims[2]+1)) +
      theme_minimal()+
      coord_cartesian(expand = FALSE) +

      theme(axis.text = element_blank(),
            axis.title = element_blank(),
            axis.ticks = element_blank())

    #Display the image
    Initial_plot <-  ggplot_to_EBImage(Initial_plot)
    print(magick::image_read(Initial_plot))

    #Ask the user if the process should proceed
    Proceed <- menu(choices = c("Proceed", "Abort"), title = paste0("Should the core extraction process proceed?"))
    if(Proceed == 2){
      stop("Function aborted. Please review mask generation parameters")
    }


    ####Generate the coordinates tibble and ask user if OK####
    #Obtain mask data (centroids and min and max X Y coords of box)
    cat("\n Computing bounding box of TMA cores")
    Look_up_table <- as_tibble(EBImage::computeFeatures.moment(Tissue_mask))

    #Get only the X and Y centroid coordinates
    Look_up_table <- dplyr::bind_cols(tibble::tibble(TMA_core = stringr::str_c("Core_", 1:nrow(Look_up_table))),
                                      Look_up_table[1:2])
    Look_up_table$Object_size <- table(Tissue_mask)[-1]
    Look_up_table <- Look_up_table %>% dplyr::filter(Object_size >= Min_size)
    #Obtain the bounding box coordinates
    Look_up_table <- dplyr::bind_cols(Look_up_table,
                                     purrr::map_dfr(1:nrow(Look_up_table), function(Mask_index){
                                       coords <- which(Tissue_mask == Mask_index, arr.ind = TRUE)
                                       bbox <- c(
                                         xmin = min(coords[, 1]),
                                         xmax = max(coords[, 1]),
                                         ymin = min(coords[, 2]),
                                         ymax = max(coords[, 2])
                                       )
                                       return(bbox)
                                     })
    )
    #Change centroid names
    names(Look_up_table)[2:3] <- c("X_centroid", "Y_centroid")

    #Add the row and col number according to the clustering
    Look_up_table <- dplyr::bind_cols(Look_up_table,
                                     assign_irregular_grid(x = Look_up_table$X_centroid, y = Look_up_table$Y_centroid, tol_x = Tolerance_column, tol_y = Tolerance_row)
    )
    #Arrange by row (important if Core_names_matrix is provided)
    Look_up_table <- Look_up_table %>% dplyr::arrange(row, col)
    Look_up_table <- Look_up_table %>% dplyr::mutate(Core_basic_label = str_c("Row_", row, "-Col_", col))

    #Plot the result
    Initial_plot <-
      ggplot() +
      annotation_raster(Image_downsized$Image %>% magick::image_flip(), xmin = 1, xmax = Image_downsized$Current_Dims[1],  ymin = 1, ymax = Image_downsized$Current_Dims[2]) +
      scale_x_continuous("", limits = c(0, Image_downsized$Current_Dims[1]+1)) +
      scale_y_continuous("", limits = c(0, ymax = Image_downsized$Current_Dims[2]+1)) +
      theme_minimal()+
      coord_cartesian(expand = FALSE) +

      geom_rect(
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
        fill = NA,
        colour = "red",
        linewidth = 1.2,
        data = Look_up_table
      ) +

      geom_text(aes(label = Core_basic_label, x = X_centroid, y = Y_centroid), size = Label_text_size, fontface = "bold", color = "red",
                data = Look_up_table) +

      theme(axis.text = element_blank(),
            axis.title = element_blank(),
            axis.ticks = element_blank())

    #Display the image
    Initial_plot <-  ggplot_to_EBImage(Initial_plot)
    print(magick::image_read(Initial_plot))

    #Ask the user if the process should proceed
    Proceed <- menu(choices = c("Proceed", "Abort"), title = paste0("Should the core extraction process proceed?"))
    if(Proceed == 2){
      stop("Function aborted. Please review mask generation parameters")
    }

    ####Process Core_names_matrix####
    if(!is.null(Core_names_matrix)){

      cat("\n Adding core names provided by user")

      if(any(duplicated(na.omit(as.vector(Core_names_matrix))))) stop("Names provided in Core_names_matrix must be unique with the exception of NA values")

      #Get the initial matrix
      Label_matrix <- as.matrix(Core_names_matrix)

      N_rows_Core_names <- nrow(Label_matrix)
      N_cols_Core_names <- ncol(Label_matrix)

      #Check that number of rows and number of columns match the columns and rows identified
      if(!any(N_rows_Core_names == max(Look_up_table$row) & N_cols_Core_names == max(Look_up_table$col),
              N_rows_Core_names == max(Look_up_table$col) & N_cols_Core_names == max(Look_up_table$row))
         ) stop(paste0("Number of rows and columns provided in the Core_names_matrix do not match the ones identified", "\n",
                       "Rows/columns identified: ", max(Look_up_table$row), "/", max(Look_up_table$col), "\n",
                       "Rows/columns provided in Core_names_matrix: ", N_rows_Core_names, "/", N_cols_Core_names))

      #If the rows and columns are transversed then transverse the matrix as well
      if(N_rows_Core_names == max(Look_up_table$col) & N_cols_Core_names == max(Look_up_table$row)) Label_matrix <- t(Label_matrix)

      #Add names according to column and matrix number
      colnames(Label_matrix) <- 1:ncol(Label_matrix)
      rownames(Label_matrix) <- 1:nrow(Label_matrix)

      #Generate a tibble of the TMA names (for plotting purposes)
      #Obtain the Row and column centroids from the basic tibble
      Average_Y_centroid <- Look_up_table %>% dplyr::group_by(row) %>% dplyr::summarise(Y_centroid = mean(Y_centroid))
      Average_X_centroid <- Look_up_table %>% dplyr::group_by(col) %>% dplyr::summarise(X_centroid = mean(X_centroid))


      #Enter the loop
      Not_ready_to_exit <- TRUE

      while(Not_ready_to_exit){
        #Generate the final tibble containing information for all possible rows and columns
        Merged_label_tibble <- tidyr::expand_grid(row = as.integer(unique(Look_up_table$row)), col = as.integer(unique(Look_up_table$col)))


        #Generate the basic core names tibble
        TMA_label_tibble <- tibble::as_tibble(Label_matrix) %>% dplyr::mutate(row = 1:nrow(Label_matrix)) %>% tidyr::pivot_longer(1:ncol(Label_matrix))
        names(TMA_label_tibble)[2] <- "col" #Change the name column by col
        TMA_label_tibble$row <- as.integer(TMA_label_tibble$row) #Guarantee they are integers
        TMA_label_tibble$col <- as.integer(TMA_label_tibble$col)#Guarantee they are integers


        #Bind the label tibble to the provided labels (by row and column position)
        Merged_label_tibble <- dplyr::left_join(Merged_label_tibble, TMA_label_tibble, dplyr::join_by(row == row, col == col))

        #Bind the tibble to the actual TMA_core IDs present in the look up table
        Merged_label_tibble <- dplyr::left_join(Merged_label_tibble,
                                                Look_up_table %>% dplyr::select(TMA_core, X_centroid, Y_centroid, row, col),
                                                dplyr::join_by(row == row, col == col))

        #If any value is empty then substitute by "Empty"
        Merged_label_tibble$value[is.na(Merged_label_tibble$value)] <- "Empty"


        #For loop to substitute X_centroids and Y centroids
        if(any(is.na(Merged_label_tibble$X_centroid))){
          for(i in which(is.na(Merged_label_tibble$X_centroid))){
            Merged_label_tibble$X_centroid[i] <- Average_X_centroid$X_centroid[Average_X_centroid$col == Merged_label_tibble$col[i]]
            Merged_label_tibble$Y_centroid[i] <- Average_Y_centroid$Y_centroid[Average_Y_centroid$row == Merged_label_tibble$row[i]]
          }
        }

        #Generate the plot
        Initial_plot <-
          ggplot() +
          annotation_raster(Image_downsized$Image %>% magick::image_flip(), xmin = 1, xmax = Image_downsized$Current_Dims[1],  ymin = 1, ymax = Image_downsized$Current_Dims[2]) +

          scale_x_continuous("", limits = c(0, Image_downsized$Current_Dims[1]+1)) +
          scale_y_continuous("", limits = c(0, ymax = Image_downsized$Current_Dims[2]+1)) +
          theme_minimal()+
          coord_cartesian(expand = FALSE) +

          geom_rect(
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = NA,
            colour = "red",
            linewidth = 1.2,
            data = Look_up_table
          ) +

          geom_text(aes(label = value, x = X_centroid, y = Y_centroid), size = Label_text_size, fontface = "bold", color = "red",
                    data = Merged_label_tibble) +

          theme(axis.text = element_blank(),
                axis.title = element_blank(),
                panel.grid = element_blank(),
                axis.ticks = element_blank(),
                margins = margin(0,0,0,0))

        #Display the image
        Initial_plot <-  ggplot_to_EBImage(Initial_plot)
        print(magick::image_read(Initial_plot))

        #Ask the user if the process should proceed if required modify the matrix
        Proceed <- menu(choices = c("Proceed", "Abort", "Flip rows", "Flip cols", "Flip both", "Transverse"), title = paste0("Should the core extraction process proceed?"))
        if(Proceed == 1) break
        if(Proceed == 2) stop("Function aborted. Please review mask generation parameters")
        if(Proceed == 3) Label_matrix <- flip_order(Label_matrix, "rows")
        if(Proceed == 4){
          Label_matrix <- flip_order(Label_matrix, "cols")
          colnames(Label_matrix) <- 1:ncol(Label_matrix) #To account for name of cols to re-generate the Merged tibble
        }
        if(Proceed == 5){
          Label_matrix <- flip_order(Label_matrix, "both")
          colnames(Label_matrix) <- 1:ncol(Label_matrix) #To account for name of cols to re-generate the Merged tibble
        }
        if(Proceed == 6) Label_matrix <- t(Label_matrix)
      }

      #Add the labels to the Look_up_table
      #Add the names using the TMA_core ID column present in both tibbles
      Look_up_table <- dplyr::left_join(Look_up_table, Merged_label_tibble %>% dplyr::select(TMA_core, value), by = "TMA_core") %>% dplyr::rename("TMA_core_name" = "value")

      }

    ####Extract cores####

    cat("\n Obtaining individual tumor cores from original image and saving results")

    #Obtain the TMA name from the Original_Image_Directory
    TMA_image_name <- stringr::str_extract(Original_Image_Directory, "(?<=^|[\\\\/])[^\\\\/.]+(?=\\.[^.]+$)")

    #Obtain the original coordinates
    Multiplication_factor <- Image_downsized$Original_Dims / Image_downsized$Current_Dims

    Look_up_table$Original_xmin <- floor(((Look_up_table$xmin - 1) * Multiplication_factor[1]) + 1)
    Look_up_table$Original_xmax <- ceiling(((Look_up_table$xmax + 1) * Multiplication_factor[1]) - 1)
    Look_up_table$Original_ymin <- floor(((Look_up_table$ymin - 1) * Multiplication_factor[2]) + 1)
    Look_up_table$Original_ymax <- ceiling(((Look_up_table$ymax + 1) * Multiplication_factor[2]) - 1)

    #If any coordinate is beyond the max X or Y correct
    Look_up_table$Original_xmax[Look_up_table$Original_xmax > Image_downsized$Original_Dims[1]] <- Image_downsized$Original_Dims[1]
    Look_up_table$Original_ymax[Look_up_table$Original_ymax > Image_downsized$Original_Dims[2]] <- Image_downsized$Original_Dims[2]


    #Plan the multicore session
    future::plan("future::multisession", workers = N_cores)
    options(future.globals.maxSize = Inf, future.rng.onMisuse = "ignore")
    furrr::furrr_options(scheduling = Inf)

    #Extract the actual cores
    furrr::future_walk(1:nrow(Look_up_table), function(tibble_row_index){

      #Generate a list for the Image name
      TMA_core_name <- list(TMA_name = TMA_image_name)
      if(!is.null(Core_names_matrix)) TMA_core_name$Core_name <- Look_up_table$TMA_core_name[tibble_row_index]
      Row_index <- Look_up_table$row[tibble_row_index]
      Col_index <- Look_up_table$col[tibble_row_index]
      TMA_core_name$ROW_COL <- stringr::str_c("Row", Row_index, "Col", Col_index)


      #Collapse the image name list
      Image_name <- stringr::str_c(TMA_core_name, collapse = "_")

      #Add the extension
      Image_name <- stringr::str_c(Image_name, ".tif")


      #Add the path
      Image_name <- stringr::str_c(Output_Directory, "/", Image_name)

      #Obtain the coordinates to be cropped
      Subset_xmin <- Look_up_table$Original_xmin[tibble_row_index]
      Subset_xmax <- Look_up_table$Original_xmax[tibble_row_index]
      Subset_ymin <- Look_up_table$Original_ymin[tibble_row_index]
      Subset_ymax <- Look_up_table$Original_ymax[tibble_row_index]

      TMA_image <-
        Smart_image_importer(
          N_cores = 1,
          Image_directory = Original_Image_Directory,
          Log10_pixel_output = 10, #An incredibly high number of pixels

          X_crop_coordinates = c(Subset_xmin, Subset_xmax),

          Y_crop_coordinates = c(Subset_ymin, Subset_ymax),

          Ordered_Channels = Ordered_Channels,
          Channels_to_keep = Channels_to_keep,
          Save_processed_images = FALSE
        )

      EBImage::writeImage(TMA_image$Image %>% magick::as_EBImage(), Image_name,
                          bits.per.sample = 16, compression = "LZW")
      cat(paste0("\n", Image_name))
    },
    .progress = TRUE)

    future::plan("future::sequential")
    gc()

    cat("DONE!")

  }
