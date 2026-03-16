#' Calculates which cellular or pixel elements are in the path from COO to target cells
#'
#' Given a Cell of Origin (COO) and its closest target cells, the function will compute the space polygon between them. Then cells or pixels inside this polygon are quantified.
#'
#'
#' @param N_cores Integer. Number of cores to parallelize your computation.
#'
#' @param DATA A dataframe or tibble containing a column named 'Phenotype' containing cell phenotype labels.
#' @param Cell_Of_Origin A character value indicating the cell phenotype label of the Cell of Origin.
#' @param Target_Cell A character value indicating the cell phenotype label of the Target cell.
#'
#' @param Barrier_Cells A character vector containing the names of the cell phenotypes that will be measured.
#' @param Image_parameter_list A list of lists containing information of images being used for pixel quantification (see details).
#' @param Pixel_distance_ratio (OPTIONAL if pixel quantification is required) A numeric value indicating the ratio between pixel size / cell coordinates. Use this argument when cell coordinates have been transformed to a distance unit like microns.
#'
#' @param Distance_sampled The maximum distance allowed to consider a target cell in proximity to a COO.
#' @param Perform_edge_correction A logical value indicating if COO close to the edge of the tissue should be removed from the analysis. If TRUE, COO that are less than the Distance_sampled to the edge will be removed.
#' @param Hull_ratio If edge correction needs to be performed, a numeric value indicating the hull ratio. Smaller values calculate more precise edge silhouettes at the cost of being more computationally demanding.
#' @param Polygon_angle A numeric value indicating the angle used to triangulate COO-target space. Higher values will result in coarser space polygons.
#' @param Target_cell_Angle_tolerance (OPTIONAL) If not NULL, a numeric value indicating the amount of angular distance required to keep two target cells. If below the tolerance value, only the closest target cell to the COO will be kept in the analysis.
#'
#'
#' @details
#' If image pixels need to be quantified the Image_parameter_list argument must be provided. Image_parameter_list must be a named list of lists.
#' Each named element must be a list containing the following items:
#' \itemize{
#' \item{Directory: Character specifying the path to the folder where images are present.}
#' \item{Channel: Either NULL (for single channel images) or a integer value indicating the channel number to be analyzed.}
#' \item{Image_rotate: Either NULL or a integer value indicating the degrees of rotation of the image. If NULL no rotation is applied.}
#' \item{Image_x_flip: A logical value indicating if X image flip should be performed.}
#' \item{Image_y_flip: A logical value indicating if X image flip should be performed.}
#' }
#'
#' The function returns the absolute cell and pixel counts within the barrier. The function also computes the ratios normalized by
#' the total number of target cells studied, the area of the barrier polygon and the total distance of the COO to the target cells studied.
#'
#' @returns A tibble containing information of barrier cells and pixels per COO analyzed.
#'
#' @seealso [Barrier_effect_analyzer()], [Barrier_graph_maker()]
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

Barrier_effect_calculator <-
  function(N_cores = 1,

           DATA = NULL,
           Cell_Of_Origin = NULL,
           Target_Cell = NULL,

           Barrier_Cells = NULL,
           Image_parameter_list = NULL,
           Pixel_distance_ratio = NULL,


           Distance_sampled = NULL,
           Perform_edge_correction = FALSE,
           Hull_ratio = NULL,
           Polygon_angle = 10,
           Target_cell_Angle_tolerance = 1
  ){

    ########Check suggested packages########
    {
      if(!requireNamespace("sf", quietly = FALSE)) stop(
        paste0("sf CRAN package is required to execute the function. Please install using the following code: ",
               expression(install.packages("sf")))
      )

      if(!requireNamespace("rtree", quietly = FALSE)) stop(
        paste0("rtree GitHub package is required to execute the function. Please install using the following code: ",
               expression(remotes::install_github("akoyabio/rtree")))
      )

      if(!is.null(Image_parameter_list)){
        if(!requireNamespace("magick", quietly = FALSE)) stop(
          paste0("magick CRAN package is required to execute the function. Please install using the following code: ",
                 expression(install.packages("magick")))
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
    }

    ########Check arguments########
    #Obtain the data
    DATA_Phenotypes <- DATA

    if(!identical(names(DATA_Phenotypes)[1:4],  c("Cell_no", "X", "Y", "Subject_Names"))) { #Check if Data is correctly formatted
      stop("DATA provided should have an adecuate format")
    }
    if(!("Phenotype" %in% names(DATA_Phenotypes))) {
      stop("DATA should contain a column named Phenotype specifying the cell types")
    }
    if(!all(c(Cell_Of_Origin %in% unique(DATA_Phenotypes$Phenotype),
              Target_Cell %in% unique(DATA_Phenotypes$Phenotype)
    )
    )) { #Check if provided cell types are in the phenotype variable
      stop(paste0("Cell of origin and Target cell provided should be one of: ", stringr::str_c(unique(DATA_Phenotypes$Phenotype), collapse = ", ")))
    }
    #Check that COO and target cells are single element vectors
    if(any(length(Cell_Of_Origin) > 1, length(Target_Cell) > 1)) stop("Only one type of cells should be provided for Cell_Of_Origin and Target_Cell")

    #Check barrier cells if required
    if(!is.null(Barrier_Cells)){
      #COO and target cells do not belong to barrier cells
      if(any(Barrier_Cells %in% c(Cell_Of_Origin, Target_Cell))) stop("Barrier_Cells cannot be the same as Cell_Of_Origin or Target_Cell")
      if(!all(Barrier_Cells %in% unique(DATA_Phenotypes$Phenotype))) stop(paste0("Barrier_Cells provided should be any of: ", stringr::str_c(unique(DATA_Phenotypes$Phenotype), collapse = ", ")))
    }

    if(!all(N_cores >= 1 & N_cores%%1 == 0)) stop("N_cores must be an integer value > 0")

    if(!all(is.numeric(Distance_sampled), Distance_sampled > 0)) stop("Distance_sampled must be a numeric value > 0")

    if(!is.logical(Perform_edge_correction)) stop("Perform_edge_correction must be a logical value")
    if(Perform_edge_correction){
      if(!all(is.numeric(Hull_ratio), Hull_ratio >= 0, Hull_ratio <= 1)) stop("Hull_ratio must be a numeric value between 0 and 1")
    }
    if(!all(is.numeric(Polygon_angle), Polygon_angle > 0, Polygon_angle <= 90)) stop("Polygon_angle must be a numeric value between 1 and 90")
    if(!any(is.null(Target_cell_Angle_tolerance),
            all(is.numeric(Target_cell_Angle_tolerance), Target_cell_Angle_tolerance > 0, Target_cell_Angle_tolerance <= 90))
    ) stop("Target_cell_Angle_tolerance must be either NULL or a numeric value between 0 and 90")

    if(!is.null(Image_parameter_list)){
      if(!is.list(Image_parameter_list)) stop("Image_parameter_list must be a list")
      if(any(names(Image_parameter_list) == "")) stop("Lists in Image_parameter_list must be named")
      if(length(unique(names(Image_parameter_list))) != length(Image_parameter_list)) stop("Names in Image_parameter_list must be unique")
      if(any(!purrr::map_lgl(Image_parameter_list, ~is.list(.)))) stop("Image_paremeter_list must contain lists")
      if(any(!purrr::map_lgl(Image_parameter_list, function(Image_parameter){
        identical(names(Image_parameter), c("Directory", "Channel", "Image_rotate", "Image_x_flip", "Image_y_flip"))
      }))) stop("All lists in Image_paremeter_list must contain the following named items: Directory, Channel, Image_rotate, Image_x_flip, Image_y_flip")
      if(!any(is.null(Pixel_distance_ratio),
              all(is.numeric(Pixel_distance_ratio), Pixel_distance_ratio > 0, length(Pixel_distance_ratio) == 1))
         ) stop("Pixel_distance_ratio must be NULL or a numeric value > 0")

      #Check the individual image parameters
      purrr::map(names(Image_parameter_list), function(Image_parameter){
        if(!length(dir(Image_parameter_list[[Image_parameter]][["Directory"]])) >= 1) stop(paste0(Image_parameter, ": No files found at the directory provided. Please check the path."))
        if(!any(is.null(Image_parameter_list[[Image_parameter]][["Channel"]]),
                all(is.numeric(Image_parameter_list[[Image_parameter]][["Channel"]]), Image_parameter_list[[Image_parameter]][["Channel"]]%%1 == 0, Image_parameter_list[[Image_parameter]][["Channel"]] >= 1))
           ) stop(paste0(Image_parameter, ": Channel must be either NULL or an integer value larger than 0"))
        if(!any(is.null(Image_parameter_list[[Image_parameter]][["Image_rotate"]]),
                all(is.numeric(Image_parameter_list[[Image_parameter]][["Image_rotate"]]),
                    Image_parameter_list[[Image_parameter]][["Image_rotate"]] >= 0,
                    Image_parameter_list[[Image_parameter]][["Image_rotate"]]<= 360,
                    Image_parameter_list[[Image_parameter]][["Image_rotate"]]%%1 == 0))
           )stop(paste0(Image_parameter, ": Image_rotate must be either NULL or a integer value between 0 and 360"))
        if(!is.logical(Image_parameter_list[[Image_parameter]][["Image_x_flip"]])) stop(paste0(Image_parameter, ": Image_x_flip must be a logical value"))
        if(!is.logical(Image_parameter_list[[Image_parameter]][["Image_y_flip"]])) stop(paste0(Image_parameter, ": Image_y_flip must be a logical value"))
      })
    }

    #Check that the user has provided either Image_parameter_list or Barrier_Cells
    if(!any(!is.null(Barrier_Cells), !is.null(Image_parameter_list))) stop("Either Barrier_Cells or Image_parameter_list should be provided")

    #########Pixel to distance ratio correction (if required)#########
    if(all(!is.null(Image_parameter_list), !is.null(Pixel_distance_ratio))){
      #Modify X Y coordinates
      DATA_Phenotypes <- DATA_Phenotypes %>% dplyr::mutate(X = X*Pixel_distance_ratio,
                                                           Y = Y*Pixel_distance_ratio)
      #Modify search space distance
      Distance_sampled <- Distance_sampled*Pixel_distance_ratio
    }

    #########Edge correction (if required)########
    if(Perform_edge_correction){
      print("Running edge correction example on a random sample")
      Sample <- DATA_Phenotypes %>% dplyr::filter(Subject_Names == sample(unique(DATA_Phenotypes$Subject_Names), size = 1))
      Cells_sf <- sf::st_as_sf(Sample , coords = c("X", "Y"))
      Edge_line <- sf::st_cast((Cells_sf %>% summarise() %>% sf::st_concave_hull(ratio = Hull_ratio) %>% summarise), "LINESTRING")
      Cells_in_Border_vector <- unlist(sf::st_is_within_distance(Cells_sf, Edge_line, sparse = F, dist = Distance_sampled))

      plot(Sample %>%
             dplyr::mutate(Removed = Cells_in_Border_vector) %>%
             ggplot(aes(x = X, y = Y, color = Cells_in_Border_vector)) +
             geom_point() +
             scale_color_manual("", labels = c("Included", "Removed"), values = c("black", "grey")) +
             theme_minimal() +
             scale_x_continuous("") +
             scale_y_continuous("") +
             theme(panel.grid = element_blank(),
                   axis.text = element_blank(),
                   legend.position = "bottom",
                   legend.text = element_text(size = 12)))

      #Ask the user if the algorihtm should proceed
      answer <- menu(c("Proceed", "Abort"), title = "Check edge correction results. Should the analysis proceed?")
      #If user decides to stop then abort function and return stop message
      if(answer == 2) stop("The function has been stopped. Please tune edge correction parameters for a better result")

      #Remove cells
      Keep_vector <- unlist(
        purrr::map(unique(DATA_Phenotypes$Subject_Names), function(x){
          #Prepare our data
          Image_tibble <- DATA_Phenotypes %>% dplyr::filter(Subject_Names == x)
          Cells_sf <- sf::st_as_sf(Image_tibble , coords = c("X", "Y"))
          Edge_line <- sf::st_cast((Cells_sf %>% summarise() %>% sf::st_concave_hull(ratio = Hull_ratio) %>% summarise), "LINESTRING")
          Cells_in_Border_vector <- unlist(sf::st_is_within_distance(Cells_sf, Edge_line, sparse = F, dist = Distance_sampled))
          #Calculate cells in border
          COO_in_Border_vector <- Image_tibble$Phenotype == Cell_Of_Origin & Cells_in_Border_vector

          #Print message to warn COO removed in analysis
          COO_in_Border_N <- sum(COO_in_Border_vector)
          COO_total_N <- sum(Image_tibble$Phenotype == Cell_Of_Origin)

          if(all(COO_in_Border_N == COO_total_N, COO_total_N > 0)){
            Colored_print(text = paste0("\n",
                                        "Sample ", as.character(x), ": ", COO_in_Border_N, " / ", COO_total_N, " ", Cell_Of_Origin, " cell/s will be removed due to edge proximity."),
                          color = "red")
          }
          if(COO_in_Border_N < COO_total_N){
            Colored_print(text = paste0("\n",
                                        "Sample ", as.character(x), ": ", COO_in_Border_N, " / ", COO_total_N, " ", Cell_Of_Origin, " cell/s will be removed due to edge proximity."),
                          color = "green")
          }

          #Return a vector with the cells to keep (either no COO or COO not in edge)
          return(!COO_in_Border_vector)
        })
      )

      #Remove the COO that are not required
      DATA_Phenotypes <- DATA_Phenotypes[Keep_vector,]
    }

    ########Sample filtering########
    #Prepare the look up table (Subject_Names, COO and Targets, photos...)
    #Filter the dataset to contain only cells of interest
    print("Removing samples without COO or target cells")
    DATA_Phenotypes <- DATA_Phenotypes %>% dplyr::filter(Phenotype %in% c(Cell_Of_Origin, Target_Cell, Barrier_Cells))

    Lookup_table <-
      DATA_Phenotypes %>% dplyr::group_by(Subject_Names) %>%
      dplyr::count(Phenotype) %>% dplyr::ungroup() %>%
      tidyr::pivot_wider(names_from = Phenotype, values_from = n)
    Lookup_table[is.na(Lookup_table)] <- 0

    #Remove samples without Cell of origin
    if(any(Lookup_table[Cell_Of_Origin] == 0)){
      if(!any(Lookup_table[Cell_Of_Origin] > 0)) stop("All samples are devoid of Cells_Of_Origin after Pre-processing, please refine edge correction parameters.")

      message(paste0("The following samples are devoid of cells of origin and will be removed: ", "\n",
                     stringr::str_c(Lookup_table$Subject_Names[Lookup_table[Cell_Of_Origin] == 0], collapse = "\n"))
      )
      Lookup_table <- Lookup_table[!Lookup_table[Cell_Of_Origin] == 0, ]
    }
    #Remove samples without Target_cells
    if(any(Lookup_table[Target_Cell] == 0)){
      if(!any(Lookup_table[Target_Cell] > 0)) stop("After removing samples without cells of origin or target cells no samples remained for analysis.")

      message(paste0("The following samples are devoid of target cells and will be removed: ", "\n",
                     stringr::str_c(Lookup_table$Subject_Names[Lookup_table[Target_Cell] == 0], collapse = "\n"))
      )
      Lookup_table <- Lookup_table[!Lookup_table[Target_Cell] == 0, ]
    }

    ########Match Subject_Names to images (if required)########
    if(!is.null(Image_parameter_list)){
      #Now pair each Subject_Name in the dataset with an image in each of the directories provided
      #Iterate over the Image_parameter_list
      Image_Parameters <-
        purrr::map_dfc(names(Image_parameter_list), function(Image_type_name){
          Parameter_list <- Image_parameter_list[[Image_type_name]]
          Image_names_FULL <- dir(Parameter_list[["Directory"]], full.names = TRUE)
          Image_names_SHORT<- dir(Parameter_list[["Directory"]], full.names = FALSE)

          #Find if any of the image in the directory is a tissue mask generated using the Pixel_thresholder and remove if necessary
          Image_names_FULL <- Image_names_FULL[!grepl("_Tissue_mask\\.tiff$", Image_names_FULL)]
          Image_names_SHORT <- Image_names_SHORT[!grepl("_Tissue_mask\\.tiff$", Image_names_SHORT)]

          #Find if image is binary or mnultilevel
          Image_type <- 'Original'
          if(any(grepl("_BinaryThresholded\\.tiff$", Image_names_SHORT))) Image_type <- "Binary"
          if(any(grepl("_MultiThresholded.tiff", Image_names_SHORT))) Image_type <- "Multi"

          #Find which images match with the Subject_Names column of the Lookup_table
          #For every Subject_Names in the Lookup_table find the closest name in the Image_names_SHORT
          Image_name_order <- purrr::map_dbl(Lookup_table$Subject_Names, function(Image_name){
            which.min(adist(Image_name, Image_names_SHORT, fixed = FALSE, ignore.case = FALSE))
          })

          #Build the dataframe
          Image_tibble <-
            tibble::tibble(Image_type = Image_type,
                           Image_name = Image_names_SHORT[Image_name_order],
                           Image_path = Image_names_FULL[Image_name_order],
                           Image_channel = Parameter_list[["Channel"]],
                           Image_rotate = Parameter_list[["Image_rotate"]],
                           Image_x_flip = Parameter_list[["Image_x_flip"]],
                           Image_y_flip = Parameter_list[["Image_y_flip"]])
          names(Image_tibble) <- stringr::str_c(Image_type_name, names(Image_tibble), sep = "_")
          return(Image_tibble)
        })

      Lookup_table <- dplyr::bind_cols(Lookup_table, Image_Parameters)

      #Check that Image names are unique and that no one is repeated
      Conflictive_image_directories <- !(Lookup_table %>% dplyr::select(dplyr::ends_with("_Image_name")) %>%
                                           purrr::map_lgl(function(Image_names) length(unique(Image_names)) == nrow(Lookup_table))
      )
      #If any is repeated then find the conflictive cases and STOP the function
      if(any(Conflictive_image_directories)){
        purrr::walk(names(Conflictive_image_directories)[Conflictive_image_directories],
                    function(Image_directory){
                      Interim <- Lookup_table %>% dplyr::select(Subject_Names, dplyr::all_of(Image_directory))
                      Conflictive_images <- unique(Interim[[Image_directory]][duplicated(Interim[[Image_directory]])])

                      Interim <- Interim[Interim[[Image_directory]] %in% Conflictive_images, ]

                      print("One or more images are matched to multiple Subject_Names. Please solve before proceeding.")
                      print(Interim[order(Interim[[Image_directory]]),], n = Inf)
                    })
        stop("Please solve conflicts before proceeding")
      }

      #If everything is matched, run a quick test to check adequate cell to image information alignment
      purrr::walk(names(Image_parameter_list), function(Image_type_name){

        #Generate the unique look up table
        Lookup_table <- Lookup_table %>% dplyr::select(Subject_Names, dplyr::starts_with(stringr::str_c(Image_type_name, "_", collapse = "")))


        #To enter the while loop set the answer to 3
        Answer <- 3
        while(Answer == 3){
          #Obtain a random image
          Sample_DATA <- Lookup_table %>% dplyr::slice_sample(n = 1)

          #Obtain parameters
          Name_of_Subject <- unlist(Sample_DATA$Subject_Names)
          Image_path <- unlist(Sample_DATA %>% dplyr::select(dplyr::ends_with("_Image_path")))
          Image_channel <- unlist(Sample_DATA %>% dplyr::select(dplyr::ends_with("_Image_channel")))
          Image_rotate <- unlist(Sample_DATA %>% dplyr::select(dplyr::ends_with("_Image_rotate")))
          Image_x_flip <- unlist(Sample_DATA %>% dplyr::select(dplyr::ends_with("_Image_x_flip")))
          Image_y_flip <- unlist(Sample_DATA %>% dplyr::select(dplyr::ends_with("_Image_y_flip")))

          #Read image
          if(is.null(Image_channel)) Image <- magick::image_read(Image_path)
          if(!is.null(Image_channel)) Image <- magick::image_read(Image_path)[Image_channel]

          #Execute rotation if required
          if(!is.null(Image_rotate)) Image <- Image %>% magick::image_rotate(degrees = Image_rotate)
          #Execute x or y flip if required
          if(Image_x_flip) Image <- Image %>% magick::image_flop()
          if(Image_y_flip) Image <- Image %>% magick::image_flip() #Here we will flip the image as required by user as annotation_raster will flip it for plotting purposes

          #obtain the max width and height
          X_max <- magick::image_info(Image)$width
          Y_max <- magick::image_info(Image)$height

          #print the ggplot
          suppressWarnings(
            plot(
            ggplot() +
              annotation_raster(Image, xmin = 1, xmax = X_max, ymin = 0, ymax = Y_max) +
              geom_point(aes(X, Y), color = "red", size = 1.2, data = DATA %>% dplyr::filter(Subject_Names == Name_of_Subject)) +
              scale_x_continuous(limits = c(1, X_max)) +
              scale_y_continuous(limits = c(1, Y_max)) +
              ggtitle(Image_type_name) +
              theme_minimal() +
              theme(plot.title = element_text(hjust = 0.5),
                    panel.grid = element_blank())
          )
          )


          #Exit the loop if required with an stop or with a proceed
          Answer <- menu(choices = c("Proceed", "Abort", "Test again"), title = "Check the alignment between images and cells. Should the analysis proceed?")
          if(Answer == 2) stop("Analysis aborted")
        }
      })
    }

    ########Compute COO-target distances########
    # Compute the target cells within distance to each COO
    print("Calculating target cells within distance to every COO")

    #Set the on exit statement
    on.exit({
      future::plan("future::sequential")
      gc()
    })

    #Set the cores
    future::plan("future::multisession", workers = N_cores)
    options(future.globals.maxSize = Inf, future.rng.onMisuse = "ignore")
    furrr::furrr_options(scheduling = Inf)

    COO_list <-
      furrr::future_map(unique(Lookup_table$Subject_Names), function(Image){
        #Obtain individual image data
        DATA <- DATA_Phenotypes %>% dplyr::filter(Subject_Names == Image)

        #Obtain datasets for COO and target cells
        DATA_COO <- DATA %>% dplyr::filter(Phenotype == Cell_Of_Origin)
        DATA_Target <- DATA %>% dplyr::filter(Phenotype == Target_Cell)

        #Generate the COO and Target objects and generate the index
        COO <- cbind(DATA_COO[[2]], DATA_COO[[3]])
        Targets <- rtree::RTree(cbind(DATA_Target[[2]], DATA_Target[[3]]))
        Index <- rtree::withinDistance(Targets, COO, Distance_sampled)

        #Add the index to the COO tibble and count the number of target cells
        DATA_COO$Targe_cell_ID <- Index
        DATA_COO$N_Targets <- purrr::map_dbl(Index, ~length(.))

        #Remove COO without target cells in the distance sampled
        N_COO_without_targets <- sum(DATA_COO$N_Targets == 0)
        Total_N_COO <- nrow(DATA_COO)
        N_COO_without_targets_message <- stringr::str_c(Image, ": ", N_COO_without_targets, " out of ", Total_N_COO, " COO without target cells within sampled distance", "(", Distance_sampled, ").", collapse = "")
        DATA_COO <- DATA_COO %>% dplyr::filter(N_Targets > 0)

        #Exchange the Target_cell_ID by the Cell_no of the actual target cells
        DATA_COO$Targe_cell_ID <- purrr::map(DATA_COO[["Targe_cell_ID"]], ~DATA_Target$Cell_no[.])

        #Return a list with the COO data and the Number of COO without targets
        return(list(DATA_COO = DATA_COO,
                    N_COO_without_targets = N_COO_without_targets,
                    Total_N_COO = Total_N_COO,
                    N_COO_without_targets_message = N_COO_without_targets_message))

      },
      .progress = TRUE)

    future::plan("future::sequential")
    gc()

    #Print message with the number of COO removed due to absence of target cells
    purrr::walk(COO_list, function(Image){
      if(Image[["N_COO_without_targets"]] > 0){
        if(Image[["Total_N_COO"]] > Image[["N_COO_without_targets"]]) Colored_print(paste0("\n", Image[["N_COO_without_targets_message"]]), color = "green")
        else Colored_print(paste0("\n", Image[["N_COO_without_targets_message"]]), color = "red")
      }
    })

    #Obtain the final COO Dataframe and remove the COO list
    DATA_COO <- purrr::map_dfr(COO_list, ~dplyr::bind_rows(.[["DATA_COO"]]))
    rm(COO_list)
    gc()

    #If no COO remain, stop the computation
    if(nrow(DATA_COO) == 0) stop("No COO remain after aplying filters, analysis aborted.")

    #Update the lookup_table to remove subject names without COO-target samples
    Lookup_table <- Lookup_table %>% dplyr::filter(Subject_Names %in% unique(DATA_COO$Subject_Names))

    ########COO-Target triangulation########
    #Now triangulate the target cells
    #First generate a tibble with all the COO - target cell combinations
    print("Triangulating COO-Target search space")
    Target_complete_list <- unlist(DATA_COO$Targe_cell_ID)
    DATA_COO_triangulation <- DATA_COO %>%
      dplyr::slice(rep(1:n(), times = N_Targets))
    DATA_COO_triangulation$Targe_cell_ID <- Target_complete_list

    #Find the target cells
    Target_cell_coords <-
      dplyr::left_join(DATA_COO_triangulation %>% dplyr::select(Targe_cell_ID),
                       DATA_Phenotypes %>% dplyr::select(Cell_no, X, Y),
                       dplyr::join_by(Targe_cell_ID == Cell_no)) %>% dplyr::select(-Targe_cell_ID)
    names(Target_cell_coords) <- c("Target_X", "Target_Y")

    #Bind Target cell coordinates to COO tibble
    DATA_COO_triangulation <-
      DATA_COO_triangulation %>% dplyr::rename("COO_X" = "X",
                                               "COO_Y" = "Y") %>%
      dplyr::bind_cols(Target_cell_coords)

    #Generate the triangulation coordinates and info using the dedicated function Triangle_generator
    Triangle_points <-
      Triangle_generator(DATA = DATA_COO_triangulation[c("COO_X", "COO_Y", "Target_X", "Target_Y")],
                         Angle = Polygon_angle)
    DATA_COO_triangulation <- dplyr::bind_cols(DATA_COO_triangulation, Triangle_points[-c(1:4)])

    #########Random TEST########
    print("Running a test on a random Cell of Origin")


    #Set answer result to 3 and enter the loop
    Test_OK <- 3
    while(Test_OK == 3){
      #Generate test with a random COO to show the user
      Random_COO <- DATA_COO_triangulation %>% dplyr::filter(Cell_no == sample(DATA_COO_triangulation$Cell_no, size = 1))


      #If required by the user remove overlapping target cells
      if(!is.null(Target_cell_Angle_tolerance)){
        #if only one target cell do not perform angle check
        if(nrow(Random_COO) == 1) {Random_COO <- Random_COO %>% dplyr::mutate(Kept_Target = TRUE)}

        #If more than one target cell check the angle differences between cells
        else{
          #Compute the angles
          Angle_check <- Random_COO %>% dplyr::select(Targe_cell_ID, height, COO_tar_angl_deg)
          Angles_norm <- Angle_check$COO_tar_angl_deg + 180
          Angle_check <- Angle_check %>% dplyr::mutate(Angles_norm = Angles_norm)
          Angle_check <- Angle_check %>% arrange(Angles_norm)
          Angles_sorted <- Angle_check$Angles_norm
          Differences <- diff(Angles_sorted)
          Differences <- c(Differences, 360 - (tail(Angles_sorted, 1) - Angles_sorted[1])) #To close the circle
          Angle_check$Angle_diff <- Differences


          if(min(Angle_check$Angle_diff) < Target_cell_Angle_tolerance){
            #NEW design of while loop
            Buzz_word <- "NOT_OK"
            while(Buzz_word == "NOT_OK"){
              #Get the orther of the points from furthest to closest
              Target_ID_order <- Angle_check$Targe_cell_ID[order(Angle_check$height, decreasing = TRUE)]

              #Enter the for loop
              for(i in Target_ID_order){
                #Obtain the challenged cell
                Challenged_cell <- which(i == Angle_check$Targe_cell_ID)

                #Find the index angles to check
                if(Challenged_cell == 1) Angles_to_check <- c(1L, nrow(Angle_check))
                if(Challenged_cell != 1) Angles_to_check <- c(Challenged_cell, Challenged_cell-1)

                #Evaluate the angles
                Angles_to_eval <- Angle_check$Angle_diff[Angles_to_check]

                #If anyone is positive remove the cell and break the loop
                if(any(Angles_to_eval < Target_cell_Angle_tolerance)){
                  Angle_check <- Angle_check[-Challenged_cell, ]

                  break #To break the loop
                }
              }

              #Once out of the loop compute again the differences if more than one row remains
              #If after removal, a single cell remains then proceed
              if(nrow(Angle_check) == 1) Buzz_word <- "PROCEED!!!"
              else{
                Angle_check <- Angle_check %>% arrange(Angles_norm)
                Angles_sorted <- Angle_check$Angles_norm
                Differences <- diff(Angles_sorted)
                Differences <- c(Differences, 360 - (tail(Angles_sorted, 1) - Angles_sorted[1])) #To close the circle
                Angle_check$Angle_diff <- Differences

                #If the min angle is above the tolerance no need to remove any target cells
                if(min(Angle_check$Angle_diff) >=  Target_cell_Angle_tolerance) Buzz_word <- "PROCEED!!!"
              }
            }
          }

          #Add a tag indicating if COO should be kept or not
          Random_COO <- Random_COO %>% dplyr::mutate(Kept_Target = Random_COO$Targe_cell_ID %in% Angle_check$Targe_cell_ID)
        }
      }
      else{Random_COO <- Random_COO %>% dplyr::mutate(Kept_Target = TRUE)}

      #Start bulding the polygons for the target cells
      Random_COO_polygons <- Random_COO %>% dplyr::filter(Kept_Target)
      #Generate the coordinate list
      Polygon_edge_list <-
        purrr::map(1:nrow(Random_COO_polygons), function(Target_cell_index){
          COO_coords <- Random_COO_polygons[Target_cell_index, c(which(names(Random_COO_polygons) == "COO_X"), which(names(Random_COO_polygons) == "COO_Y"))]
          Vertex_1_coords <- Random_COO_polygons[Target_cell_index, c(which(names(Random_COO_polygons) == "Vertex_1_X"), which(names(Random_COO_polygons) == "Vertex_1_Y"))]
          Vertex_2_coords <- Random_COO_polygons[Target_cell_index, c(which(names(Random_COO_polygons) == "Vertex_2_X"), which(names(Random_COO_polygons) == "Vertex_2_Y"))]

          #return a list
          return(
            list(matrix(c(unlist(COO_coords), unlist(Vertex_1_coords), unlist(Vertex_2_coords), unlist(COO_coords)),
                        byrow = TRUE,
                        ncol = 2)
            )
          )
        })
      #Generate the polygons
      Polygons <- lapply(Polygon_edge_list, function(coords) sf::st_polygon(coords))

      #Unify the polygons in an st_sf object
      sf_polygons <- sf::st_sf(
        id = 1:length(Polygons),
        geometry = sf::st_sfc(Polygons)
      )

      Unified_polygons <- sf::st_union(sf_polygons)

      suppressMessages(
        plot(ggplot() +
             geom_point(aes(X, Y), color = "grey", size = 2, data = DATA_Phenotypes %>% dplyr::filter(Subject_Names == unique(Random_COO$Subject_Names)[1])) +
             geom_point(aes(COO_X, COO_Y), color = "blue",size = 2, data = Random_COO %>% dplyr::select(COO_X, COO_Y) %>% dplyr::distinct()) +
             geom_point(aes(Target_X, Target_Y, color = Kept_Target), size = 2, data = Random_COO) +
             geom_sf(linewidth = 0.5, fill = NA, color = "black", data = Unified_polygons) +
             scale_x_continuous(limits = c(min(Random_COO$Target_X, Random_COO$COO_X), max(Random_COO$Target_X, Random_COO$COO_X))) +
             scale_y_continuous(limits = c(min(Random_COO$Target_Y, Random_COO$COO_Y), max(Random_COO$Target_Y, Random_COO$COO_Y))) +
             theme_minimal()
      )
      )


      Test_OK <- menu(choices = c("Proceed", "Abort", "Test again"), title = "Check the barrier polygon generated.Should the analysis proceed?")
      if(Test_OK == 2) stop("Analysis aborted")
    }

    #########Perform Polygon computation for every cell#####
    print("Computing barrier polygons and cells in barrier")

    #Generate a list with reduced datasets to look for barrier cells if required
    if(!is.null(Barrier_Cells)){
      Barrier_cells_list <- purrr::map(unique(Lookup_table$Subject_Names), function(Image) DATA_Phenotypes %>% dplyr::filter(Subject_Names == Image, Phenotype %in% Barrier_Cells))
      names(Barrier_cells_list) <- unique(Lookup_table$Subject_Names)
    }

    #Set the on exit statement
    on.exit({
      future::plan("future::sequential")
      gc()
    })

    #Set the cores
    future::plan("future::multisession", workers = N_cores)
    options(future.globals.maxSize = Inf, future.rng.onMisuse = "ignore")
    furrr::furrr_options(scheduling = Inf)

    #Execute the algorithm
    Results <- furrr::future_map_dfr(unique(DATA_COO_triangulation$Cell_no), function(COO_ID){
      #Obtain COO data
      COO_cell_DATA <- DATA_COO_triangulation %>% dplyr::filter(Cell_no == COO_ID)

      #If required by the user remove overlapping target cells
      if(!is.null(Target_cell_Angle_tolerance)){

        #if only one target cell do not perform angle check
        if(nrow(COO_cell_DATA) == 1) {COO_cell_DATA <- COO_cell_DATA %>% dplyr::mutate(Kept_Target = TRUE)}

        #If more than one target cell check the angle differences between cells
        else{
          #Compute the angles
          Angle_check <- COO_cell_DATA %>% dplyr::select(Targe_cell_ID, height, COO_tar_angl_deg)
          Angles_norm <- Angle_check$COO_tar_angl_deg + 180
          Angle_check <- Angle_check %>% dplyr::mutate(Angles_norm = Angles_norm)
          Angle_check <- Angle_check %>% arrange(Angles_norm)
          Angles_sorted <- Angle_check$Angles_norm
          Differences <- diff(Angles_sorted)
          Differences <- c(Differences, 360 - (tail(Angles_sorted, 1) - Angles_sorted[1])) #To close the circle
          Angle_check$Angle_diff <- Differences

          #If the min angle is smallest than tolerance then start clean up
          if(min(Angle_check$Angle_diff) < Target_cell_Angle_tolerance){
            #NEW design of while loop
            Buzz_word <- "NOT_OK"
            while(Buzz_word == "NOT_OK"){
              #Get the orther of the points from furthest to closest
              Target_ID_order <- Angle_check$Targe_cell_ID[order(Angle_check$height, decreasing = TRUE)]

              #Enter the for loop
              for(i in Target_ID_order){
                #Obtain the challenged cell
                Challenged_cell <- which(i == Angle_check$Targe_cell_ID)

                #Find the index angles to check
                if(Challenged_cell == 1) Angles_to_check <- c(1L, nrow(Angle_check))
                if(Challenged_cell != 1) Angles_to_check <- c(Challenged_cell, Challenged_cell-1)

                #Evaluate the angles
                Angles_to_eval <- Angle_check$Angle_diff[Angles_to_check]

                #If anyone is positive remove the cell and break the loop
                if(any(Angles_to_eval < Target_cell_Angle_tolerance)){
                  Angle_check <- Angle_check[-Challenged_cell, ]

                  break #To break the loop
                }
              }

              #Once out of the loop compute again the differences if more than one row remains
              #If after removal, a single cell remains then proceed
              if(nrow(Angle_check) == 1) Buzz_word <- "PROCEED!!!"
              else{
                Angle_check <- Angle_check %>% arrange(Angles_norm)
                Angles_sorted <- Angle_check$Angles_norm
                Differences <- diff(Angles_sorted)
                Differences <- c(Differences, 360 - (tail(Angles_sorted, 1) - Angles_sorted[1])) #To close the circle
                Angle_check$Angle_diff <- Differences

                #If the min angle is above the tolerance no need to remove any target cells
                if(min(Angle_check$Angle_diff) >=  Target_cell_Angle_tolerance) Buzz_word <- "PROCEED!!!"
              }
            }
          }

          #Add a tag indicating if COO should be kept or not
          COO_cell_DATA <- COO_cell_DATA %>% dplyr::mutate(Kept_Target = COO_cell_DATA$Targe_cell_ID %in% Angle_check$Targe_cell_ID)
        }
      }
      else{COO_cell_DATA <- COO_cell_DATA %>% dplyr::mutate(Kept_Target = TRUE)}

      #Compute the polygons
      COO_cell_DATA_polygons <- COO_cell_DATA %>% dplyr::filter(Kept_Target)
      #Generate the coordinate list
      Polygon_edge_list <-
        purrr::map(1:nrow(COO_cell_DATA_polygons), function(Target_cell_index){
          COO_coords <- COO_cell_DATA_polygons[Target_cell_index, c(which(names(COO_cell_DATA_polygons) == "COO_X"), which(names(COO_cell_DATA_polygons) == "COO_Y"))] #Obtain coords of COO
          Vertex_1_coords <- COO_cell_DATA_polygons[Target_cell_index, c(which(names(COO_cell_DATA_polygons) == "Vertex_1_X"), which(names(COO_cell_DATA_polygons) == "Vertex_1_Y"))] #Obtain coords for vertex1
          Vertex_2_coords <- COO_cell_DATA_polygons[Target_cell_index, c(which(names(COO_cell_DATA_polygons) == "Vertex_2_X"), which(names(COO_cell_DATA_polygons) == "Vertex_2_Y"))] #Obtain coords for vertex2

          #return a list containing the polygon edge coordinates
          return(
            list(matrix(c(unlist(COO_coords), unlist(Vertex_1_coords), unlist(Vertex_2_coords), unlist(COO_coords)), #COO, vertex1, vertex2, COO again (to close polygon)
                        byrow = TRUE,
                        ncol = 2)
            )
          )
        })
      #Generate the polygons
      Polygons <- lapply(Polygon_edge_list, function(coords) sf::st_polygon(coords))

      #Unify the polygons in an st_sf object
      sf_polygons <- sf::st_sf(
        id = 1:length(Polygons),
        geometry = sf::st_sfc(Polygons)
      )
      Unified_polygons <- sf::st_union(sf_polygons)

      #Build the final tibble
      Final_tibble <- COO_cell_DATA %>% dplyr::select(1:4, Phenotype, N_Targets) %>% dplyr::distinct() %>% dplyr::rename("Original_N_Targets" = "N_Targets")
      Final_tibble$Analyzed_N_Targets <- nrow(COO_cell_DATA_polygons)
      Final_tibble$Cumulative_distance <- sum(COO_cell_DATA_polygons$height)
      Final_tibble$Analyzed_Area <- sf::st_area(Unified_polygons)
      Final_tibble$Polygon_objects <- list(list(All_polygons = sf_polygons,
                                           Unified_polygons = Unified_polygons))

      #Compute the Barrier cells within the polygon if required
      if(!is.null(Barrier_Cells)){
        Barrier_cells_df <- Barrier_cells_list[[unique(COO_cell_DATA$Subject_Names)[1]]]
        Barrier_points_sf <- sf::st_as_sf(Barrier_cells_df, coords = c("X", "Y"))
        Barrier_index <- unlist(sf::st_within(Barrier_points_sf, Unified_polygons, sparse = FALSE))
        Barrier_cells_df <- Barrier_cells_df[Barrier_index, ]

        Barrier_cell_count <-
          suppressMessages(
            purrr::map_dfc(Barrier_Cells, function(Barrier_cell_name){
              nrow(Barrier_cells_df %>% dplyr::filter(Phenotype == Barrier_cell_name))
            })
          )
        names(Barrier_cell_count) <- Barrier_Cells
        Final_tibble <- dplyr::bind_cols(Final_tibble, Barrier_cell_count)
      }

      return(Final_tibble)
    },
    .progress = TRUE)

    future::plan("future::sequential")
    gc()

    #########Compute positive pixels within area (if required)#########
    if(!is.null(Image_parameter_list)){
      Pixel_score_list <-
        #1st Iterate for the type of image in the image_parameter_list
        purrr::map(names(Image_parameter_list), function(Image_type_name){

          #Create a message and print it
          Image_number <- which(Image_type_name == names(Image_parameter_list))
          Total_image_types <- length(names(Image_parameter_list))
          cat(paste0("\n", "Computing pixel counts for image set ", Image_number, " out of ", Total_image_types,  ". Computing pixels in barrier for ", Image_type_name, ".", "\n"))

          #Get the information
          Image_information_table <- Lookup_table %>% dplyr::select(Subject_Names, dplyr::starts_with(stringr::str_c(Image_type_name, "_", collapse = ""))) %>%
            dplyr::select(Subject_Names, dplyr::ends_with("_Image_path"))
          Image_channel <- Image_parameter_list[[Image_type_name]][["Channel"]]
          Image_rotate <- Image_parameter_list[[Image_type_name]][["Image_rotate"]]
          Image_x_flip <- Image_parameter_list[[Image_type_name]][["Image_x_flip"]]
          Image_y_flip <- Image_parameter_list[[Image_type_name]][["Image_y_flip"]]

          #2nd Iterate for every ROW in the Image_information_table (that contains individual Subject_Names and image paths)
          Individual_Subject_Names_tibble <-
            purrr::map_dfr(1:nrow(Image_information_table), function(Index_row){

              #Obtain the name of the subject name being evaluated and the path to this image
              Subject_Name_Evaluated <- Image_information_table$Subject_Names[Index_row]
              Image_path <- Image_information_table[[2]][Index_row]

              #Create a message and print it
              cat(paste0("\n", Image_type_name, ": Analyzing sample ", Index_row, "/", nrow(Image_information_table), " - ", Subject_Name_Evaluated, "\n"))

              #Get the COO data for the required subject name
              Image_COO_data <- Results %>% dplyr::filter(Subject_Names == Subject_Name_Evaluated)

              #Obtain the image and perform pre-processing
              #Read image
              if(is.null(Image_channel)) Image <- magick::image_read(Image_path)
              if(!is.null(Image_channel)) Image <- magick::image_read(Image_path)[Image_channel]

              #Execute rotation if required
              if(!is.null(Image_rotate)) Image <- Image %>% magick::image_rotate(degrees = Image_rotate)
              #Execute x or y flip if required
              if(Image_x_flip) Image <- Image %>% magick::image_flop()
              if(!Image_y_flip) Image <- Image %>% magick::image_flip() #Since the function is prepared to work with native Y axis keep it as it is if the user decides to flip it

              #Turn to an EBImage object
              Image <- magick::as_EBImage(Image)
              #Colapse it to a tibble
              Image_tibble <- as_tibble(expand.grid(X = 1:dim(Image)[[1]], Y = 1:dim(Image)[[2]]))
              Image_tibble$Value <- as.vector(Image)

              #remove the image and run gc
              rm(Image)
              gc()

              #Keep non-0 values
              Image_tibble <- Image_tibble %>% dplyr::filter(Value != 0)

              #IF the image contains no positive pixels
              if(nrow(Image_tibble) == 0){
                Individual_COO_tibble <- tibble(Cell_no = Image_COO_data$Cell_no,
                                                Value = 0)
                names(Individual_COO_tibble)[2] <- stringr::str_c(Image_type_name, "_PixelValueSum", collapse = "")
              }

              #If there is at least 1 positive pixel then compute if within polygons
              else{
                #Set the cores for the individual COO iteration
                future::plan("future::multisession", workers = N_cores)
                options(future.globals.maxSize = Inf, future.rng.onMisuse = "ignore")
                furrr::furrr_options(scheduling = Inf)

                #3rd iterate across the COO
                Individual_COO_tibble <-
                  furrr::future_map_dfr(1:nrow(Image_COO_data), function(Row_index){
                    COO_polygon <- Image_COO_data[Row_index,]
                    COO_polygon <- COO_polygon[["Polygon_objects"]][[1]]$Unified_polygons

                    #Find the polygon min and max for X and Y coordinates and filter the Image_tibble
                    COO_polygon_coords <- sf::st_coordinates(COO_polygon)
                    Polygon_X_min <- min(COO_polygon_coords[, "X"])
                    Polygon_X_max <- max(COO_polygon_coords[, "X"])
                    Polygon_Y_min <- min(COO_polygon_coords[, "Y"])
                    Polygon_Y_max <- max(COO_polygon_coords[, "Y"])
                    Candidate_pixels <- Image_tibble %>% dplyr::filter(X >= Polygon_X_min, X <= Polygon_X_max,
                                                                       Y >= Polygon_Y_min, Y <= Polygon_Y_max)

                    #If no positive pixels in candidate pixels skip the sf point process
                    if(nrow(Candidate_pixels) == 0){
                      Pixel_in_barrier_sum <- 0
                    }
                    else{
                      #Compute the pixels that are within the polygon
                      Pixel_points_sf <- sf::st_as_sf(Candidate_pixels, coords = c("X", "Y"))
                      Pixel_Barrier_index <- unlist(sf::st_within(Pixel_points_sf, COO_polygon, sparse = FALSE))
                      #Add the resulting value
                      Pixel_in_barrier_sum <- sum(unlist(Candidate_pixels[Pixel_Barrier_index, "Value"]))
                    }

                    #Prepare the exiting tibble
                    Exiting_tibble <- tibble(Cell_no = Image_COO_data$Cell_no[Row_index],
                                             Value = Pixel_in_barrier_sum)

                    names(Exiting_tibble)[2] <- stringr::str_c(Image_type_name, "_PixelValueSum", collapse = "")
                    return(Exiting_tibble)

                  },
                  .progress = TRUE)
                future::plan("future::sequential")
                gc()
              }

              #Return the COO from each Subject names
              return(Individual_COO_tibble)
            })

          #Return the dataframe with the results from all the subject names
          return(Individual_Subject_Names_tibble)

        })

      #Bind together the subject_names lists
      if(length(Pixel_score_list) > 1) Pixel_results <- purrr::reduce(Pixel_score_list, function(Tibble_A, Tibble_B) dplyr::left_join(Tibble_A, Tibble_B, by = "Cell_no"))
      else{ Pixel_results <- Pixel_score_list[[1]] }

      #Bind the pixel results to the general results tibble
      Results <- dplyr::left_join(Results, Pixel_results, by = "Cell_no")
    }


    #########Return the final result#########

    #Compute the ratio of metrics to Target cells accumulated distance and polygon area
    Results_N_targets <- Results[-c(1:10)] / Results$Analyzed_N_Targets
    names(Results_N_targets) <- stringr::str_c(names(Results_N_targets), "_TargetCell_Ratio", sep = "")

    Results_Distance <- Results[-c(1:10)] / Results$Cumulative_distance
    names(Results_Distance) <- stringr::str_c(names(Results_Distance), "_Distance_Ratio", sep = "")

    Results_Area <- Results[-c(1:10)] / Results$Analyzed_Area
    names(Results_Area) <- stringr::str_c(names(Results_Area), "_Area_Ratio", sep = "")


    #Bind to the results and return
    Results <- dplyr::bind_cols(Results, Results_N_targets, Results_Distance, Results_Area)

    #If pixel_distance ratio has been applied, print a final message
    if(all(!is.null(Image_parameter_list), !is.null(Pixel_distance_ratio))){
      Colored_print(paste0("\n", "WARNING!!\nCell of Origin and barrier polygon coordinates have been modified according to user provided Pixel_distance_ratio"), color = "red")
    }

    #Return the final tibble
    return(Results)
  }





