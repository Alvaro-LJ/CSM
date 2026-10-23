#' Generate a summary image of thumbnails showing the area surrounding a cell with a given category across images
#'
#' After quantifying cell phenotypes in every image using [Phenotype_quantifier()] function, the user can plot the results using this function.
#'
#' @param DATA A dataframe or tibble containing a column classifying cells in categories (eg: phenotypes, neighborhoods, etc).
#' @param Radius An numeric value indicating the distance from each selected cell to include in the crops.
#' @param N_crops An integer indicating the number of crops to have on each column/row.
#' @param Variable_name A character value indicating the column name containing the cell labels.
#' @param Feature_to_analyze A character vector indicating the labels within the Variable_name that will be analyzed.
#' @param Ordered_Channels Character vector specifying image channels.
#' @param Channels_to_keep Character vector specifying image channels to be kept.
#' @param Channel_for_test A character vector indicating the channel name to be displayed in the pre-processing test.
#' @param Directory Character string specifying the path to the folder where images to be used in tile generation are.
#' @param Output_directory A string indicating the path to the folder in which to store the resulting thumbnails.
#' @param Pixel_distance_ratio (OPTIONAL) A numeric value indicating the ratio between pixel size / cell coordinates. Use this argument when cell coordinates have been transformed to a distance unit like microns.
#'
#' @returns A mosaic of thumbails for each level present on the given category column.
#'
#' @examples
#' \dontrun{
#' #Generate input and output directory
#' Image_Dir <- tempfile(pattern = "tempdir1_Input")
#' dir.create(Image_Dir, recursive = TRUE)
#' Output_Dir <- tempfile(pattern = "tempdir1_Output")
#' dir.create(Output_Dir, recursive = TRUE)
#'
#'#Save images in Input directory
#'purrr::map(1:2,
#'           function(Image){
#'           EBImage::writeImage(CSM_MiniMultiTiff_test[[Image]],
#'                               file.path(Image_Dir, names(CSM_MiniMultiTiff_test)[Image]))
#'         })
#'
#'
#' #Export the thumbnails to the output directory
#'Thumbnail_generator(DATA = CSM_Phenotypecell_test,  
#'                  N_cores = 1,
#'                  Radius = 75, 
#'                  N_crops = 25, 
#'                  Variable_name = 'Phenotype',
#'                  Feature_to_analyze = unique(CSM_Phenotypecell_test$Phenotype)[-1], 
#'                  Directory = Image_Dir,
#'                  Output_directory = 'OUTPUT_TESTS',
#'                  Ordered_Channels = c("DAPI", "PDL1", "GZMB", "PD1", "CK-EPCAM", "CD8a", "FOXP3"), 
#'                  Channels_to_keep = c("DAPI", "PDL1", "GZMB", "PD1", "CK-EPCAM", "CD8a", "FOXP3"),
#'                  Channel_for_test = "DAPI",
#'                  RGB_Color_images = FALSE,
#'                  Image_x_flip = FALSE,
#'                  Image_y_flip = FALSE, 
#'                  Pixel_distance_ratio = 1
#')
#'
#'
#' #Get the crops
#' CD8_GZMBneg <-  magick::image_read(dir(Output_Dir, full.names = TRUE)[1])
#' CD8_GZMBpos <-  magick::image_read(dir(Output_Dir, full.names = TRUE)[2])
#' TUMOR <- magick::image_read(dir(Output_Dir, full.names = TRUE)[3])
#'
#' #Set the merging params
#' Merge_parameter_list <- 
#'   list(
#'       DAPI = list(channel_index = 1,
#'                   color = "blue",
#'                   gamma = 1,
#'                   min = 0,
#'                   max = 80,
#'                   equalize = FALSE),
#'
#'       CK = list(channel_index = 5,
#'                 color = "white",
#'                 gamma = 1.1,
#'                 min = 5,
#'                 max = 20,
#'                 equalize = FALSE),
#'  
#'       CD8 = list(channel_index = 6,
#'                  color = "#f0703a",
#'                  gamma = 1,
#'                  min = 10,
#'                  max = 30,
#'                  equalize = FALSE),
#'
#'       GZMB = list(channel_index = 3,
#'                   color = "cyan",
#'                   gamma = 1,
#'                   min = 5,
#'                   max = 20,
#'                   equalize = FALSE)
#'
#' )
#' 
#' #Merge channels
#' Image_channel_color_merger(Image = CD8_GZMBneg, 
#'                            Parameter_list = Merge_parameter_list)
#' 
#' Image_channel_color_merger(Image = CD8_GZMBpos, 
#'                            Parameter_list = Merge_parameter_list)
#' 
#' Image_channel_color_merger(Image = TUMOR, 
#'                           Parameter_list = Merge_parameter_list)
#'
#' }
#' 
#' @export

Thumbnail_generator <- function(
    N_cores = 1,
    DATA,
    Radius,
    N_crops = 10,
    Variable_name,
    Feature_to_analyze,
    Ordered_Channels = NULL,
    Channels_to_keep = NULL,
    Channel_for_test,
    RGB_Color_images = FALSE,
    Image_x_flip = FALSE,
    Image_y_flip = FALSE,
    Directory,
    Output_directory,
    Pixel_distance_ratio = NULL
) {
  #######Check suggested packages#######
  if(!requireNamespace("magick", quietly = FALSE)) stop(
    paste0("magick CRAN package is required to execute the function. Please install using the following code: ",
           expression(install.packages("magick")))
  )
  
  if(!requireNamespace("ijtiff", quietly = FALSE)) stop(
    paste0("ijtiff CRAN package is required to execute the function. Please install using the following code: ",
           expression(install.packages("ijtiff")))
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
  
  #####CHECK ARGUMENTS#####
  if(!all(N_cores >= 1 & N_cores%%1 == 0)) stop("N_cores must be an integer value > 0")
  if(!identical(names(DATA)[c(1:4)], c("Cell_no", "X", "Y", "Subject_Names"))) stop("DATA has not been adequately formatted with the Data_arrange_function")
  if(!all(length(Variable_name) == 1, Variable_name %in% names(DATA))) stop(paste0(Variable_name, " not present in DATA"))
  if(!all(Feature_to_analyze %in% unique(DATA[[Variable_name]]))){
    Absent_features <- Feature_to_analyze[!Feature_to_analyze %in% unique(DATA[[Variable_name]])]
    stop("The following cell labels are not present in ", Variable_name, 
         '\n',
         paste(Absent_features, sep = ', '))
  }
  if(!all(Radius >= 1, is.numeric(Radius), length(Radius) == 1)) stop("Radius must be an single numeric value >= 1")
  if(!all(N_crops >= 1 & N_crops%%1 == 0)) stop("N_crops must be a single integer value > 0")
  if(!is.logical(RGB_Color_images)) stop("RGB_Color_images must be a logical value")
  #If dealing with RGB images assign channels to RGB
  if(RGB_Color_images){
    Ordered_Channels <- c("R", "G", "B")
    Channels_to_keep <- c("R", "G", "B")
  }
  
  if(!all(Channels_to_keep %in% Ordered_Channels)){
    Problematic_Channels <- Channels_to_keep[!Channels_to_keep %in% Ordered_Channels]
    stop(paste0("The following Channels_to_keep are not present in Ordered_Channels: ", stringr::str_c(Problematic_Channels, collapse = ", ")))
  }
  if(!all(length(Channel_for_test) == 1, Channel_for_test %in% Ordered_Channels)) stop(paste0(Channel_for_test, " not present in ", Ordered_Channels, " or multiple values provided"))
  if(!all(length(Image_x_flip), is.logical(Image_x_flip))) stop("Image_x_flip must be a logical value")
  if(!all(length(Image_y_flip), is.logical(Image_y_flip))) stop("Image_y_flip must be a logical value")
  if(!length(dir(Directory)) >= 1) stop("Directory does not contain any files")
  if(length(dir(Output_directory)) != 0) stop("Output directory should be an empty folder")
  if(!is.null(Pixel_distance_ratio)) {
    if(!all(is.numeric(Pixel_distance_ratio), Pixel_distance_ratio > 0)) stop("Pixel_distance_ratio must be NULL or a numeric value > 0")
  }
   
  #####PRE-Process cell data####
  #If pixel_distance ratio has been provided apply it here
  DATA <- DATA
  if(!is.null(Pixel_distance_ratio)){
    DATA$X <- DATA$X * Pixel_distance_ratio
    DATA$Y <- DATA$Y * Pixel_distance_ratio
  }
  
  #Perform image flipping as required
  if(any(Image_x_flip, Image_y_flip)){
    DATA <-
      purrr::map_dfr(unique(DATA$Subject_Names), function(Subject){
        Subject_data <- DATA %>% dplyr::filter(Subject_Names == Subject)
        
        if(Image_x_flip) Subject_data[c("X", "Y")] <- Flip_xy(MATRIX = as.matrix(Subject_data[c("X", "Y")]), mode = "x")
        if(Image_y_flip) Subject_data[c("X", "Y")] <- Flip_xy(MATRIX = as.matrix(Subject_data[c("X", "Y")]), mode = "y")
        
        return(Subject_data)
      })
  }
  
  #####Run alignment test######
  cat(paste0('\n', 'Running an Image/Cell data alignment test'))
  #Get image names
  Full_image_names <- dir(Directory, full.names = TRUE)
  Image_names <- dir(Directory, full.names = FALSE)
  
  Full_image_names_for_sampling <- dir(Directory, full.names = TRUE)
  Image_names_names_for_sampling <- dir(Directory, full.names = FALSE)
  
  Answer_value <- 3
  while(Answer_value == 3){
    if(length(Full_image_names_for_sampling) < 1){
      Colored_print("All images have been tested. Restarting the sampling process")
      Full_image_names_for_sampling <- dir(Directory, full.names = TRUE)
      Image_names_names_for_sampling <- dir(Directory, full.names = FALSE)
    }
    
    #Sample a random Full_image
    Sample_index <- sample(1:length(Image_names_names_for_sampling), size = 1)
    Current_name <- Image_names_names_for_sampling[Sample_index]
    Current_path <- Full_image_names_for_sampling[Sample_index]
    
    #Remove the index from the list (to avoid selecting it again)
    Full_image_names_for_sampling <- Full_image_names_for_sampling[-Sample_index]
    Image_names_names_for_sampling <- Image_names_names_for_sampling[-Sample_index]
    
    #Get the Subject_Names from data that is closer to the Image_name
    Similarity_scores <- 
      as.vector(adist(Current_name,
                      unique(DATA$Subject_Names),
                      fixed = TRUE,
                      ignore.case = TRUE))
    Current_Subject_Names <- unique(DATA$Subject_Names)[which.min(Similarity_scores)]

    #Obtain the cell data
    Current_cell_data <- DATA %>% dplyr::select(1:4) %>% dplyr::filter(Subject_Names == Current_Subject_Names)
    
    #If the cell data is too large simplify it
    if(nrow(Current_cell_data) > (10^5/2)){
      Colored_print(paste0("Current_Subject_Names", " has over 50K cells. Subsetting the dataset for plotting purposes"), color = 'green')
      Current_cell_data <- Current_cell_data %>% dplyr::slice_sample(n = (10^5/2))
    } 
    
    #Get the image
    Image <- Smart_image_importer(N_cores = N_cores, 
                                  Image_directory = Current_path, 
                                  Log10_pixel_output = 5.69897, 
                                  Ordered_Channels = Ordered_Channels, 
                                  Channels_to_keep = Channel_for_test, 
                                  Save_processed_images = FALSE)
    #Get the widht height and Image
    Image_width <- Image$Original_Dims[1]
    Image_height <- Image$Original_Dims[2]
    Image <- Image$Image %>% magick::image_flip()
    
    #Plot the image
    cat(paste0("\n", "Generating a sample overlay image using ", Current_name))
    plot(
      ggplot() +
        annotation_raster(Image, xmin = 1, xmax = Image_width, ymin = 1, ymax = Image_height) +
        geom_point(aes(x = X, y = Y), color = "red", data = Current_cell_data) +
        theme_minimal() +
        guides(color = "none") +
        scale_x_continuous("", labels = NULL, limits = c(1, Image_width)) +
        scale_y_continuous("", labels = NULL, limits = c(1, Image_height)) +
        ggtitle(Current_name) +
        theme(aspect.ratio = Image_height/Image_width,
              panel.grid = element_blank(),
              panel.background = element_rect(fill = "black"),
              plot.title = element_text(hjust = 0.5))
    )
    
    #Generate a menu to proceed with compuation
    answer <- menu(c("Proceed", "Abort", "Test again"), title = paste0("\n", "Check parameters provided and sample image generated. Should the analysis proceed?"))
    #If user decides to stop then abort function and return stop message
    if(answer == 2) stop("The function has been stopped. Please tune parameters and try again")
    if(answer == 1) Answer_value <- "Proceed"
  }

  #####Extract cells that will be cropped######
  cat(paste0('\n', 'Selecting random cell crops', '\n'))
  #Get only desired candidate cells
  Cell_candidates <- DATA[DATA[[Variable_name]] %in% Feature_to_analyze,]
  
  Cell_candidates <- 
    purrr::map(Feature_to_analyze, function(Variable){
      #Get the cells belonging to the variable
      Current_candidates <- Cell_candidates[Cell_candidates[[Variable_name]] == Variable, ]
      #If there are fewer cells than crops required then reduce the number of crops
      if(N_crops > nrow(Current_candidates)) N_crops <- nrow(Current_candidates)
      #Randomly select the candidates
      Current_candidates <- Current_candidates %>% dplyr::slice_sample(n = N_crops)
    
      #Get only desired variables
      Current_candidates <- Current_candidates %>% dplyr::select(1:4, dplyr::all_of(Variable_name))
      
      #Align the Subject_names with the image names
      #Generate an alignment tibble with  the subject names
      Alignment_tibble <- tibble::tibble(Subject_Names = unique(Current_candidates$Subject_Names))

      #Compute alignment scores
      Alignment_scores <-
        purrr::map(Alignment_tibble$Subject_Names, function(Subject) {
          Scores <- adist(Image_names,
                          Subject,
                          fixed = TRUE,
                          ignore.case = TRUE)
          Image_match <- Image_names[which.min(Scores)]
          Path <- Full_image_names[which.min(Scores)]
          Scores <- min(Scores)
          
          return(list(
            Score = Scores,
            Image = Image_match,
            Path = Path
          ))
        })
      #Place the results of the alignment in the tibble
      Alignment_tibble$Image <- purrr::map_chr(Alignment_scores, function(Item) Item[['Image']])
      Alignment_tibble$Path <- purrr::map_chr(Alignment_scores, function(Item) Item[['Path']])
      Alignment_tibble$Score <- purrr::map_dbl(Alignment_scores, function(Item) Item[['Score']])
     
      #Check for duplicates
      if(any(duplicated(Alignment_tibble$Image))){
        #If duplicates are present keep only the ones with the smallest score
        for(i in unique(Alignment_tibble$Image[duplicated(Alignment_tibble$Image)])){
          Conflictive_Images <- Alignment_tibble %>% dplyr::filter(Image == i)
          Eliminated_images <- Conflictive_Images$Subject_Names[-which.min(Conflictive_Images$Score)]
          Alignment_tibble <- Alignment_tibble %>% dplyr::filter(!Subject_Names %in% Eliminated_images)
          
          #Print a warning message
          Colored_print(paste0("Issue with randomly selected ", Variable, " cell candidates.",
                               "\n",
                               "For the following Subject Names the function was unable to find their matching image in the directory provided:",
                               stringr::str_c(Eliminated_images, collapse = "\n")),
                        color = 'red')
        }
      }
      
      #Bind the tibbles and remove NA values
      Current_candidates <- dplyr::left_join(Current_candidates, Alignment_tibble[-4], by = "Subject_Names")
      Current_candidates <- na.omit(Current_candidates)
      
      #Find the cropping coords (if any below or equal to 0 upgrade to 1)
      Current_candidates$X_Position_Min <- floor(Current_candidates$X - ceiling(Radius/2))
      Current_candidates$X_Position_Min[Current_candidates$X_Position_Min <= 0] <- 1
      Current_candidates$Y_Position_Min <- floor(Current_candidates$Y - ceiling(Radius/2))
      Current_candidates$Y_Position_Min[Current_candidates$Y_Position_Min <= 0] <- 1
      
      Current_candidates$X_Position_Max <- Current_candidates$X_Position_Min + Radius - 1
      Current_candidates$Y_Position_Max <- Current_candidates$Y_Position_Min + Radius - 1
      
      return(Current_candidates)
    })
  
  names(Cell_candidates) <- Feature_to_analyze
  
  #Get the metadata of the unique images to obtain width and height
  Unique_image_paths <- unique(unlist(purrr::map(Cell_candidates, function(Cell_type_canditates) Cell_type_canditates[['Path']])))
  Image_dim_tibble <- tibble::tibble(Path = Unique_image_paths)
  Dim_results <-
    purrr::map(Image_dim_tibble$Path, function(Path_current_image){
      #Get the metadata
      metadata <- RBioFormats::read.metadata(Path_current_image,  filter.metadata = TRUE)
      width <- metadata$coreMetadata$sizeX
      height <- metadata$coreMetadata$sizeY
      #Return width and height
      return(list(Width = width, Height = height))
    })
  
  Image_dim_tibble$X_Image_Max <- purrr::map_dbl(Dim_results, ~.[['Width']])
  Image_dim_tibble$Y_Image_Max <- purrr::map_dbl(Dim_results, ~.[['Height']])
  
  #Add the obtained info to the test
  Cell_candidates <- 
    purrr::map(Cell_candidates, function(Cell_type_canditates){
      return(dplyr::left_join(Cell_type_canditates, Image_dim_tibble, by = 'Path'))
    })
  names(Cell_candidates) <- Feature_to_analyze
  
  ######GENERATE THE THUMBNAILS#####
  cat(paste0('\n', 'Generating thumbnails images', '\n'))
  
  on.exit({
    future::plan("future::sequential")
    gc()
  }
  )
  #Plan the multissesion
  future::plan("future::multisession", workers = N_cores)
  options(future.globals.maxSize = Inf, future.rng.onMisuse = "ignore")
  furrr::furrr_options(scheduling = Inf)
  

  purrr::walk(1:length(Cell_candidates), function(Thumbnail_index){
    #Print a message
    cat(paste0('\n', 'Generating thumbnails image for ', names(Cell_candidates)[Thumbnail_index], '\n'))
    
    #Obtain the current images information and arrange by the subject Names
    Current_images <- Cell_candidates[[Thumbnail_index]] %>% dplyr::arrange(Subject_Names)
    
    #Get the number of rows and columns
    nrow <- floor(sqrt(nrow(Current_images)))
    ncol <- ceiling(nrow(Current_images) / nrow)
    while (nrow * ncol < nrow(Current_images)) {
      nrow <- nrow + 1
      ncol <- ceiling(nrow(Current_images) / nrow)
    }
    
    #Obtain the images
    Image_list <- furrr::future_map(1:nrow(Current_images), function(Image_index){
      #Obtain the current image and the individual characteristics
      Current_image <- Current_images[Image_index, ]
      Image_path <- Current_image[['Path']]
      X_min <- Current_image[['X_Position_Min']]
      X_max <- Current_image[['X_Position_Max']]
      Y_min <- Current_image[['Y_Position_Min']]
      Y_max <- Current_image[['Y_Position_Max']]
      Image_width <- Current_image[['X_Image_Max']]
      Image_height <- Current_image[['Y_Image_Max']]
      
      #Set the padding to false
      X_padding <- FALSE
      Y_padding <- FALSE
      X_padding_pixels <- 0
      Y_padding_pixels <- 0
      
      #If subsetting is larger than available turn to image max and set padding to true
      if(X_max > Image_width){
        X_padding_pixels <- X_max - Image_width
        X_max <- Image_width
        X_padding <- TRUE
      }  
      if(Y_max > Image_height){
        Y_padding_pixels <- Y_max - Image_height
        Y_max <- Image_height
        Y_padding <- TRUE
      } 
      
      #Obtain the image
      Image <- 
        Smart_image_importer(N_cores = 1, 
                             Image_directory = Image_path, 
                             Log10_pixel_output = NULL, Percentage_downsize = 1, 
                             
                             X_crop_coordinates = c(X_min, X_max),
                             Y_crop_coordinates = c(Y_min, Y_max), 
                             
                             Ordered_Channels = Ordered_Channels, Channels_to_keep = Channels_to_keep, 
                             
                             Save_processed_images = FALSE)$Image %>% magick::as_EBImage()
      
      #Turn third dimension to 1 if only a single channel
      if(any(is.null(dim(Image)[3]), is.na(dim(Image)[3]))) dim(Image)[3] <- 1
      
      #Change color mode if necessary
      if(RGB_Color_images) EBImage::colorMode(Image) <- "Color"
      if(!RGB_Color_images) EBImage::colorMode(Image) <- "Grayscale"
      
      #If padding required, perform padding
      if(any(X_padding, Y_padding)){
        Image <- Pad_image(img = Image, 
                           Right_padding = X_padding_pixels, 
                           Bottom_padding = Y_padding_pixels)
      }
      
      #Return the final product
      return(Image)
    })
    
    #If the length of the Image list is smaller than the col*row then add fake images
    if((nrow*ncol) > length(Image_list)){
      #Fake image for color
      if(RGB_Color_images){
        Fake_image <- 
          EBImage::Image(
            array(
              1, #white color
              dim = c(
                Radius,
                Radius,
                length(Channels_to_keep)
              )
            ),
            colormode = "Color"
          )
      }
      
      #Fake image for Grayscale
      if(!RGB_Color_images){
        Fake_image <- 
          EBImage::Image(
            array(
              0, #Black color
              dim = c(
                Radius,
                Radius,
                length(Channels_to_keep)
              )
            ),
            colormode = "Grayscale"
          )
      }
      
      #Generate a fake list
      Fake_list <- list()
      for(i in 1:((nrow*ncol)-length(Image_list))){
        Fake_list[[i]] <- Fake_image
      }
      
      #Append the results to the Image list
      Image_list <- append(Image_list, Fake_list)
      
    }
    
    #Generate the final image
    #Provide dimensions to the list
    dim(Image_list) <- c(nrow, ncol)
    
    #Generate the rows
    Image_panels <- purrr::map(1:nrow, function(Row_index){
      return(EBImage::abind(Image_list[Row_index, ], along = 1))
    })
    
    #Bind the rows together
    Image_panels <- EBImage::abind(Image_panels, along = 2)
    
    #Save the results in the required directory
    Image_name <- paste0(Output_directory, "/", "Thumbnail-", names(Cell_candidates)[Thumbnail_index], ".tif")
    ijtiff::write_tif(img = Image_panels, path = Image_name, bits_per_sample = 32, msg = FALSE, compression = "LZW")
    
    #Print message
    cat(paste0("\n", "Successfully saved ", names(Cell_candidates)[Thumbnail_index], " thumbnail collage(s) in: ", Output_directory, "\n"))
    
  })
    
  future::plan("future::sequential")
  gc()

}