#' Identify image issues based on correlation hubs
#'
#'
#' @param N_cores Integer. Number of cores to parallelize your computation.
#' 
#' @param Directory Character string specifying the path to the folder where images are present.
#' @param Ordered_Channels Character vector specifying image channels in their exact order.
#' @param Pixel_downsize (OPTIONAL) A numeric value between 0 and 1 indicating the percentage of pixel downsize. If NULL no pixel downsize will be performed.
#' @param Channel_for_test A character vector indicating the channel name to be displayed in the pre-processing test.
#' @param Channels_to_keep Character vector indicating channels used in tissue mask generation.
#' @param Threshold_type_tissueMask Type of threshold to performe tissue mask. Either 'Otsu', 'Arbitrary' or 'Absolute'.
#' @param Threshold_value_tissueMask Numeric value used if Arbitrary is the threshold type of choice.
#' @param Blurr_tissueMask Logical value indicating if image blurring be performed before tissue mask generation.
#' @param Sigma_tissueMask Numeric value indicating the sigma value to perform Gaussian blurring.
#' 
#' @param DATA A tibble containing the cell features.
#' 
#' @param List_correlation_hubs A list containing vectors of channel names where positive correlation is expected.
#' @param Spearman_threshold A single numeric value indicating the positive threshold to consider the correlation value to be adequate (0.05 by default).
#' 
#' @details
#' The function is based on the prior knowledge of protein/biomolecule spatial distribution. If two proteins are expected to be expressed in the same cell/tissue structure they should show
#' positive correlation. If the correlation is not found, this may warrant the review of specific images or channels. The function is based on computing the Spearman correlation index.
#' In order to ameliorate the computational burden, the images can be downsized with the Pixel_downsize argument.
#' 
#' 
#' @returns A tibble containing flagged image channel and samples according to correlation metrics
#'
#' @seealso [Data_QC_Check_function()]
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
#'#Generate list of correlation hubs
#'Correlation_hubs <- list(Nuclear_markers = c("DAPI", "FOXP3"),
#'                         Tumor_markers = c("CK-EPCAM", "PDL1"),
#'                         Lymphocyte_markers = c("CD8a", "GZMB", "PD1"))
#'
#'
#'#Run the test
#'Channel_Sample_diagnostics(
#'   N_cores = 1,
#'  
#'   Directory = Input_Dir,
#'   Ordered_Channels = c("DAPI", "PDL1", "GZMB", "PD1", "CK-EPCAM", "CD8a", "FOXP3", "RANDOM_1", "RANDOM_2", "RANDOM_3", "RANDOM_4"),
#'   Pixel_downsize = 0.2,
#'   Channel_for_test = "DAPI",
#'   Channels_to_keep = c("DAPI", "PDL1", "GZMB", "PD1", "CK-EPCAM", "CD8a", "FOXP3"),
#'   Threshold_type_tissueMask = "Arbitrary",
#'   Threshold_value_tissueMask = 0.05,
#'   Blurr_tissueMask = FALSE,
#'   Sigma_tissueMask = NULL,
#'  
#'   DATA = NULL,
#'  
#'   List_correlation_hubs = Correlation_hubs,
#'   Spearman_threshold = 0.05
#')
#'}
#'
#'
#' @export

Channel_Sample_diagnostics <-
  function(N_cores = 1,
           
           Directory = NULL,
           Ordered_Channels = NULL,
           Pixel_downsize = NULL,
           Channel_for_test = NULL,
           Channels_to_keep = NULL,
           Threshold_type_tissueMask = NULL,
           Threshold_value_tissueMask = NULL,
           Blurr_tissueMask = FALSE,
           Sigma_tissueMask = NULL,
           
           DATA = NULL,
           
           List_correlation_hubs,
           Spearman_threshold = 0.05){
    
    #######CHECK SUGGESTED PACKAGES#######
    if(!is.null(Directory)){
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
      
      #Require patchwork namespace
      requireNamespace("patchwork", quietly = TRUE)
    }
    
    #######CHECK ARGUMENTS#######
    if(all(!is.null(DATA), !is.null(Directory))) stop("Only image data (Directory) or cell data (DATA) should be provided")
    if(all(is.null(DATA), is.null(Directory))) stop("At least image data (Directory) or cell data (DATA) should be provided")
    if(!all(N_cores >= 1 & N_cores%%1 == 0)) stop("N_cores must be an integer value > 0")
    if(!all(length(Spearman_threshold) == 1, Spearman_threshold > 0, Spearman_threshold < 1)) stop("Spearman_threshold must be a single numeric value between 0 and 1")
    
    #Check the correlation hub list
    if(!is.list(List_correlation_hubs)) stop("List_correlation_hubs should be a list")
    if(any(names(List_correlation_hubs) == "")) stop("List_correlation_hubs elements must be named")
    if(length(unique(names(List_correlation_hubs))) != length(List_correlation_hubs)) stop("Names in List_correlation_hubs must be unique")
    if(any(purrr::map_lgl(List_correlation_hubs, ~is.list(.)))) stop("List_correlation_hubs must contain vectors")
    
    #Check arguments if DATA is provided
    if(!is.null(DATA)){
      if(!identical(c("Cell_no", "X", "Y", "Subject_Names"), names(DATA)[c(1:4)])) {
        stop("Your data does not contain adequate format (Cell_no, X, Y, Subject_Names). Please format using the Data_arrange_function.")
      }
      #Check that hub member names are in DATA column names
      purrr::walk(1:length(List_correlation_hubs), function(Hub_index){
        Hub_members <- List_correlation_hubs[[Hub_index]]
        if(!all(Hub_members %in% names(DATA))) stop(
          paste0("Issue with ", names(List_correlation_hubs)[Hub_index], ": one or more members of the hub are not present in DATA")
        )
      })
      
    }
    
    #Check arguments if images are provided
    if(!is.null(Directory)){
      Argument_checker <- c(Empty_directory = length(dir(Directory)) >= 1,
                            Pixel_downsize_OK = any(is.null(Pixel_downsize), all(is.numeric(Pixel_downsize), Pixel_downsize > 0, Pixel_downsize < 1)),
                            Channels_OK = all(Channels_to_keep %in% Ordered_Channels),
                            Channel_for_test_OK = all(length(Channel_for_test) == 1, Channel_for_test %in% Ordered_Channels),
                            Threshold_type_tissueMask_OK = Threshold_type_tissueMask %in% c("Otsu", "Arbitrary", "Absolute"),
                            Threshold_value_tissueMask_OK = if(Threshold_type_tissueMask == "Arbitrary"){
                              all(is.numeric(Threshold_value_tissueMask), Threshold_value_tissueMask >=0, Threshold_value_tissueMask <= 1)
                            }else(TRUE),
                            Blurr_tissueMask_OK = is.logical(Blurr_tissueMask),
                            Sigma_tissueMask_OK = if(Blurr_tissueMask){
                              all(is.numeric(Sigma_tissueMask), Sigma_tissueMask > 0)
                            }else(TRUE)
      )
      
      Stop_messages <- c(Empty_directory = "No files found at the directory provided. Please check out the path.",
                         Pixel_downsize_OK = "Pixel_downsize must be eigher NULL or a numeric value between 0 and 1",
                         Channels_OK = stringr::str_c(
                           "The following channels are not present the channel names provided: ",
                           stringr::str_c(Channels_to_keep[!(Channels_to_keep %in% Ordered_Channels)], collapse = ", "),
                           sep = ""),
                         Channel_for_test_OK = "Channel_for_test must be a single character cotained in Ordered_Channels",
                         Threshold_value_tissueMask_OK = "Threshold_value_tissueMask must be a single numeric value between 0 and 1",
                         Blurr_tissueMask_OK = "Blurr_tissueMask must be a logical value",
                         Sigma_tissueMask_OK = "Sigma_tissueMask must be a positive numeric value > 0"
      )
      #Check arguments and stop if necessary
      if(!all(Argument_checker)){
        stop(cat(Stop_messages[!Argument_checker],
                 fill = sum(!Argument_checker)))
      }
      
      #check names in correlation hub list are present in channel names
      purrr::walk(1:length(List_correlation_hubs), function(Hub_index){
        Hub_members <- List_correlation_hubs[[Hub_index]]
        if(!all(Hub_members %in% Ordered_Channels)) stop(
          paste0("Issue with ", names(List_correlation_hubs)[Hub_index], ": one or more members of the hub are not present in Ordered_Channels")
        )
      })
    }
    
    #######IF DIRECTORY PROVIDED#######
    if(!is.null(Directory)){
      print("Processing Images")
      
      #Get image paths and names
      Image_paths <- dir(Directory, full.names = TRUE)
      Image_names <- dir(Directory, full.names = FALSE)
      
      #Get channel indexes to build tissue mask
      Channel_to_keep_index <- match(Channels_to_keep, Ordered_Channels)
      
      #Run a test of Image, Image downsize and tissue mask
      Answer <- 3
      while(Answer == 3){
        print("Running a tissue mask random test")
        
        Random_image_index <- sample(1:length(Image_paths), size = 1)
        Image_Original <- magick::image_read(Image_paths[Random_image_index])

        #Downsize if required
        if(!is.null(Pixel_downsize)){
         Image_processed <- Image_Original %>% magick::image_resize(magick::geometry_size_percent(width = Pixel_downsize*100))
        }
        else Image_processed <- Image_Original
        
        #Generate tissue mask
        Tissue_mask <- Tissue_mask_generator(Image = Image_processed[Channel_to_keep_index] %>% magick::as_EBImage(),
                                             Threshold_type = Threshold_type_tissueMask,
                                             Threshold_value = Threshold_value_tissueMask,
                                             Blurr = Blurr_tissueMask,
                                             Sigma = Sigma_tissueMask)
        
        
        #Print the result
        Image_width <- magick::image_info(Image_Original[1])$width
        Image_height <- magick::image_info(Image_Original[1])$height
        Channel_display_index <- match(Channel_for_test, Ordered_Channels)

        Original_image_plot <- ggplot() + 
          annotation_raster(Image_Original[Channel_display_index], xmin = 1, ymin = 1, xmax = Image_width, ymax = Image_height) +
          scale_x_continuous("", labels = NULL, limits = c(1, Image_width)) +
          scale_y_continuous("", labels = NULL, limits = c(1, Image_height)) +
          theme_minimal() +
          ggtitle(str_c(Image_names[Random_image_index], Channel_for_test, sep = " - "),
                  "Original image")+
          theme(panel.grid = element_blank(),
                axis.ticks = element_blank(),
                aspect.ratio = Image_height/Image_width)
        
        Pre_processed_plot <- ggplot() + 
          annotation_raster(Image_processed[Channel_display_index], xmin = 1, ymin = 1, xmax = Image_width, ymax = Image_height) +
          scale_x_continuous("", labels = NULL, limits = c(1, Image_width)) +
          scale_y_continuous("", labels = NULL, limits = c(1, Image_height)) +
          theme_minimal() +
          ggtitle(str_c(Image_names[Random_image_index], Channel_for_test, sep = " - "),
                  "Downsized image")+
          theme(panel.grid = element_blank(),
                axis.ticks = element_blank(),
                aspect.ratio = Image_height/Image_width)
        
        Tissue_mask <- ggplot() + 
          annotation_raster(Tissue_mask, xmin = 1, ymin = 1, xmax = Image_width, ymax = Image_height) +
          scale_x_continuous("", labels = NULL, limits = c(1, Image_width)) +
          scale_y_continuous("", labels = NULL, limits = c(1, Image_height)) +
          theme_minimal() +
          ggtitle(Image_names[Random_image_index],
                  "Tissue_mask")+
          theme(panel.grid = element_blank(),
                axis.ticks = element_blank(),
                aspect.ratio = Image_height/Image_width)
        
        plot(Original_image_plot + Pre_processed_plot + Tissue_mask + patchwork::plot_layout(nrow = 1, ncol = 3))
        
        #Exit the loop if required with an stop or with a proceed
        Answer <- menu(choices = c("Proceed", "Abort", "Test again"), title = "Check tissue mask generated. Should the analysis proceed?")
        if(Answer == 2) stop("Analysis aborted")
      }
       
      
      #GENERATE A MASTER TIBBLE with all pixel information by image
      on.exit({
        future::plan("future::sequential")
        gc()
      })
      
      print("Obtaining pixel values from images in foreground and background pixels")
      
      #Multicore settings
      future::plan("future::multisession", workers = N_cores)
      options(future.globals.maxSize = Inf, future.rng.onMisuse = "ignore")
      furrr::furrr_options(scheduling = Inf)
      
      #Generate a tibble with the pixel information
      INFO_tibble <-
        furrr::future_map_dfr(1:length(Image_paths), function(Image_index){
          Image_Original <- magick::image_read(Image_paths[Image_index])
          
          #Downsize if required
          if(!is.null(Pixel_downsize)){
            Image_processed <- Image_Original %>% magick::image_resize(magick::geometry_size_percent(width = Pixel_downsize*100))
          }
          else Image_processed <- Image_Original
          
          #Generate the tissue mask
          Tissue_mask <- Tissue_mask_generator(Image = Image_processed[Channel_to_keep_index] %>% magick::as_EBImage(),
                                               Threshold_type = Threshold_type_tissueMask,
                                               Threshold_value = Threshold_value_tissueMask,
                                               Blurr = Blurr_tissueMask,
                                               Sigma = Sigma_tissueMask)
          
          #Generate the tibble
          Pixel_tibble <- tibble::tibble(Subject_Names = Image_names[Image_index],
                                         Foreground_pixel = as.vector(Tissue_mask))
          
          #Add the pixel values for every required channel and bind it to the pixel tibble
          Used_channels <- unique(unlist(List_correlation_hubs))
          Image_processed <- Image_processed %>% magick::as_EBImage()
          
          Pixel_channel <- suppressMessages(
            purrr::map_dfc(Used_channels, function(channel_name){
              as.vector(Image_processed[,,match(channel_name, Ordered_Channels)])
            })
          ) 
          names(Pixel_channel) <- Used_channels
          Pixel_tibble <- dplyr::bind_cols(Pixel_tibble, Pixel_channel)
          
          #RETURN THE FINAL RESULT
          return(Pixel_tibble)
          
          
        },.progress = TRUE)
      
      #Return to single core
      future::plan("future::sequential")
      gc()
    }
    
    
    #####DEALING WITH TILES####
    print("Detecting if dealing with tiled images...")
    
    #Detect if (TileX-) is present in image names
    Pattern_identified <- stringr::str_detect(unique(INFO_tibble$Subject_Names), "^.*\\(Tile.*\\)\\.tiff$")
    
    #If tiles have been identified then ask the user if pixels should be re-grouped by original image
    if(any(Pattern_identified)){
      Answer <- menu(choices = c("Proceed", "Don't do anything"), title = "\nTiled images have been identified. Pixel information will be analyzed taking into account the orignal images. Should the analysis proceed?")
      if(Answer == 1) INFO_tibble$Subject_Names <- gsub("\\(Tile[^)]*\\)", "", INFO_tibble$Subject_Names)
    }
    
    #If no Tile patterns found then print message and proceed
    else print("No tiled images identified")
    
    #######ANALYZE THE CORRELATION HUBS#######
    #First for images
    if(!is.null(Directory)){
      
      #If no foreground pixels stop the computation
      if(sum(INFO_tibble$Foreground_pixel) < 2) stop("Not enough foreground pixels to compute correlations")
      
      #Iterate over correlation hubs
      RESULT <- 
        purrr::map_dfr(1:length(List_correlation_hubs), function(Hub_index){
          print(paste0("Analyzing pixel correlation for ", names(List_correlation_hubs)[Hub_index]))
          
          #Get Hub members
          Hub_members <- List_correlation_hubs[[Hub_index]]
          
          #Get pixel tibble with specific members
          HUB_tibble <- INFO_tibble %>% dplyr::select(Subject_Names, dplyr::all_of(c("Foreground_pixel", Hub_members)))
          HUB_tibble_foreground <- HUB_tibble %>% dplyr::filter(Foreground_pixel)
          HUB_tibble_background <- HUB_tibble %>% dplyr::filter(!Foreground_pixel)
          
          #Create the correlation tibble
          Correlation_strategy_tibble <- suppressWarnings(as_tibble(expand.grid.unique(Hub_members, Hub_members, include.equals = FALSE)))
          names(Correlation_strategy_tibble) <- c("Var1", "Var2")
          
          
          ####GLOBAL CORRELATION IMAGES####
          #Compute the global correlation for foreground pixels
          Global_correlation_foreground <- map_dbl(1:nrow(Correlation_strategy_tibble), function(Correlation_tibble_row){
            #Get marker 1 and marker 2 names
            Marker_1 <- Correlation_strategy_tibble[[Correlation_tibble_row, 1]]
            Marker_2 <- Correlation_strategy_tibble[[Correlation_tibble_row, 2]]
            
            #Compute the spearman correlation between these two markers
            return(suppressWarnings(cor.test(HUB_tibble_foreground[[Marker_1]], HUB_tibble_foreground[[Marker_2]], method = "spearman")$estimate))   
          })
          Global_cor_tibble <- dplyr::bind_cols(Correlation_strategy_tibble, tibble(Foreground = Global_correlation_foreground))
          
          ####LOCAL CORRELATION IMAGES####
          #Multicore settings
          future::plan("future::multisession", workers = N_cores)
          options(future.globals.maxSize = Inf, future.rng.onMisuse = "ignore")
          furrr::furrr_options(scheduling = Inf)
          
          LOCAL_Correlation_List <- 
            furrr::future_map(unique(HUB_tibble$Subject_Names), function(Image_name){
              #Generate the subject tibbles
              HUB_tibble_subject <- HUB_tibble %>% dplyr::filter(Subject_Names == Image_name)
              HUB_tibble_subject_foreground <- HUB_tibble_subject %>% dplyr::filter(Foreground_pixel)
              
              #If not enough foreground pixels substitute for NA
              if(sum(HUB_tibble_subject$Foreground_pixel) < 2) Local_correlation_foreground <- NA
              else{
                Local_correlation_foreground <- map_dbl(1:nrow(Correlation_strategy_tibble), function(Correlation_tibble_row){
                  #Get marker 1 and marker 2 names
                  Marker_1 <- Correlation_strategy_tibble[[Correlation_tibble_row, 1]]
                  Marker_2 <- Correlation_strategy_tibble[[Correlation_tibble_row, 2]]
                  
                  #Compute the spearman correlation between these two markers
                  return(suppressWarnings(cor.test(HUB_tibble_subject_foreground[[Marker_1]], HUB_tibble_subject_foreground[[Marker_2]], method = "spearman")$estimate))   
                })
              }
              
              Local_cor_tibble <- dplyr::bind_cols(tibble(Subject_Names = Image_name),
                                                   Correlation_strategy_tibble, 
                                                   tibble(Foreground = Local_correlation_foreground))
            })
          
          #Return to single core
          future::plan("future::sequential")
          gc()
          names(LOCAL_Correlation_List) <- unique(HUB_tibble$Subject_Names)
          
          ####ANALYZE CORRELATIONS#####
          #For global analyze those correlations below the threshold, if any channel appears more than 2 times (two problematic correlations) it will be flagged
          Global_cor_tibble_problematic <- Global_cor_tibble %>% dplyr::filter(Foreground < Spearman_threshold)
          
          if(nrow(Global_cor_tibble_problematic) == 0){
            Hub_results <- "No channel issues detected"
            Flagged_channels <- NA
          }
          else{
            Markers_flagged <- (
              tibble(Markers = c(Global_cor_tibble_problematic[["Var1"]], Global_cor_tibble_problematic[["Var2"]])) %>% 
                dplyr::count(Markers) %>% dplyr::filter(n >= 2)
            )[[1]]
            
            Hub_results <- "Problem with hub"
            if(length(Markers_flagged) < 1) Flagged_channels <- NA
            else Flagged_channels <- Markers_flagged
          }
          
          #For sample analysis iterate across the samples and locate samples with at least one flagged marker
          Flagged_images <-
            purrr::map_lgl(LOCAL_Correlation_List, function(Image){
              any(Image$Foreground < Spearman_threshold)
            })
          
          if(all(!Flagged_images)) {
            Images_results <- "No issues with images"
            Flagged_images <- NA
          }
          else{
            Images_results <- "Problem with images"
            Flagged_images <- names(LOCAL_Correlation_List)[Flagged_images]
          }
          
          
          #GENERATE PLOTS
          #Global plot
          Global_hub_correlation <- 
            Global_cor_tibble %>% 
            ggplot(aes(x = Var1, y = fct_reorder(Var2, Var1), fill = Foreground)) +
            geom_tile() +
            geom_text(aes(label = round(Foreground, digits = 2))) +
            scale_x_discrete("") +
            scale_y_discrete("") +
            scale_fill_gradient2("Spearman rho", 
                                 limits = c(-1, 1), 
                                 breaks = c(-1, Spearman_threshold, 1),
                                 midpoint = Spearman_threshold,
                                 low = "#f23a3a", 
                                 mid = "white", 
                                 high = "#3af256"
            ) +
            theme_minimal() +
            guides(fill = "none") +
            ggtitle(names(List_correlation_hubs)[Hub_index], "Overall pixel correlation") +
            theme(panel.grid = element_blank(),
                  axis.text = element_text(color = "black"))
          
          suppressWarnings(plot(Global_hub_correlation))
          
          #Individual image plot
          Image_correlation_tibble <- 
            purrr::map_dfr(LOCAL_Correlation_List, dplyr::bind_rows) %>% 
            dplyr::mutate(Comparisons = str_c(Var1, Var2, sep = " - "))
          
          Global_cor_tibble <- Global_cor_tibble %>% mutate(Comparisons = str_c(Var1, Var2, sep = " - "))
          
          Shadding_tibble <- as_tibble(expand.grid(Comparisons = Global_cor_tibble$Comparisons,
                                                   Foreground = seq(from = -1, to = 1, by = 0.05)))
          
          
          Individual_image_correlation <-
            ggplot() +
            geom_rect(aes(xmin = 0.5, xmax = length(unique(Comparisons))+0.5, ymin = Foreground-0.05, ymax = Foreground+0.05, fill = Foreground), data = Shadding_tibble) +
            geom_hline(yintercept = Spearman_threshold, color = "black") +
            geom_point(aes(x = Comparisons, y = Foreground), position = position_jitter(height = 0, width = 0.1), data = Image_correlation_tibble) +
            geom_point(aes(x = Comparisons, y = Foreground), shape = 8, size = 2, color = "black", data = Global_cor_tibble) +
            scale_fill_gradient2("Spearman rho", 
                                 limits = c(-1, 1), 
                                 breaks = c(-1, Spearman_threshold, 1),
                                 midpoint = Spearman_threshold,
                                 low = "#f23a3a", 
                                 mid = "white", 
                                 high = "#3af256") +
            scale_x_discrete("")+
            scale_y_continuous("Spearman rho")+
            theme_minimal() +
            guides(fill = "none") +
            ggtitle(names(List_correlation_hubs)[Hub_index], "Individual image pixel correlation (asterisc represents overall correlation)") +
            theme(panel.grid = element_blank(),
                  axis.text = element_text(color = "black"))
          
          plot(Individual_image_correlation)
          
          
          #Prepare the summary results
          Summary_tibble_result <- tibble(Correlation_hub = names(List_correlation_hubs)[Hub_index],
                                          Global_issues = Hub_results,
                                          Flagged_channels = list(Flagged_channels),
                                          Individual_images = Images_results,
                                          Flagged_images = list(Flagged_images)
          )
          return(Summary_tibble_result)
        })
      
      return(RESULT)
      
    }
    
    #For data
    if(!is.null(DATA)){
      RESULT <- 
        purrr::map_dfr(1:length(List_correlation_hubs), function(Hub_index){
          
          print(paste0("Analyzing pixel correlation for ", names(List_correlation_hubs)[Hub_index]))
          
          #Get Hub members
          Hub_members <- List_correlation_hubs[[Hub_index]]
          
          #Get Cell information for specific members
          HUB_tibble <- DATA %>% dplyr::select(Subject_Names, dplyr::all_of(Hub_members))
          
          #Create the correlation tibble
          Correlation_strategy_tibble <- suppressWarnings(as_tibble(expand.grid.unique(Hub_members, Hub_members, include.equals = FALSE)))
          names(Correlation_strategy_tibble) <- c("Var1", "Var2")
          
          ####GLOBAL CORRELATION CELLS####
          #Compute the global correlation for all CELLS
          Global_correlation_foreground <- map_dbl(1:nrow(Correlation_strategy_tibble), function(Correlation_tibble_row){
            #Get marker 1 and marker 2 names
            Marker_1 <- Correlation_strategy_tibble[[Correlation_tibble_row, 1]]
            Marker_2 <- Correlation_strategy_tibble[[Correlation_tibble_row, 2]]
            
            #Compute the spearman correlation between these two markers
            return(suppressWarnings(cor.test(HUB_tibble[[Marker_1]], HUB_tibble[[Marker_2]], method = "spearman")$estimate))   
          })
          Global_cor_tibble <- dplyr::bind_cols(Correlation_strategy_tibble, tibble(Spearman_rho_GLOBAL = Global_correlation_foreground))
          
          ####LOCAL CORRELATION CELLS####
          #Multicore settings
          future::plan("future::multisession", workers = N_cores)
          options(future.globals.maxSize = Inf, future.rng.onMisuse = "ignore")
          furrr::furrr_options(scheduling = Inf)
          
          LOCAL_Correlation_List <- 
            furrr::future_map(unique(HUB_tibble$Subject_Names), function(Image_name){
              #Generate the subject tibbles
              HUB_tibble_subject <- HUB_tibble %>% dplyr::filter(Subject_Names == Image_name)
              
              #If not enough foreground pixels substitute for NA
              if(nrow(HUB_tibble_subject) < 2) Local_correlation <- NA
              else{
                Local_correlation <- map_dbl(1:nrow(Correlation_strategy_tibble), function(Correlation_tibble_row){
                  #Get marker 1 and marker 2 names
                  Marker_1 <- Correlation_strategy_tibble[[Correlation_tibble_row, 1]]
                  Marker_2 <- Correlation_strategy_tibble[[Correlation_tibble_row, 2]]
                  
                  #Compute the spearman correlation between these two markers
                  return(suppressWarnings(cor.test(HUB_tibble_subject[[Marker_1]], HUB_tibble_subject[[Marker_2]], method = "spearman")$estimate))   
                })
              }
              
              Local_cor_tibble <- dplyr::bind_cols(tibble(Subject_Names = Image_name),
                                                   Correlation_strategy_tibble, 
                                                   tibble(Spearman_rho_LOCAL = Local_correlation))
            })
          
          #Return to single core
          future::plan("future::sequential")
          gc()
          names(LOCAL_Correlation_List) <- unique(HUB_tibble$Subject_Names)
          
          
          ####ANALYZE CORRELATIONS#####
          #For global analyze those correlations below the threshold, if any channel appears more than 2 times (two problematic correlations) it will be flagged
          Global_cor_tibble_problematic <- Global_cor_tibble %>% dplyr::filter(Spearman_rho_GLOBAL < Spearman_threshold)
          
          if(nrow(Global_cor_tibble_problematic) == 0){
            Hub_results <- "No channel issues detected"
            Flagged_channels <- NA
          }
          else{
            Markers_flagged <- (
              tibble(Markers = c(Global_cor_tibble_problematic[["Var1"]], Global_cor_tibble_problematic[["Var2"]])) %>% 
                dplyr::count(Markers) %>% dplyr::filter(n >= 2)
            )[[1]]
            
            Hub_results <- "Problem with hub"
            if(length(Markers_flagged) < 1) Flagged_channels <- NA
            else Flagged_channels <- Markers_flagged
          }
          
          #For sample analysis iterate across the samples and locate samples with at least one flagged marker
          Flagged_images <-
            purrr::map_lgl(LOCAL_Correlation_List, function(Sample){
              any(Sample$Spearman_rho_LOCAL < Spearman_threshold)
            })
          
          if(all(!Flagged_images)) {
            Images_results <- "No issues with images"
            Flagged_images <- NA
          }
          else{
            Images_results <- "Problem with images"
            
            
            Flagged_images <- names(LOCAL_Correlation_List)[Flagged_images]
          }
          
          
          #GENERATE PLOTS
          #Global plot
          Global_hub_correlation <- 
            Global_cor_tibble %>% 
            ggplot(aes(x = Var1, y = fct_reorder(Var2, Var1), fill = Spearman_rho_GLOBAL)) +
            geom_tile() +
            geom_text(aes(label = round(Spearman_rho_GLOBAL, digits = 2))) +
            scale_x_discrete("") +
            scale_y_discrete("") +
            scale_fill_gradient2("Spearman rho", 
                                 limits = c(-1, 1), 
                                 breaks = c(-1, Spearman_threshold, 1),
                                 midpoint = Spearman_threshold,
                                 low = "#f23a3a", 
                                 mid = "white", 
                                 high = "#3af256"
            ) +
            theme_minimal() +
            guides(fill = "none") +
            ggtitle(names(List_correlation_hubs)[Hub_index], "Overall cell marker expression correlation") +
            theme(panel.grid = element_blank(),
                  axis.text = element_text(color = "black"))
          
          suppressWarnings(plot(Global_hub_correlation))
          
          #Individual image plot
          Image_correlation_tibble <- 
            purrr::map_dfr(LOCAL_Correlation_List, dplyr::bind_rows) %>% 
            dplyr::mutate(Comparisons = str_c(Var1, Var2, sep = " - "))
          
          Global_cor_tibble <- Global_cor_tibble %>% mutate(Comparisons = str_c(Var1, Var2, sep = " - "))
          
          Shadding_tibble <- as_tibble(expand.grid(Comparisons = Global_cor_tibble$Comparisons,
                                                   Spearman_rho_LOCAL = seq(from = -1, to = 1, by = 0.05)))
          
          
          Individual_image_correlation <-
            ggplot() +
            geom_rect(aes(xmin = 0.5, xmax = length(unique(Comparisons))+0.5, ymin = Spearman_rho_LOCAL-0.05, ymax = Spearman_rho_LOCAL+0.05, fill = Spearman_rho_LOCAL), data = Shadding_tibble) +
            geom_hline(yintercept = Spearman_threshold, color = "black") +
            geom_point(aes(x = Comparisons, y = Spearman_rho_LOCAL), position = position_jitter(height = 0, width = 0.1), data = Image_correlation_tibble) +
            geom_point(aes(x = Comparisons, y = Spearman_rho_GLOBAL), shape = 8, size = 2, color = "black", data = Global_cor_tibble) +
            scale_fill_gradient2("Spearman rho", 
                                 limits = c(-1, 1), 
                                 breaks = c(-1, Spearman_threshold, 1),
                                 midpoint = Spearman_threshold,
                                 low = "#f23a3a", 
                                 mid = "white", 
                                 high = "#3af256") +
            scale_x_discrete("")+
            scale_y_continuous("Spearman rho")+
            theme_minimal() +
            guides(fill = "none") +
            ggtitle(names(List_correlation_hubs)[Hub_index], "Individual sample cell marker expression correlation (asterisc represents overall correlation)") +
            theme(panel.grid = element_blank(),
                  axis.text = element_text(color = "black"))
          
          plot(Individual_image_correlation)
          
          
          #Prepare the summary results
          Summary_tibble_result <- tibble(Correlation_hub = names(List_correlation_hubs)[Hub_index],
                                          Global_issues = Hub_results,
                                          Flagged_channels = list(Flagged_channels),
                                          Individual_images = Images_results,
                                          Flagged_images = list(Flagged_images)
          )
          return(Summary_tibble_result)
        })
      
      return(RESULT)
    }
  }



