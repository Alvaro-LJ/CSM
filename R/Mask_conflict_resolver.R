#' Solves conflicts between overlapping image tiles
#'
#' The function is able to resolve conflicts between tiled images that have overlapping edges, reducing edge bias.
#'
#' @param N_cores Integer. Number of cores to parallelize your computation. 
#' @param Directory Character specifying the path to the folder where cell segmentation masks are present. They must be tiled images with overlap.
#' @param Output_directory A character value indicating the path to the empty folder where the refined cell mask images should be written.
#' @param Resolution_approach A character value indicating how conflict between overlapping tiles should be solved. One of the following: Remove_whole_cells, Erode_pixels (see details).  
#' @param Priority_type The type of metric to decide which cells or pixels should prevail after conflict identification. Either Dist_to_edge, Cell_size or both (see details). If both criteria are provided, the first one
#' will be used first, and if there is still a draw, the second one will be used.
#' @param Min_cell_size If Resolution_approach is Erode_pixels, an integer value that indicates the minimum pixels a cell must have to be kept in the analysis. Cells below threshold will be removed from the cell mask image. 
#' 
#' @param Mask_refinement_path (Optional) If NULL, a temporary directory will be created to temporarily save erosion masks (see details). Alternatively, the user can provided the path to an empty folder were erosion masks should be writen.
#'
#' @returns Writes the refined cell mask images in the output directory
#'
#' @seealso [Image_tile_deconstruction_function()]
#' 
#' @details
#' The function first identifies the amount of tile overlap. Afterwards identifies the conflict areas between tiles and solves them using to possible priority scores: 
#' \itemize{
#' \item{Distance to edge (Dist_to_edge): If two cells or pixels collide, the pixel or cell that is closest to it's image edge is removed.}
#' \item{Cell size (Cell_size): If two cells or pixels collide, the one belonging to the smallest cell is removed.}
#' }
#' 
#' Once the losing cell or pixel has been identified then one of the following actions will be taken
#' \itemize{
#' \item{Removing the whole cell (Remove_whole_cells): The loser cell is removed.}
#' \item{Remove the pixel (Erode_pixels): The loser pixel is removed. Afterwards if the cell size is below the Min_cell_size threshold, the whole cell is removed.}
#' }
#' 
#'
#' @examples
#' \dontrun{
#' 
#' Mask_conflict_resolver(
#'   Directory = Path/to/Cell_mask/folder,
#'   Output_directory = Path/to/output/folder
#' )
#' }
#'
#' @export

Mask_conflict_resolver <- 
  function(N_cores = 1,
           Directory,
           Output_directory,
           
           Resolution_approach = "Remove_whole_cells",
           Priority_type = c("Dist_to_edge", "Cell_size"),
           Min_cell_size = NULL,
           
           Mask_refinement_path = NULL){
    
    #####Check suggested packages#####
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
    
    #####Check arguments#####
    if(!all(N_cores >= 1, N_cores%%1 == 0)) stop("N_cores must be an integer value > 0")
    if(!dir.exists(Directory)) stop(paste0("Issue with ", Directory, ". Please check"))
    if(length(list.files(Directory)) == 0) stop(paste0("No files found at ", Directory, ". Please check"))
    if(!dir.exists(Output_directory)) stop(paste0("Issue with ", Output_directory, ". Please check"))
    if(length(list.files(Output_directory) != 0)) stop("Output_directory must be an empty directory")
    if(!all(length(Resolution_approach) == 1, Resolution_approach %in% c("Remove_whole_cells", "Erode_pixels"))) stop("Resolution approach must be one of the following: Remove_whole_cells, Erode_pixels")
    if(!all(Priority_type %in% c("Dist_to_edge", "Cell_size"), length(Priority_type) <= 2)) stop("Priority_type must be one or both of the following: Dist_to_edge, Cell_size")
    if(Resolution_approach == "Erode_pixels" & !is.null(Min_cell_size)){
      if(!all(Min_cell_size%%1 == 0, Min_cell_size > 0, is.numeric(Min_cell_size))) stop("Min_cell_size must be a integer value > 0")
    }
    if(!is.null(Mask_refinement_path)){
      if(!dir.exists(Mask_refinement_path)) stop(paste0("Issue with ", Mask_refinement_path, ". Please check"))
      if(length(list.files(Mask_refinement_path)) > 0) stop("Mask_refinement_path must be an empty directory")
    }
    
    #####Generate the temporary directory if not provided####
    if(is.null(Mask_refinement_path)){
      cat("\n Mask_refinement_path not provided. A temporary directory will be created to perform computations")
      Final_Mask_refinement_path <- tempfile(pattern = "temp_dir")
      dir.create(Final_Mask_refinement_path, recursive = TRUE)
      
      on.exit({
        unlink(Final_Mask_refinement_path, recursive = TRUE)
      })
    } else Final_Mask_refinement_path <- Mask_refinement_path
    
    #####Generate look up table from image names#####
    Mask_image_names <- dir(Directory, full.names = FALSE)
    Mask_image_paths <- dir(Directory, full.names = TRUE)
    
    #Check that all images come from tiles
    Image_from_tile <- stringr::str_detect(Mask_image_names, "\\(Tile")
    if(!all(Image_from_tile)){
      Problematic_images <- Mask_image_names[!Image_from_tile]
      stop(paste0("The following images are not tiles or their file names have been severely modified. Please check before running the function: ",
                  stringr::str_c(Problematic_images, collapse = ", ")))
    }
    
    #Start generating the look_up table
    #First two elements are mask names and paths
    Look_up_table <- tibble::tibble(Mask_names = Mask_image_names,
                                    Mask_paths = Mask_image_paths)
    
    #Obtain the original image mask without the tile info
    Look_up_table$Original_image_name <- gsub("\\([^)]*\\)", "", Look_up_table$Mask_names)
    
    #Isolate the tile info
    Look_up_table$Tile_info <- regmatches(Look_up_table$Mask_names, regexpr("(?<=\\()[^)]*(?=\\))", Look_up_table$Mask_names, perl = TRUE))
    
    #Obtain the coordinates from the tile information
    Tile_info_list <- stringr::str_split(Look_up_table$Tile_info, "-")
    Tile_info_list <- stringr::str_split(Look_up_table$Tile_info, "-")
    Look_up_table$Tile_ID <- purrr::map_chr(Tile_info_list, ~.[[1]])
    Look_up_table$X_Position_Min <-purrr::map_int(Tile_info_list, ~as.integer(.[[2]]))
    Look_up_table$X_Position_Max <-purrr::map_int(Tile_info_list, ~as.integer(.[[3]]))
    Look_up_table$Y_Position_Min <-purrr::map_int(Tile_info_list, ~as.integer(.[[4]]))
    Look_up_table$Y_Position_Max <-purrr::map_int(Tile_info_list, ~as.integer(.[[5]]))
    
    #Compute the tile overlap
    Overlap_value_tibble <-
      purrr::map_dfr(unique(Look_up_table$Original_image_name), 
                     function(Original_image_name_index){
                       #Obtain the isolated images
                       Interim <- Look_up_table %>% dplyr::filter(Original_image_name == Original_image_name_index)
                       
                       #FIRST X POSITIONS
                       #If there is only a single column, then set to NA
                       if(length(sort(unique(Interim$X_Position_Min))) == 1){
                         X_overlap <- NA
                       } else{ #if more than one column then proceed
                         Min_x_position <- sort(unique(Interim$X_Position_Min))
                         Min_x_position <- Min_x_position[-1]
                         Max_x_position <- sort(unique(Interim$X_Position_Max))
                         Max_x_position <- Max_x_position[-length(Max_x_position)]
                         X_overlap <- unique(Max_x_position-Min_x_position) + 1
                         if(length(X_overlap) > 1) stop(paste0("Issue with ", Original_image_name_index, ". Multiple tile overlap values identified. Please check."))
                       }
                       
                       
                       #SECOND Y POSITIONS
                       #If there is only a single row then set to NA
                       if(length(sort(unique(Interim$Y_Position_Min))) == 1){
                         Y_overlap <- NA
                       } else{ #If more than one row then proceed
                         Min_y_position <- sort(unique(Interim$Y_Position_Min))
                         Min_y_position <- Min_y_position[-1]
                         Max_y_position <- sort(unique(Interim$Y_Position_Max))
                         Max_y_position <- Max_y_position[-length(Max_y_position)]
                         Y_overlap <- unique(Max_y_position-Min_y_position) + 1
                         if(length(Y_overlap) > 1) stop(paste0("Issue with ", Original_image_name_index, ". Multiple tile overlap values identified. Please check."))
                       }
                       
                       #Check that they are both equal if both are available
                       Overlap_vector <- c(X_overlap, Y_overlap)

                       #If both values are complete and NO NA values are identified then check that X and Y overlap match
                       if(all(!is.na(Overlap_vector))){
                         #If X and Y overlap do no match, then stop
                         if(X_overlap != Y_overlap) stop(paste0("Issue with ", Original_image_name_index, ". Multiple tile overlap values identified. Please check."))
                       }
                       
                       #If both values are NA, then return and error
                       if(all(is.na(Overlap_vector))) stop(paste0("Issue with ", Original_image_name_index, ". No tile overlap identified. Please check."))
                      
                       #Return the final product
                       return(
                         tibble::tibble(Original_image_name = Original_image_name_index,
                                        Overlap = unique(Overlap_vector[!is.na(Overlap_vector)]))
                       )
                     })
    

    #If there are multiple overlap values then stop the function and return the tibble with the overlap value by image
    if(length(unique(Overlap_value_tibble$Overlap)) > 1){
      tryCatch(
        stop(),
        error = function(e) {
          message("Multiple overlap values identified for the images in the Directory provided")
          return(Overlap_value_tibble)
        }
      )
    }
    
    #If no error, then obtain the overlap value
    Overlap_value_identified <- unique(Overlap_value_tibble$Overlap)
    
    #If the overlap value is not higher than 0 then return an error
    if(Overlap_value_identified < 1) stop(paste0("\n", 
                                                 "An overlap value between tiles of ",
                                                 Overlap_value_identified,
                                                 " pixels has been detected. Impossible to compute conflict between tiles."))
    
    #Else print a nice message and proceed
    cat(paste0("\n", 
               "An overlap value between tiles of ",
               Overlap_value_identified,
               " pixels has been detected. Computing between-tile conflict space...",
               "\n"))
    
    #####Find the conflict space for every image#####
    Overlap_strategy_list <-
      purrr::map(unique(Look_up_table$Original_image_name),
                 function(Original_image_index){
   
                   #Filter the look_up_table
                   Interim <- Look_up_table %>% dplyr::filter(Original_image_name == Original_image_index)
                   
                   #Find conflict in the X space (columns)
                   #If single column images then Columns_overlap tibble is NA
                   if(length(sort(unique(Interim$X_Position_Min))) == 1){
                     Columns_overlap <- tibble::tibble(FROM_X = 1, TO_X = max(Interim$X_Position_Max))
                   } else{ #If multi row compute as usual
                     Min_x_position <- sort(unique(Interim$X_Position_Min))
                     Min_x_position <- Min_x_position[-1]
                     Max_x_position <- sort(unique(Interim$X_Position_Max))
                     Max_x_position <- Max_x_position[-length(Max_x_position)]
                     Columns_overlap <- tibble::tibble(FROM_X = Min_x_position,
                                                       TO_X = Max_x_position)
                   }
                   
                   
                   #Find conflict in the Y space (rows)
                   #If single row images then Rows_overlap tibble is NA
                   if(length(sort(unique(Interim$Y_Position_Min))) == 1){
                     Rows_overlap <- tibble::tibble(FROM_Y = 1, TO_Y = max(Interim$Y_Position_Max))
                   } else{
                     Min_y_position <- sort(unique(Interim$Y_Position_Min))
                     Min_y_position <- Min_y_position[-1]
                     Max_y_position <- sort(unique(Interim$Y_Position_Max))
                     Max_y_position <- Max_y_position[-length(Max_y_position)]
                     Rows_overlap <- tibble::tibble(FROM_Y = Min_y_position,
                                                    TO_Y = Max_y_position)
                   }
                   
                   #If there is any overlap tibble with one single value then there is no need to combine results
                   if(any(length(sort(unique(Interim$X_Position_Min))) == 1, length(sort(unique(Interim$Y_Position_Min))) == 1)){
                     Overlap_crossroads_tibble <- dplyr::bind_cols(Columns_overlap, Rows_overlap)
                   }
                   #Alternatively expand both tibbles to find all the crossroads
                   else{
                     Overlap_crossroads_tibble <- tidyr::expand_grid(Columns_overlap, Rows_overlap)
                   }
                   
                   Overlap_crossroads_tibble$Crossroad <- stringr::str_c("Crossroad_", 1:nrow(Overlap_crossroads_tibble))
                   
                   #For every row in the tibble find the conflicted images in every crossroad
                   Images_in_crossroads <- 
                     purrr::map(seq_along(1:nrow(Overlap_crossroads_tibble)),
                              function(Crossroads_index){
                                #Obtain the from to coordinates
                                FROM_X <- Overlap_crossroads_tibble$FROM_X[[Crossroads_index]]
                                TO_X <- Overlap_crossroads_tibble$TO_X[[Crossroads_index]]
                                FROM_Y <- Overlap_crossroads_tibble$FROM_Y[[Crossroads_index]]
                                TO_Y <- Overlap_crossroads_tibble$TO_Y[[Crossroads_index]]
                                
                                #Select Images corresponding to 1,2,3,4 (topleft, topright, bottomleft, bottomright)
                                Image_1 <- Interim %>% dplyr::filter(X_Position_Max == TO_X, Y_Position_Max == TO_Y) %>%
                                  dplyr::select(Mask_names, Mask_paths, X_Position_Min, X_Position_Max, Y_Position_Min, Y_Position_Max)
                                Image_2 <- Interim %>% dplyr::filter(X_Position_Min == FROM_X, Y_Position_Max == TO_Y) %>%
                                  dplyr::select(Mask_names, Mask_paths, X_Position_Min, X_Position_Max, Y_Position_Min, Y_Position_Max)
                                Image_3 <- Interim %>% dplyr::filter(X_Position_Min == FROM_X, Y_Position_Min == FROM_Y) %>%
                                  dplyr::select(Mask_names, Mask_paths, X_Position_Min, X_Position_Max, Y_Position_Min, Y_Position_Max) 
                                Image_4 <- Interim %>% dplyr::filter(X_Position_Max == TO_X, Y_Position_Min == FROM_Y) %>%
                                  dplyr::select(Mask_names, Mask_paths, X_Position_Min, X_Position_Max, Y_Position_Min, Y_Position_Max) 
                                #Generate the tibble
                                Image_tibble <- dplyr::bind_rows(Image_1, Image_2, Image_3, Image_4) %>% 
                                  dplyr::mutate(Image = stringr::str_c("Image_", 1:4))
                                
                                #Remove duplicates
                                Image_tibble <- Image_tibble[!duplicated(Image_tibble$Mask_names), ]
                              })
                   
                   #Modify the overlap crossroads tibble according to images in crossroads
                   for(index in seq_along(1:length(Images_in_crossroads))){
                     #Isolate current conflict area
                     Current_conflict_area <- Images_in_crossroads[[index]]
                     
                     #If conflict area has 4 images then it's a real crossroads
                     if(nrow(Current_conflict_area) == 4) next
                     
                     #Else check if is a row based or column based 2 image conflict area
                     #First if its a column based
                     if(all(identical(Current_conflict_area$X_Position_Min[1], Current_conflict_area$X_Position_Min[2]),
                            identical(Current_conflict_area$X_Position_Max[1], Current_conflict_area$X_Position_Max[2]))){
                       
                       Overlap_crossroads_tibble$Crossroad[index] <- stringr::str_c("Column_", index)
                     } else{
                       Overlap_crossroads_tibble$Crossroad[index] <- stringr::str_c("Row_", index)
                     }
                   }
                   #Modify the crossroad tibble to include if the conflict zone belongs to the first (upper) row or first(right) column
                   Min_x <- min(Overlap_crossroads_tibble$FROM_X)
                   Max_x <- max(Overlap_crossroads_tibble$FROM_X)
                   Min_y <- min(Overlap_crossroads_tibble$FROM_Y)
                   Max_y <- max(Overlap_crossroads_tibble$FROM_Y)
                   Overlap_crossroads_tibble <- Overlap_crossroads_tibble %>% 
                     dplyr::mutate(First_Column = dplyr::case_when(FROM_X == Min_x ~ TRUE, 
                                                                   FROM_X > Min_x ~ FALSE),
                                   First_Row = dplyr::case_when(FROM_Y == Min_y ~ TRUE, 
                                                                FROM_Y > Min_y ~ FALSE),
                                   Last_Column = dplyr::case_when(FROM_X == Max_x ~ TRUE, 
                                                                  FROM_X < Max_x ~ FALSE),
                                   Last_Row = dplyr::case_when(FROM_Y == Max_y ~ TRUE, 
                                                               FROM_Y < Max_y ~ FALSE))
                   
                   #Now modify the names of Images_in_crossroads accordingly
                   names(Images_in_crossroads) <- Overlap_crossroads_tibble$Crossroad
                   
                   #Generate the final list
                   Final_list <- list(Overlap_pixel_value = Overlap_value_identified,
                                      Conflict_zones = Overlap_crossroads_tibble,
                                      Images_involved = Images_in_crossroads)
                   
                   return(Final_list)
                  
                
                 })
    
    #Add names
    names(Overlap_strategy_list) <- unique(Look_up_table$Original_image_name)
    
    #Generate a tibble containing all the conflict zones
    Overlap_strategy_tibble <- 
      purrr::map_dfr(1:length(Overlap_strategy_list), function(Overlap_strategy_list_index){
        #Get only the conflict zones
        Conflict_zones <- Overlap_strategy_list[[Overlap_strategy_list_index]][["Conflict_zones"]]
        #Add the name of the whole image
        Conflict_zones$Whole_image_name <- names(Overlap_strategy_list)[Overlap_strategy_list_index]
        #Generate a unique ID
        Conflict_zones$Unique_ID <- stringr::str_c(Conflict_zones$Whole_image_name, Conflict_zones$Crossroad, sep = "_")
        #Return the product
        return(Conflict_zones)
      })
    Overlap_strategy_tibble$LOCK <- FALSE

    #####COMPUTE EROSION MASKS AND SAVE THEM IN THE Final_Mask_refinement_path#####
    #Print some messages before proceeding
    if(Resolution_approach == "Remove_whole_cells"){
      cat(paste0("\n",
                 "Conflicts between cell segmentation masks will be resolved by removing whole cells"))
      
    } else{
      cat(paste0("\n",
                 "Conflicts between cell masks will be resolved by eroding pixels from segmentation masks"))
    }
    
    cat(paste0("\n",
               "Primary criteria will be: ", Priority_type[1]))
    if(length(Priority_type) == 2) cat(paste0("\n",
                                              "Secondary criteria will be: ", Priority_type[2]))
    
    cat(paste0("\n",
               "Computing erosion masks for every tiled image",
               "\n"))
    
    #Prepare for parallel computing
    #Set the on exit statement
    on.exit({
      future::plan("future::sequential")
      gc()
    }, add = TRUE)
    
    #Set the cores
    future::plan("future::multisession", workers = N_cores)
    options(future.globals.maxSize = Inf, future.rng.onMisuse = "ignore")
    furrr::furrr_options(scheduling = Inf)
    
    #Compute the total conflicts to be solved
    Total_conflicts <- nrow(Overlap_strategy_tibble)
    
    #Start iteration over conflict zones
    while(nrow(Overlap_strategy_tibble) > 0){
      #Print the percentage done
      Percentage_done <- (1 - (nrow(Overlap_strategy_tibble)/Total_conflicts))*100
      cat("\r", sprintf("%d%%", floor(Percentage_done)))
      
      #Generate the current conflict list round that will be filled by the for loop
      Current_round_conflict_list <- list()
      
      #Fill the current round conflict list
      for(i in 1:N_cores){
        #If no UNLOCKED conflicts then skip
        if(nrow(Overlap_strategy_tibble %>% dplyr::filter(!LOCK)) == 0) next
        
        #Else select a random crossroad from the UNLOCKED samples
        Crossroad_to_be_analyzed <- sample((Overlap_strategy_tibble %>% dplyr::filter(!LOCK))$Unique_ID, size = 1)
        
        #Get the random crossroad data
        Crossroad_to_be_analyzed_data <- Overlap_strategy_tibble %>% dplyr::filter(Unique_ID == Crossroad_to_be_analyzed)
        Crossroad_to_be_analyzed_Image_Name <- Crossroad_to_be_analyzed_data$Whole_image_name[1]
        Crossroad_to_be_analyzed_Conflict_Zone_name <- Crossroad_to_be_analyzed_data$Crossroad[1]
        
        #Get the data for the Conflicted_image_mask_modifier function
        Current_whole_image_data <- Overlap_strategy_list[[Crossroad_to_be_analyzed_Image_Name]]
        Image_param_list <- Current_whole_image_data$Images_involved[[Crossroad_to_be_analyzed_Conflict_Zone_name]]
        Conflict_type <- Crossroad_to_be_analyzed_Conflict_Zone_name
        First_Column <- Crossroad_to_be_analyzed_data$First_Column[1]
        First_Row <- Crossroad_to_be_analyzed_data$First_Row[1]
        Last_Column <- Crossroad_to_be_analyzed_data$Last_Column[1]
        Last_Row <- Crossroad_to_be_analyzed_data$Last_Row[1]
        Overlap_value <- Current_whole_image_data$Overlap_pixel_value
        
        #Place it all in the current round conflict list
        Current_round_conflict_list[[i]] <- list(Image_param_list = Image_param_list,
                                                 Conflict_type = Conflict_type, 
                                                 First_Column = First_Column,
                                                 First_Row = First_Row,
                                                 Last_Column = Last_Column,
                                                 Last_Row = Last_Row,
                                                 Overlap_value = Overlap_value
        )
        
        #Remove the current sample from the work tibble
        Overlap_strategy_tibble <- Overlap_strategy_tibble %>% dplyr::filter(Unique_ID != Crossroad_to_be_analyzed)
        
        #Get the names that will be locked in the current round
        Current_tiles_involved <- Image_param_list$Mask_names
        #Check the crossroads from the current image that need to be locked
        Conflict_index <- purrr::map_lgl(Current_whole_image_data[["Images_involved"]], function(Current_whole_image_data){
          return(any(Current_tiles_involved %in% Current_whole_image_data$Mask_names))
        })
        
        #If no conflicts then proceed to the next round
        if(!any(Conflict_index)) next
        
        #else get the names of the crossroads
        Crossroads_to_be_locked <- names(Current_whole_image_data[["Images_involved"]])[Conflict_index]
        #Generate the unique IDs
        Crossroads_to_be_locked <- stringr::str_c(Crossroad_to_be_analyzed_Image_Name, Crossroads_to_be_locked, sep = "_")
        
        #Set the lock to all the required crossroads
        Overlap_strategy_tibble$LOCK[Overlap_strategy_tibble$Unique_ID %in% Crossroads_to_be_locked] <- TRUE
      }
      
      #Unlock again all the crossroads
      Overlap_strategy_tibble$LOCK <- FALSE
      
      #Perform the computation and repeat as necessary
      furrr::future_walk(Current_round_conflict_list,
                         function(Conflict_index){
                           #Execute the mask modifier
                           Conflicted_image_mask_modifier(Image_param_list = Conflict_index[["Image_param_list"]],
                                                          Conflict_type = Conflict_index[["Conflict_type"]],
                                                          First_Column = Conflict_index[["First_Column"]],
                                                          First_Row = Conflict_index[["First_Row"]],
                                                          Last_Column = Conflict_index[["Last_Column"]],
                                                          Last_Row = Conflict_index[["Last_Row"]],
                                                          Overlap_value = Conflict_index[["Overlap_value"]], 
                                                          
                                                          Resolution_approach = Resolution_approach, 
                                                          Priority_type = Priority_type,
                                                          Min_cell_size = Min_cell_size,
                                                          Mask_refinement_path = Final_Mask_refinement_path)
                           
                         })
      
    }
    
    cat(paste0("\n",
               "Erosion masks computed...")
    )
    
    
    #####Iterate over the mask refinement outcome to erode initial cell masks#####
    if(length(list.files(Final_Mask_refinement_path)) == 0) return("No conflict between tiles has been found!")
    Erosion_mask_present <- list.files(Final_Mask_refinement_path)
    
    cat(paste0("\n",
               length(Erosion_mask_present), " segmentation masks require adjustments. Proceeding with the segmentation mask refinement..."))

    furrr::future_walk(1:length(Mask_image_names),
                function(Mask_index){
                  
                  #Get name and path
                  Mask_name <- Mask_image_names[Mask_index]
                  Mask_path_original <- Mask_image_paths[Mask_index]
                  
                  #If mask name is present in the erosion directory, then it needs erosion
                  if(Mask_name %in% Erosion_mask_present){
                    #Obtain the segmentation mask and the erosion mask
                    Mask <- EBImage::readImage(Mask_path_original)
                    Erosion_mask <- EBImage::readImage(paste0(Final_Mask_refinement_path, "/", Mask_name)) > 0
                    #Apply refinement
                    Mask[Erosion_mask] <- 0
                    
                    #Write results
                    EBImage::writeImage(Mask, paste0(Output_directory, "/", Mask_name), bits.per.sample = 16, compression = "LZW")
                    
                  } else{
                    #If no erosion needed just copy the file
                    file.copy(from = Mask_path_original,
                              to  = paste0(Output_directory, "/", Mask_name))
                    
                  }
                })
    future::plan("future::sequential")
    gc()
    
    #If the Mask refinement path was not provided remove the temporary directory
    if(is.null(Mask_refinement_path)){
      unlink(Final_Mask_refinement_path, recursive = TRUE)
    }
    
    cat(paste0("\n", "DONE!!"))
  }

