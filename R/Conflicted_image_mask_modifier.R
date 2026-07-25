#' Resolves conflicts between cell masks and generates erosion masks that are saved in an output directory
#'
#' Intended for internal use only. 
#'
#' @param Image_param_list A tibble containing information regarding images involved.
#' @param Conflict_type A character value indicating the type of conflict.
#' @param First_Column A logical value indicating if conflict area is located in the left-side column.
#' @param First_Row A logical value indicating if the conflict area is located in the top-side row.
#' @param Last_Column A logical value indicating if conflict area is located in the right-side column.
#' @param Last_Row A logical value indicating if the conflict area is located in the bottom-side row.
#' @param Overlap_value A integer value indicating the overlap between tiles.
#' 
#' @inheritParams Mask_conflict_resolver Resolution_approach Priority_type Min_cell_size Mask_refinement_path 
#'
#' @details
#' Used in [Mask_conflict_resolver()]
#'
#' @returns It will save erosion masks in Mask_refinement_path directory. No other outputs are generated
#'
#' @keywords Internal

Conflicted_image_mask_modifier <-
  function(Image_param_list = NULL,
           Conflict_type = NULL,
           First_Column = NULL,
           First_Row = NULL,
           Last_Column = NULL,
           Last_Row = NULL,
           Overlap_value = NULL,
           
           Resolution_approach = NULL,
           Priority_type = NULL,
           Min_cell_size = NULL,
           
           Mask_refinement_path = NULL
  ){
    
    ######Load the images in a list######
    Mask_list <- purrr::map(Image_param_list$Mask_paths, function(Image_path) EBImage::readImage(Image_path))
    names(Mask_list) <- Image_param_list$Image
    
    #Image_1 => TOP_LEFT
    #Image_2 => TOP_RIGHT
    #Image_3 => BOTTOM_RIGHT
    #Image_4 => BOTTOM_LEFT
    
    ######See if any mask can be deleted if empty######
    Mask_to_keep <- purrr::map_lgl(Mask_list, ~max(.) > 0)
    if(sum(Mask_to_keep) < 2) return() #if only one image then don't do anything 
    Mask_list <- Mask_list[Mask_to_keep]
    
    Final_names_mask_list <- names(Mask_list) #Keep the final names (they will be stripped by all functions)

    ######Modify Mask list to contain the mask, the dimensions, the cell tibble and the pixel tibble######
    Mask_list <- purrr::map(seq_along(1:length(Mask_list)), 
                            function(Mask_index){
                              
                              #Get Mask
                              Mask <- Mask_list[[Mask_index]]
                              
                              #Obtain the name
                              Mask_name <- names(Mask_list)[Mask_index]
                              
                              #Get dimensions
                              Dimensions_Mask <- dim(Mask)
                              
                              #Modify mask values to integer values
                              Object_id_tibble <- tibble::tibble(value = unique(sort(as.vector(Mask))))
                              Object_id_tibble <- Object_id_tibble[Object_id_tibble$value != 0,]
                              Object_id_tibble$Object_id <- 1:nrow(Object_id_tibble)
                              Match <- match(Mask, Object_id_tibble$value)
                              Mask[!is.na(Match)] <- Object_id_tibble$Object_id[Match[!is.na(Match)]]
                              
                              #Start building the conflicted cell database getting all cell sizes
                              Mask_cellData <- tibble::as_tibble(table(Mask)) %>% dplyr::filter(Mask != 0) %>%
                                dplyr::rename("CellID" = "Mask", "CellSize" = "n") %>% dplyr::mutate_if(is.character, as.integer) %>%
                                dplyr::mutate(CellID = stringr::str_c(Mask_name, "_", CellID))
                              
                              #Now turn each mask to a tibble containing only cell mask pixels
                              Mask_tibble <- tibble::tibble(X_original = as.integer(rep(1:Dimensions_Mask[1], times = Dimensions_Mask[2])),
                                                            Y_original = as.integer(rep(1:Dimensions_Mask[2], each = Dimensions_Mask[1])),
                                                            Value_Mask = as.integer(Mask)) %>%
                                dplyr::filter(Value_Mask != 0) %>% #Remove 0 value values
                                dplyr::mutate(Value_Mask = stringr::str_c(Mask_name, "_", Value_Mask))
                              
                              #Return this in a new list
                              return(list(Mask = Mask,
                                          Dimensions_Mask = Dimensions_Mask,
                                          Mask_cellData = Mask_cellData,
                                          Mask_tibble = Mask_tibble))
                            })
    names(Mask_list) <- Final_names_mask_list
    
    ######Divide each pixel tibble for every mask into one or three areas according to pixel identity and conflict type (horizontal vertical or crossroad)######
    Mask_list_conflict_zones <- purrr::map(seq_along(1:length(Mask_list)), 
                                           function(Mask_index){
                                             #Get the Image type
                                             Image_type <- names(Mask_list)[Mask_index]
                                             
                                             #Get the mask tibble
                                             Mask_tibble <- Mask_list[[Mask_index]][["Mask_tibble"]]
                                             
                                             #Get the mask dimension
                                             Mask_dims <- Mask_list[[Mask_index]][["Dimensions_Mask"]]
                                             
                                             #Get the original pixel position of the mask
                                             Index_in_tibble <- match(Image_type, Image_param_list$Image)
                                             
                                             X_Position_Min <- Image_param_list$X_Position_Min[Index_in_tibble]
                                             X_Position_Max <- Image_param_list$X_Position_Max[Index_in_tibble]
                                             Y_Position_Min <- Image_param_list$Y_Position_Min[Index_in_tibble]
                                             Y_Position_Max <- Image_param_list$Y_Position_Max[Index_in_tibble]
                                             
                                             #IF Row
                                             if(stringr::str_detect(Conflict_type, "Row_")){
                                               if(Image_type == "Image_1"){
                                                 #Horizontal conflict area
                                                 Final_tibble <- 
                                                   Mask_tibble %>%
                                                   dplyr::filter(X_original > (Mask_dims[1] - Overlap_value)) %>%
                                                   dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                 Y_final = Y_original + Y_Position_Min - 1) %>%
                                                   dplyr::mutate(Dist_to_edge = abs(X_original - Mask_dims[1])) %>%
                                                   dplyr::mutate(Zone = "Horizontal",
                                                                 Image_type = "Image_1")
                                                 
                                               }
                                               
                                               if(Image_type != "Image_1"){
                                                 #Horizontal conflict area
                                                 Final_tibble <- 
                                                   Mask_tibble %>%
                                                   dplyr::filter(X_original <= Overlap_value) %>%
                                                   dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                 Y_final = Y_original + Y_Position_Min - 1) %>%
                                                   dplyr::mutate(Dist_to_edge = X_original - 1) %>%
                                                   dplyr::mutate(Zone = "Horizontal",
                                                                 Image_type = "Image_2")
                                               }
                                             }
                                             
                                             #IF Column
                                             if(stringr::str_detect(Conflict_type, "Column_")){
                                               if(Image_type == "Image_1"){
                                                 #Vertical conflict area
                                                 Final_tibble <- 
                                                   Mask_tibble %>%
                                                   dplyr::filter(Y_original > (Mask_dims[2] - Overlap_value)) %>%
                                                   dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                 Y_final = Y_original + Y_Position_Min - 1) %>%
                                                   dplyr::mutate(Dist_to_edge = abs(Y_original - Mask_dims[2])) %>%
                                                   dplyr::mutate(Zone = "Vertical",
                                                                 Image_type = "Image_1")
                                               }
                                               
                                               if(Image_type != "Image_1"){
                                                 #Vertical conflict area
                                                 Final_tibble <- 
                                                   Mask_tibble %>%
                                                   dplyr::filter(Y_original <= Overlap_value) %>%
                                                   dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                 Y_final = Y_original + Y_Position_Min - 1) %>%
                                                   dplyr::mutate(Dist_to_edge = abs(Y_original - 1)) %>%
                                                   dplyr::mutate(Zone = "Vertical",
                                                                 Image_type = "Image_3")
                                               }
                                             }
                                             
                                             #IF crossroads
                                             if(stringr::str_detect(Conflict_type, "Crossroad_")){
                                               #Obtain the three tibbles depending on the Image type and the row/col position of the crossroad
                                               
                                               #Single crossroad
                                               if(First_Column & First_Row & Last_Column & Last_Row){
                                                 if(Image_type == "Image_1"){
                                                   #Horizontal conflict area
                                                   Tibble_horizontal <- 
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original > (Mask_dims[1] - Overlap_value),
                                                                   Y_original <= (Mask_dims[2] - Overlap_value)) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = abs(X_original - Mask_dims[1])) %>%
                                                     dplyr::mutate(Zone = "Horizontal",
                                                                   Image_type = "Image_1")
                                                   
                                                   #Vertical conflict area
                                                   Tibble_vertical <- 
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original <= (Mask_dims[1] - Overlap_value),
                                                                   Y_original > (Mask_dims[2] - Overlap_value)) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = abs(Y_original - Mask_dims[2])) %>%
                                                     dplyr::mutate(Zone = "Vertical",
                                                                   Image_type = "Image_1")
                                                   
                                                   #Shared area
                                                   Tibble_shared <-
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original > (Mask_dims[1] - Overlap_value),
                                                                   Y_original > (Mask_dims[2] - Overlap_value)) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = pmin(abs(X_original - Mask_dims[1]),
                                                                                       abs(Y_original - Mask_dims[2]),
                                                                                       na.rm = TRUE)) %>%
                                                     dplyr::mutate(Zone = "Shared",
                                                                   Image_type = "Image_1")
                                                 }
                                                 
                                                 if(Image_type == "Image_2"){
                                                   #Horizontal conflict area
                                                   Tibble_horizontal <- 
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original <= Overlap_value,
                                                                   Y_original <= (Mask_dims[2] - Overlap_value)) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = X_original - 1) %>%
                                                     dplyr::mutate(Zone = "Horizontal",
                                                                   Image_type = "Image_2")
                                                   
                                                   #Vertical conflict area
                                                   Tibble_vertical <- 
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original > Overlap_value,
                                                                   Y_original > (Mask_dims[2] - Overlap_value)) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = abs(Y_original - Mask_dims[2])) %>%
                                                     dplyr::mutate(Zone = "Vertical",
                                                                   Image_type = "Image_2")
                                                   
                                                   #Shared area
                                                   Tibble_shared <-
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original <= Overlap_value,
                                                                   Y_original > (Mask_dims[2] - Overlap_value)) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = pmin(abs(X_original - 1), 
                                                                                       abs(Y_original - Mask_dims[2]),
                                                                                       na.rm = TRUE)) %>%
                                                     dplyr::mutate(Zone = "Shared",
                                                                   Image_type = "Image_2")
                                                 }
                                                 
                                                 if(Image_type == "Image_3"){
                                                   #Horizontal conflict area
                                                   Tibble_horizontal <- 
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original <= Overlap_value,
                                                                   Y_original >  Overlap_value) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = abs(X_original - 1)) %>%
                                                     dplyr::mutate(Zone = "Horizontal",
                                                                   Image_type = "Image_3")
                                                   
                                                   #Vertical conflict area
                                                   Tibble_vertical <- 
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original > Overlap_value,
                                                                   Y_original <= Overlap_value) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = abs(Y_original - 1)) %>%
                                                     dplyr::mutate(Zone = "Vertical",
                                                                   Image_type = "Image_3")
                                                   
                                                   #Shared area
                                                   Tibble_shared <-
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original <= Overlap_value,
                                                                   Y_original <= Overlap_value) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = pmin(abs(X_original - 1), 
                                                                                       abs(Y_original - 1),
                                                                                       na.rm = TRUE)) %>%
                                                     dplyr::mutate(Zone = "Shared",
                                                                   Image_type = "Image_3")
                                                 }
                                                 
                                                 if(Image_type == "Image_4"){
                                                   #Horizontal conflict area
                                                   Tibble_horizontal <- 
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original > (Mask_dims[1] - Overlap_value),
                                                                   Y_original >  Overlap_value) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = abs(X_original - Mask_dims[1])) %>%
                                                     dplyr::mutate(Zone = "Horizontal",
                                                                   Image_type = "Image_4")
                                                   
                                                   #Vertical conflict area
                                                   Tibble_vertical <- 
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original <= (Mask_dims[1] - Overlap_value),
                                                                   Y_original <= Overlap_value) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = abs(Y_original - 1)) %>%
                                                     dplyr::mutate(Zone = "Vertical",
                                                                   Image_type = "Image_4")
                                                   
                                                   #Shared area
                                                   Tibble_shared <-
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original > (Mask_dims[1] - Overlap_value),
                                                                   Y_original <= Overlap_value) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = pmin(abs(X_original - Mask_dims[1]), 
                                                                                       abs(Y_original - 1), 
                                                                                       na.rm = TRUE)) %>%
                                                     dplyr::mutate(Zone = "Shared",
                                                                   Image_type = "Image_4")
                                                 }
                                               }
                                               
                                               #Top Left 
                                               if(First_Column & First_Row & !Last_Column & !Last_Row){
                                                 if(Image_type == "Image_1"){
                                                   #Horizontal conflict area
                                                   Tibble_horizontal <- 
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original > (Mask_dims[1] - Overlap_value),
                                                                   Y_original <= (Mask_dims[2] - Overlap_value)) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = abs(X_original - Mask_dims[1])) %>%
                                                     dplyr::mutate(Zone = "Horizontal",
                                                                   Image_type = "Image_1")
                                                   
                                                   #Vertical conflict area
                                                   Tibble_vertical <- 
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original <= (Mask_dims[1] - Overlap_value),
                                                                   Y_original > (Mask_dims[2] - Overlap_value)) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = abs(Y_original - Mask_dims[2])) %>%
                                                     dplyr::mutate(Zone = "Vertical",
                                                                   Image_type = "Image_1")
                                                   
                                                   #Shared area
                                                   Tibble_shared <-
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original > (Mask_dims[1] - Overlap_value),
                                                                   Y_original > (Mask_dims[2] - Overlap_value)) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = pmin(abs(X_original - Mask_dims[1]),
                                                                                       abs(Y_original - Mask_dims[2]),
                                                                                       na.rm = TRUE)) %>%
                                                     dplyr::mutate(Zone = "Shared",
                                                                   Image_type = "Image_1")
                                                 }
                                                 
                                                 if(Image_type == "Image_2"){
                                                   #Horizontal conflict area
                                                   Tibble_horizontal <- 
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original <= Overlap_value,
                                                                   Y_original <= (Mask_dims[2] - Overlap_value)) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = X_original - 1) %>%
                                                     dplyr::mutate(Zone = "Horizontal",
                                                                   Image_type = "Image_2")
                                                   
                                                   #Vertical conflict area
                                                   Tibble_vertical <- 
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original > Overlap_value, X_original <= (Mask_dims[1] - Overlap_value), #To avoid shared area with next crossroad
                                                                   Y_original > (Mask_dims[2] - Overlap_value)) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = abs(Y_original - Mask_dims[2])) %>%
                                                     dplyr::mutate(Zone = "Vertical",
                                                                   Image_type = "Image_2")
                                                   
                                                   #Shared area
                                                   Tibble_shared <-
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original <= Overlap_value,
                                                                   Y_original > (Mask_dims[2] - Overlap_value)) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = pmin(abs(X_original - 1), 
                                                                                       abs(Y_original - Mask_dims[2]),
                                                                                       na.rm = TRUE)) %>%
                                                     dplyr::mutate(Zone = "Shared",
                                                                   Image_type = "Image_2")
                                                 }
                                                 
                                                 if(Image_type == "Image_3"){
                                                   #Horizontal conflict area
                                                   Tibble_horizontal <- 
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original <= Overlap_value,
                                                                   Y_original >  Overlap_value, Y_original <= (Mask_dims[2] - Overlap_value)) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = abs(X_original - 1)) %>%
                                                     dplyr::mutate(Zone = "Horizontal",
                                                                   Image_type = "Image_3")
                                                   
                                                   #Vertical conflict area
                                                   Tibble_vertical <- 
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original > Overlap_value, X_original <= (Mask_dims[1] - Overlap_value), #To avoid shared area with next crossroad 
                                                                   Y_original <= Overlap_value) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = abs(Y_original - 1)) %>%
                                                     dplyr::mutate(Zone = "Vertical",
                                                                   Image_type = "Image_3")
                                                   
                                                   #Shared area
                                                   Tibble_shared <-
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original <= Overlap_value,
                                                                   Y_original <= Overlap_value) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = pmin(abs(X_original - 1), 
                                                                                       abs(Y_original - 1),
                                                                                       na.rm = TRUE)) %>%
                                                     dplyr::mutate(Zone = "Shared",
                                                                   Image_type = "Image_3")
                                                 }
                                                 
                                                 if(Image_type == "Image_4"){
                                                   #Horizontal conflict area
                                                   Tibble_horizontal <- 
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original > (Mask_dims[1] - Overlap_value),
                                                                   Y_original >  Overlap_value, Y_original <= (Mask_dims[2] - Overlap_value)) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = abs(X_original - Mask_dims[1])) %>%
                                                     dplyr::mutate(Zone = "Horizontal",
                                                                   Image_type = "Image_4")
                                                   
                                                   #Vertical conflict area
                                                   Tibble_vertical <- 
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original <= (Mask_dims[1] - Overlap_value),
                                                                   Y_original <= Overlap_value) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = abs(Y_original - 1)) %>%
                                                     dplyr::mutate(Zone = "Vertical",
                                                                   Image_type = "Image_4")
                                                   
                                                   #Shared area
                                                   Tibble_shared <-
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original > (Mask_dims[1] - Overlap_value),
                                                                   Y_original <= Overlap_value) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = pmin(abs(X_original - Mask_dims[1]), 
                                                                                       abs(Y_original - 1), 
                                                                                       na.rm = TRUE)) %>%
                                                     dplyr::mutate(Zone = "Shared",
                                                                   Image_type = "Image_4")
                                                 }
                                               }
                                               
                                               #top middle 
                                               if(!First_Column & First_Row & !Last_Column & !Last_Row){
                                                 if(Image_type == "Image_1"){
                                                   #Horizontal conflict area
                                                   Tibble_horizontal <- 
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original > (Mask_dims[1] - Overlap_value),
                                                                   Y_original <= (Mask_dims[2] - Overlap_value)) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = abs(X_original - Mask_dims[1])) %>%
                                                     dplyr::mutate(Zone = "Horizontal",
                                                                   Image_type = "Image_1")
                                                   
                                                   #Vertical conflict area (DOES NOT EXIST)
                                                   Tibble_vertical <- tibble::tibble()
                                                   
                                                   #Shared area
                                                   Tibble_shared <-
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original > (Mask_dims[1] - Overlap_value),
                                                                   Y_original > (Mask_dims[2] - Overlap_value)) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = pmin(abs(X_original - Mask_dims[1]),
                                                                                       abs(Y_original - Mask_dims[2]),
                                                                                       na.rm = TRUE)) %>%
                                                     dplyr::mutate(Zone = "Shared",
                                                                   Image_type = "Image_1")
                                                 }
                                                 
                                                 if(Image_type == "Image_2"){
                                                   #Horizontal conflict area
                                                   Tibble_horizontal <- 
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original <= Overlap_value,
                                                                   Y_original <= (Mask_dims[2] - Overlap_value)) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = X_original - 1) %>%
                                                     dplyr::mutate(Zone = "Horizontal",
                                                                   Image_type = "Image_2")
                                                   
                                                   #Vertical conflict area
                                                   Tibble_vertical <- 
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original > Overlap_value, X_original <= (Mask_dims[1] - Overlap_value),
                                                                   Y_original > (Mask_dims[2] - Overlap_value)) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = abs(Y_original - Mask_dims[2])) %>%
                                                     dplyr::mutate(Zone = "Vertical",
                                                                   Image_type = "Image_2")
                                                   
                                                   #Shared area
                                                   Tibble_shared <-
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original <= Overlap_value,
                                                                   Y_original > (Mask_dims[2] - Overlap_value)) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = pmin(abs(X_original - 1), 
                                                                                       abs(Y_original - Mask_dims[2]),
                                                                                       na.rm = TRUE)) %>%
                                                     dplyr::mutate(Zone = "Shared",
                                                                   Image_type = "Image_2")
                                                 }
                                                 
                                                 if(Image_type == "Image_3"){
                                                   #Horizontal conflict area
                                                   Tibble_horizontal <- 
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original <= Overlap_value,
                                                                   Y_original >  Overlap_value, Y_original <= (Mask_dims[2] - Overlap_value)) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = abs(X_original - 1)) %>%
                                                     dplyr::mutate(Zone = "Horizontal",
                                                                   Image_type = "Image_3")
                                                   
                                                   #Vertical conflict area
                                                   Tibble_vertical <- 
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original > Overlap_value, X_original <= (Mask_dims[1] - Overlap_value),
                                                                   Y_original <= Overlap_value) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = abs(Y_original - 1)) %>%
                                                     dplyr::mutate(Zone = "Vertical",
                                                                   Image_type = "Image_3")
                                                   
                                                   #Shared area
                                                   Tibble_shared <-
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original <= Overlap_value,
                                                                   Y_original <= Overlap_value) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = pmin(abs(X_original - 1), 
                                                                                       abs(Y_original - 1),
                                                                                       na.rm = TRUE)) %>%
                                                     dplyr::mutate(Zone = "Shared",
                                                                   Image_type = "Image_3")
                                                 }
                                                 
                                                 if(Image_type == "Image_4"){
                                                   #Horizontal conflict area
                                                   Tibble_horizontal <- 
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original > (Mask_dims[1] - Overlap_value),
                                                                   Y_original >  Overlap_value, Y_original <= (Mask_dims[2] - Overlap_value)) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = abs(X_original - Mask_dims[1])) %>%
                                                     dplyr::mutate(Zone = "Horizontal",
                                                                   Image_type = "Image_4")
                                                   
                                                   #Vertical conflict area (DOES NOT EXIST)
                                                   Tibble_vertical <- tibble::tibble()
                                                   
                                                   #Shared area
                                                   Tibble_shared <-
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original > (Mask_dims[1] - Overlap_value),
                                                                   Y_original <= Overlap_value) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = pmin(abs(X_original - Mask_dims[1]), 
                                                                                       abs(Y_original - 1), 
                                                                                       na.rm = TRUE)) %>%
                                                     dplyr::mutate(Zone = "Shared",
                                                                   Image_type = "Image_4")
                                                 }
                                               }
                                               
                                               #top right 
                                               if(!First_Column & First_Row & Last_Column & !Last_Row){
                                                 if(Image_type == "Image_1"){
                                                   #Horizontal conflict area
                                                   Tibble_horizontal <- 
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original > (Mask_dims[1] - Overlap_value),
                                                                   Y_original <= (Mask_dims[2] - Overlap_value)) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = abs(X_original - Mask_dims[1])) %>%
                                                     dplyr::mutate(Zone = "Horizontal",
                                                                   Image_type = "Image_1")
                                                   
                                                   #Vertical conflict area (DOES NOT EXIST)
                                                   Tibble_vertical <- tibble::tibble()
                                                   
                                                   #Shared area
                                                   Tibble_shared <-
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original > (Mask_dims[1] - Overlap_value),
                                                                   Y_original > (Mask_dims[2] - Overlap_value)) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = pmin(abs(X_original - Mask_dims[1]),
                                                                                       abs(Y_original - Mask_dims[2]),
                                                                                       na.rm = TRUE)) %>%
                                                     dplyr::mutate(Zone = "Shared",
                                                                   Image_type = "Image_1")
                                                 }
                                                 
                                                 if(Image_type == "Image_2"){
                                                   #Horizontal conflict area
                                                   Tibble_horizontal <- 
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original <= Overlap_value,
                                                                   Y_original <= (Mask_dims[2] - Overlap_value)) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = X_original - 1) %>%
                                                     dplyr::mutate(Zone = "Horizontal",
                                                                   Image_type = "Image_2")
                                                   
                                                   #Vertical conflict area
                                                   Tibble_vertical <- 
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original > Overlap_value,
                                                                   Y_original > (Mask_dims[2] - Overlap_value)) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = abs(Y_original - Mask_dims[2])) %>%
                                                     dplyr::mutate(Zone = "Vertical",
                                                                   Image_type = "Image_2")
                                                   
                                                   #Shared area
                                                   Tibble_shared <-
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original <= Overlap_value,
                                                                   Y_original > (Mask_dims[2] - Overlap_value)) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = pmin(abs(X_original - 1), 
                                                                                       abs(Y_original - Mask_dims[2]),
                                                                                       na.rm = TRUE)) %>%
                                                     dplyr::mutate(Zone = "Shared",
                                                                   Image_type = "Image_2")
                                                 }
                                                 
                                                 if(Image_type == "Image_3"){
                                                   #Horizontal conflict area
                                                   Tibble_horizontal <- 
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original <= Overlap_value,
                                                                   Y_original >  Overlap_value, Y_original <= (Mask_dims[2] - Overlap_value)) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = abs(X_original - 1)) %>%
                                                     dplyr::mutate(Zone = "Horizontal",
                                                                   Image_type = "Image_3")
                                                   
                                                   #Vertical conflict area
                                                   Tibble_vertical <- 
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original > Overlap_value,
                                                                   Y_original <= Overlap_value) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = abs(Y_original - 1)) %>%
                                                     dplyr::mutate(Zone = "Vertical",
                                                                   Image_type = "Image_3")
                                                   
                                                   #Shared area
                                                   Tibble_shared <-
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original <= Overlap_value,
                                                                   Y_original <= Overlap_value) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = pmin(abs(X_original - 1), 
                                                                                       abs(Y_original - 1),
                                                                                       na.rm = TRUE)) %>%
                                                     dplyr::mutate(Zone = "Shared",
                                                                   Image_type = "Image_3")
                                                 }
                                                 
                                                 if(Image_type == "Image_4"){
                                                   #Horizontal conflict area
                                                   Tibble_horizontal <- 
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original > (Mask_dims[1] - Overlap_value),
                                                                   Y_original >  Overlap_value, Y_original <= (Mask_dims[2] - Overlap_value)) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = abs(X_original - Mask_dims[1])) %>%
                                                     dplyr::mutate(Zone = "Horizontal",
                                                                   Image_type = "Image_4")
                                                   
                                                   #Vertical conflict area (DOES NOT EXIST)
                                                   Tibble_vertical <- tibble::tibble()

                                                   #Shared area
                                                   Tibble_shared <-
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original > (Mask_dims[1] - Overlap_value),
                                                                   Y_original <= Overlap_value) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = pmin(abs(X_original - Mask_dims[1]), 
                                                                                       abs(Y_original - 1), 
                                                                                       na.rm = TRUE)) %>%
                                                     dplyr::mutate(Zone = "Shared",
                                                                   Image_type = "Image_4")
                                                 }
                                               }
                                               
                                               #Bottom left  
                                               if(First_Column & !First_Row & !Last_Column & Last_Row){
                                                 if(Image_type == "Image_1"){
                                                   #Horizontal conflict area (DOES NOT EXIST)
                                                   Tibble_horizontal <- tibble::tibble()
                                                   
                                                   #Vertical conflict area
                                                   Tibble_vertical <- 
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original <= (Mask_dims[1] - Overlap_value),
                                                                   Y_original > (Mask_dims[2] - Overlap_value)) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = abs(Y_original - Mask_dims[2])) %>%
                                                     dplyr::mutate(Zone = "Vertical",
                                                                   Image_type = "Image_1")
                                                   
                                                   #Shared area
                                                   Tibble_shared <-
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original > (Mask_dims[1] - Overlap_value),
                                                                   Y_original > (Mask_dims[2] - Overlap_value)) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = pmin(abs(X_original - Mask_dims[1]),
                                                                                       abs(Y_original - Mask_dims[2]),
                                                                                       na.rm = TRUE)) %>%
                                                     dplyr::mutate(Zone = "Shared",
                                                                   Image_type = "Image_1")
                                                 }
                                                 
                                                 if(Image_type == "Image_2"){
                                                   #Horizontal conflict area (DOES NOT EXIST)
                                                   Tibble_horizontal <- tibble::tibble()
                                                   
                                                   #Vertical conflict area
                                                   Tibble_vertical <- 
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original > Overlap_value, X_original <= (Mask_dims[1] - Overlap_value),
                                                                   Y_original > (Mask_dims[2] - Overlap_value)) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = abs(Y_original - Mask_dims[2])) %>%
                                                     dplyr::mutate(Zone = "Vertical",
                                                                   Image_type = "Image_2")
                                                   
                                                   #Shared area
                                                   Tibble_shared <-
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original <= Overlap_value,
                                                                   Y_original > (Mask_dims[2] - Overlap_value)) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = pmin(abs(X_original - 1), 
                                                                                       abs(Y_original - Mask_dims[2]),
                                                                                       na.rm = TRUE)) %>%
                                                     dplyr::mutate(Zone = "Shared",
                                                                   Image_type = "Image_2")
                                                 }
                                                 
                                                 if(Image_type == "Image_3"){
                                                   #Horizontal conflict area
                                                   Tibble_horizontal <- 
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original <= Overlap_value,
                                                                   Y_original >  Overlap_value) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = abs(X_original - 1)) %>%
                                                     dplyr::mutate(Zone = "Horizontal",
                                                                   Image_type = "Image_3")
                                                   
                                                   #Vertical conflict area
                                                   Tibble_vertical <- 
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original > Overlap_value, X_original <= (Mask_dims[1] - Overlap_value),
                                                                   Y_original <= Overlap_value) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = abs(Y_original - 1)) %>%
                                                     dplyr::mutate(Zone = "Vertical",
                                                                   Image_type = "Image_3")
                                                   
                                                   #Shared area
                                                   Tibble_shared <-
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original <= Overlap_value,
                                                                   Y_original <= Overlap_value) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = pmin(abs(X_original - 1), 
                                                                                       abs(Y_original - 1),
                                                                                       na.rm = TRUE)) %>%
                                                     dplyr::mutate(Zone = "Shared",
                                                                   Image_type = "Image_3")
                                                 }
                                                 
                                                 if(Image_type == "Image_4"){
                                                   #Horizontal conflict area
                                                   Tibble_horizontal <- 
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original > (Mask_dims[1] - Overlap_value),
                                                                   Y_original >  Overlap_value) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = abs(X_original - Mask_dims[1])) %>%
                                                     dplyr::mutate(Zone = "Horizontal",
                                                                   Image_type = "Image_4")
                                                   
                                                   #Vertical conflict area
                                                   Tibble_vertical <- 
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original <= (Mask_dims[1] - Overlap_value),
                                                                   Y_original <= Overlap_value) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = abs(Y_original - 1)) %>%
                                                     dplyr::mutate(Zone = "Vertical",
                                                                   Image_type = "Image_4")
                                                   
                                                   #Shared area
                                                   Tibble_shared <-
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original > (Mask_dims[1] - Overlap_value),
                                                                   Y_original <= Overlap_value) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = pmin(abs(X_original - Mask_dims[1]), 
                                                                                       abs(Y_original - 1), 
                                                                                       na.rm = TRUE)) %>%
                                                     dplyr::mutate(Zone = "Shared",
                                                                   Image_type = "Image_4")
                                                 }
                                               }
                                               
                                               #Bottom middle
                                               if(!First_Column & !First_Row & !Last_Column & Last_Row){
                                                 if(Image_type == "Image_1"){
                                                   #Horizontal conflict area (DOES NOT EXIST)
                                                   Tibble_horizontal <- tibble::tibble()
                                                   
                                                   #Vertical conflict area (DOES NOT EXIST)
                                                   Tibble_vertical <- tibble::tibble()
                                                   
                                                   #Shared area
                                                   Tibble_shared <-
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original > (Mask_dims[1] - Overlap_value),
                                                                   Y_original > (Mask_dims[2] - Overlap_value)) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = pmin(abs(X_original - Mask_dims[1]),
                                                                                       abs(Y_original - Mask_dims[2]),
                                                                                       na.rm = TRUE)) %>%
                                                     dplyr::mutate(Zone = "Shared",
                                                                   Image_type = "Image_1")
                                                 }
                                                 
                                                 if(Image_type == "Image_2"){
                                                   #Horizontal conflict area (DOES NOT EXIST)
                                                   Tibble_horizontal <- tibble::tibble()
                                                   
                                                   #Vertical conflict area
                                                   Tibble_vertical <- 
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original > Overlap_value, X_original <= (Mask_dims[1] - Overlap_value),
                                                                   Y_original > (Mask_dims[2] - Overlap_value)) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = abs(Y_original - Mask_dims[2])) %>%
                                                     dplyr::mutate(Zone = "Vertical",
                                                                   Image_type = "Image_2")
                                                   
                                                   #Shared area
                                                   Tibble_shared <-
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original <= Overlap_value,
                                                                   Y_original > (Mask_dims[2] - Overlap_value)) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = pmin(abs(X_original - 1), 
                                                                                       abs(Y_original - Mask_dims[2]),
                                                                                       na.rm = TRUE)) %>%
                                                     dplyr::mutate(Zone = "Shared",
                                                                   Image_type = "Image_2")
                                                 }
                                                 
                                                 if(Image_type == "Image_3"){
                                                   #Horizontal conflict area
                                                   Tibble_horizontal <- 
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original <= Overlap_value,
                                                                   Y_original >  Overlap_value) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = abs(X_original - 1)) %>%
                                                     dplyr::mutate(Zone = "Horizontal",
                                                                   Image_type = "Image_3")
                                                   
                                                   #Vertical conflict area
                                                   Tibble_vertical <- 
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original > Overlap_value, X_original <= (Mask_dims[1] - Overlap_value),
                                                                   Y_original <= Overlap_value) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = abs(Y_original - 1)) %>%
                                                     dplyr::mutate(Zone = "Vertical",
                                                                   Image_type = "Image_3")
                                                   
                                                   #Shared area
                                                   Tibble_shared <-
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original <= Overlap_value,
                                                                   Y_original <= Overlap_value) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = pmin(abs(X_original - 1), 
                                                                                       abs(Y_original - 1),
                                                                                       na.rm = TRUE)) %>%
                                                     dplyr::mutate(Zone = "Shared",
                                                                   Image_type = "Image_3")
                                                 }
                                                 
                                                 if(Image_type == "Image_4"){
                                                   #Horizontal conflict area
                                                   Tibble_horizontal <- 
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original > (Mask_dims[1] - Overlap_value),
                                                                   Y_original >  Overlap_value) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = abs(X_original - Mask_dims[1])) %>%
                                                     dplyr::mutate(Zone = "Horizontal",
                                                                   Image_type = "Image_4")
                                                   
                                                   #Vertical conflict area (DOES NOT EXIST)
                                                   Tibble_vertical <- tibble::tibble()
                                                   
                                                   #Shared area
                                                   Tibble_shared <-
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original > (Mask_dims[1] - Overlap_value),
                                                                   Y_original <= Overlap_value) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = pmin(abs(X_original - Mask_dims[1]), 
                                                                                       abs(Y_original - 1), 
                                                                                       na.rm = TRUE)) %>%
                                                     dplyr::mutate(Zone = "Shared",
                                                                   Image_type = "Image_4")
                                                 }
                                               }
                                               
                                               #Bottom right 
                                               if(!First_Column & !First_Row & Last_Column & Last_Row){
                                                 if(Image_type == "Image_1"){
                                                   #Horizontal conflict area (DOES NOT EXIST)
                                                   Tibble_horizontal <- tibble::tibble()
                                                   
                                                   #Vertical conflict area (DOES NOT EXIST)
                                                   Tibble_vertical <- tibble::tibble()
                                                   
                                                   #Shared area
                                                   Tibble_shared <-
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original > (Mask_dims[1] - Overlap_value),
                                                                   Y_original > (Mask_dims[2] - Overlap_value)) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = pmin(abs(X_original - Mask_dims[1]),
                                                                                       abs(Y_original - Mask_dims[2]),
                                                                                       na.rm = TRUE)) %>%
                                                     dplyr::mutate(Zone = "Shared",
                                                                   Image_type = "Image_1")
                                                 }
                                                 
                                                 if(Image_type == "Image_2"){
                                                   #Horizontal conflict area (DOES NOT EXIST)
                                                   Tibble_horizontal <- tibble::tibble()
                                                   
                                                   #Vertical conflict area
                                                   Tibble_vertical <- 
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original > Overlap_value,
                                                                   Y_original > (Mask_dims[2] - Overlap_value)) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = abs(Y_original - Mask_dims[2])) %>%
                                                     dplyr::mutate(Zone = "Vertical",
                                                                   Image_type = "Image_2")
                                                   
                                                   #Shared area
                                                   Tibble_shared <-
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original <= Overlap_value,
                                                                   Y_original > (Mask_dims[2] - Overlap_value)) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = pmin(abs(X_original - 1), 
                                                                                       abs(Y_original - Mask_dims[2]),
                                                                                       na.rm = TRUE)) %>%
                                                     dplyr::mutate(Zone = "Shared",
                                                                   Image_type = "Image_2")
                                                 }
                                                 
                                                 if(Image_type == "Image_3"){
                                                   #Horizontal conflict area
                                                   Tibble_horizontal <- 
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original <= Overlap_value,
                                                                   Y_original >  Overlap_value) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = abs(X_original - 1)) %>%
                                                     dplyr::mutate(Zone = "Horizontal",
                                                                   Image_type = "Image_3")
                                                   
                                                   #Vertical conflict area
                                                   Tibble_vertical <- 
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original > Overlap_value,
                                                                   Y_original <= Overlap_value) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = abs(Y_original - 1)) %>%
                                                     dplyr::mutate(Zone = "Vertical",
                                                                   Image_type = "Image_3")
                                                   
                                                   #Shared area
                                                   Tibble_shared <-
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original <= Overlap_value,
                                                                   Y_original <= Overlap_value) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = pmin(abs(X_original - 1), 
                                                                                       abs(Y_original - 1),
                                                                                       na.rm = TRUE)) %>%
                                                     dplyr::mutate(Zone = "Shared",
                                                                   Image_type = "Image_3")
                                                 }
                                                 
                                                 if(Image_type == "Image_4"){
                                                   #Horizontal conflict area
                                                   Tibble_horizontal <- 
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original > (Mask_dims[1] - Overlap_value),
                                                                   Y_original >  Overlap_value) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = abs(X_original - Mask_dims[1])) %>%
                                                     dplyr::mutate(Zone = "Horizontal",
                                                                   Image_type = "Image_4")
                                                   
                                                   #Vertical conflict area (DOES NOT EXIST)
                                                   Tibble_vertical <- tibble::tibble()
                                                   
                                                   #Shared area
                                                   Tibble_shared <-
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original > (Mask_dims[1] - Overlap_value),
                                                                   Y_original <= Overlap_value) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = pmin(abs(X_original - Mask_dims[1]), 
                                                                                       abs(Y_original - 1), 
                                                                                       na.rm = TRUE)) %>%
                                                     dplyr::mutate(Zone = "Shared",
                                                                   Image_type = "Image_4")
                                                 }
                                               }
                                               
                                               #Left middle
                                               if(First_Column & !First_Row & !Last_Column & !Last_Row){
                                                 if(Image_type == "Image_1"){
                                                   #Horizontal conflict area (DOES NOT EXIST)
                                                   Tibble_horizontal <- tibble::tibble()
                                                   
                                                   #Vertical conflict area
                                                   Tibble_vertical <- 
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original <= (Mask_dims[1] - Overlap_value),
                                                                   Y_original > (Mask_dims[2] - Overlap_value)) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = abs(Y_original - Mask_dims[2])) %>%
                                                     dplyr::mutate(Zone = "Vertical",
                                                                   Image_type = "Image_1")
                                                   
                                                   #Shared area
                                                   Tibble_shared <-
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original > (Mask_dims[1] - Overlap_value),
                                                                   Y_original > (Mask_dims[2] - Overlap_value)) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = pmin(abs(X_original - Mask_dims[1]),
                                                                                       abs(Y_original - Mask_dims[2]),
                                                                                       na.rm = TRUE)) %>%
                                                     dplyr::mutate(Zone = "Shared",
                                                                   Image_type = "Image_1")
                                                 }
                                                 
                                                 if(Image_type == "Image_2"){
                                                   #Horizontal conflict area (DOES NOT EXIST)
                                                   Tibble_horizontal <- tibble::tibble()
                                                   
                                                   #Vertical conflict area
                                                   Tibble_vertical <- 
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original > Overlap_value, X_original <= (Mask_dims[1] - Overlap_value),
                                                                   Y_original > (Mask_dims[2] - Overlap_value)) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = abs(Y_original - Mask_dims[2])) %>%
                                                     dplyr::mutate(Zone = "Vertical",
                                                                   Image_type = "Image_2")
                                                   
                                                   #Shared area
                                                   Tibble_shared <-
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original <= Overlap_value,
                                                                   Y_original > (Mask_dims[2] - Overlap_value)) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = pmin(abs(X_original - 1), 
                                                                                       abs(Y_original - Mask_dims[2]),
                                                                                       na.rm = TRUE)) %>%
                                                     dplyr::mutate(Zone = "Shared",
                                                                   Image_type = "Image_2")
                                                 }
                                                 
                                                 if(Image_type == "Image_3"){
                                                   #Horizontal conflict area
                                                   Tibble_horizontal <- 
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original <= Overlap_value,
                                                                   Y_original >  Overlap_value, Y_original <= (Mask_dims[2] - Overlap_value)) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = abs(X_original - 1)) %>%
                                                     dplyr::mutate(Zone = "Horizontal",
                                                                   Image_type = "Image_3")
                                                   
                                                   #Vertical conflict area
                                                   Tibble_vertical <- 
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original > Overlap_value, X_original <= (Mask_dims[1] - Overlap_value),
                                                                   Y_original <= Overlap_value) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = abs(Y_original - 1)) %>%
                                                     dplyr::mutate(Zone = "Vertical",
                                                                   Image_type = "Image_3")
                                                   
                                                   #Shared area
                                                   Tibble_shared <-
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original <= Overlap_value,
                                                                   Y_original <= Overlap_value) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = pmin(abs(X_original - 1), 
                                                                                       abs(Y_original - 1),
                                                                                       na.rm = TRUE)) %>%
                                                     dplyr::mutate(Zone = "Shared",
                                                                   Image_type = "Image_3")
                                                 }
                                                 
                                                 if(Image_type == "Image_4"){
                                                   #Horizontal conflict area
                                                   Tibble_horizontal <- 
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original > (Mask_dims[1] - Overlap_value),
                                                                   Y_original >  Overlap_value, Y_original <= (Mask_dims[2] - Overlap_value)) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = abs(X_original - Mask_dims[1])) %>%
                                                     dplyr::mutate(Zone = "Horizontal",
                                                                   Image_type = "Image_4")
                                                   
                                                   #Vertical conflict area
                                                   Tibble_vertical <- 
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original <= (Mask_dims[1] - Overlap_value),
                                                                   Y_original <= Overlap_value) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = abs(Y_original - 1)) %>%
                                                     dplyr::mutate(Zone = "Vertical",
                                                                   Image_type = "Image_4")
                                                   
                                                   #Shared area
                                                   Tibble_shared <-
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original > (Mask_dims[1] - Overlap_value),
                                                                   Y_original <= Overlap_value) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = pmin(abs(X_original - Mask_dims[1]), 
                                                                                       abs(Y_original - 1), 
                                                                                       na.rm = TRUE)) %>%
                                                     dplyr::mutate(Zone = "Shared",
                                                                   Image_type = "Image_4")
                                                 }
                                               }
                                               
                                               #Right middle
                                               if(!First_Column & !First_Row & Last_Column & !Last_Row){
                                                 if(Image_type == "Image_1"){
                                                   #Horizontal conflict area (DOES NOT EXIST)
                                                   Tibble_horizontal <- tibble::tibble()
                                                   
                                                   #Vertical conflict area (DOES NOT EXIST)
                                                   Tibble_vertical <- tibble::tibble()
                                                   
                                                   #Shared area
                                                   Tibble_shared <-
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original > (Mask_dims[1] - Overlap_value),
                                                                   Y_original > (Mask_dims[2] - Overlap_value)) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = pmin(abs(X_original - Mask_dims[1]),
                                                                                       abs(Y_original - Mask_dims[2]),
                                                                                       na.rm = TRUE)) %>%
                                                     dplyr::mutate(Zone = "Shared",
                                                                   Image_type = "Image_1")
                                                 }
                                                 
                                                 if(Image_type == "Image_2"){
                                                   #Horizontal conflict area (DOES NOT EXIST)
                                                   Tibble_horizontal <- tibble::tibble()
                                                   
                                                   #Vertical conflict area
                                                   Tibble_vertical <- 
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original > Overlap_value,
                                                                   Y_original > (Mask_dims[2] - Overlap_value)) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = abs(Y_original - Mask_dims[2])) %>%
                                                     dplyr::mutate(Zone = "Vertical",
                                                                   Image_type = "Image_2")
                                                   
                                                   #Shared area
                                                   Tibble_shared <-
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original <= Overlap_value,
                                                                   Y_original > (Mask_dims[2] - Overlap_value)) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = pmin(abs(X_original - 1), 
                                                                                       abs(Y_original - Mask_dims[2]),
                                                                                       na.rm = TRUE)) %>%
                                                     dplyr::mutate(Zone = "Shared",
                                                                   Image_type = "Image_2")
                                                 }
                                                 
                                                 if(Image_type == "Image_3"){
                                                   #Horizontal conflict area
                                                   Tibble_horizontal <- 
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original <= Overlap_value,
                                                                   Y_original >  Overlap_value, Y_original <= (Mask_dims[2] - Overlap_value)) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = abs(X_original - 1)) %>%
                                                     dplyr::mutate(Zone = "Horizontal",
                                                                   Image_type = "Image_3")
                                                   
                                                   #Vertical conflict area
                                                   Tibble_vertical <- 
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original > Overlap_value,
                                                                   Y_original <= Overlap_value) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = abs(Y_original - 1)) %>%
                                                     dplyr::mutate(Zone = "Vertical",
                                                                   Image_type = "Image_3")
                                                   
                                                   #Shared area
                                                   Tibble_shared <-
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original <= Overlap_value,
                                                                   Y_original <= Overlap_value) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = pmin(abs(X_original - 1), 
                                                                                       abs(Y_original - 1),
                                                                                       na.rm = TRUE)) %>%
                                                     dplyr::mutate(Zone = "Shared",
                                                                   Image_type = "Image_3")
                                                 }
                                                 
                                                 if(Image_type == "Image_4"){
                                                   #Horizontal conflict area
                                                   Tibble_horizontal <- 
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original > (Mask_dims[1] - Overlap_value),
                                                                   Y_original >  Overlap_value, Y_original <= (Mask_dims[2] - Overlap_value)) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = abs(X_original - Mask_dims[1])) %>%
                                                     dplyr::mutate(Zone = "Horizontal",
                                                                   Image_type = "Image_4")
                                                   
                                                   #Vertical conflict area (DOES NOT EXIST)
                                                   Tibble_vertical <- tibble::tibble()
                                                   
                                                   #Shared area
                                                   Tibble_shared <-
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original > (Mask_dims[1] - Overlap_value),
                                                                   Y_original <= Overlap_value) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = pmin(abs(X_original - Mask_dims[1]), 
                                                                                       abs(Y_original - 1), 
                                                                                       na.rm = TRUE)) %>%
                                                     dplyr::mutate(Zone = "Shared",
                                                                   Image_type = "Image_4")
                                                 }
                                               }
                                               
                                               #Surrounded
                                               if(!First_Column & !First_Row & !Last_Column & !Last_Row){
                                                 if(Image_type == "Image_1"){
                                                   #Horizontal conflict area (DOES NOT EXIST)
                                                   Tibble_horizontal <- tibble::tibble()
                                                   
                                                   #Vertical conflict area (DOES NOT EXIST)
                                                   Tibble_vertical <- tibble::tibble()
                                                   
                                                   #Shared area
                                                   Tibble_shared <-
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original > (Mask_dims[1] - Overlap_value),
                                                                   Y_original > (Mask_dims[2] - Overlap_value)) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = pmin(abs(X_original - Mask_dims[1]),
                                                                                       abs(Y_original - Mask_dims[2]),
                                                                                       na.rm = TRUE)) %>%
                                                     dplyr::mutate(Zone = "Shared",
                                                                   Image_type = "Image_1")
                                                 }
                                                 
                                                 if(Image_type == "Image_2"){
                                                   #Horizontal conflict area (DOES NOT EXIST)
                                                   Tibble_horizontal <- tibble::tibble()
                                                   
                                                   #Vertical conflict area
                                                   Tibble_vertical <- 
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original > Overlap_value, X_original <= (Mask_dims[1] - Overlap_value),
                                                                   Y_original > (Mask_dims[2] - Overlap_value)) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = abs(Y_original - Mask_dims[2])) %>%
                                                     dplyr::mutate(Zone = "Vertical",
                                                                   Image_type = "Image_2")
                                                   
                                                   #Shared area
                                                   Tibble_shared <-
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original <= Overlap_value,
                                                                   Y_original > (Mask_dims[2] - Overlap_value)) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = pmin(abs(X_original - 1), 
                                                                                       abs(Y_original - Mask_dims[2]),
                                                                                       na.rm = TRUE)) %>%
                                                     dplyr::mutate(Zone = "Shared",
                                                                   Image_type = "Image_2")
                                                 }
                                                 
                                                 if(Image_type == "Image_3"){
                                                   #Horizontal conflict area
                                                   Tibble_horizontal <- 
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original <= Overlap_value,
                                                                   Y_original >  Overlap_value, Y_original <= (Mask_dims[2] - Overlap_value)) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = abs(X_original - 1)) %>%
                                                     dplyr::mutate(Zone = "Horizontal",
                                                                   Image_type = "Image_3")
                                                   
                                                   #Vertical conflict area
                                                   Tibble_vertical <- 
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original > Overlap_value, X_original <= (Mask_dims[1] - Overlap_value),
                                                                   Y_original <= Overlap_value) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = abs(Y_original - 1)) %>%
                                                     dplyr::mutate(Zone = "Vertical",
                                                                   Image_type = "Image_3")
                                                   
                                                   #Shared area
                                                   Tibble_shared <-
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original <= Overlap_value,
                                                                   Y_original <= Overlap_value) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = pmin(abs(X_original - 1), 
                                                                                       abs(Y_original - 1),
                                                                                       na.rm = TRUE)) %>%
                                                     dplyr::mutate(Zone = "Shared",
                                                                   Image_type = "Image_3")
                                                 }
                                                 
                                                 if(Image_type == "Image_4"){
                                                   #Horizontal conflict area
                                                   Tibble_horizontal <- 
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original > (Mask_dims[1] - Overlap_value),
                                                                   Y_original >  Overlap_value, Y_original <= (Mask_dims[2] - Overlap_value)) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = abs(X_original - Mask_dims[1])) %>%
                                                     dplyr::mutate(Zone = "Horizontal",
                                                                   Image_type = "Image_4")
                                                   
                                                   #Vertical conflict area (DOES NOT EXIST)
                                                   Tibble_vertical <- tibble::tibble()
                                                   
                                                   #Shared area
                                                   Tibble_shared <-
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original > (Mask_dims[1] - Overlap_value),
                                                                   Y_original <= Overlap_value) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = pmin(abs(X_original - Mask_dims[1]), 
                                                                                       abs(Y_original - 1), 
                                                                                       na.rm = TRUE)) %>%
                                                     dplyr::mutate(Zone = "Shared",
                                                                   Image_type = "Image_4")
                                                 }
                                               }
                                               
                                               
                                               #Special case ROW ONLY crossroads
                                               #Special ROW only left
                                               if(First_Column & First_Row & !Last_Column & Last_Row){
                                                 if(Image_type == "Image_1"){
                                                   #Horizontal conflict area
                                                   Tibble_horizontal <- 
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original > (Mask_dims[1] - Overlap_value),
                                                                   Y_original <= (Mask_dims[2] - Overlap_value)) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = abs(X_original - Mask_dims[1])) %>%
                                                     dplyr::mutate(Zone = "Horizontal",
                                                                   Image_type = "Image_1")
                                                   
                                                   #Vertical conflict area
                                                   Tibble_vertical <- 
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original <= (Mask_dims[1] - Overlap_value),
                                                                   Y_original > (Mask_dims[2] - Overlap_value)) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = abs(Y_original - Mask_dims[2])) %>%
                                                     dplyr::mutate(Zone = "Vertical",
                                                                   Image_type = "Image_1")
                                                   
                                                   #Shared area
                                                   Tibble_shared <-
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original > (Mask_dims[1] - Overlap_value),
                                                                   Y_original > (Mask_dims[2] - Overlap_value)) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = pmin(abs(X_original - Mask_dims[1]),
                                                                                       abs(Y_original - Mask_dims[2]),
                                                                                       na.rm = TRUE)) %>%
                                                     dplyr::mutate(Zone = "Shared",
                                                                   Image_type = "Image_1")
                                                 }
                                                 
                                                 if(Image_type == "Image_2"){
                                                   #Horizontal conflict area
                                                   Tibble_horizontal <- 
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original <= Overlap_value,
                                                                   Y_original <= (Mask_dims[2] - Overlap_value)) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = X_original - 1) %>%
                                                     dplyr::mutate(Zone = "Horizontal",
                                                                   Image_type = "Image_2")
                                                   
                                                   #Vertical conflict area
                                                   Tibble_vertical <- 
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original > Overlap_value, X_original <= (Mask_dims[1] - Overlap_value),
                                                                   Y_original > (Mask_dims[2] - Overlap_value)) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = abs(Y_original - Mask_dims[2])) %>%
                                                     dplyr::mutate(Zone = "Vertical",
                                                                   Image_type = "Image_2")
                                                   
                                                   #Shared area
                                                   Tibble_shared <-
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original <= Overlap_value,
                                                                   Y_original > (Mask_dims[2] - Overlap_value)) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = pmin(abs(X_original - 1), 
                                                                                       abs(Y_original - Mask_dims[2]),
                                                                                       na.rm = TRUE)) %>%
                                                     dplyr::mutate(Zone = "Shared",
                                                                   Image_type = "Image_2")
                                                 }
                                                 
                                                 if(Image_type == "Image_3"){
                                                   #Horizontal conflict area
                                                   Tibble_horizontal <- 
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original <= Overlap_value,
                                                                   Y_original >  Overlap_value) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = abs(X_original - 1)) %>%
                                                     dplyr::mutate(Zone = "Horizontal",
                                                                   Image_type = "Image_3")
                                                   
                                                   #Vertical conflict area
                                                   Tibble_vertical <- 
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original > Overlap_value, X_original <= (Mask_dims[1] - Overlap_value),
                                                                   Y_original <= Overlap_value) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = abs(Y_original - 1)) %>%
                                                     dplyr::mutate(Zone = "Vertical",
                                                                   Image_type = "Image_3")
                                                   
                                                   #Shared area
                                                   Tibble_shared <-
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original <= Overlap_value,
                                                                   Y_original <= Overlap_value) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = pmin(abs(X_original - 1), 
                                                                                       abs(Y_original - 1),
                                                                                       na.rm = TRUE)) %>%
                                                     dplyr::mutate(Zone = "Shared",
                                                                   Image_type = "Image_3")
                                                 }
                                                 
                                                 if(Image_type == "Image_4"){
                                                   #Horizontal conflict area
                                                   Tibble_horizontal <- 
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original > (Mask_dims[1] - Overlap_value),
                                                                   Y_original >  Overlap_value) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = abs(X_original - Mask_dims[1])) %>%
                                                     dplyr::mutate(Zone = "Horizontal",
                                                                   Image_type = "Image_4")
                                                   
                                                   #Vertical conflict area
                                                   Tibble_vertical <- 
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original <= (Mask_dims[1] - Overlap_value),
                                                                   Y_original <= Overlap_value) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = abs(Y_original - 1)) %>%
                                                     dplyr::mutate(Zone = "Vertical",
                                                                   Image_type = "Image_4")
                                                   
                                                   #Shared area
                                                   Tibble_shared <-
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original > (Mask_dims[1] - Overlap_value),
                                                                   Y_original <= Overlap_value) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = pmin(abs(X_original - Mask_dims[1]), 
                                                                                       abs(Y_original - 1), 
                                                                                       na.rm = TRUE)) %>%
                                                     dplyr::mutate(Zone = "Shared",
                                                                   Image_type = "Image_4")
                                                 }
                                                 
                                               }
                                               
                                               #Special ROW only right
                                               if(!First_Column & First_Row & Last_Column & Last_Row){
                                                 if(Image_type == "Image_1"){
                                                   #Horizontal conflict area
                                                   Tibble_horizontal <- 
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original > (Mask_dims[1] - Overlap_value),
                                                                   Y_original <= (Mask_dims[2] - Overlap_value)) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = abs(X_original - Mask_dims[1])) %>%
                                                     dplyr::mutate(Zone = "Horizontal",
                                                                   Image_type = "Image_1")
                                                   
                                                   #Vertical conflict area (DOES NOT EXIST)
                                                   Tibble_vertical <- tibble::tibble()
                                                   
                                                   #Shared area
                                                   Tibble_shared <-
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original > (Mask_dims[1] - Overlap_value),
                                                                   Y_original > (Mask_dims[2] - Overlap_value)) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = pmin(abs(X_original - Mask_dims[1]),
                                                                                       abs(Y_original - Mask_dims[2]),
                                                                                       na.rm = TRUE)) %>%
                                                     dplyr::mutate(Zone = "Shared",
                                                                   Image_type = "Image_1")
                                                 }
                                                 
                                                 if(Image_type == "Image_2"){
                                                   #Horizontal conflict area
                                                   Tibble_horizontal <- 
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original <= Overlap_value,
                                                                   Y_original <= (Mask_dims[2] - Overlap_value)) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = X_original - 1) %>%
                                                     dplyr::mutate(Zone = "Horizontal",
                                                                   Image_type = "Image_2")
                                                   
                                                   #Vertical conflict area
                                                   Tibble_vertical <- 
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original > Overlap_value,
                                                                   Y_original > (Mask_dims[2] - Overlap_value)) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = abs(Y_original - Mask_dims[2])) %>%
                                                     dplyr::mutate(Zone = "Vertical",
                                                                   Image_type = "Image_2")
                                                   
                                                   #Shared area
                                                   Tibble_shared <-
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original <= Overlap_value,
                                                                   Y_original > (Mask_dims[2] - Overlap_value)) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = pmin(abs(X_original - 1), 
                                                                                       abs(Y_original - Mask_dims[2]),
                                                                                       na.rm = TRUE)) %>%
                                                     dplyr::mutate(Zone = "Shared",
                                                                   Image_type = "Image_2")
                                                 }
                                                 
                                                 if(Image_type == "Image_3"){
                                                   #Horizontal conflict area
                                                   Tibble_horizontal <- 
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original <= Overlap_value,
                                                                   Y_original >  Overlap_value) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = abs(X_original - 1)) %>%
                                                     dplyr::mutate(Zone = "Horizontal",
                                                                   Image_type = "Image_3")
                                                   
                                                   #Vertical conflict area
                                                   Tibble_vertical <- 
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original > Overlap_value,
                                                                   Y_original <= Overlap_value) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = abs(Y_original - 1)) %>%
                                                     dplyr::mutate(Zone = "Vertical",
                                                                   Image_type = "Image_3")
                                                   
                                                   #Shared area
                                                   Tibble_shared <-
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original <= Overlap_value,
                                                                   Y_original <= Overlap_value) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = pmin(abs(X_original - 1), 
                                                                                       abs(Y_original - 1),
                                                                                       na.rm = TRUE)) %>%
                                                     dplyr::mutate(Zone = "Shared",
                                                                   Image_type = "Image_3")
                                                 }
                                                 
                                                 if(Image_type == "Image_4"){
                                                   #Horizontal conflict area
                                                   Tibble_horizontal <- 
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original > (Mask_dims[1] - Overlap_value),
                                                                   Y_original >  Overlap_value) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = abs(X_original - Mask_dims[1])) %>%
                                                     dplyr::mutate(Zone = "Horizontal",
                                                                   Image_type = "Image_4")
                                                   
                                                   #Vertical conflict area (DOES NOT EXIST)
                                                   Tibble_vertical <- tibble::tibble()
                                                   
                                                   #Shared area
                                                   Tibble_shared <-
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original > (Mask_dims[1] - Overlap_value),
                                                                   Y_original <= Overlap_value) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = pmin(abs(X_original - Mask_dims[1]), 
                                                                                       abs(Y_original - 1), 
                                                                                       na.rm = TRUE)) %>%
                                                     dplyr::mutate(Zone = "Shared",
                                                                   Image_type = "Image_4")
                                                 }
                                               }
                                               
                                               #Special ROW only middle
                                               if(!First_Column & First_Row & !Last_Column & Last_Row){
                                                 if(Image_type == "Image_1"){
                                                   #Horizontal conflict area
                                                   Tibble_horizontal <- 
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original > (Mask_dims[1] - Overlap_value),
                                                                   Y_original <= (Mask_dims[2] - Overlap_value)) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = abs(X_original - Mask_dims[1])) %>%
                                                     dplyr::mutate(Zone = "Horizontal",
                                                                   Image_type = "Image_1")
                                                   
                                                   #Vertical conflict area (DOES NOT EXIST)
                                                   Tibble_vertical <- tibble::tibble()
                                                   
                                                   #Shared area
                                                   Tibble_shared <-
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original > (Mask_dims[1] - Overlap_value),
                                                                   Y_original > (Mask_dims[2] - Overlap_value)) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = pmin(abs(X_original - Mask_dims[1]),
                                                                                       abs(Y_original - Mask_dims[2]),
                                                                                       na.rm = TRUE)) %>%
                                                     dplyr::mutate(Zone = "Shared",
                                                                   Image_type = "Image_1")
                                                 }
                                                 
                                                 if(Image_type == "Image_2"){
                                                   #Horizontal conflict area
                                                   Tibble_horizontal <- 
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original <= Overlap_value,
                                                                   Y_original <= (Mask_dims[2] - Overlap_value)) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = X_original - 1) %>%
                                                     dplyr::mutate(Zone = "Horizontal",
                                                                   Image_type = "Image_2")
                                                   
                                                   #Vertical conflict area
                                                   Tibble_vertical <- 
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original > Overlap_value, X_original <= (Mask_dims[1] - Overlap_value),
                                                                   Y_original > (Mask_dims[2] - Overlap_value)) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = abs(Y_original - Mask_dims[2])) %>%
                                                     dplyr::mutate(Zone = "Vertical",
                                                                   Image_type = "Image_2")
                                                   
                                                   #Shared area
                                                   Tibble_shared <-
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original <= Overlap_value,
                                                                   Y_original > (Mask_dims[2] - Overlap_value)) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = pmin(abs(X_original - 1), 
                                                                                       abs(Y_original - Mask_dims[2]),
                                                                                       na.rm = TRUE)) %>%
                                                     dplyr::mutate(Zone = "Shared",
                                                                   Image_type = "Image_2")
                                                 }
                                                 
                                                 if(Image_type == "Image_3"){
                                                   #Horizontal conflict area
                                                   Tibble_horizontal <- 
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original <= Overlap_value,
                                                                   Y_original >  Overlap_value) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = abs(X_original - 1)) %>%
                                                     dplyr::mutate(Zone = "Horizontal",
                                                                   Image_type = "Image_3")
                                                   
                                                   #Vertical conflict area
                                                   Tibble_vertical <- 
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original > Overlap_value, X_original <= (Mask_dims[1] - Overlap_value),
                                                                   Y_original <= Overlap_value) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = abs(Y_original - 1)) %>%
                                                     dplyr::mutate(Zone = "Vertical",
                                                                   Image_type = "Image_3")
                                                   
                                                   #Shared area
                                                   Tibble_shared <-
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original <= Overlap_value,
                                                                   Y_original <= Overlap_value) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = pmin(abs(X_original - 1), 
                                                                                       abs(Y_original - 1),
                                                                                       na.rm = TRUE)) %>%
                                                     dplyr::mutate(Zone = "Shared",
                                                                   Image_type = "Image_3")
                                                 }
                                                 
                                                 if(Image_type == "Image_4"){
                                                   #Horizontal conflict area
                                                   Tibble_horizontal <- 
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original > (Mask_dims[1] - Overlap_value),
                                                                   Y_original >  Overlap_value) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = abs(X_original - Mask_dims[1])) %>%
                                                     dplyr::mutate(Zone = "Horizontal",
                                                                   Image_type = "Image_4")
                                                   
                                                   #Vertical conflict area (DOES NOT EXIS)
                                                   Tibble_vertical <- tibble::tibble()
                                                   
                                                   #Shared area
                                                   Tibble_shared <-
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original > (Mask_dims[1] - Overlap_value),
                                                                   Y_original <= Overlap_value) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = pmin(abs(X_original - Mask_dims[1]), 
                                                                                       abs(Y_original - 1), 
                                                                                       na.rm = TRUE)) %>%
                                                     dplyr::mutate(Zone = "Shared",
                                                                   Image_type = "Image_4")
                                                 }
                                               }
                                               
                                               #Special case COLUMN ONLY crossroads
                                               #Special COLUMN only top
                                               if(First_Column & First_Row & Last_Column & !Last_Row){
                                                 if(Image_type == "Image_1"){
                                                   #Horizontal conflict area
                                                   Tibble_horizontal <- 
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original > (Mask_dims[1] - Overlap_value),
                                                                   Y_original <= (Mask_dims[2] - Overlap_value)) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = abs(X_original - Mask_dims[1])) %>%
                                                     dplyr::mutate(Zone = "Horizontal",
                                                                   Image_type = "Image_1")
                                                   
                                                   #Vertical conflict area
                                                   Tibble_vertical <- 
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original <= (Mask_dims[1] - Overlap_value),
                                                                   Y_original > (Mask_dims[2] - Overlap_value)) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = abs(Y_original - Mask_dims[2])) %>%
                                                     dplyr::mutate(Zone = "Vertical",
                                                                   Image_type = "Image_1")
                                                   
                                                   #Shared area
                                                   Tibble_shared <-
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original > (Mask_dims[1] - Overlap_value),
                                                                   Y_original > (Mask_dims[2] - Overlap_value)) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = pmin(abs(X_original - Mask_dims[1]),
                                                                                       abs(Y_original - Mask_dims[2]),
                                                                                       na.rm = TRUE)) %>%
                                                     dplyr::mutate(Zone = "Shared",
                                                                   Image_type = "Image_1")
                                                 }
                                                 
                                                 if(Image_type == "Image_2"){
                                                   #Horizontal conflict area
                                                   Tibble_horizontal <- 
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original <= Overlap_value,
                                                                   Y_original <= (Mask_dims[2] - Overlap_value)) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = X_original - 1) %>%
                                                     dplyr::mutate(Zone = "Horizontal",
                                                                   Image_type = "Image_2")
                                                   
                                                   #Vertical conflict area
                                                   Tibble_vertical <- 
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original > Overlap_value,
                                                                   Y_original > (Mask_dims[2] - Overlap_value)) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = abs(Y_original - Mask_dims[2])) %>%
                                                     dplyr::mutate(Zone = "Vertical",
                                                                   Image_type = "Image_2")
                                                   
                                                   #Shared area
                                                   Tibble_shared <-
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original <= Overlap_value,
                                                                   Y_original > (Mask_dims[2] - Overlap_value)) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = pmin(abs(X_original - 1), 
                                                                                       abs(Y_original - Mask_dims[2]),
                                                                                       na.rm = TRUE)) %>%
                                                     dplyr::mutate(Zone = "Shared",
                                                                   Image_type = "Image_2")
                                                 }
                                                 
                                                 if(Image_type == "Image_3"){
                                                   #Horizontal conflict area
                                                   Tibble_horizontal <- 
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original <= Overlap_value,
                                                                   Y_original >  Overlap_value, Y_original <= (Mask_dims[2] - Overlap_value)) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = abs(X_original - 1)) %>%
                                                     dplyr::mutate(Zone = "Horizontal",
                                                                   Image_type = "Image_3")
                                                   
                                                   #Vertical conflict area
                                                   Tibble_vertical <- 
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original > Overlap_value,
                                                                   Y_original <= Overlap_value) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = abs(Y_original - 1)) %>%
                                                     dplyr::mutate(Zone = "Vertical",
                                                                   Image_type = "Image_3")
                                                   
                                                   #Shared area
                                                   Tibble_shared <-
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original <= Overlap_value,
                                                                   Y_original <= Overlap_value) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = pmin(abs(X_original - 1), 
                                                                                       abs(Y_original - 1),
                                                                                       na.rm = TRUE)) %>%
                                                     dplyr::mutate(Zone = "Shared",
                                                                   Image_type = "Image_3")
                                                 }
                                                 
                                                 if(Image_type == "Image_4"){
                                                   #Horizontal conflict area
                                                   Tibble_horizontal <- 
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original > (Mask_dims[1] - Overlap_value),
                                                                   Y_original >  Overlap_value, Y_original <= (Mask_dims[2] - Overlap_value)) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = abs(X_original - Mask_dims[1])) %>%
                                                     dplyr::mutate(Zone = "Horizontal",
                                                                   Image_type = "Image_4")
                                                   
                                                   #Vertical conflict area
                                                   Tibble_vertical <- 
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original <= (Mask_dims[1] - Overlap_value),
                                                                   Y_original <= Overlap_value) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = abs(Y_original - 1)) %>%
                                                     dplyr::mutate(Zone = "Vertical",
                                                                   Image_type = "Image_4")
                                                   
                                                   #Shared area
                                                   Tibble_shared <-
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original > (Mask_dims[1] - Overlap_value),
                                                                   Y_original <= Overlap_value) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = pmin(abs(X_original - Mask_dims[1]), 
                                                                                       abs(Y_original - 1), 
                                                                                       na.rm = TRUE)) %>%
                                                     dplyr::mutate(Zone = "Shared",
                                                                   Image_type = "Image_4")
                                                 }
                                               }
                                               
                                               #Special COLUMN only bottom
                                               if(First_Column & !First_Row & Last_Column & Last_Row){
                                                 if(Image_type == "Image_1"){
                                                   #Horizontal conflict area (DOES NOT EXIST)
                                                   Tibble_horizontal <- tibble::tibble()
                                                   
                                                   #Vertical conflict area
                                                   Tibble_vertical <- 
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original <= (Mask_dims[1] - Overlap_value),
                                                                   Y_original > (Mask_dims[2] - Overlap_value)) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = abs(Y_original - Mask_dims[2])) %>%
                                                     dplyr::mutate(Zone = "Vertical",
                                                                   Image_type = "Image_1")
                                                   
                                                   #Shared area
                                                   Tibble_shared <-
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original > (Mask_dims[1] - Overlap_value),
                                                                   Y_original > (Mask_dims[2] - Overlap_value)) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = pmin(abs(X_original - Mask_dims[1]),
                                                                                       abs(Y_original - Mask_dims[2]),
                                                                                       na.rm = TRUE)) %>%
                                                     dplyr::mutate(Zone = "Shared",
                                                                   Image_type = "Image_1")
                                                 }
                                                 
                                                 if(Image_type == "Image_2"){
                                                   #Horizontal conflict area (DOES NOT EXIST)
                                                   Tibble_horizontal <- tibble::tibble()
                                                   
                                                   #Vertical conflict area
                                                   Tibble_vertical <- 
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original > Overlap_value,
                                                                   Y_original > (Mask_dims[2] - Overlap_value)) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = abs(Y_original - Mask_dims[2])) %>%
                                                     dplyr::mutate(Zone = "Vertical",
                                                                   Image_type = "Image_2")
                                                   
                                                   #Shared area
                                                   Tibble_shared <-
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original <= Overlap_value,
                                                                   Y_original > (Mask_dims[2] - Overlap_value)) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = pmin(abs(X_original - 1), 
                                                                                       abs(Y_original - Mask_dims[2]),
                                                                                       na.rm = TRUE)) %>%
                                                     dplyr::mutate(Zone = "Shared",
                                                                   Image_type = "Image_2")
                                                 }
                                                 
                                                 if(Image_type == "Image_3"){
                                                   #Horizontal conflict area
                                                   Tibble_horizontal <- 
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original <= Overlap_value,
                                                                   Y_original >  Overlap_value) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = abs(X_original - 1)) %>%
                                                     dplyr::mutate(Zone = "Horizontal",
                                                                   Image_type = "Image_3")
                                                   
                                                   #Vertical conflict area
                                                   Tibble_vertical <- 
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original > Overlap_value,
                                                                   Y_original <= Overlap_value) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = abs(Y_original - 1)) %>%
                                                     dplyr::mutate(Zone = "Vertical",
                                                                   Image_type = "Image_3")
                                                   
                                                   #Shared area
                                                   Tibble_shared <-
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original <= Overlap_value,
                                                                   Y_original <= Overlap_value) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = pmin(abs(X_original - 1), 
                                                                                       abs(Y_original - 1),
                                                                                       na.rm = TRUE)) %>%
                                                     dplyr::mutate(Zone = "Shared",
                                                                   Image_type = "Image_3")
                                                 }
                                                 
                                                 if(Image_type == "Image_4"){
                                                   #Horizontal conflict area
                                                   Tibble_horizontal <- 
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original > (Mask_dims[1] - Overlap_value),
                                                                   Y_original >  Overlap_value) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = abs(X_original - Mask_dims[1])) %>%
                                                     dplyr::mutate(Zone = "Horizontal",
                                                                   Image_type = "Image_4")
                                                   
                                                   #Vertical conflict area
                                                   Tibble_vertical <- 
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original <= (Mask_dims[1] - Overlap_value),
                                                                   Y_original <= Overlap_value) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = abs(Y_original - 1)) %>%
                                                     dplyr::mutate(Zone = "Vertical",
                                                                   Image_type = "Image_4")
                                                   
                                                   #Shared area
                                                   Tibble_shared <-
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original > (Mask_dims[1] - Overlap_value),
                                                                   Y_original <= Overlap_value) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = pmin(abs(X_original - Mask_dims[1]), 
                                                                                       abs(Y_original - 1), 
                                                                                       na.rm = TRUE)) %>%
                                                     dplyr::mutate(Zone = "Shared",
                                                                   Image_type = "Image_4")
                                                 }
                                               }
                                              
                                               #Special COLUMN only middle
                                               if(First_Column & !First_Row & Last_Column & !Last_Row){
                                                 if(Image_type == "Image_1"){
                                                   #Horizontal conflict area (DOES NOT EXIST)
                                                   Tibble_horizontal <- tibble::tibble()
                                                   
                                                   #Vertical conflict area
                                                   Tibble_vertical <- 
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original <= (Mask_dims[1] - Overlap_value),
                                                                   Y_original > (Mask_dims[2] - Overlap_value)) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = abs(Y_original - Mask_dims[2])) %>%
                                                     dplyr::mutate(Zone = "Vertical",
                                                                   Image_type = "Image_1")
                                                   
                                                   #Shared area
                                                   Tibble_shared <-
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original > (Mask_dims[1] - Overlap_value),
                                                                   Y_original > (Mask_dims[2] - Overlap_value)) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = pmin(abs(X_original - Mask_dims[1]),
                                                                                       abs(Y_original - Mask_dims[2]),
                                                                                       na.rm = TRUE)) %>%
                                                     dplyr::mutate(Zone = "Shared",
                                                                   Image_type = "Image_1")
                                                 }
                                                 
                                                 if(Image_type == "Image_2"){
                                                   #Horizontal conflict area (DOES NOT EXIST)
                                                   Tibble_horizontal <- tibble::tibble()
                                                   
                                                   #Vertical conflict area
                                                   Tibble_vertical <- 
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original > Overlap_value,
                                                                   Y_original > (Mask_dims[2] - Overlap_value)) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = abs(Y_original - Mask_dims[2])) %>%
                                                     dplyr::mutate(Zone = "Vertical",
                                                                   Image_type = "Image_2")
                                                   
                                                   #Shared area
                                                   Tibble_shared <-
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original <= Overlap_value,
                                                                   Y_original > (Mask_dims[2] - Overlap_value)) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = pmin(abs(X_original - 1), 
                                                                                       abs(Y_original - Mask_dims[2]),
                                                                                       na.rm = TRUE)) %>%
                                                     dplyr::mutate(Zone = "Shared",
                                                                   Image_type = "Image_2")
                                                 }
                                                 
                                                 if(Image_type == "Image_3"){
                                                   #Horizontal conflict area
                                                   Tibble_horizontal <- 
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original <= Overlap_value,
                                                                   Y_original >  Overlap_value, Y_original <= (Mask_dims[2] - Overlap_value)) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = abs(X_original - 1)) %>%
                                                     dplyr::mutate(Zone = "Horizontal",
                                                                   Image_type = "Image_3")
                                                   
                                                   #Vertical conflict area
                                                   Tibble_vertical <- 
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original > Overlap_value,
                                                                   Y_original <= Overlap_value) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = abs(Y_original - 1)) %>%
                                                     dplyr::mutate(Zone = "Vertical",
                                                                   Image_type = "Image_3")
                                                   
                                                   #Shared area
                                                   Tibble_shared <-
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original <= Overlap_value,
                                                                   Y_original <= Overlap_value) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = pmin(abs(X_original - 1), 
                                                                                       abs(Y_original - 1),
                                                                                       na.rm = TRUE)) %>%
                                                     dplyr::mutate(Zone = "Shared",
                                                                   Image_type = "Image_3")
                                                 }
                                                 
                                                 if(Image_type == "Image_4"){
                                                   #Horizontal conflict area
                                                   Tibble_horizontal <- 
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original > (Mask_dims[1] - Overlap_value),
                                                                   Y_original >  Overlap_value, Y_original <= (Mask_dims[2] - Overlap_value)) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = abs(X_original - Mask_dims[1])) %>%
                                                     dplyr::mutate(Zone = "Horizontal",
                                                                   Image_type = "Image_4")
                                                   
                                                   #Vertical conflict area
                                                   Tibble_vertical <- 
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original <= (Mask_dims[1] - Overlap_value),
                                                                   Y_original <= Overlap_value) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = abs(Y_original - 1)) %>%
                                                     dplyr::mutate(Zone = "Vertical",
                                                                   Image_type = "Image_4")
                                                   
                                                   #Shared area
                                                   Tibble_shared <-
                                                     Mask_tibble %>%
                                                     dplyr::filter(X_original > (Mask_dims[1] - Overlap_value),
                                                                   Y_original <= Overlap_value) %>%
                                                     dplyr::mutate(X_final = X_original + X_Position_Min - 1,
                                                                   Y_final = Y_original + Y_Position_Min - 1) %>%
                                                     dplyr::mutate(Dist_to_edge = pmin(abs(X_original - Mask_dims[1]), 
                                                                                       abs(Y_original - 1), 
                                                                                       na.rm = TRUE)) %>%
                                                     dplyr::mutate(Zone = "Shared",
                                                                   Image_type = "Image_4")
                                                 }
                                               }
                                               
                                               #Merge all the tibbles
                                               Final_tibble <- dplyr::bind_rows(Tibble_horizontal, Tibble_vertical, Tibble_shared)
                                             }
                                             
                                             #Return the final 3 tibbles combined
                                             return(Final_tibble)
                                           })
    
    names(Mask_list_conflict_zones) <- Final_names_mask_list

    ######Modify the cell tibble in the Mask list to contain only cells in conflict and add distance to edge######
    #Do it with a for loop to modify them in place
    for(Index in 1:length(Mask_list_conflict_zones)){
      Distance_by_cell <- Mask_list_conflict_zones[[Index]] %>% dplyr::group_by(Value_Mask) %>% 
        dplyr::summarise(Dist_to_edge = min(Dist_to_edge)) %>% dplyr::ungroup()
      
      Mask_list[[Index]][["Mask_cellData"]] <- Mask_list[[Index]][["Mask_cellData"]] %>% dplyr::filter(CellID %in% Distance_by_cell$Value_Mask) %>%
        dplyr::left_join(Distance_by_cell, dplyr::join_by(CellID == Value_Mask))
      
      
    }
    ######Compute the pixels that are disputed and generate Horizontal/vertical and Shared conflictive pixel tibbles######
    #If ROW or COLUMN compute the Pixels in conflict directly
    if(!stringr::str_detect(Conflict_type, "Crossroad_")){
      
      #The shared will always be an empty tibble (it cannot exist)
      Pixels_in_conflict_Shared <- tibble::tibble()
      
      #Any of the contenders is missing return a blank tibble
      if(any(nrow(Mask_list_conflict_zones[[1]]) == 0, nrow(Mask_list_conflict_zones[[2]]) == 0)){
        Pixels_in_conflict_Horizontal_Vertical <- tibble::tibble()
      } else{
        
        #Obtain both masks and join them by any of the two
        Mask_A <- Mask_list_conflict_zones[[1]] %>% dplyr::select(Value_Mask,  X_final, Y_final) %>%
          dplyr::rename("Value_Mask_A" = "Value_Mask")
        Mask_B <- Mask_list_conflict_zones[[2]] %>% dplyr::select(Value_Mask,  X_final, Y_final) %>%
          dplyr::rename("Value_Mask_B" = "Value_Mask")
        
        Pixels_in_conflict_Horizontal_Vertical <- dplyr::inner_join(Mask_A, Mask_B, dplyr::join_by(X_final == X_final, Y_final == Y_final))
      }
      }
    
    #If CROSSROADS then compute the pixels in conflict shared in a special manner
    if(stringr::str_detect(Conflict_type, "Crossroad_")){
      #Now we need to find conflict pixels for every pair of conflict zones and for the shared area
      Conflict_area_list <- list(Horizontal_A = list(Contenders = c("Image_1", "Image_2"), Type = "Horizontal"),
                                 Horizontal_B = list(Contenders = c("Image_3", "Image_4"), Type = "Horizontal"),
                                 Vertical_A = list(Contenders = c("Image_1", "Image_4"), Type = "Vertical"),
                                 Vertical_B = list(Contenders = c("Image_2", "Image_3"), Type = "Vertical"))
      
      #Start with horizontal and vertical
      Pixels_in_conflict_Horizontal_Vertical <- 
        purrr::map_dfr(Conflict_area_list,
                       function(Conflict_match){
                         Contenders <- Conflict_match[["Contenders"]]
                         Type <- Conflict_match[["Type"]]
                         
                         #If both contenders are not present then return a simple tibble
                         if(!all(Contenders %in% names(Mask_list_conflict_zones))) return(tibble::tibble())
                         
                         #Obtain both masks
                         Mask_A <- Mask_list_conflict_zones[[Contenders[1]]] %>% dplyr::filter(Zone == Type) %>% dplyr::select(Value_Mask,  X_final, Y_final) %>%
                           dplyr::rename("Value_Mask_A" = "Value_Mask")
                         Mask_B <- Mask_list_conflict_zones[[Contenders[2]]] %>% dplyr::filter(Zone == Type) %>% dplyr::select(Value_Mask,  X_final, Y_final) %>%
                           dplyr::rename("Value_Mask_B" = "Value_Mask")
                         
                         #Compute the inner join and return
                         return(dplyr::inner_join(Mask_A, Mask_B, dplyr::join_by(X_final == X_final, Y_final == Y_final)))
                         
                       })
      
      
      #Then proceed with shared
      Combinations_shared <- expand.grid.unique(1:length(Mask_list_conflict_zones), 1:length(Mask_list_conflict_zones))
      colnames(Combinations_shared) <- c("A", "B")
      Combinations_shared <- tibble::as_tibble(Combinations_shared)
      
      Pixels_in_conflict_Shared <- 
        purrr::pmap_dfr(Combinations_shared, function(A, B){
          #Obtain both masks
          Mask_A <- Mask_list_conflict_zones[[A]] %>% dplyr::filter(Zone == "Shared") %>% dplyr::select(Value_Mask,  X_final, Y_final) %>%
            dplyr::rename("Value_Mask_A" = "Value_Mask")
          Mask_B <- Mask_list_conflict_zones[[B]] %>% dplyr::filter(Zone == "Shared") %>% dplyr::select(Value_Mask,  X_final, Y_final) %>%
            dplyr::rename("Value_Mask_B" = "Value_Mask")
          
          #If any of the masks has no rows then return empty tibble
          if(any(nrow(Mask_A) == 0, nrow(Mask_B) == 0)) return(tibble::tibble())
          
          #Else
          return(dplyr::inner_join(Mask_A, Mask_B, dplyr::join_by(X_final == X_final, Y_final == Y_final)))
        })
    }
    

    #####If no pixels in conflict the function will return blank
    if(nrow(Pixels_in_conflict_Horizontal_Vertical) + nrow(Pixels_in_conflict_Shared) == 0) return()
    
    ######Solve conflicts by removing cells######
    if(Resolution_approach == "Remove_whole_cells"){
      
      ######Generate all the cell conflicts identified######
      Conflictive_cells <- dplyr::bind_rows(Pixels_in_conflict_Horizontal_Vertical, Pixels_in_conflict_Shared) %>%
        dplyr::select(Value_Mask_A, Value_Mask_B) %>% dplyr::distinct()
      
      #####Enter the loop to compute loser cells####
      Loser_cells_final <- vector(mode = "character", length = 0L)
      Whole_cell_data <- purrr::map_dfr(Mask_list,  ~.[["Mask_cellData"]])
      Loop_iter <- 0
      
      while(nrow(Conflictive_cells) > 0 | Loop_iter == 3){ #While there are cells in conflict or a max of 3 iters (only 3 images can be potentially removed)
        #####Generate a graph to find matches between cells in dispute
        Conflictive_cells_list <- igraph::graph_from_data_frame(Conflictive_cells, directed = FALSE)
        Conflictive_cells_list <- igraph::components(Conflictive_cells_list)
        Conflictive_cells_list <- split(names(Conflictive_cells_list$membership), Conflictive_cells_list$membership)
        
        ######Dispute the match between cell types####
        ######If more than one image the absolute loser losses all its cells
        Loser_cells <-
          unlist(purrr::map(Conflictive_cells_list, function(Conflict){
            Conflict_tibble <- Conflictive_cells %>% dplyr::filter(Value_Mask_A %in% Conflict | Value_Mask_B %in% Conflict)
            
            
            
            Conflict_tibble <- dplyr::left_join(Conflict_tibble, Whole_cell_data %>% dplyr::rename("CellSize_A" = "CellSize",
                                                                                                   "Dist_to_edge_A" = "Dist_to_edge"),
                                                dplyr::join_by(Value_Mask_A == CellID))
            Conflict_tibble <- dplyr::left_join(Conflict_tibble, Whole_cell_data %>% dplyr::rename("CellSize_B" = "CellSize",
                                                                                                   "Dist_to_edge_B" = "Dist_to_edge"),
                                                dplyr::join_by(Value_Mask_B == CellID))
            
            #Apply priority score based on edge
            if(Priority_type[1] == "Dist_to_edge"){
              Conflict_tibble <- Conflict_tibble %>% dplyr::mutate(Mask_A_wins = dplyr::case_when(Dist_to_edge_A > Dist_to_edge_B ~ TRUE,
                                                                                                  Dist_to_edge_A < Dist_to_edge_B ~ FALSE,
                                                                                                  Dist_to_edge_A == Dist_to_edge_B ~ NA))
              #If second priority apply it for NA values
              if(all(length(Priority_type) == 2, Priority_type[2] == "Cell_size")){
                Conflict_tibble$Mask_A_wins[all(is.na(Conflict_tibble$Mask_A_wins),
                                                Conflict_tibble$CellSize_A[is.na(Conflict_tibble$Mask_A_wins)] > Conflict_tibble$CellSize_B[is.na(Conflict_tibble$Mask_A_wins)]
                )] <-  TRUE
                
                Conflict_tibble$Mask_A_wins[all(is.na(Conflict_tibble$Mask_A_wins),
                                                Conflict_tibble$CellSize_A[is.na(Conflict_tibble$Mask_A_wins)] < Conflict_tibble$CellSize_B[is.na(Conflict_tibble$Mask_A_wins)]
                )] <- FALSE
              }
              
              #If no second priority or still NA values then choose random
              if(anyNA(Conflict_tibble$Mask_A_wins)){
                Conflict_tibble$Mask_A_wins[is.na(Conflict_tibble$Mask_A_wins)] <- sample(c(TRUE, FALSE), size = sum(is.na(Conflict_tibble$Mask_A_wins)), replace = TRUE)
              }
            }
            
            #Apply priority based on cell size
            if(Priority_type[1] == "Cell_size"){
              Conflict_tibble <- Conflict_tibble %>% dplyr::mutate(Mask_A_wins = dplyr::case_when(CellSize_A > CellSize_B ~ TRUE,
                                                                                                  CellSize_A < CellSize_B ~ FALSE,
                                                                                                  CellSize_A == CellSize_B ~ NA))
              #If second priority apply it for NA values
              if(all(length(Priority_type) == 2, Priority_type[2] == "Dist_to_edge")){
                Conflict_tibble$Mask_A_wins[all(is.na(Conflict_tibble$Mask_A_wins),
                                                Conflict_tibble$Dist_to_edge_A[is.na(Conflict_tibble$Mask_A_wins)] > Conflict_tibble$Dist_to_edge_B[is.na(Conflict_tibble$Mask_A_wins)]
                )] <-  TRUE
                
                Conflict_tibble$Mask_A_wins[all(is.na(Conflict_tibble$Mask_A_wins),
                                                Conflict_tibble$Dist_to_edge_A[is.na(Conflict_tibble$Mask_A_wins)] < Conflict_tibble$Dist_to_edge_B[is.na(Conflict_tibble$Mask_A_wins)]
                )] <- FALSE
              }
              
              #If no second priority or still NA values then choose random
              if(anyNA(Conflict_tibble$Mask_A_wins)){
                Conflict_tibble$Mask_A_wins[is.na(Conflict_tibble$Mask_A_wins)] <- sample(c(TRUE, FALSE), size = sum(is.na(Conflict_tibble$Mask_A_wins)), replace = TRUE)
              }
            }
            
            #Count the votes 
            Conflict_tibble <- Conflict_tibble %>% dplyr::mutate(Mask_vote = dplyr::case_when(Mask_A_wins ~ Value_Mask_A,
                                                                                              !Mask_A_wins ~ Value_Mask_B)) %>%
              dplyr::mutate(Mask_vote = substr(Mask_vote, start = 0, stop = 7))
            
            Votes_tibble <- Conflict_tibble %>% dplyr::count(Mask_vote)
            
            Votes_tibble <- dplyr::left_join(tibble::tibble(Mask_vote = unique(c(substr(Conflict_tibble$Value_Mask_A, start = 0, stop = 7), 
                                                                                 substr(Conflict_tibble$Value_Mask_B, start = 0, stop = 7)))),
                                             Votes_tibble,
                                             by = "Mask_vote")
            
            Votes_tibble$n[is.na(Votes_tibble$n)] <- 0
            Votes_tibble <- Votes_tibble %>% dplyr::arrange(n)
            
            #Select the loser mask and all the cells belonging to that match
            Loser_mask <- Votes_tibble[[1,1]]
            Loser_cells <- 
              unique(c(Conflict_tibble$Value_Mask_A, Conflict_tibble$Value_Mask_B)[
                stringr::str_detect(string = c(Conflict_tibble$Value_Mask_A, Conflict_tibble$Value_Mask_B), pattern = Loser_mask)
              ])
            
            #Return the loser cells
            return(Loser_cells)
            
          }))
        
        #Add these loser cells to the loser cells final
        Loser_cells_final <- unique(c(Loser_cells_final, Loser_cells))
        #Remove loser cells from conflictive cells
        Conflictive_cells <- Conflictive_cells %>% dplyr::filter(!Value_Mask_A %in% Loser_cells_final,
                                                                 !Value_Mask_B %in% Loser_cells_final)
        
        #Update loser cells list
        Loop_iter <- Loop_iter + 1
      }
      
      ######With the loser cells compute the final erosion masks####
      #with the loser cells generate the tibble
      Whole_cells_to_remove <- tibble::tibble(Image_type = substr(Loser_cells_final, start = 0, stop = 7),
                                              CellID = substr(Loser_cells_final, start = 9, stop = nchar(Loser_cells_final))) %>% dplyr::arrange(Image_type)
      
      Final_erosion_mask_list <- purrr::map(1:length(Mask_list),
                                            function(Mask_index){
                                              #Get the image type
                                              Image_type <- names(Mask_list)[Mask_index]
                                              
                                              #Get the cells from that image that should be removed
                                              Cells_to_remove <- Whole_cells_to_remove$CellID[Whole_cells_to_remove$Image_type == Image_type]
                                              
                                              if(length(Cells_to_remove) > 0){
                                                Erosion_mask <- Mask_list[[Mask_index]][["Mask"]]
                                                Erosion_mask[]<- as.numeric(!is.na(match(Erosion_mask[], Cells_to_remove)))
                                              } else return(NULL) #If no cells need to be removed return NULL
                                              
                                              return(Erosion_mask > 0)
                                            })
      
      names(Final_erosion_mask_list) <- Final_names_mask_list
    }
    
    ######Solve conflicts by eroding pixels######
    if(Resolution_approach == "Erode_pixels"){
      ##Obtain whole cell data and whole pixel data
      Whole_cell_data <- purrr::map_dfr(Mask_list,  ~.[["Mask_cellData"]])
      Whole_pixel_data <- purrr::map_dfr(Mask_list_conflict_zones,  ~.)
      
      #Work first with Horizontal_vertical pixels
      #If horizontal vertical is absent don't do anything else compute as usual
      if(nrow(Pixels_in_conflict_Horizontal_Vertical) == 0) Pixels_to_remove_horizontal_vertical <- tibble::tibble() else{
        
        #Get the values of the closest pixels and the size of the cell where the pixel belongs
        Conflict_tibble <- dplyr::left_join(Pixels_in_conflict_Horizontal_Vertical, 
                                            Whole_pixel_data %>% dplyr::rename("Value_Mask_A" = "Value_Mask",
                                                                               "Dist_to_edge_A" = "Dist_to_edge") %>%
                                              dplyr::select(Value_Mask_A, Dist_to_edge_A, X_final, Y_final),
                                            dplyr::join_by(Value_Mask_A == Value_Mask_A,
                                                           X_final == X_final,
                                                           Y_final == Y_final))
        Conflict_tibble <- dplyr::left_join(Conflict_tibble, 
                                            Whole_pixel_data %>% dplyr::rename("Value_Mask_B" = "Value_Mask",
                                                                               "Dist_to_edge_B" = "Dist_to_edge") %>%
                                              dplyr::select(Value_Mask_B, Dist_to_edge_B, X_final, Y_final),
                                            dplyr::join_by(Value_Mask_B == Value_Mask_B,
                                                           X_final == X_final,
                                                           Y_final == Y_final))
        
        Conflict_tibble <- dplyr::left_join(Conflict_tibble, Whole_cell_data %>% 
                                              dplyr::rename("CellSize_A" = "CellSize") %>%
                                              dplyr::select(CellID, CellSize_A),
                                            dplyr::join_by(Value_Mask_A == CellID))
        Conflict_tibble <- dplyr::left_join(Conflict_tibble, Whole_cell_data %>% 
                                              dplyr::rename("CellSize_B" = "CellSize") %>%
                                              dplyr::select(CellID, CellSize_B),
                                            dplyr::join_by(Value_Mask_B == CellID))
        
        #Compute results based on user preferences
        #Apply priority score based on edge
        if(Priority_type[1] == "Dist_to_edge"){
          Conflict_tibble <- Conflict_tibble %>% dplyr::mutate(Mask_A_wins = dplyr::case_when(Dist_to_edge_A > Dist_to_edge_B ~ TRUE,
                                                                                              Dist_to_edge_A < Dist_to_edge_B ~ FALSE,
                                                                                              Dist_to_edge_A == Dist_to_edge_B ~ NA))
          #If second priority apply it for NA values
          if(all(length(Priority_type) == 2, Priority_type[2] == "Cell_size")){
            Conflict_tibble$Mask_A_wins[all(is.na(Conflict_tibble$Mask_A_wins),
                                            Conflict_tibble$CellSize_A[is.na(Conflict_tibble$Mask_A_wins)] > Conflict_tibble$CellSize_B[is.na(Conflict_tibble$Mask_A_wins)]
            )] <-  TRUE
            
            Conflict_tibble$Mask_A_wins[all(is.na(Conflict_tibble$Mask_A_wins),
                                            Conflict_tibble$CellSize_A[is.na(Conflict_tibble$Mask_A_wins)] < Conflict_tibble$CellSize_B[is.na(Conflict_tibble$Mask_A_wins)]
            )] <- FALSE
          }
          
          #If no second priority or still NA values then choose random
          if(anyNA(Conflict_tibble$Mask_A_wins)){
            Conflict_tibble$Mask_A_wins[is.na(Conflict_tibble$Mask_A_wins)] <- sample(c(TRUE, FALSE), size = sum(is.na(Conflict_tibble$Mask_A_wins)), replace = TRUE)
          }
        }
        
        #Apply priority based on cell size
        if(Priority_type[1] == "Cell_size"){
          Conflict_tibble <- Conflict_tibble %>% dplyr::mutate(Mask_A_wins = dplyr::case_when(CellSize_A > CellSize_B ~ TRUE,
                                                                                              CellSize_A < CellSize_B ~ FALSE,
                                                                                              CellSize_A == CellSize_B ~ NA))
          #If second priority apply it for NA values
          if(all(length(Priority_type) == 2, Priority_type[2] == "Dist_to_edge")){
            Conflict_tibble$Mask_A_wins[all(is.na(Conflict_tibble$Mask_A_wins),
                                            Conflict_tibble$Dist_to_edge_A[is.na(Conflict_tibble$Mask_A_wins)] > Conflict_tibble$Dist_to_edge_B[is.na(Conflict_tibble$Mask_A_wins)]
            )] <-  TRUE
            
            Conflict_tibble$Mask_A_wins[all(is.na(Conflict_tibble$Mask_A_wins),
                                            Conflict_tibble$Dist_to_edge_A[is.na(Conflict_tibble$Mask_A_wins)] < Conflict_tibble$Dist_to_edge_B[is.na(Conflict_tibble$Mask_A_wins)]
            )] <- FALSE
          }
          
          #If no second priority or still NA values then choose random
          if(anyNA(Conflict_tibble$Mask_A_wins)){
            Conflict_tibble$Mask_A_wins[is.na(Conflict_tibble$Mask_A_wins)] <- sample(c(TRUE, FALSE), size = sum(is.na(Conflict_tibble$Mask_A_wins)), replace = TRUE)
          }
        }
        
        #Return the loser pixel with its cell identity
        Pixels_to_remove_horizontal_vertical <- Conflict_tibble %>% dplyr::mutate(CellID = dplyr::case_when(Mask_A_wins ~ Value_Mask_B,
                                                                                                            !Mask_A_wins ~ Value_Mask_A)) %>%
          dplyr::select(CellID, X_final, Y_final) %>% dplyr::mutate(Image = substr(CellID, start = 0, stop = 7))
       
      }
      
      #For pixels in shared area we will generate a graph according to X Y values
      if(nrow(Pixels_in_conflict_Shared) == 0) Pixels_to_remove_shared <- tibble::tibble() else{
        #Obtain the X and Y positions and the cells involved and place them in a list inside the tibble
        Conflict_tibble <- Pixels_in_conflict_Shared %>%
          dplyr::group_by(X_final, Y_final) %>%
          summarise(
            CellID = list(unique(c(Value_Mask_A, Value_Mask_B))),
            .groups = "drop"
          )
        
        #Prepare a function for NA if ties
        which.max.na_ties <- function(x) {
          if (length(x) == 0 || all(is.na(x))) return(NA_integer_)
          
          i <- which(x == max(x, na.rm = TRUE))
          if (length(i) == 1) i else NA_integer_
        }
        
        Pixels_to_remove_shared <-
          purrr::map_dfr(1:nrow(Conflict_tibble), function(Index_row){
            #Arrange the pixels in conflict data
            Interim_tibble <- tibble::tibble(X_final = as.integer(Conflict_tibble$X_final[Index_row]),
                                             Y_final = as.integer(Conflict_tibble$Y_final[Index_row]),
                                             CellID = unlist(Conflict_tibble$CellID[Index_row]))
            #Obtain pixel data to perform winner selection
            Interim_tibble <- dplyr::left_join(Interim_tibble, Whole_cell_data %>% dplyr::select(CellID, CellSize), dplyr::join_by(CellID == CellID))
            Interim_tibble <- dplyr::left_join(Interim_tibble, Whole_pixel_data %>% dplyr::filter(X_final == Interim_tibble$X_final[1], Y_final == Interim_tibble$Y_final[1]) %>%
                                                 dplyr::select(Value_Mask, Dist_to_edge), 
                                               dplyr::join_by(CellID == Value_Mask))
            #Apply priority score based on edge
            if(Priority_type[1] == "Dist_to_edge"){
              Winner_pixel_index <- which.max.na_ties(Interim_tibble$Dist_to_edge)
              #If NA try alternative or choose random winner
              if(is.na(Winner_pixel_index)){
                if(all(length(Priority_type) == 2, Priority_type[2] == "Cell_size")) Winner_pixel_index <- which.max.na_ties(Interim_tibble$CellSize) else{
                  Winner_pixel_index <- sample(1:nrow(Interim_tibble), size = 1)
                } 
              }
              Loser_pixels <- Interim_tibble[-Winner_pixel_index,]

            }
            #Apply priority score based on cell size
            if(Priority_type[1] == "Cell_size"){
              Winner_pixel_index <- which.max.na_ties(Interim_tibble$CellSize)
              #If NA try alternative or choose random winner
              if(is.na(Winner_pixel_index)){
                if(all(length(Priority_type) == 2, Priority_type[2] == "Dist_to_edge")) Winner_pixel_index <- which.max.na_ties(Interim_tibble$Dist_to_edge) else{
                  Winner_pixel_index <- sample(1:nrow(Interim_tibble), size = 1)
                } 
              }
              Loser_pixels <- Interim_tibble[-Winner_pixel_index,]
            }
            
            return(Loser_pixels)
          }, .progress = TRUE)
        
        Pixels_to_remove_shared <- Pixels_to_remove_shared %>% dplyr::select(CellID, X_final, Y_final) %>% dplyr::mutate(Image = substr(CellID, start = 0, stop = 7))
      
      }
      
      #Merge both tibbles containing pixels to remove
      Pixels_to_remove <- dplyr::bind_rows(Pixels_to_remove_horizontal_vertical, Pixels_to_remove_shared)
      
      #Proceed to obtain all data from pixels to remove (Image, X_original and Y_original)
      Pixels_to_remove <- dplyr::left_join(Pixels_to_remove, 
                                           Whole_pixel_data %>% 
                                             dplyr::filter(Value_Mask %in% Pixels_to_remove$CellID) %>%
                                             dplyr::select(X_original, Y_original, X_final, Y_final, Value_Mask, Zone, Image_type),
                                           dplyr::join_by(CellID == Value_Mask,
                                                          X_final == X_final,
                                                          Y_final == Y_final)) %>% dplyr::select(-Image)
      
      
      #Prepare the erosion mask list
      Final_erosion_mask_list <- 
        purrr::map(1:length(Mask_list),
                   function(Mask_index){
                     #Get the name of the image
                     Image_type_index <- names(Mask_list)[Mask_index]
                     
                     #Get the mask dims
                     Dimensions_Mask <- Mask_list[[Mask_index]][["Dimensions_Mask"]]
                     
                     #Get the pixels to remove
                     Pixels_to_remove <- Pixels_to_remove %>% dplyr::filter(Image_type == Image_type_index) %>%
                       dplyr::select(X_original, Y_original) %>% dplyr::mutate(Value_Mask = 1)
                     
                     #If no pixels need to be removed return null
                     if(nrow(Pixels_to_remove) < 1) return(NULL)
                     
                     #Generate a tibble that will then be reconstituted to an image
                     Erosion_mask <- tibble::tibble(X_original = as.integer(rep(1:Dimensions_Mask[1], times = Dimensions_Mask[2])),
                                                    Y_original = as.integer(rep(1:Dimensions_Mask[2], each = Dimensions_Mask[1])),
                                                    Value_Mask = as.integer(0))
                     
                     #Substitute values
                     Erosion_mask <- dplyr::rows_update(Erosion_mask, Pixels_to_remove, by = c("X_original", "Y_original"))
                     #Rebuild image
                     Erosion_mask <-  matrix(Erosion_mask$Value_Mask, 
                                             nrow = Dimensions_Mask[1],
                                             ncol = Dimensions_Mask[2],
                                             byrow = FALSE)
                     Erosion_mask <- EBImage::as.Image(Erosion_mask)
                     return(Erosion_mask > 0)
                   })
      names(Final_erosion_mask_list) <- Final_names_mask_list
      
      
      #If minimum cell size requirements provided then proceed to generate masks of whole cells to remove
      if(!is.null(Min_cell_size)){
        #Compute the size of the erosion for every cell
        Pixel_size_erosion <- Pixels_to_remove %>% dplyr::count(CellID)
        
        #Add the full size of the cell and compute the final size after erosion
        Pixel_size_erosion <- dplyr::left_join(Pixel_size_erosion, 
                                               Whole_cell_data %>% dplyr::select(CellID, CellSize),
                                               dplyr::join_by(CellID == CellID)) %>%
          dplyr::mutate(Final_size = CellSize - n) %>%
          dplyr::mutate(Remove_whole_cell = dplyr::case_when(Final_size >= Min_cell_size ~ FALSE,
                                                             Final_size < Min_cell_size &  Final_size != 0 ~ TRUE,
                                                             Final_size == 0 ~ FALSE #If 0 no need the erosion will already remove the whole cell
                                                             )) 
        
        Pixel_size_erosion <- Pixel_size_erosion %>% dplyr::filter(Remove_whole_cell) %>% 
          dplyr::mutate(Image_type = substr(CellID, start = 0, stop = 7),
                        CellID = as.integer(substr(CellID, start = 9, stop = nchar(CellID))))
        
        Whole_cell_remove_masks <- 
          purrr::map(1:length(Mask_list), 
                     function(Mask_index){
                       #Get the name of the image
                       Image_type <- names(Mask_list)[Mask_index]
                       
                       #Get the cell mask ID to be removed
                       Cells_to_remove_mask <- Pixel_size_erosion$CellID[Pixel_size_erosion$Image_type == Image_type]
                       
                       #If no cell to remove then return NULL
                       if(length(Cells_to_remove_mask) < 1) return(NULL)
                       
                       #Get the original mask
                       Mask <- Mask_list[[Mask_index]][["Mask"]]
                       
                       #prepare the cell mask of whole cells
                       Mask[] <- as.numeric(!is.na(match(Mask[], Cells_to_remove_mask)))
                       Mask <- Mask > 0
                       return(Mask)
                     })
        
        names(Whole_cell_remove_masks) <- Final_names_mask_list
        
        #Bind both results
        Final_erosion_mask_list <- 
          purrr::map2(.x = Final_erosion_mask_list, .y = Whole_cell_remove_masks,
                      function(.x, .y){
                        #If both are null return NULL
                        if(all(is.null(.x), is.null(.y))) return(NULL)
                        #If no whole cells need to be removed return pixel erosion mask
                        if(is.null(.y)) return(.x)
                        
                        #Else compute .x or .y
                        return(.x | .y)
                      })
        names(Final_erosion_mask_list) <- Final_names_mask_list
        
      }

    }
    
    #####Remove null image mask####
    Final_erosion_mask_list <- Final_erosion_mask_list[!purrr::map_lgl(Final_erosion_mask_list, ~is.null(.))]

    #####Substite the names of the list by the final names that will be used#####
    names(Final_erosion_mask_list) <- Image_param_list$Mask_names[match(names(Final_erosion_mask_list), Image_param_list$Image)]
    
    
    #####WRITE RESULTING MASKS####
    purrr::walk(1:length(Final_erosion_mask_list),
                function(Mask_index){
                  #Get the mask name
                  Mask_name <- names(Final_erosion_mask_list)[Mask_index]
                  
                  #Get mask
                  Mask <- Final_erosion_mask_list[[Mask_index]]
                  
                  #If already present import the mask and find intersect
                  if(Mask_name %in% dir(Mask_refinement_path)){
                    Already_mask <- EBImage::readImage(paste0(Mask_refinement_path, "/", Mask_name)) > 0
                    Mask <- Mask | Already_mask
                  }
                  
                  #Write mask
                  EBImage::writeImage(Mask, paste0(Mask_refinement_path, "/", Mask_name), bits.per.sample = 16, compression = "LZW")
                  
                })
    
    #Perform GC
    gc()
    
  }

