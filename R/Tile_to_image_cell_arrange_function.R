#' Arrange cell coordinates according to position in whole image using tile information
#'
#' `Tile_to_image_cell_arrange_function()` modifies cell X and Y coordinates to reflect cell position in the entire image. For cell segmentation pipelines that
#' are based on tiled images (those processed using the [Image_tile_deconstruction_function()] function), the function will switch the tile-based X Y coordinates to image-based coordinates.
#' The function performs the operation based on the tile X Y coordinates info contained in the Subject_Names variable. If names have been modified the function could return unexpected results. In addition,
#' the Subject_Names variable of the dataset is also modified to reflect entire image information. If tiles present overlap, duplicated cells are removed. Further, cells close to the tile edge can
#' also be discarded according to user preferences using the 'Dist_to_edge' argument.
#'
#' @param DATA Dataframe or tibble with cell features. This data should have been processed at some point using [Data_arrange_function()].
#' @param Dist_to_edge (OPTIONAL) For overlapping tiles, cells closer to tile edge than Dist_to_edge will be removed.
#' @param Tile_overlap (OPTIONAL) A numeric value indicating the tile overlap. If provided, automatic tile overlap detection will be aborted. If a value of 0 is provided, no overlap correction will be performed at all.
#'
#' @returns A tibble containing cell features with modified X and Y coordinates and Subject_Names information.
#'
#' @seealso [Image_tile_deconstruction_function()], [Image_from_tile_rebuilder()]
#'
#' @examples
#' \dontrun{
#' #Returns an error since Subject_Names in dataframe does not contain adequate tile information
#' Tile_to_image_cell_arrange_function(
#'  DATA = CSM_Arrangedcellfeaturedata_test,
#'  Dist_to_edge = 0
#' )
#'
#' }
#'
#'
#'
#' @export
Tile_to_image_cell_arrange_function <-
  function(DATA,
           Dist_to_edge = NULL,
           Tile_overlap = NULL){
    #Check arguments
    if(!identical(names(DATA)[1:4],  c("Cell_no", "X", "Y", "Subject_Names")))  stop("DATA provided should have been processed using the DATA_arrange_function")

    #Check that all images come from tiles
    Image_from_tile <- stringr::str_detect(unique(DATA$Subject_Names), "\\(Tile")
    if(!all(Image_from_tile)){
      Problematic_images <- unique(DATA$Subject_Names)[!Image_from_tile]
      stop(paste0("The followin images are not tiles or their file names have been severely modified. Please check before running the function: ",
                  stringr::str_c(Problematic_images, collapse = ", ")))
    }

    #Now lets generate the look-up table for every cell
    Look_up_table <- tibble(Old_Cell_no = DATA$Cell_no) #Old data
    Look_up_table$Old_X <- DATA$X #Old data
    Look_up_table$Old_Y <- DATA$Y #Old data
    Look_up_table$Old_subject_Names <- sub("(.*)\\.[^.]+$", "\\1", DATA$Subject_Names) #Old data

    Look_up_table$New_Subject_Names <- gsub("\\([^)]*\\)", "", DATA$Subject_Names)
    Look_up_table$Tile_info <- regmatches(DATA$Subject_Names, regexpr("(?<=\\()[^)]*(?=\\))", DATA$Subject_Names, perl = TRUE))
    Tile_info_list <- stringr::str_split(Look_up_table$Tile_info, "-")
    Look_up_table$Tile_ID <-purrr::map_chr(Tile_info_list, ~.[[1]])
    Look_up_table$X_Position_Min <-purrr::map_int(Tile_info_list, ~as.integer(.[[2]]))
    Look_up_table$X_Position_Max <-purrr::map_int(Tile_info_list, ~as.integer(.[[3]]))
    Look_up_table$Y_Position_Min <-purrr::map_int(Tile_info_list, ~as.integer(.[[4]]))
    Look_up_table$Y_Position_Max <-purrr::map_int(Tile_info_list, ~as.integer(.[[5]]))

    #Basically we now add the X_position_min and Y_position_min to every cell according to the tile
    Look_up_table$New_X <- Look_up_table$Old_X + Look_up_table$X_Position_Min - 1
    Look_up_table$New_Y <- Look_up_table$Old_Y + Look_up_table$Y_Position_Min - 1


    #If user has not provided an overlap value, check if there is overlap between tiles by analyzing the first two tiles between the first image
    if(is.null(Tile_overlap)){
      #Calculate if there is overlap and the amount of overlap by image
      Overlap_info_list <-
        purrr::map(unique(Look_up_table$New_Subject_Names), function(Image){

          #Obtain individual image tiles
          Tiles <- Look_up_table %>% dplyr::filter(New_Subject_Names == Image)

          #If single tile return FALSE (No overlap)
          if(nrow(Tiles) == 1) return(list(Overlap = FALSE,
                                           Overlap_size = NA)
          )

          #Calculate the min_x max_x and min_y max_y
          min_x <- min(Tiles$X_Position_Min)
          max_x <- max(Tiles$X_Position_Max)
          Size_x <- max_x - min_x

          min_y <- min(Tiles$Y_Position_Min)
          max_y <- max(Tiles$Y_Position_Max)
          Size_y <- max_y - min_y

          #Calculate the expected tile number based on tile size
          Expected_size_x <- Tiles$X_Position_Max[1] - Tiles$X_Position_Min[1] + 1
          Expected_size_y <- Tiles$Y_Position_Max[1] - Tiles$Y_Position_Min[1] + 1

          Expected_tiles_n_x <- ceiling(Size_x/Expected_size_x)
          Expected_tiles_n_y <- ceiling(Size_y/Expected_size_y)

          #Claculate the observed number of tiles
          Observed_tiles_n_x <- length(unique(Tiles$X_Position_Min))
          Observed_tiles_n_y <- length(unique(Tiles$Y_Position_Min))

          #If more than expected then there is overlap
          if(any(Observed_tiles_n_x > Expected_tiles_n_x,
                 Observed_tiles_n_y > Expected_tiles_n_y)){

            #Obtain the X positions and calculate overlap
            Unique_x_min <- sort(unique(Tiles$X_Position_Min))[-1]
            Unique_x_max <- sort(unique(Tiles$X_Position_Max))[-length(unique(Tiles$X_Position_Max))]

            Overlap_X_value <- min(Unique_x_max - Unique_x_min) + 1

            Unique_y_min <- sort(unique(Tiles$Y_Position_Min))[-1]
            Unique_y_max <- sort(unique(Tiles$Y_Position_Max))[-length(unique(Tiles$Y_Position_Max))]

            Overlap_Y_value <- min(Unique_y_max - Unique_y_min) + 1

            #Return the final result
            return(list(Overlap = TRUE,
                        Overlap_size = max(c(Overlap_X_value, Overlap_Y_value))))
          }
          else{
            return(list(Overlap = FALSE,
                        Overlap_size = NA)
            )
          }
        }
        )

      #Check that all overlaps are the same size if not, throw error
      if(any(purrr::map_lgl(Overlap_info_list, function(Image) Image[["Overlap"]]))){
        Overlap_sizes <- na.omit(map_dbl(Overlap_info_list, function(Image) Image[["Overlap_size"]]))
        if(length(unique(Overlap_sizes)) > 1) stop("Overlap size deteceted in the dataset varies from image to image, please review tiling process, or specify Tile_overlap argument.")

        Overlap_size <- unique(Overlap_sizes)

        Overlap_presence <- TRUE
      }

      #If no Overlap then set presence to FALSE
      else{
        Overlap_presence <- FALSE
      }
    }

    else{
      if(!all(Tile_overlap >= 0, is.numeric(Tile_overlap), Tile_overlap%%1 == 0)) stop("Tile_overlap must be a integer value >= 0")
      if(Tile_overlap == 0){
        message("Tile overlap set to 0, aborting automatic tile overlap evaluation. Arranging cell data without tile overlap correction...")
        Overlap_presence <- FALSE
      }

      else{
        message("Tile overlap value manually provided. Aborting automatic tile overlap evaluation. Arranging cell data with tile overlap correction...")
        Overlap_presence <- TRUE
        Overlap_size <- Tile_overlap
      }
    }

    #If there is overlap presence, calculate the overlap size and calculate the conflict zone
    if(Overlap_presence){

      #Print a message saying that an overlap has been found and will be resolved
      message(paste0("An overlap of ", Overlap_size, " pixels has been identified. ", "Conflicts between overlapping cell masks will be resolved. During this process some cells may be eliminated from DATA"))

      if(!all(Dist_to_edge >= 0, Dist_to_edge%%1 == 0, Dist_to_edge < Overlap_size)) stop(paste0("Dist_to_edge must be a numeric integer larger than or equal to 0 and smaller than ", Overlap_size))

      Look_up_table <-
        purrr::map_dfr(unique(Look_up_table$New_Subject_Names), function(Image_name){
          #Select the individual images
          Image_data <- Look_up_table %>% dplyr::filter(New_Subject_Names == Image_name)

          #Calculate the image edges (these do belong to the conflict area)
          Min_Image_X <- min(Image_data$X_Position_Min)
          Min_Image_Y <- min(Image_data$Y_Position_Min)
          Max_Image_X <- max(Image_data$X_Position_Max)
          Max_Image_Y <- max(Image_data$Y_Position_Max)

          #Calculate the row and col bands were cells in conflict are present
          X_conflict <- unique(Image_data$X_Position_Min)
          X_conflict <- X_conflict[-which(X_conflict == Min_Image_X)]
          X_conflict_tibble <- tibble(From_x = X_conflict, To_x = X_conflict + Overlap_size -1)

          Y_conflict <- unique(Image_data$Y_Position_Min)
          Y_conflict <- Y_conflict[-which(Y_conflict == Min_Image_Y)]
          Y_conflict_tibble <- tibble(From_y = Y_conflict, To_y = Y_conflict + Overlap_size -1)

          Cells_in_conflict_area_X <- purrr::pmap(X_conflict_tibble, function(From_x, To_x) Image_data$New_X >= From_x & Image_data$New_X <= To_x)
          Cells_in_conflict_area_Y <- purrr::pmap(Y_conflict_tibble, function(From_y, To_y) Image_data$New_Y >= From_y & Image_data$New_Y <= To_y)

          #reduce to a single vector
          Cell_in_conflict_vector <- purrr::reduce(c(Cells_in_conflict_area_X, Cells_in_conflict_area_Y), function(A, B) A | B)

          #First we get the cells not in conflict and in conflict
          Cells_out_of_conflict <- Image_data %>% dplyr::filter(!unlist(Cell_in_conflict_vector))
          Cells_in_conflict <- Image_data %>% dplyr::filter(unlist(Cell_in_conflict_vector))

          #Now we remove cells that are in the conflict zone and very close to the edges according to their tiles
          Cells_in_conflict <-
            purrr::map_dfr(unique(Cells_in_conflict$Tile_ID), function(Tile){
              Interim <- Cells_in_conflict %>% dplyr::filter(Tile_ID == Tile)
              Interim %>% dplyr::filter(abs(New_X - X_Position_Min) >= Dist_to_edge,
                                        abs(New_X - X_Position_Max) >= Dist_to_edge,
                                        abs(New_Y - Y_Position_Min) >= Dist_to_edge,
                                        abs(New_Y - Y_Position_Max) >= Dist_to_edge)
            })

          #Remove duplicated cells
          Cells_in_conflict$Duplicated <- unlist(Cells_in_conflict %>% dplyr::select(New_X, New_Y) %>%
                                                   dplyr::mutate(Is_duplicated = duplicated(.)) %>%
                                                   dplyr::select(Is_duplicated))
          Cells_in_conflict <- Cells_in_conflict %>% dplyr::filter(!Duplicated) %>% dplyr::select(-Duplicated)

          Final_tibble <- dplyr::bind_rows(Cells_out_of_conflict, Cells_in_conflict) %>% dplyr::arrange(New_X, New_Y)
          return(Final_tibble)
        }, .progress = list(clear = F,
                            name = "Correcting edge bias",
                            show_after = 1,
                            type = "iterator")
        )
    }

    #Now generate the new data
    New_DATA <- Look_up_table %>% dplyr::select(New_X, New_Y, New_Subject_Names, Old_Cell_no)
    New_DATA <-dplyr::left_join(New_DATA, DATA %>% dplyr::select(-c(2:4)), by = dplyr::join_by(Old_Cell_no == Cell_no)) %>% dplyr::select(-Old_Cell_no)
    names(New_DATA)[c(1:3)] <- c("X", "Y", "Subject_Names")

    New_DATA <- New_DATA %>% dplyr::mutate(Cell_no = stringr::str_c("CELL", as.character(unlist(purrr::map(
      purrr::map_dbl(unique(New_DATA$Subject_Names), function(x) {
        nrow(New_DATA %>% dplyr::filter(Subject_Names == x))
      }), function(x) {
        1:x
      }, .progress = list(clear = F,
                          name = "Adding new Cell ID",
                          show_after = 1,
                          type = "iterator")))), sep ="_"))
    #Include the Subject_Name in the Cell_no
    New_DATA$Cell_no <- stringr::str_c(New_DATA$Cell_no, New_DATA$Subject_Names, sep = "__")
    New_DATA <- New_DATA[c(ncol(New_DATA), 1:(ncol(New_DATA)-1))]
    return(New_DATA)
  }
