#' Calculates the average texture features of a single cell type
#'
#' The function calculates the average texture features of a single cell type for every image. Image must have been previously tiled using [Image_tiling_processing_function()].
#' The tiled image is first transformed into a raster like object, then the grey level co-ocurrence matrix is calculated and texture features are obtained and averaged.
#'
#' @param Tiled_images A list containing tiled images obtained using [Image_tiling_processing_function()].
#' @param Variable_name A character value indicating the column in the tiled dataset where the cell label to analyzed is located ("Phenotype" by default).
#' @param Phenotype_included A character value indicating which cell phenotype (or neighborhood type) to analyze.
#'
#' @details
#' Metrics are calculated using the glcm R package.
#'
#' @returns A tibble containing averaged texture features by image.
#'
#' @seealso [Tiled_images_graphicator()]
#'
#' @examples
#' \dontrun{
#' #Divide cells into tiles-------------------------------------
#' Tiled_Images <-
#' Image_tiling_processing_function(
#'    N_cores = 2,
#'    DATA = CSM_Phenotypecell_test,
#'    Tile_width = 125,
#'    Tile_height = 125,
#'    Variables_to_keep = "Phenotype"
#')
#'
#'#Calculate average texture features for a given cell phenotype---
#'Texture_features_calculator(
#'     Tiled_images = Tiled_Images,
#'     Variable_name = "Phenotype",
#'     Phenotype_included = "CD8_GZMBneg"
#')
#' }
#'
#' @export

Texture_features_calculator <-
  function(Tiled_images,
           Variable_name = "Phenotype",
           Phenotype_included) {

    ####Check suggested packages####
    {
      if(!requireNamespace("raster", quietly = FALSE)) stop(
        paste0("raster CRAN package is required to execute the function. Please install using the following code: ",
               expression(install.packages("raster")))
      )
      if(!requireNamespace("glcm", quietly = FALSE)) stop(
        paste0("glcm CRAN package is required to execute the function. Please install using the following code: ",
               expression(install.packages("glcm")))
      )
      if(!requireNamespace("tabularaster", quietly = FALSE)) stop(
        paste0("tabularaster CRAN package is required to execute the function. Please install using the following code: ",
               expression(install.packages("tabularaster")))
      )
    }
    
    #####ADAPT THE FUNCTION TO WORK WITH ANY DATA SOURCE####
    #We will adapt the function to work with any source of tiled images we will change the name of the Tiled_image variable name to `Phenotype` and proceed as usual
    #Add checks to test if Tiled_images has been adequately created
    if(!is.list(Tiled_images)) stop("Tiled_images should be a list created with the Image_tiling_processing_function or Tiled_Image_Clustering_function")
    
    List_content_tibble <- unique(purrr::map_lgl(names(Tiled_images), ~tibble::is_tibble(Tiled_images[[.]])))
    if(length(List_content_tibble) > 1) stop("Tiled_images should be a list created with the Image_tiling_processing_function or Tiled_Image_Clustering_function")
    
    #If is tibble then they are derived from Tiled_Image_Clustering_function. Check specific variables
    if(List_content_tibble){
      #If variable name not found stop and print an error
      if(!all(purrr::map_lgl(Tiled_images, ~(Variable_name %in% names(.))))) stop(paste0(Variable_name, " not found in the names of one or more of Tiled_images elements"))
      
      #If no target found in any image then stop the computation
      if(!Phenotype_included %in% unique(unlist(purrr::map(Tiled_images, ~.[[Variable_name]])))) stop(paste0(Phenotype_included,
                                                                                                             " not found in any of the Tiled_images provided"))
      
      #Else proceed by generating a list with two elements one containing all the tiles, and the second containing information from the valid tiles
      Interim <- 
        map(Tiled_images, function(Image){
          #Select desired variables and change the name to `Phenotype`
          Image <- Image %>% dplyr::select(Subject_Names, tile_X_centroid, tile_Y_centroid, tile_id, tile_xmin, tile_xmax, tile_ymin, tile_ymax, dplyr::all_of(Variable_name))
          names(Image)[9] <- "Phenotype"
          
          #Now build all the possible tiles
          Full_tile_tibble <- tidyr::expand_grid(tile_xmin = unique(Image$tile_xmin), tile_ymin = unique(Image$tile_ymin))
          #Compute the tile width and height
          Tile_width <- unique(Image$tile_xmax - Image$tile_xmin)
          Tile_height <- unique(Image$tile_ymax - Image$tile_ymin)
          
          Full_tile_tibble <- Full_tile_tibble %>% dplyr::mutate(tile_xmax = tile_xmin + Tile_width,
                                                                 tile_ymax = tile_ymin + Tile_height)
          #Join with existing data
          Full_tile_tibble <- 
            dplyr::left_join(Full_tile_tibble, Image, dplyr::join_by(tile_xmin, tile_xmax, tile_ymin, tile_ymax))  
          
          Full_tile_tibble$tile_id[is.na(Full_tile_tibble$tile_id)] <- "ZGenerated_tile"
          Full_tile_tibble$tile_X_centroid[is.na(Full_tile_tibble$tile_X_centroid)] <- Full_tile_tibble$tile_xmin[is.na(Full_tile_tibble$tile_X_centroid)] + (Tile_width/2)
          Full_tile_tibble$tile_Y_centroid[is.na(Full_tile_tibble$tile_Y_centroid)] <- Full_tile_tibble$tile_ymin[is.na(Full_tile_tibble$tile_Y_centroid)] + (Tile_height/2)
          Full_tile_tibble <- Full_tile_tibble %>% dplyr::select(tile_X_centroid, tile_Y_centroid, tile_id, tile_xmin, tile_xmax, tile_ymin, tile_ymax)
          
          #Return a list with the results
          return(list(Full_tile_tibble,
                      Image))
        })
      names(Interim) <- names(Tiled_images)
      
      Tiled_images <- Interim
    }
    
    #If the content is a list check that there name of the Variable_name is `Phenotype`
    if(!List_content_tibble){
      Interim <- 
        map(Tiled_images, function(Image){
          Tiles <- Image[[1]]
          Cells <- Image[[2]] %>% dplyr::rename("Phenotype" = dplyr::all_of(Variable_name))
          
          return(list(Tiles,
                      Cells))
        })
      names(Interim) <- names(Tiled_images)
      Tiled_images <- Interim
    }

    #Check that the phenotypes included are present in the data
    if(!all(Phenotype_included %in% unique(unlist(purrr::map(Tiled_images, function(df) df[[2]]$Phenotype))))) {
      stop(paste0("Phenotypes included must be any of: ", stringr::str_c(unique(unlist(purrr::map(Tiled_images, function(df) df[[2]]$Phenotype))), collapse = ", ")))
    }
    
    ####COMPUTE THE ACTUAL METRICS####
    #First import our data
    Tiled_images <- Tiled_images
    #Now we need to calculate the number of cells of interest by tile and transform the tile pattern into a pseudo image, then we calculate the glcm and the texture metrics
    Texture_metric_results <-
      purrr::map_dfr(1:length(Tiled_images), function(Image) {
        #Calculate the cell count by tile
        Cell_count_by_tile <- Tiled_images[[Image]][[2]] %>% dplyr::filter(Phenotype == Phenotype_included) %>%
          group_by(tile_id) %>% dplyr::count(Phenotype)

        #Bind the tile matrix with each cell count
        Result <- dplyr::left_join(Tiled_images[[Image]][[1]], Cell_count_by_tile, by = "tile_id") %>% dplyr::select(-Phenotype)
        Result[is.na(Result)] <- 0

        #Start generating the pseudo image
        x_length <- length(unique(Result$tile_xmin))#Calculte the X pixel length
        y_length <- length(unique(Result$tile_ymin))#Calculate the Y pixel length
        Interim <- Result %>% dplyr::arrange(tile_xmin, desc(tile_ymin))#arrange the tibble in an adequate format
        
        #Build the pseudo image according to the cell number by the cell count
        Pseudo_image <- matrix(Interim[[8]], nrow = y_length, ncol = x_length)#Calculate the matrix
        colnames(Pseudo_image) <- unique(Interim$tile_xmin)
        rownames(Pseudo_image) <- unique(Interim$tile_ymin)

        #If our matrix contain diverse pixel types then execute the functions
        if(length(unique(as.vector(Pseudo_image))) > 1){
          #Now transform our matrix into a raster object and calculate the texture features based on the gray-level-coocurrence-matrix
          Raster_values <- glcm::glcm(raster::raster(Pseudo_image), n_grey = length(unique(as.vector(Pseudo_image))), window = c(3, 3), shift = c(1, 1), statistics =
                                        c("mean", "variance", "homogeneity", "contrast", "dissimilarity", "entropy",
                                          "second_moment", "correlation"), na_opt = "center", na_val = NA, scale_factor = 1, asinteger = FALSE)

          #We average the results of the different texture metrics
          Averaged_results <-purrr::map_dbl(names(Raster_values), function(characteristics) {
            Interim <- tabularaster::as_tibble(raster::subset(Raster_values, characteristics, drop = F), cell = F)
            mean(Interim[[1]], na.rm = T)#Calculate mean values for all texture metrics
          })
          names(Averaged_results) <-stringr::str_c("Mean_", names(Raster_values))
        }

        #If only a unique value present in the sample print message and return a NA tibble
        else{
          message(paste0(names(Tiled_images)[Image], ": the GLCM based texture features could not be calculated for the selected cell type"))
          Averaged_results <- tibble(Mean_glcm_mean = NA,
                                     Mean_glcm_variance = NA,
                                     Mean_glcm_homogeneity = NA,
                                     Mean_glcm_contrast = NA,
                                     Mean_glcm_dissimilarity = NA,
                                     Mean_glcm_entropy = NA,
                                     Mean_glcm_second_moment = NA,
                                     Mean_glcm_correlation = NA)
        }

        return(Averaged_results)
      },
      .progress = list(clear = F,
                       name = "Calculating texture features",
                       show_after = 2,
                       type = "iterator")
      )


    #Bind the results to the Image name and return the final result
    return(dplyr::bind_cols(tibble(Subject_names = names(Tiled_images)),
                            Texture_metric_results)
           )

  }

