#' Calculates the average texture features of a single cell type.
#'
#' The function calculates the average texture features of a single cell type for every image. Image must have been previously tiled using [Image_tiling_processing_function()].
#' The tiled image is first transformed into a raster like object, then the grey level co-ocurrence matrix is calculated and texture features are obtained and averaged.
#'
#' @param Tiled_images A list containing tiled images obtained using [Image_tiling_processing_function()].
#' @param Phenotype_included A character value indicating which cell phenotype to analyze.
#'
#' @details
#' Metrics are calculated using the glcm R package.
#'
#' @returns A tibble containing averaged texture features by image.
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
#'     Phenotype_included = "CD8_GZMBneg"
#')
#' }
#'
#' @export

Texture_features_calculator <-
  function(Tiled_images = NULL,
           Phenotype_included = NULL) {

    #Check suggested packages
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

    #Check that the phenotypes included are present in the data
    if(!all(Phenotype_included %in% unique(unlist(purrr::map(Tiled_images, function(df) df[[2]]$Phenotype))))) {
      stop(paste0("Phenotypes included must be any of: ", stringr::str_c(unique(unlist(purrr::map(Tiled_images, function(df) df[[2]]$Phenotype))), collapse = ", ")))
    }

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
        x_length <- length(unique(Result$tile_X_centroid))#Calculte the X pixel length
        y_length <- length(unique(Result$tile_Y_centroid))#Calculate the Y pixel length
        Interim <- Result %>% dplyr::arrange(tile_X_centroid, desc(tile_Y_centroid))#arrange the tibble in an adequate format

        #Build the pseudo image according to the cell number by the cell count
        Pseudo_image <- matrix(Interim[[8]], nrow = y_length, ncol = x_length)#Calculate the matrix
        colnames(Pseudo_image) <- unique(Interim$tile_X_centroid)
        rownames(Pseudo_image) <- unique(Interim$tile_Y_centroid)

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
