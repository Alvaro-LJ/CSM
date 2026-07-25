#' Computes cell features based on image and cell segmentation mask
#'
#' Intended for internal use only
#'
#' @param Image An EBImage object containing the Image.
#' @param Mask An EBImage object containing the cell segmentation mask. Each element in the mask must be an integer value.
#' 
#' @param Ordered_Channels Character vector specifying Image channel names in their exact order.
#' 
#' @param Compute_basic_features A logical value indicating if basic features (mean, sd, quantiles) should be computed.
#' @param quantiles_to_calculate Numeric vector specifying the expression quantiles calculated for every cell and channel.
#' @param Compute_texture_features A logical value indicating if texture features should be computed for all the markers measured.
#' @param Texture_pixel_distance If Compute_texture_features is TRUE, de pixel distance taken into account to compute the GLCM and texture features.
#' 
#' @param SubCellular_compartment_list A named list containing the cell masks for subcellular expression computation.
#'
#' @details
#' Used in [Cell_feature_extractor()]
#'
#'
#' @returns A tibble containing cell data.
#' @keywords Internal


Cell_image_mask_profiler <- function(
    Image = NULL,
    Mask = NULL,
    
    
    Ordered_Channels = NULL,
    
    
    Compute_basic_features = TRUE,
    quantiles_to_calculate = c(0.25, 0.5, 0.75), 
    Compute_texture_features = FALSE,
    Texture_pixel_distance = 1,
    
    SubCellular_compartment_list = NULL
){
  
  #Compute mask-only based tibbles for all cells
  Position_tibble <- tibble::as_tibble(EBImage::computeFeatures.moment(Mask))
  Shape_tibble <- tibble::as_tibble(EBImage::computeFeatures.shape(Mask))
  Cell_id <- tibble::tibble(Cell_id = stringr::str_c("Cell_", rownames(Position_tibble)))
  
  #Generate the preliminary final tibble
  Final_tibble <- dplyr::bind_cols(Cell_id, Position_tibble, Shape_tibble)
  
  #If basic features are required, compute the basic features for the whole cell
  if(Compute_basic_features){
    #Compute channel expression values for every channel in the Image
    Expression_results <- 
      purrr::map_dfc(seq_along(1:length(Ordered_Channels)), 
                     function(Image_index){
                       #Get channel name
                       Channel_name <- Ordered_Channels[Image_index]
                       
                       
                       #Compute the result tibble and change names
                       Result_tibble <- tibble::as_tibble(EBImage::computeFeatures.basic(x = Mask,
                                                                                         ref = Image[,,Image_index],
                                                                                         basic.quantiles = quantiles_to_calculate))
                       names(Result_tibble) <- stringr::str_c(Channel_name, 
                                                              "_", 
                                                              substr(names(Result_tibble), start = 3, stop = nchar(names(Result_tibble))))
                       
                       #Return the result
                       return(Result_tibble)
                     })
    
    Final_tibble <- dplyr::bind_cols(Final_tibble, Expression_results)
  }
  
  #If texture features are required, compute texture features
  if(Compute_texture_features){
    Texture_results <- 
      purrr::map_dfc(seq_along(1:length(Ordered_Channels)), 
                     function(Image_index){
                       #Get channel name
                       Channel_name <- Ordered_Channels[Image_index]
                       
                       
                       #Compute the result tibble and change names
                       Result_tibble <- tibble::as_tibble(EBImage::computeFeatures.haralick(x = Mask,
                                                                                            ref = Image[,,Image_index],
                                                                                            haralick.scales = Texture_pixel_distance))
                       names(Result_tibble) <- stringr::str_c(Channel_name, 
                                                              "_", 
                                                              substr(names(Result_tibble), start = 1, stop = nchar(names(Result_tibble))))
                       
                       #Return the result
                       return(Result_tibble)
                     })
    Final_tibble <- dplyr::bind_cols(Final_tibble, Texture_results)
    
  }
  
  #If cell subcellular mask provided then proceed and compute subcellular elements
  if(!is.null(SubCellular_compartment_list)){
    Subcellular_result_list <- 
      purrr::map(seq_along(1:length(SubCellular_compartment_list)),
                 function(Subcellular_region_index){
                   #Get the subcellular name
                   Subcellular_name <- names(SubCellular_compartment_list)[Subcellular_region_index]
                   Subcellular_mask <- SubCellular_compartment_list[[Subcellular_region_index]]

                   
                   #Obtain the basic features of the subcellular region
                   Position_tibble <- tibble::as_tibble(EBImage::computeFeatures.moment(Subcellular_mask))
                   Shape_tibble <- tibble::as_tibble(EBImage::computeFeatures.shape(Subcellular_mask))
                   Cell_id <- tibble::tibble(Cell_id = stringr::str_c("Cell_", rownames(Position_tibble)))
                   
                   #Generate the basic tibble
                   Subcellular_tibble <- dplyr::bind_cols(Cell_id, Position_tibble, Shape_tibble)
                   
                   #Generate additional data if required
                   #If basic features are required, compute the basic features for the whole cell
                   if(Compute_basic_features){
                     #Compute channel expression values for every channel in the Image
                     Expression_results <- 
                       purrr::map_dfc(seq_along(1:length(Ordered_Channels)), 
                                      function(Image_index){
                                        #Get channel name
                                        Channel_name <- Ordered_Channels[Image_index]
                                        
                                        
                                        #Compute the result tibble and change names
                                        Result_tibble <- tibble::as_tibble(EBImage::computeFeatures.basic(x = Subcellular_mask,
                                                                                                          ref = Image[,,Image_index],
                                                                                                          basic.quantiles = quantiles_to_calculate))
                                        names(Result_tibble) <- stringr::str_c(Channel_name, 
                                                                               "_", 
                                                                               substr(names(Result_tibble), start = 3, stop = nchar(names(Result_tibble))))
                                        
                                        #Return the result
                                        return(Result_tibble)
                                      })
                     
                     Subcellular_tibble <- dplyr::bind_cols(Subcellular_tibble, Expression_results)
                   }
                   
                   #Change the names adding the subcellular region
                   names(Subcellular_tibble)[-1] <- stringr::str_c(Subcellular_name, "_", names(Subcellular_tibble)[-1])
                   
                   #return the final tibble
                   return(Subcellular_tibble)
                 }
      )
    
    #Reduce the list by performing full join
    
    Subcellular_tibble <- purrr::reduce(Subcellular_result_list, function(x, y) dplyr::full_join(x, y, by = "Cell_id"))
    
    Final_tibble <- dplyr::left_join(Final_tibble, Subcellular_tibble, by = "Cell_id")
  }
  

  #Return the final product
  return(Final_tibble)
}

