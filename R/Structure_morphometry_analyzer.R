#' Generates by-sample summary metrics for structure morphometry analysis
#'
#' Given a list containing structure morphometry metrics (generated using the [Structure_morphometry_calculator()]) the function will compute by-sample summary metrics.
#' 
#' 
#' @param DATA_Morphometry A list containing structure morphometry data (generated using the [Structure_morphometry_calculator()]).
#' 
#' 
#' @returns A list containing a tibble with the summary results.
#'
#' @seealso [Image_tiling_processing_function()], [Tiled_Image_Clustering_function()], [Structure_morphometry_calculator()], [Structure_morphometry_graph_maker()]
#'
#' @examples
#' \dontrun{
#' #Divide cells into discrete tiles
#'Tiled_Images <-
#' Image_tiling_processing_function(
#'  N_cores = 1,
#'  DATA = CSM_Phenotypecell_test,
#'  Tile_width = 25,
#'  Tile_height = 25,
#'  Variables_to_keep = "Phenotype")
#'
#' #Try to identify morphological structures for OTHER cells according to tile distribution
#'Structures <- 
#' Structure_morphometry_calculator(N_cores = 1,
#'                                 
#'                                 Tiled_Images = Tiled_Images,
#'                                 Variable_name = "Phenotype",
#'                                 Feature_to_analyze = "OTHER",
#'                                 
#'                                 Min_tile_number = 1, 
#'                                 Smooth_amount = 0.01, 
#'                                 Tolerance_value = 1, 
#'                                 Neighborhood_distance = 1,
#'                                 Minimum_cell_no_per_tile = 2)
#'
#' #Analyze the results
#' Structure_morphometry_analyzer(DATA = Structures)
#'
#' #Graph the results
#' Structure_morphometry_graph_maker(Tiled_Images = Tiled_Images,
#'                                   Variable_name = "Phenotype",
#'                                   DATA_Morphometry = Structures,
#'                                   Color_by = "Structure_ID",
#'                                   Image_name = "ABCM22001_B09_MiniCrop.tif")
#'}
#'
#'
#' @export

Structure_morphometry_analyzer <-
  function(
    DATA_Morphometry
  ){
    
    ####Check arguments provided####
    #Check that DATA_Morphometry is a list
    if(!is.list(DATA_Morphometry)) stop("DATA_Morphometry must be generated using the Structure_morphometry_calculator function")
    
    #Check content of DATA_Morphometry is adecuate
    Conflictive_elements <- 
      !purrr::map_lgl(DATA_Morphometry, function(Image){
        identical(names(Image), c("Morphology_tibble", "Structure_mask"))
      })
    
    if(any(Conflictive_elements)) stop("DATA_Morphometry must be generated using the Structure_morphometry_calculator function")
    
    #Remove images without adequate analysis
    
    ####Compute the actual summary metrics##
    Metrics_to_analyze <- c("s.area", "s.perimeter", "s.radius.mean", "s.radius.sd", "s.radius.min", "s.radius.max",
                            "m.majoraxis", "m.eccentricity", "m.theta_AVERAGE")
    
    #Iterate over relevant metrics (the final result will be a list with a single element per metric)
    FINAL_LIST <-
      purrr::map(Metrics_to_analyze, function(Metric_name){
        
        #Then iterate across samples
        purrr::map_dfr(1:length(DATA_Morphometry), function(Image_index){
          Subject_Names <- names(DATA_Morphometry)[Image_index]
          N_structures <- nrow(DATA_Morphometry[[Image_index]][[1]])
          
          Variable <- DATA_Morphometry[[Image_index]][[1]][[Metric_name]]
          Average_variable <- mean(Variable)
          P25_variable <- quantile(Variable, 0.25)
          Median_variable <- quantile(Variable, 0.5)
          P75_variable <- quantile(Variable, 0.75)
          Minimum_variable <- min(Variable)
          Maximum_variable <- max(Variable)
          SD_variable <- sd(Variable)
          
          return(tibble(Subject_Names = Subject_Names,
                        N_structures = N_structures,
                        Min = Minimum_variable,
                        P25 = P25_variable,
                        Median = Median_variable,
                        P75 = P75_variable,
                        Max = Maximum_variable,
                        Average = Average_variable,
                        SD = SD_variable))
        })
      })
    
    names(FINAL_LIST) <- Metrics_to_analyze
    
    #Return the final product
    return(FINAL_LIST)
  }

