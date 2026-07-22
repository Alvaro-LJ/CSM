#' Identifies noisy cells 
#'
#' Intended for internal use only
#'
#' @details
#' Used in [Clustering_Phenotyper()], [UTAG_Neighborhood_identifier()]
#'
#' @returns A list containing four items: NOISE_VECTOR containing logical values indicating if cells are noise or not, DATA_NOISE a tibble containing noisy cells,
#' DATA a tibble containing non-noisy cells, MARKERS the markers of non-noisy cells.
#' 
#' @keywords Internal

CSM_Denoising_function <- 
  function(
    Original_data = NULL,
    
    Denoising_strategy = NULL, 
    
    Percentile = NULL, 
    N_Standard_Deviations = NULL, 
    Selected_threshold = NULL, 
    
    Perform_Dimension_reduction = NULL,
    DATA_Reduction = NULL,
    Min_cell_no = NULL, 
    Distance_radius = NULL
    
    ){
    
    ####ADD UNIQUE CELL ID####
    DATA <- Original_data %>% dplyr::mutate(Unique_ID = 1:nrow(Original_data))
    DATA <- DATA[c(ncol(DATA), 1:(ncol(DATA)-1))]
    
    ####QUANTILE FILTERING####
    if("Quantile" %in% Denoising_strategy) {
      print("Filtering by quantiles")
      #Check arguments
      if(Percentile < 0.01 | Percentile > 0.99) stop("Percentile must be between 0.01 and 0.99")
      
      FILTER_Quantile <- purrr::map_dfc(DATA[-(1:5)], function(x){
        x <= quantile(x, Percentile)
      })
      
      #Generate the variable to specify if the cell is noise or not
      NOISE_Quantile <- unlist(apply(FILTER_Quantile, MARGIN = 1, function(x) sum(x) == ncol(FILTER_Quantile)))
    }
    
    ####STANDARD DEVIATION FILTERING####
    if("Standard_Deviation" %in% Denoising_strategy){
      print("Filtering by SD")
      #Check arguments
      if(!is.numeric(N_Standard_Deviations)) stop("N_Standard_Deviations must be a numeric value")
      
      FILTER_SD <-purrr::map_dfc(DATA[-(1:5)], function(x){
        x <= (mean(x) - (N_Standard_Deviations*sd(x)))
      })
      
      #Generate the variable to specify if the cell is noise or not
      NOISE_SD <- unlist(apply(FILTER_SD, MARGIN = 1, function(x) sum(x) == ncol(FILTER_SD)))
    }
    
    ####THRESHOLD FILTERING####
    if("Threshold" %in% Denoising_strategy){
      print("Filtering by user selected threshold")
      #Check arguments
      if(!is.numeric(Selected_threshold)) stop("Selected_threshold must be a numeric value")
      
      FILTER_Threshold <- purrr::map_dfc(DATA[-(1:5)], function(x){
        x <= Selected_threshold
      })
      
      #Generate the variable to specify if the cell is noise or not
      NOISE_Threshold <- unlist(apply(FILTER_Threshold, MARGIN = 1, function(x) sum(x) == ncol(FILTER_Threshold)))
    }
    
    ####OTSU FILTERING####
    if("Otsu" %in% Denoising_strategy){
      print("Filtering by Otsu threshold")
      FILTER_OTSU <- purrr::map2_df(.x = DATA[-c(1:5)],
                                    .y =purrr::map_dbl(DATA[-c(1:5)], function(z){
                                      EBImage::otsu(array(z, dim = c(1, length(z))), range = c(min(z), max(z)), levels = length(unique(z)))
                                    }),
                                    function(.x, .y) .x <= .y)
      
      #Generate the variable to specify if the cell is noise or not
      NOISE_OTSU <- unlist(apply(FILTER_OTSU, MARGIN = 1, function(x) sum(x) == ncol(FILTER_OTSU)))
    }
    
    ####DIMENSION REDUCTION FOLLOWED BY DB-SCAN####
    if("DimRed_DBscan" %in% Denoising_strategy){
      #Requires previous dimension reduction
      if(!Perform_Dimension_reduction) stop("DBscan clustering requires Dimension reduction. Please set Perform_Dimension_reduction to TRUE")
      #Check other arguments
      if(!all(is.numeric(Min_cell_no), Min_cell_no%%1 == 0, Min_cell_no > 1)) stop("Min_cell_no must be an integer value > 1")
      if(!all(is.numeric(Distance_radius), Distance_radius > 0)) stop("Distance_radius must be a numeric value > 0")
      
      print("Filtering by DBSCAN clustering")
      
      #Proceed with algorithm
      DB_results <- dbscan::dbscan(DATA_Reduction[c("DIMENSION_1", "DIMENSION_2")], eps = Distance_radius, minPts = Min_cell_no, borderPoints = FALSE)
      
      #whole plot for small samples
      if(length(DB_results$cluster) <= 100000){
        plot(
          tibble(Dim_1 = DATA_Reduction[["DIMENSION_1"]], Dim_2 = DATA_Reduction[["DIMENSION_2"]], Cluster = DB_results$cluster) %>%
            dplyr::mutate(Cluster = case_when(Cluster == 0 ~ "Noise",
                                              TRUE ~ "Approved")) %>%
            ggplot(aes(x = Dim_1, y = Dim_2, color = Cluster)) + geom_point(size = 0.8) +
            scale_color_manual("", values = c("black", "grey"))+
            cowplot::theme_cowplot()+
            theme(panel.grid = element_blank())
        )
      }
      
      #Subsample plot for large dataset
      if(length(DB_results$cluster) > 100000){
        message(">100K observations to generate plots. A random subset containing 10% of the dataset will be selected for Dimension reduction plots")
        plot(
          tibble(Dim_1 = DATA_Reduction[["DIMENSION_1"]], Dim_2 = DATA_Reduction[["DIMENSION_2"]], Cluster = DB_results$cluster) %>%
            dplyr::mutate(Cluster = case_when(Cluster == 0 ~ "Noise",
                                              TRUE ~ "Approved")) %>%
            dplyr::slice_sample(n = 100000) %>%
            ggplot(aes(x = Dim_1, y = Dim_2, color = Cluster)) + geom_point(size = 1.5) +
            scale_color_manual("", values = c("black", "grey"))+
            cowplot::theme_cowplot()+
            theme(panel.grid = element_blank())
        )
      }
      
      DB_OK <- menu(choices = c("Proceed", "Abort"), title = "Are the results of the filtering OK?")
      if(DB_OK == 2) stop("Please refine Distance_radius and Min_cell_no parameters and retry")
      
      #Generate the variable to specify if the cell is noise or not
      NOISE_DBSCAN <- DB_results$cluster == 0
    }
    
    ####MERGE ALL IN A SINGLE VECTOR#####
    NOISE_VECTOR_LIST <- list()
    if("Quantile" %in% Denoising_strategy) NOISE_VECTOR_LIST$NOISE_Quantile <- NOISE_Quantile
    if("Standard_Deviation" %in% Denoising_strategy) NOISE_VECTOR_LIST$NOISE_SD <- NOISE_SD
    if("Threshold" %in% Denoising_strategy) NOISE_VECTOR_LIST$NOISE_Threshold <- NOISE_Threshold
    if("Otsu" %in% Denoising_strategy) NOISE_VECTOR_LIST$NOISE_OTSU <- NOISE_OTSU
    if("DimRed_DBscan" %in% Denoising_strategy) NOISE_VECTOR_LIST$NOISE_DBSCAN <- NOISE_DBSCAN
    
    #If single method then is the first element
    if(length(NOISE_VECTOR_LIST) == 1) NOISE_VECTOR <- NOISE_VECTOR_LIST[[1]]
    #If multiple methods then reduce
    else{
      NOISE_VECTOR <- as.logical(purrr::reduce(NOISE_VECTOR_LIST, function(Method_1, Method_2) Method_1 | Method_2))
    }

    ####GENERATE THE FINAL OBJECTS TO BE RETURNED####
    #Generate two tibbles, one with noisy cells and other (DATA) with the actual cells
    DATA_NOISE <- DATA %>% dplyr::filter(NOISE_VECTOR) %>% dplyr::mutate(Cluster = 1)
    DATA <- DATA %>% dplyr::filter(!NOISE_VECTOR)
    MARKERS <- DATA %>% dplyr::select(-Unique_ID, -Cell_no, -X, -Y, -Subject_Names)
    
    return(
      list(NOISE_VECTOR = NOISE_VECTOR,
           DATA_NOISE = DATA_NOISE,
           DATA = DATA,
           MARKERS = MARKERS)
    )
    
  }
