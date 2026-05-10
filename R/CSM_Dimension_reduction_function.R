CSM_Dimension_reduction_function <- 
  function(
    Original_data = NULL,
    
    Dimension_reduction_strategy = NULL,
    Dimension_reduction_prop = NULL
    ){
    
    ####PRINCIPAL COMPONENT ANALYSIS (PCA) ####
    if(Dimension_reduction_strategy == "PCA"){
      if(Dimension_reduction_prop != 1) stop("PCA must be performed using Dimension_reduction_prop = 1")
      print("Generating PCA projections")
      #Scale and turn into matrix
      DATA_matrix <- Original_data %>% dplyr::select(-c(1:4)) %>% scale() %>% as.matrix()
      Result_PCA <- svd::propack.svd(DATA_matrix, neig = 2)$u
      DATA_Reduction <- tibble(Cell_no = Original_data$Cell_no, DIMENSION_1 = unlist(Result_PCA[,1]), DIMENSION_2 = unlist(Result_PCA[,2]))
    }
    
    ####t-DISTRIBUTED STOCHASTIC NEIGHBORHOOD EMBEDDING (TSNE)####
    if(Dimension_reduction_strategy == "TSNE"){
      if(Dimension_reduction_prop == 1) {
        print("Generating TSNE projections")
        if(nrow(Original_data) > 50000) print("Warning! Data set contains more than 50K observations. tSNE embedding can take a long time")
        #scale and turn into matrix
        DATA_matrix <- Original_data %>% dplyr::select(-c(1:4)) %>% scale() %>% as.matrix()
        Result_TSNE <- snifter::fitsne(DATA_matrix,
                                       simplified = TRUE,
                                       n_components = 2L,
                                       n_jobs = 1L,
                                       perplexity = 30,
                                       n_iter = 500L,
                                       initialization = "pca",
                                       pca = FALSE,
                                       neighbors = "auto",
                                       negative_gradient_method = "fft",
                                       learning_rate = "auto",
                                       early_exaggeration = 12,
                                       early_exaggeration_iter = 250L,
                                       exaggeration = NULL,
                                       dof = 1,
                                       theta = 0.5,
                                       n_interpolation_points = 3L,
                                       min_num_intervals = 50L,
                                       ints_in_interval = 1,
                                       metric = "euclidean",
                                       metric_params = NULL,
                                       initial_momentum = 0.5,
                                       final_momentum = 0.8,
                                       max_grad_norm = NULL,
                                       random_state = NULL,
                                       verbose = FALSE)
        DATA_Reduction <- dplyr::bind_cols(Original_data["Cell_no"], DIMENSION_1 = unlist(Result_TSNE[,1]), DIMENSION_2 = unlist(Result_TSNE[,2]))
      }
      
      if(Dimension_reduction_prop != 1) {
        print("Generating TSNE projections")
        DATA_matrix <- Original_data %>% dplyr::group_by(Subject_Names) %>% dplyr::slice_sample(prop = Dimension_reduction_prop) %>% dplyr::ungroup() %>%
          dplyr::select(-c(1:4)) %>% scale() %>% as.matrix()
        if(nrow(DATA_matrix) > 50000) print("Warning! Data set contains more than 50K observations. tSNE embedding can take a long time")
        #scale and turn into matrix
        Result_TSNE <- snifter::fitsne(DATA_matrix,
                                       simplified = FALSE,
                                       n_components = 2L,
                                       n_jobs = 1L,
                                       perplexity = 30,
                                       n_iter = 500L,
                                       initialization = "pca",
                                       pca = FALSE,
                                       neighbors = "auto",
                                       negative_gradient_method = "fft",
                                       learning_rate = "auto",
                                       early_exaggeration = 12,
                                       early_exaggeration_iter = 250L,
                                       exaggeration = NULL,
                                       dof = 1,
                                       theta = 0.5,
                                       n_interpolation_points = 3L,
                                       min_num_intervals = 50L,
                                       ints_in_interval = 1,
                                       metric = "euclidean",
                                       metric_params = NULL,
                                       initial_momentum = 0.5,
                                       final_momentum = 0.8,
                                       max_grad_norm = NULL,
                                       random_state = NULL,
                                       verbose = FALSE)
        Coords <- snifter::project(Result_TSNE,
                                   new = Original_data  %>% dplyr::select(-c(1:4)) %>% scale() %>% as.matrix(),
                                   old = DATA_matrix)
        DATA_Reduction <-dplyr::bind_cols(Original_data["Cell_no"], DIMENSION_1 = unlist(Coords[,1]), DIMENSION_2 = unlist(Coords[,2]))
      }
    }
    
    ####UNIFORM MANIFOLD APPROXIMATION AND PROJECTION (UMAP)####
    if(Dimension_reduction_strategy == "UMAP"){
      if(Dimension_reduction_prop == 1) {
        print("Generating UMAP projections")
        if(nrow(Original_data) > 50000) print("Warning! Data set contains more than 50K observations. UMAP embedding can take some time")
        #scale and turn into matrix
        DATA_matrix <- Original_data %>% dplyr::select(-c(1:4)) %>% scale() %>% as.matrix()
        Result_UMAP <- uwot::tumap(DATA_matrix, n_components = 2L)
        DATA_Reduction <- dplyr::bind_cols(Original_data["Cell_no"], DIMENSION_1 = unlist(Result_UMAP[,1]), DIMENSION_2 = unlist(Result_UMAP[,2]))
      }
      
      if(Dimension_reduction_prop != 1) {
        print("Generating UMAP projections")
        DATA_matrix <- Original_data %>% dplyr::group_by(Subject_Names) %>% dplyr::slice_sample(prop = Dimension_reduction_prop) %>% dplyr::ungroup() %>%
          dplyr::select(-c(1:4)) %>% scale() %>% as.matrix()
        if(nrow(DATA_matrix) > 50000) print("Warning! Data set contains more than 50K observations. UMAP embedding can take some time")
        #scale and turn into matrix
        Result_UMAP <- uwot::tumap(DATA_matrix, n_components = 2L, ret_model = TRUE)
        Coords <- uwot::umap_transform(X = Original_data  %>% dplyr::select(-c(1:4)) %>% scale() %>% as.matrix(),
                                       model = Result_UMAP)
        DATA_Reduction <- dplyr::bind_cols(Original_data["Cell_no"], DIMENSION_1 = unlist(Coords[,1]), DIMENSION_2 = unlist(Coords[,2]))
      }
    }
    
    #### RETURN THE REDUCTION DATA ####
    return(DATA_Reduction)
  }
