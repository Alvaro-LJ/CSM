CSM_Clustering_function <- 
  function(
    Original_data,
    MARKERS,
    
    Strategy,
    Force_N_Clusters,
    
    Max_N_clusters_Consensus = NULL, 
    Consensus_reps = NULL, 
    Consensus_p_Items = NULL, 
    Consensus_Cluster_Alg = NULL,
    Consensus_Distance = NULL, 
    Consensus_Name = NULL, 
    
    Max_SOM_clusters = NULL, 
    
    Graph_type = NULL,
    Graph_Distance_method = NULL,
    Nearest_neighbors_for_graph = NULL, 
    Graph_Method = NULL,
    Graph_Resolution = NULL, 
    N_steps = NULL, 
    
    N_K_centroids = NULL, 
    Max_N_clusters_Meta = NULL, 
    Consensus_reps_Meta = NULL, 
    Consensus_p_Items_Meta = NULL,
    Consensus_Name_Meta = NULL, 
    
    Batch_size = NULL,
    Max_N_clusters_Batch = NULL, 
    N_initiations = NULL, 
    Max_iterations = NULL, 
    
    Quality_metric = NULL, 
    Max_N_clusters_GMM = NULL, 
    Max_iterations_km = NULL, 
    Max_iterations_em = NULL,
    GMM_Distance = NULL, 
    
    Samples_CLARA = NULL, 
    Sample_per_CLARA = NULL, 
    Max_N_clusters_CLARA = NULL, 
    Distance_CLARA = NULL,
    N_cores = NULL 
  ){
    
    ####CONSENSUS CLUSTERING####
    if(Strategy == "Consensus_Clustering"){
      #Perform consensus clustering
      Clustering_result <- try(ConsensusClusterPlus::ConsensusClusterPlus(t(as.matrix((MARKERS %>% scale()))),
                                                                          maxK = Max_N_clusters_Consensus,
                                                                          reps = Consensus_reps,
                                                                          pItem = Consensus_p_Items,
                                                                          pFeature = 1,
                                                                          title = Consensus_Name,
                                                                          clusterAlg = Consensus_Cluster_Alg,
                                                                          distance = Consensus_Distance,
                                                                          plot = "png",
                                                                          verbose = T)
      )
      #Test if consensus clustering returned an error
      if(berryFunctions::is.error(Clustering_result)) {
        stop("Data is too large for Consensus Clustering. Please try another strategy")
      }
      else {
        #Make the user decide the number of neighborhoods according to results
        N_Clusters <- menu(choices = as.character(1:Max_N_clusters_Consensus), title = paste0("Check the results at: ", getwd(), ". Then decide the appropiate number of clusters"))
        DATA_Clusters <- Original_data %>% dplyr::mutate(Cluster = Clustering_result[[as.double(N_Clusters)]][["consensusClass"]])
      }
    }
    
    ####SELF ORGANIZING MAPS####
    else if(Strategy == "SOM"){
      print("Executing Self Organizing Map algorithm")
      #CLUSTER NUMBER STIMATION REQUIRED
      if(!Force_N_Clusters){
        #Transform data into a scaled matrix and perform Self Organizing Map
        Clustering_result <- try(FlowSOM::FlowSOM(MARKERS %>% scale() %>% as.matrix(),
                                                  scale = F,
                                                  colsToUse = c(1:ncol(MARKERS)),
                                                  maxMeta = Max_SOM_clusters,#To find optimal meta clusters
                                                  silent = F,
                                                  seed = 21)
        )
      }
      
      #CLUSTER STIMATION NOT REQUIRED
      if(Force_N_Clusters){
        #Transform data into a scaled matrix and perform Self Organizing Map
        Clustering_result <- try(FlowSOM::FlowSOM(MARKERS %>% scale() %>% as.matrix(),
                                                  scale = F,
                                                  colsToUse = c(1:ncol(MARKERS)),
                                                  nClus = Max_SOM_clusters,#To specify the exact number of clusters
                                                  silent = F,
                                                  seed = 21)
        )
      }
      
      #Test if SOM returned an error
      if(berryFunctions::is.error(Clustering_result)) {
        stop("Data is too large or too small for Self-Organizing Maps. Please try another strategy")
      }
      else{
        #Assign phenotypes to each cell
        DATA_Clusters <- Original_data %>% dplyr::mutate(Cluster = FlowSOM::GetMetaclusters(Clustering_result))
      }
      
    }
    
    ####GRAPH-BASED CLUSTERING####
    else if(Strategy == "Graph_Based"){
      print("Start Graph building process")
      
      #Generate graphs according to user defined preferences
      
      
      #####NEED TO TEST FOR CLUSTERING TILING FUNCTION!!!!!######
      if(Graph_type == "Complete"){
        print("Generating the complete graph")
        #Generate scaled MARKERS
        Scaled_Markers <- MARKERS %>% scale()
        
        #Calculate distance matrix and then calculate the graph
        #We define the number and ID of edges of the graph
        Graph_tibble <- as_tibble(expand.grid.unique(1:nrow(Scaled_Markers), 1:nrow(Scaled_Markers)))
        names(Graph_tibble) <- c("from", "to")
        Graph_tibble <- Graph_tibble %>% dplyr::mutate(ID = stringr::str_c(from, to, sep = "_"))
        
        #We determine the distance between nodes that will be the features of the edges
        DISTANCE_MATRIX <- as_tibble(as.matrix(dist(Scaled_Markers, method = Graph_Distance_method)))
        DISTANCE_MATRIX <- DISTANCE_MATRIX %>% dplyr::mutate(from = as.character(1:nrow(DISTANCE_MATRIX)))
        DISTANCE_MATRIX <- DISTANCE_MATRIX[c(ncol(DISTANCE_MATRIX), 2:(ncol(DISTANCE_MATRIX)-1))] %>% tidyr::pivot_longer(-1, names_to = "to", values_to = "weight") %>%
          dplyr::mutate(ID = stringr::str_c(from, to, sep = "_")) %>% dplyr::select(-from, -to)
        
        #We bind the edges to their features (distance) and we build the graph
        GRAPH_DF <- dplyr::left_join(Graph_tibble, DISTANCE_MATRIX, by = "ID") %>% dplyr::select(-ID)
        GRAPH_DF <- GRAPH_DF %>% dplyr::mutate(weight = 1/weight)
        Node_ID <- tibble(Name = as.character(1:nrow(Scaled_Markers)))
        Graph_for_clustering <- igraph::graph_from_data_frame(GRAPH_DF, directed = F, vertices = Node_ID)
      }

      if(Graph_type == "SNN"){
        print("Generating the SNN graph")
        #Transform data into a nearest neighbor graph
        Graph_for_clustering <- try(bluster::makeSNNGraph(MARKERS %>% scale() %>% as.matrix(),
                                                          k = Nearest_neighbors_for_graph)
        )
      }
      
      #Test if Graph construction process returned an error
      if(berryFunctions::is.error(Graph_for_clustering)) {
        stop("Data is too large to build a graph. Please try another strategy or refine graph generation parameters")
      }
      
      else{
        print("Performing Graph-based clustering")
        #Go for graph clustering
        #Cluster the graph with louvain or leiden clustering
        if(Graph_Method == "Louvain") {
          DATA_Clusters <- Original_data %>% dplyr::mutate(Cluster = igraph::cluster_louvain(Graph_for_clustering,
                                                                                             weights = NULL,
                                                                                             resolution = Graph_Resolution)$membership)
          
        }
        
        else if(Graph_Method == "Leiden") {
          DATA_Clusters <- Original_data %>%dplyr::mutate(Cluster = igraph::cluster_leiden(Graph_for_clustering,
                                                                                           objective_function = "modularity",
                                                                                           weights = NULL,
                                                                                           resolution = Graph_Resolution,
                                                                                           beta = 0.01,
                                                                                           initial_membership = NULL,
                                                                                           n_iterations = 100,
                                                                                           vertex_weights = NULL)$membership)
          
        }
        
        else if(Graph_Method == "Greedy"){
          DATA_Clusters <- Original_data %>% dplyr::mutate(Cluster = igraph::cluster_fast_greedy(Graph_for_clustering)$membership)
        }
        
        else if(Graph_Method == "WalkTrap"){
          DATA_Clusters <- Original_data %>% dplyr::mutate(Cluster = igraph::cluster_walktrap(Graph_for_clustering,
                                                                                              steps = N_steps,
                                                                                              membership = T)$membership)
        }
        
        else if (Graph_Method == "Spinglass") {
          DATA_Clusters <- Original_data %>% dplyr::mutate(Cluster = igraph::cluster_spinglass(Graph_for_clustering,
                                                                                               weights = NULL,
                                                                                               vertex = NULL,
                                                                                               spins = 25,
                                                                                               parupdate = FALSE,
                                                                                               start.temp = 1,
                                                                                               stop.temp = 0.01,
                                                                                               cool.fact = 0.99,
                                                                                               update.rule = c("config", "random", "simple"),
                                                                                               gamma = 1,
                                                                                               implementation = c("orig", "neg"),
                                                                                               gamma.minus = 1)$membership)
          
        }
        
        else if(Graph_Method == "Leading_Eigen"){
          DATA_Clusters <- Original_data %>% dplyr::mutate(Cluster = igraph::cluster_leading_eigen(Graph_for_clustering,
                                                                                                   membership = T)$membership)
        }
        
        else if(Graph_Method == "Edge_Betweenness"){
          DATA_Clusters <- Original_data %>% dplyr::mutate(Cluster = igraph::cluster_edge_betweenness(Graph_for_clustering,
                                                                                                      weights = NULL,
                                                                                                      directed = FALSE,
                                                                                                      edge.betweenness = FALSE,
                                                                                                      merges = FALSE,
                                                                                                      bridges = FALSE,
                                                                                                      modularity = FALSE,
                                                                                                      membership = TRUE)$membership)
          
        }
      }
      
    }
    
    ####K-MEANS META CLUSTERINGS####
    else if(Strategy == "K_Means_Meta_clustering"){
      #First we need to perform K means Clustering
      print("Performing initial K-means algorithm")
      cl <- try(stats::kmeans(MARKERS %>% scale() %>% as.matrix(), #Scale it and turn it into a matrix
                              centers = N_K_centroids, #Number of centroids to be calculated
                              iter.max = 50,
                              nstart = 10))
      
      #Stop function if K means returned an error
      if(berryFunctions::is.error(cl)) {
        stop("Data is too large for K means clustering. Please try another strategy")
      }
      #Proceed if no error was returned
      else{
        #Assign this K means cluster to each observation
        DATA_filter_Markers <- Original_data %>% dplyr::mutate(K_means_Cl = cl$cluster)
        
        #Prepare data for Meta-Clustering
        #Create a tibble with the K means centroids and the format it for Consensus clustering
        K_medoids <- tibble::as_tibble(cl$centers) %>% dplyr::mutate(K_means_Cl = 1:nrow(as_tibble(cl$centers)))
        tK_medoids <- K_medoids %>% dplyr::select(-K_means_Cl) %>% as.matrix() %>% t()
        
        #Perform Consensus clustering with hierarchical clustering
        print("Perorming Consensus Clustering")
        Clustering_result <- try(ConsensusClusterPlus::ConsensusClusterPlus(tK_medoids,
                                                                            maxK = Max_N_clusters_Meta,
                                                                            reps = Consensus_reps_Meta,
                                                                            pItem = Consensus_p_Items_Meta,
                                                                            pFeature = 1,
                                                                            title = Consensus_Name_Meta,
                                                                            distance = "euclidean",
                                                                            clusterAlg = "pam",
                                                                            plot = "png",
                                                                            verbose = T)
        )
        
        #Test if consensus clustering returned an error
        if(berryFunctions::is.error(Clustering_result)) {
          stop("Data is too large for Meta Clustering. Please try another strategy or select a smaller N_K_centroids value")
        }
        else {
          #Make the user decide the number of neighborhoods according to results
          N_Clusters <- menu(choices = as.character(1:Max_N_clusters_Meta), title = paste0("Check the results at: ", getwd(), ". Then decide the appropiate number of clusters"))
          
          #Bind the final Phenotype to the K medoids tibble
          K_medoids <- K_medoids %>% dplyr::mutate(Cluster = Clustering_result[[as.double(N_Clusters)]][["consensusClass"]])
          K_medoids_for_join <- K_medoids %>% dplyr::select(K_means_Cl, Cluster)
          
          #Bind The DATA and the K_meoids to obtain the final matrix
          DATA_Clusters <- dplyr::left_join(DATA_filter_Markers, K_medoids_for_join, by = "K_means_Cl") %>% dplyr::select(-K_means_Cl)
        }
      }
      
    }
    
    ####BATCHED K-MEANS####
    else if(Strategy == "Batch_K_means"){
      #First we calculate a metric to decide the number of total phenotypes
      #Specify the params
      if(Batch_size > nrow(MARKERS)) stop(paste0("Batch_size cannot be larger than ", nrow(MARKERS)))
      
      #Compute scaled markers
      Scaled_Markers <- MARKERS %>% scale()
      
      #If cluster estimation process needs to be performed
      if(!Force_N_Clusters){
        params_mbkm <- list(batch_size = Batch_size,
                            init_fraction = 1,
                            early_stop_iter = 10)
        #Run the specified test using each of the number of clusters
        print("Starting Cluster number stimation process")
        Optimal <- try(ClusterR::Optimal_Clusters_KMeans(Scaled_Markers,
                                                         max_clusters = Max_N_clusters_Batch,
                                                         num_init = N_initiations,
                                                         max_iters = Max_iterations,
                                                         initializer = "kmeans++",
                                                         criterion = "Adjusted_Rsquared",
                                                         plot_clusters = T,
                                                         mini_batch_params = params_mbkm,
                                                         verbose = T)
        )
        #Test if optimal number of clusters returned an error
        if(berryFunctions::is.error(Optimal)) {
          stop("Could not calculate best cluster number for the data provided. Please try another strategy")
        }
        
        #Proceed if all OK
        else{
          #Make the user decide the total number of clusters to be used in the final analysis
          N_Clusters <- menu(choices = as.character(1:Max_N_clusters_Batch),
                             title = paste0("Look at the plot generated, Then decide the appropiate number of clusters"))
          
          #Calculate the desired number of clusters with batch k menas
          print("Performing Batched K means algorithm")
          Clustering_result <- ClusterR::MiniBatchKmeans(Scaled_Markers,
                                                         clusters = as.double(N_Clusters),
                                                         batch_size = Batch_size,
                                                         num_init = N_initiations,
                                                         max_iters = Max_iterations,
                                                         init_fraction = 1,
                                                         initializer = "kmeans++",
                                                         early_stop_iter = 10,
                                                         verbose = T,
                                                         tol = 1e-07, #The required improvement rate to continue with the iterations (the lower the more iterations will be required)
                                                         CENTROIDS = NULL,
                                                         seed = 21)
      }
      
      }
      
      #If cluster estimation is NOT required
      if(Force_N_Clusters){
        #Calculate the desired number of clusters with batch k menas
        print("Performing Batched K means algorithm")
        Clustering_result <- ClusterR::MiniBatchKmeans(Scaled_Markers,
                                                       clusters = as.double(Max_N_clusters_Batch),
                                                       batch_size = Batch_size,
                                                       num_init = N_initiations,
                                                       max_iters = Max_iterations,
                                                       init_fraction = 1,
                                                       initializer = "kmeans++",
                                                       early_stop_iter = 10,
                                                       verbose = T,
                                                       tol = 1e-07, #The required improvement rate to continue with the iterations (the lower the more iterations will be required)
                                                       CENTROIDS = NULL,
                                                       seed = 21)
      }
      
      #Assign the cluster to each observation of MARKER
      pr_mb <- predict(object = Clustering_result, fuzzy = F, newdata = MARKERS %>% scale())
      pr_mb <- as_tibble(pr_mb)
      names(pr_mb) <- "Cluster"
      
      #Generate the data phenotypes tibble
      DATA_Clusters <- dplyr::bind_cols(Original_data, pr_mb)
      }
      
    ####GAUSSIAN MIXTURE MODELS####
    else if(Strategy == "GMM"){
      
      #Compute scaled markers
      Scaled_Markers <- MARKERS %>% scale()
      
      #If cluster estimation process required
      if(!Force_N_Clusters){
        #First we calculate a metric to decide the number of total phenotypes
        #Run the specified test using each of the number of clusters
        print("Starting Cluster number stimation process")
        Optimal <- try(ClusterR::Optimal_Clusters_GMM(Scaled_Markers,
                                                      criterion = Quality_metric,
                                                      max_clusters = Max_N_clusters_GMM,
                                                      dist_mode = GMM_Distance,
                                                      seed_mode = "random_subset",
                                                      km_iter = Max_iterations_km,
                                                      em_iter = Max_iterations_em,
                                                      verbose = TRUE,
                                                      var_floor = 1e-10,
                                                      plot_data = TRUE)
                       
                       
        )
        #Test if optimal number of clusters returned an error
        if(berryFunctions::is.error(Optimal)) {
          stop("Could not calculate best cluster number for the data provided. Please try another strategy")
        }
        #Proceed if all OK
        else{
          #Make the user decide the total number of clusters to be used in the final analysis
          N_Clusters <- menu(choices = as.character(1:Max_N_clusters_GMM),
                             title = paste0("Look at the plot generated, Then decide the appropiate number of clusters"))
          
          #Calculate the desired number of clusters with batch k menas
          print("Calculating Gaussian Mixed Model")
          Clustering_result <- ClusterR::GMM(Scaled_Markers,
                                             gaussian_comps = as.double(N_Clusters),
                                             dist_mode = GMM_Distance,
                                             seed_mode = "random_subset",
                                             km_iter = Max_iterations_km,
                                             em_iter = Max_iterations_em,
                                             verbose = TRUE,
                                             var_floor = 1e-10,
                                             full_covariance_matrices = FALSE
          )
        }
      }
      
      #If cluster
      if(Force_N_Clusters){
        #Calculate the desired number of clusters with batch k menas
        print("Calculating Gaussian Mixed Model")
        Clustering_result <- ClusterR::GMM(Scaled_Markers,
                                           gaussian_comps = as.double(Max_N_clusters_GMM),
                                           dist_mode = GMM_Distance,
                                           seed_mode = "random_subset",
                                           km_iter = Max_iterations_km,
                                           em_iter = Max_iterations_em,
                                           verbose = TRUE,
                                           var_floor = 1e-10,
                                           full_covariance_matrices = FALSE
        )
      }
      
      #Assign the cluster to each observation of MARKER
      pr_mb <- predict(object = Clustering_result, fuzzy = F, newdata = MARKERS %>% scale())
      pr_mb <- as_tibble(pr_mb)
      names(pr_mb) <- "Cluster"
      
      #Generate the data phenotypes tibble
      DATA_Clusters <- dplyr::bind_cols(Original_data, pr_mb)
    }
    
    ####CLARA CLUSTERING####
    else if(Strategy == "CLARA_clustering"){
      #Compute scaled markers
      Scaled_Markers <- MARKERS %>% scale()
      
      #CLUSTER ESTIMATION PROCESS REQUIRED
      if(!Force_N_Clusters){
        #First we calculate a metric to decide the number of total phenotypes
        print("Starting Cluster number stimation process")
        Optimal <-  try(ClusterR::Optimal_Clusters_Medoids(Scaled_Markers,
                                                           max_clusters = Max_N_clusters_CLARA,
                                                           distance_metric = Distance_CLARA,
                                                           criterion = "silhouette" ,
                                                           clara_samples = Samples_CLARA,
                                                           clara_sample_size = Sample_per_CLARA,
                                                           swap_phase = F,
                                                           threads = N_cores,
                                                           verbose = T,
                                                           plot_clusters = T)
        )
        #Test if optimal number of clusters returned an error
        if(berryFunctions::is.error(Optimal)) {
          stop("Could not calculate best cluster number for the data provided. Please try another strategy")
        }
        #Continue if everything OK
        else{
          #Make the user decide the total number of clusters to be used in the final analysis
          N_Clusters <- menu(choices = as.character(1:Max_N_clusters_CLARA),
                             title = paste0("Based on the plots generated and you previous choice, decide the appropiate number of final clusters"))
          print("Performing CLARA (Clustering Large Applications)")
          Clustering_result <- ClusterR::Clara_Medoids(Scaled_Markers,
                                                       clusters = as.double(N_Clusters),
                                                       samples = Samples_CLARA,
                                                       sample_size = Sample_per_CLARA,
                                                       distance_metric = Distance_CLARA,
                                                       threads = N_cores,
                                                       swap_phase = F,
                                                       fuzzy = FALSE,
                                                       verbose = T,
                                                       seed = 21)
          
        }
      }
      
      #CLUSTER ESTIMATION PROCESS NOT REQUIRED
      if(Force_N_Clusters){
        print("Performing CLARA (Clustering Large Applications)")
        Clustering_result <- ClusterR::Clara_Medoids(Scaled_Markers,
                                                     clusters = as.double(Max_N_clusters_CLARA),
                                                     samples = Samples_CLARA,
                                                     sample_size = Sample_per_CLARA,
                                                     distance_metric = Distance_CLARA,
                                                     threads = N_cores,
                                                     swap_phase = F,
                                                     fuzzy = FALSE,
                                                     verbose = T,
                                                     seed = 21)
      }
      
      #Assign the cluster to each observation of MARKER
      pr_mb <- predict(object = Clustering_result, fuzzy = F, newdata = MARKERS %>% scale())
      pr_mb <- as_tibble(pr_mb)
      names(pr_mb) <- "Cluster"
      
      #Generate the data phenotypes tibble
      DATA_Clusters <- dplyr::bind_cols(Original_data, pr_mb)
    }
    
    ####RETURN THE FINAL CLUSTERING RESULT####
    return(DATA_Clusters)
  }





