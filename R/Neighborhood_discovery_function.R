#' Performs neighborhood identification
#'
#' The function assigns every cell in a dataset to a neighborhood through a clustering process. Closest neighbor cell feature matrix calculated using [Tailored_Closest_neighbor_calculator()] must be used.
#'
#' @param DATA A dataframe or tibble containing closest neighbor information calculated using [Tailored_Closest_neighbor_calculator()].
#'
#' @param Allowed_max_Dist A numeric value indicating the allowed max distance to neighbors. Cells above threshold will be removed from analysis.
#' @param Allowed_avg_Dist A numeric value indicating the allowed average distance to neighbors. Cells above threshold will be removed from analysis.
#' @param Allowed_median_Dist A numeric value indicating the allowed median distance to neighbors. Cells above threshold will be removed from analysis.
#'
#' @param Perform_Dimension_reduction Logical value. Should Dimension Reduction be performed (see details).
#' @param Dimension_reduction Dimension reduction method. One of the following: PCA, TSNE, UMAP.
#' @param Dimension_reduction_prop A numeric value between 0 and 1 to indicate the percentage of the cells to be used in dimension computation (applicable for TSNE and UMAP).
#' @param Cluster_on_Reduced A logical value indicating if clustering should be performed on new dimensions.
#'
#' @param Strategy One of the following Consensus_Clustering, SOM, Graph_Based, K_Means_Meta_clustering, Batch_K_means, GMM or CLARA_clustering (see details).
#'
#' @param Max_N_neighborhoods If Strategy is Consensus_Clustering: Number of maximum neighborhoods that can be identified.
#' @param Consensus_reps If Strategy is Consensus_Clustering: Number of iterations to converge.
#' @param Consensus_p_Items If Strategy is Consensus_Clustering: Percentage of cells that you desire to sample in each iteration.
#' @param Consensus_Cluster_Alg If Strategy is Consensus_Clustering: Clustering algorithm to be used (’hc’ hierarchical (hclust), ’pam’ for paritioning around medoids, ’km’ for k-means).
#' @param Consensus_Distance If Strategy is Consensus_Clustering: Distance metric to be used (pearson(1 - Pearson correlation), spearman(1 - Spearman correlation), euclidean, binary, maximum, canberra, minkowski.
#' @param Consensus_Name If Strategy is Consensus_Clustering: Name of the folder that is going to be created in order to place the resulting graphs.
#'
#' @param Max_SOM_neighborhoods If Strategy is SOM: umber of maximum neighborhoods that can be identified.
#'
#' @param Nearest_neighbors_for_graph If strategy is Graph_Based: The number of closest neighbors to calculate the graph.
#' @param Graph_Method If strategy is Graph_Based: One of Louvain, Leiden, Greedy, WalkTrap, Spinglass, Leading_Eigen or Edge_Betweenness.
#' @param Graph_Resolution If strategy is Graph_Based: Used for Louvain and Leiden. 1 is default. The smaller the value, the larger the clusters will be.
#' @param N_steps If strategy is Graph_Based: Number of steps given in the WalkTrap algorithm.
#'
#' @param N_K_centroids If strategy is K_Means_Meta_clustering: Number of centroids to perform K means.
#' @param Max_N_neighborhoods_Meta If strategy is K_Means_Meta_clustering: Number of maximum neighborhoods that can be identified.
#' @param Consensus_reps_Meta If strategy is K_Means_Meta_clustering: Number of iterations to converge.
#' @param Consensus_p_Items_Meta If strategy is K_Means_Meta_clustering: Percentage of cells that you desire to sample in each iteration.
#' @param Consensus_Name_Meta If strategy is K_Means_Meta_clustering: Name of the folder that is going to be created in order to place the resulting graphs.
#'
#' @param Batch_size If strategy is Batch_K_means: Number of cells to be included in each random batch.
#' @param Max_N_neighborhoods_Batch If strategy is Batch_K_means: Number of maximum neighborhoods that can be identified.
#' @param N_initiations If strategy is Batch_K_means: Number of times the algorithm is going to be tried to find the best clustering result.
#' @param Max_iterations If strategy is Batch_K_means: Max number of iterations in each try.
#'
#' @param Quality_metric If strategy is GMM:T he quality measure used to test the number of clusters ("AIC" or "BIC").
#' @param Max_N_neighborhoods_GMM If strategy is GMM: Number of maximum neighborhoods that can be identified.
#' @param Max_iterations_km If strategy is GMM: Number of max iterations in the K means clustering performed.
#' @param Max_iterations_em If strategy is GMM: Number of max iterations in the Expectation Maximization algorithm.
#' @param GMM_Distance If strategy is GMM: Distance metric used in the model ("eucl_dist" or "maha_dist").
#'
#' @param Samples_CLARA If strategy is CLARA_clustering: Number of samples the CLARA algorithm is going to use to be calculated.
#' @param Sample_per_CLARA If strategy is CLARA_clustering: Percentage (from 0 to 1) of the total cells that are going to be allocated to each sample.
#' @param Max_N_neighborhoods_CLARA If strategy is CLARA_clustering: Number of maximum neighborhoods that can be identified.
#' @param Distance_CLARA If strategy is CLARA_clustering: Distance metric used in the model (euclidean, manhattan, chebyshev, canberra, braycurtis, pearson_correlation, simple_matching_coefficient, minkowski, hamming, jaccard_coefficient, Rao_coefficient, mahalanobis, cosine)
#' @param N_cores If strategy is CLARA_clustering: Number of cores to parallelize your computation
#'
#' @details
#' Dimension reduction can be performed using PCA (svd::propack.svd function), t-SNE (snifter::fitsne function) and UMAP (uwot::tumap function). For t-SNE and UMAP a model can be build using a subset of data and then predicting coordinates for all the cells. This can be more computationally efficient.
#'
#' Consensus clustering is performed using the ConsensusClusterPlus::ConsensusClusterPlus function.
#'
#' Self Organizing Maps clustering is performed using the FlowSOM::FlowSOM function.
#'
#' For graph based clustering Nearest neighbors graphs are built using bluster::makeSNNGraph and clustered using functions included in the igraph package.
#'
#' K_Means_Meta_clustering first summarizes cell feature matrix observations using K means algorithm and the performs Consensus Clustering. Afterwards results are generalized to all cells.
#'
#' Batch K-means, Gaussian Mixture Models and Clustering Large Applications are all based on the ClusterR package.
#'
#' @seealso [Tailored_Closest_neighbor_calculator()], [DATA_neighborhoods_renamer()], [Neighborhood_Quantifier()], [Neighborhood_voting_function()], [Tiled_neighborhoods_graphicator()]
#'
#' @returns Returns a tibble with cell features and a column named 'Neighborhood_assignment' containing cell neighborhood.
#' If dimension reduction has been performed returns a list with the cell feature dataset as above and a tibble containing dimension reduction coordinates.
#'
#' @export

Neighborhood_discovery_function <-
  function(DATA = NULL,

           Allowed_max_Dist = NULL,
           Allowed_avg_Dist = NULL,
           Allowed_median_Dist = NULL,

           #Dimension reduction
           Perform_Dimension_reduction = NULL,
           Dimension_reduction = NULL,
           Dimension_reduction_prop = NULL,
           Cluster_on_Reduced = NULL,

           #Strategy for clustering
           Strategy = NULL,

           #Parameters for Consensus Clustering
           Max_N_neighborhoods = NULL,
           Consensus_reps = NULL,
           Consensus_p_Items = NULL,
           Consensus_Cluster_Alg = NULL,
           Consensus_Distance = NULL,
           Consensus_Name = NULL,

           #Parameters for Self-Organizing Maps
           Max_SOM_neighborhoods = NULL, #Maximum number of clusters (neighborhoods) to try in the algorithm

           #Parameters for Graph methods
           Graph_type = NULL,
           Graph_Method = NULL,
           Nearest_neighbors_for_graph = NULL,
           Graph_Resolution = NULL,
           Graph_Distance_method = NULL,
           N_steps = NULL,

           #Parameters for K means Meta Clustering
           N_K_centroids = NULL, #Number of centroids to perform K means
           Max_N_neighborhoods_Meta = NULL, #Number of maximum clusters (neighborhoods) that you desire to find
           Consensus_reps_Meta = NULL, #Number of iterations of the algorithm to try to converge
           Consensus_p_Items_Meta = NULL, #Percentage of cells that you desire to sample in each iteration
           Consensus_Name_Meta = NULL, #Name of the folder that is going to be created in order to place the resulting graphs

           #Parameters for Batched K means
           Batch_size = NULL, #The number of cells to be included in each random batch
           Max_N_neighborhoods_Batch = NULL, #Number of maximum clusters (neighborhoods) that you desire to find
           Percentage_centroid_initiation = NULL,
           N_initiations = NULL, #Number of times the algorithm is going to be tried to find the best clustering result
           Max_iterations = NULL, #Max number of iterations in each try

           #Parameters for Gaussian Mixture Model
           Quality_metric = NULL, #The quality measure used to test the number of clusters ("AIC" or "BIC")
           Max_N_neighborhoods_GMM = NULL, #Number of maximum clusters (phenotypes) that you desire to find
           Max_iterations_km = NULL, #Number of max iterations in the K means clustering performed
           Max_iterations_em = NULL, #Number of max iterations in the Expectation Maximization algorithm
           GMM_Distance = NULL, #Distance metric to use in the model ("eucl_dist" or "maha_dist")

           #Parameters for CLARA clustering
           Samples_CLARA = NULL, #Number of samples the CLARA algorithm is going to use to be calculated
           Sample_per_CLARA = NULL, #Percentage (from 0 to 1) of the total cells that are going to be allocated to each sample
           Max_N_neighborhoods_CLARA = NULL, #Number of maximum clusters (neighborhoods) that you desire to find
           Distance_CLARA = NULL, #euclidean, manhattan, chebyshev, canberra, braycurtis, pearson_correlation,
           #simple_matching_coefficient, minkowski, hamming, jaccard_coefficient, Rao_coefficient, mahalanobis, cosine
           N_cores = NULL #Number of cores to parallelize your computation
  ) {
    #Check arguments by generating a argument check vector and message vector
    Argument_checker <- c(DATA_OK = all(c("min_DIST", "max_DIST", "avg_DIST", "median_DIST") %in% names(DATA)),
                          Distance_OK = all(c(is.numeric(c(Allowed_max_Dist, Allowed_avg_Dist, Allowed_median_Dist)),
                                              c(Allowed_max_Dist, Allowed_avg_Dist, Allowed_median_Dist) > 0)),
                          Strategy_OK = Strategy %in% c("Consensus_Clustering", "SOM", "Graph_Based", "K_Means_Meta_clustering", "Batch_K_means", "GMM", "CLARA_clustering")
    )
    Stop_messages <- c(DATA_OK = "DATA must be obtained from Closest_neighbor_calculator or Tailored_Closest_neighbor_calculator functions",
                       Distance_OK = "Allowed distances must be numeric and > 0",
                       Strategy_OK = "Strategy must be one of the following Consensus_Clustering, SOM, Graph_Based, K_Means_Meta_clustering, Batch_K_means, GMM, CLARA_clustering")
    #Check arguments and stop if necessary
    if(!all(Argument_checker)){
      stop(cat(Stop_messages[!Argument_checker],
               fill = sum(!Argument_checker)))
    }
    if(!is.logical(Perform_Dimension_reduction)) stop("Perform_Dimension_reduction must be a logical value")
    if(Perform_Dimension_reduction){
      if(!Dimension_reduction %in% c("UMAP", "TSNE", "PCA")) stop("Dimension_reduction must be one of the following: UMAP, TSNE, PCA")
      if(!all(is.numeric(Dimension_reduction_prop), Dimension_reduction_prop > 0, Dimension_reduction_prop <= 1)) stop("Dimension_reduction_prop must be a numeric value between 0 and 1")
    }
    if(!is.logical(Cluster_on_Reduced)) stop("Cluster_on_Reduced must be a logical value")
    if(Cluster_on_Reduced){
      if(!Perform_Dimension_reduction) stop("If Clustering needst o be performed on Dimension reduced data please set Perform_Dimension_reduction to TRUE")
    }

    #Check specific arguments and suggested packages
    #Check arguments for Consensus Clustering
    if(Strategy == "Consensus_Clustering"){
      #Check suggested packages
      if(!requireNamespace("ConsensusClusterPlus", quietly = TRUE)) stop(
        paste0("ConsensusClusterPlus Bioconductor package is required to execute the function. Please install using the following code: ",
               expression({
                 if (!require("BiocManager", quietly = TRUE))
                   install.packages("BiocManager")

                 BiocManager::install("ConsensusClusterPlus")
               })
        )
      )
      #Check arguments by generating a argument check vector and message vector
      Argument_checker <- c(Max_N_neighborhoods_OK = (Max_N_neighborhoods >= 2 & Max_N_neighborhoods%%1 == 0),
                            Consensus_reps_OK = (Consensus_reps >= 1 & Consensus_reps%%1 == 0),
                            Consensus_p_Items_OK = (Consensus_p_Items > 0 & Consensus_p_Items <= 1),
                            Consensus_Cluster_Alg_OK = Consensus_Cluster_Alg %in% c("hc", "pam", "km"),
                            Consensus_Distance_OK = Consensus_Distance %in% c("pearson", "spearman", "euclidean", "binary", "maximum", "canberra", "minkowski"),
                            Consensus_Name_OK = is.character(as.character(Consensus_Name))
      )
      Stop_messages <- c(Max_N_neighborhoods_OK = "Max_N_neighborhoods must be an integer value > 1",
                         Consensus_reps_OK = "Consensus_reps_OK must be an integer value > 0",
                         Consensus_p_Items_OK = "Consensus_p_Items must be a numeric value > 0 and lower than 1",
                         Consensus_Cluster_Alg_OK = "Consensus_Cluster_Alg must be one of the following: hc, pam, km",
                         Consensus_Distance_OK = "Consensus_Distance must be one the following: pearson, spearman, euclidean, binary, maximum, canberra, minkowski",
                         Consensus_Name_OK = "Consensus_Name must ve a character value")
      #Check arguments and stop if necessary
      if(!all(Argument_checker)){
        stop(cat(Stop_messages[!Argument_checker],
                 fill = sum(!Argument_checker)))
      }
    }
    #Check arguments for Self Organizingpurrr::maps
    if(Strategy == "SOM"){
      #Check suggested packages
      if(!requireNamespace("FlowSOM", quietly = TRUE)) stop(
        paste0("FlowSOM Bioconductor package is required to execute the function. Please install using the following code: ",
               expression({
                 if (!require("BiocManager", quietly = TRUE))
                   install.packages("BiocManager")

                 BiocManager::install("FlowSOM")
               })
        )
      )
      #Check arguments
      if(!(Max_SOM_neighborhoods > 1 & Max_SOM_neighborhoods%%1 == 0)) stop("Max_SOM_neighborhoods must be an integer value > 1")
    }
    #Check arguments for Graph-Based clustering
    if(Strategy == "Graph_Based"){
      #Check suggested packages
      if(!requireNamespace("bluster", quietly = TRUE)) stop(
        paste0("bluster Bioconductor package is required to execute the function. Please install using the following code: ",
               expression({
                 if (!require("BiocManager", quietly = TRUE))
                   install.packages("BiocManager")

                 BiocManager::install("bluster")
               })
        )
      )
      if(!requireNamespace("igraph", quietly = FALSE)) stop(
        paste0("igraph CRAN package is required to execute the function. Please install using the following code: ",
               expression(install.packages("igraph")))
      )
      #Check arguments by generating a argument check vector and message vector
      Argument_checker <- c(Graph_type_OK = Graph_type %in% c("Complete", "SNN"),
                            Graph_Distance_method_OK = (Graph_Distance_method %in% c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski")),
                            Graph_Method_OK = Graph_Method %in% c("Louvain", "Leiden", "Optimal", "Greedy", "WalkTrap", "Spinglass", "Leading_Eigen", "Edge_Betweenness"),
                            Graph_Resolution_OK = all(is.numeric(Graph_Resolution), Graph_Resolution > 0),
                            N_steps_OK = any(is.null(N_steps), (N_steps >=1 & N_steps%%1 == 0))
      )
      Stop_messages <- c(Graph_type_OK = "Graph_type should be one of the following: Complete, SNN",
                         Nearest_neighbors_for_graph = "Nearest_neighbors_for_graph must be an integer value > 0",
                         Graph_Method = "Graph_Method must be one of the following: Louvain, Leiden, Greedy, WalkTrap, Spinglass, Leading_Eigen, Edge_Betweenness",
                         Graph_Resolution = "Graph_Resolution must be a numeric value > 0",
                         N_steps = "N_steps must be a integer value > 0")
      #Check arguments and stop if necessary
      if(!all(Argument_checker)){
        stop(cat(Stop_messages[!Argument_checker],
                 fill = sum(!Argument_checker)))
      }
      #Check specific argument of SNN graphs
      if(Graph_type == "SNN"){
        if(!all(Nearest_neighbors_for_graph%%1 == 0, Nearest_neighbors_for_graph > 0)) stop("Nearest_neighbors_for_graph should be a integer value > 0")
      }
    }
    #Check arguments for K means meta clustering
    if(Strategy == "K_Means_Meta_clustering"){
      #Check suggested packages
      if(!requireNamespace("ConsensusClusterPlus", quietly = TRUE)) stop(
        paste0("ConsensusClusterPlus Bioconductor package is required to execute the function. Please install using the following code: ",
               expression({
                 if (!require("BiocManager", quietly = TRUE))
                   install.packages("BiocManager")

                 BiocManager::install("ConsensusClusterPlus")
               })
        )
      )
      #Check arguments
      Argument_checker <- c(Max_N_neighborhoods_Meta_OK = (Max_N_neighborhoods_Meta >= 2 & Max_N_neighborhoods_Meta%%1 == 0),
                            Consensus_reps_Meta_OK = (Consensus_reps_Meta >= 1 & Consensus_reps_Meta%%1 == 0),
                            Consensus_p_Items_Meta_OK = (Consensus_p_Items_Meta > 0 & Consensus_p_Items_Meta <= 1),
                            Consensus_Name_Meta_OK = is.character(as.character(Consensus_Name_Meta))
      )
      Stop_messages <- c(Max_N_neighborhoods_Meta_OK = "Max_N_neighborhoods_Meta must be an integer value > 1",
                         Consensus_reps_Meta_OK = "Consensus_reps_Meta must be an integer value > 0",
                         Consensus_p_Items_Meta_OK = "Consensus_p_Items_Meta must be a numeric value > 0 and lower than 1",
                         Consensus_Name_Meta_OK = "Consensus_Name_Meta must ve a character value"
      )
      #Check arguments and stop if necessary
      if(!all(Argument_checker)){
        stop(cat(Stop_messages[!Argument_checker],
                 fill = sum(!Argument_checker)))
      }
    }
    #Check arguments for Batched K means
    if(Strategy == "Batch_K_means"){
      #Check suggested packages
      if(!requireNamespace("ClusterR", quietly = FALSE)) stop(
        paste0("ClusterR CRAN package is required to execute the function. Please install using the following code: ",
               expression(install.packages("ClusterR")))
      )
      #Check arguments
      Argument_checker <- c(Max_N_neighborhoods_Batch_OK = (Max_N_neighborhoods_Batch >= 2 & Max_N_neighborhoods_Batch%%1 == 0),
                            N_initiations_OK = (N_initiations >= 1 & N_initiations%%1 == 0),
                            Max_iterations_OK = (Max_iterations%%1 == 0)
      )
      Stop_messages <- c(Max_N_phenotypes_Batch_OK = "Max_N_neighborhoods_Batch must be an integer value > 1",
                         N_initiations_OK = "N_initiations must be an integer value > 0",
                         Max_iterations_OK = "Max_iterations must be an integer value > 0"
      )
      #Check arguments and stop if necessary
      if(!all(Argument_checker)){
        stop(cat(Stop_messages[!Argument_checker],
                 fill = sum(!Argument_checker)))
      }
    }
    #Check arguments for Gaussian mixture models
    if(Strategy == "GMM"){
      #Check suggested packages
      if(!requireNamespace("ClusterR", quietly = FALSE)) stop(
        paste0("ClusterR CRAN package is required to execute the function. Please install using the following code: ",
               expression(install.packages("ClusterR")))
      )
      #Check arguments
      Argument_checker <- c(Quality_metric_OK = Quality_metric %in% c("AIC", "BIC"),
                            Max_N_neighborhoods_GMM_OK = (Max_N_neighborhoods_GMM >= 2 & Max_N_neighborhoods_GMM%%1 == 0),
                            Max_iterations_km_OK = (Max_iterations_km >= 1 & Max_iterations_km%%1 == 0),
                            Max_iterations_em_OK = (Max_iterations_em >= 1 & Max_iterations_em%%1 == 0),
                            GMM_Distance_OK = GMM_Distance %in% c("eucl_dist", "maha_dist")
      )
      Stop_messages <- c(Quality_metric_OK = "Quality_metric must be one of the following: AIC, BIC",
                         Max_N_phenotypes_GMM_OK = "Max_N_phenotypes must be an integer value > 1",
                         Max_iterations_km_OK = "Max_iterations_km must be an integer value > 1",
                         Max_iterations_em_OK = "Max_iterations_em must be an integer value > 1",
                         GMM_Distance_OK = "GMM_Distance must be one of the following: eucl_dist, maha_dist"
      )
      #Check arguments and stop if necessary
      if(!all(Argument_checker)){
        stop(cat(Stop_messages[!Argument_checker],
                 fill = sum(!Argument_checker)))
      }
    }
    #Check arguments for CLARA clustering
    if(Strategy == "CLARA_clustering"){
      #Check suggested packages
      if(!requireNamespace("ClusterR", quietly = FALSE)) stop(
        paste0("ClusterR CRAN package is required to execute the function. Please install using the following code: ",
               expression(install.packages("ClusterR")))
      )
      #Check arguments
      Argument_checker <- c(Samples_CLARA_OK = (Samples_CLARA >= 1 & Samples_CLARA%%1 == 0),
                            Sample_per_CLARA_OK = (Sample_per_CLARA > 0 & Sample_per_CLARA <= 1),
                            Max_N_neighborhoods_CLARA_OK = (Max_N_neighborhoods_CLARA >= 2 & Max_N_neighborhoods_CLARA%%1 == 0),
                            Distance_CLARA_OK = Distance_CLARA %in% c("euclidean", "manhattan", "chebyshev", "canberra", "braycurtis",
                                                                      "pearson_correlation", "simple_matching_coefficient", "minkowski",
                                                                      "hamming", "jaccard_coefficient", "Rao_coefficient", "mahalanobis", "cosine"),
                            N_cores_OK = (N_cores >= 1 & N_cores%%1 == 0)
      )
      Stop_messages <- c(Samples_CLARA_OK = "Samples_CLARA must be an integer value > 0",
                         Sample_per_CLARA_OK = "Sample_per_CLARA must be a numeric value between 0 and 1",
                         Max_N_neighborhoods_CLARA_OK = "Max_N_neighborhoods_CLARA must be an integer value > 1",
                         Distance_CLARA_OK = "Distance_CLARA must be one of the following: euclidean, manhattan, chebyshev, canberra, braycurtis, pearson_correlation, simple_matching_coefficient, minkowski, hamming, jaccard_coefficient, Rao_coefficient, mahalanobis, cosine",
                         N_cores_OK = "N_cores must be an integer value > 0"
      )
      #Check arguments and stop if necessary
      if(!all(Argument_checker)){
        stop(cat(Stop_messages[!Argument_checker],
                 fill = sum(!Argument_checker)))
      }
    }
    #Check packages for dimension reduction
    if(Perform_Dimension_reduction){
      if(Dimension_reduction == "PCA"){
        if(!requireNamespace("svd", quietly = FALSE)) stop(
          paste0("svd CRAN package is required to execute the function. Please install using the following code: ",
                 expression(install.packages("svd")))
        )
      }
      if(Dimension_reduction == "TSNE"){
        if(!requireNamespace("snifter", quietly = TRUE)) stop(
          paste0("snifter Bioconductor package is required to execute the function. Please install using the following code: ",
                 expression({
                   if (!require("BiocManager", quietly = TRUE))
                     install.packages("BiocManager")

                   BiocManager::install("snifter")
                 })
          )
        )
      }
      if(Dimension_reduction == "UMAP"){
        if(!requireNamespace("uwot", quietly = FALSE)) stop(
          paste0("uwot CRAN package is required to execute the function. Please install using the following code: ",
                 expression(install.packages("uwot")))
        )
      }
    }
    #Check complex heatmap package
    if(!requireNamespace("ComplexHeatmap", quietly = TRUE)) stop(
      paste0("ComplexHeatmap Bioconductor package is required to execute the function. Please install using the following code: ",
             expression({
               if (!require("BiocManager", quietly = TRUE))
                 install.packages("BiocManager")

               BiocManager::install("ComplexHeatmap")
             })
      )
    )

    #Start executing the code
    print("Preparing data for analysis")
    #Import closest neighbor data
    DATA_Neighbors <- DATA

    #Filter out cells not whitin a neighborhood
    DATA_Neighbors <- DATA_Neighbors %>% dplyr::filter(max_DIST <= Allowed_max_Dist, avg_DIST <= Allowed_avg_Dist, median_DIST <= Allowed_median_Dist)

    #We select only neighbor columns
    Neighbor_patterns <- DATA_Neighbors %>% dplyr::select(-c(1:5), -contains("N_neighbors"), -c((ncol(DATA_Neighbors)-4):(ncol(DATA_Neighbors))))
    Neighbor_patterns_scaled <- Neighbor_patterns %>% scale()

    #Perform dimension reduction if required
    if(Perform_Dimension_reduction){
      #First PCA
      if(Dimension_reduction == "PCA"){
        if(Dimension_reduction_prop != 1) stop("PCA must be performed using Dimension_reduction_prop = 1")
        print("Generating PCA projections")
        #Scale and turn into matrix
        DATA_matrix <- Neighbor_patterns_scaled %>% as.matrix()
        Result_PCA <- svd::propack.svd(DATA_matrix, neig = 2)$u
        DATA_Reduction <- tibble(DIMENSION_1 = unlist(Result_PCA[,1]), DIMENSION_2 = unlist(Result_PCA[,2]))
      }

      #Second TSNE
      if(Dimension_reduction == "TSNE"){
        if(Dimension_reduction_prop == 1) {
          print("Generating TSNE projections")
          if(nrow(DATA) > 50000) print("Warning! Data set contains more than 50K observations. tSNE embedding can take a long time")
          #scale and turn into matrix
          DATA_matrix <- Neighbor_patterns_scaled %>% as.matrix()
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
          DATA_Reduction <-dplyr::bind_cols(DIMENSION_1 = unlist(Result_TSNE[,1]), DIMENSION_2 = unlist(Result_TSNE[,2]))
        }

        if(Dimension_reduction_prop != 1) {
          print("Generating TSNE projections")
          DATA_matrix <- Neighbor_patterns_scaled %>% as_tibble %>% dplyr::sample_frac(size = Dimension_reduction_prop) %>% as.matrix()
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
                                     new = Neighbor_patterns_scaled %>% as.matrix(),
                                     old = DATA_matrix)
          DATA_Reduction <-dplyr::bind_cols(DIMENSION_1 = unlist(Coords[,1]), DIMENSION_2 = unlist(Coords[,2]))
        }
      }

      #Third UMAP
      if(Dimension_reduction == "UMAP"){
        if(Dimension_reduction_prop == 1) {
          print("Generating UMAP projections")
          if(nrow(DATA) > 50000) print("Warning! Data set contains more than 50K observations. UMAP embedding can take some time")
          #scale and turn into matrix
          DATA_matrix <- Neighbor_patterns_scaled %>% as.matrix()
          Result_UMAP <- uwot::tumap(DATA_matrix, n_components = 2L)
          DATA_Reduction <-dplyr::bind_cols(DIMENSION_1 = unlist(Result_UMAP[,1]), DIMENSION_2 = unlist(Result_UMAP[,2]))
        }

        if(Dimension_reduction_prop != 1) {
          print("Generating UMAP projections")
          DATA_matrix <- Neighbor_patterns_scaled %>% as_tibble() %>% dplyr::sample_frac(size = Dimension_reduction_prop) %>% as.matrix()
          if(nrow(DATA_matrix) > 50000) print("Warning! Data set contains more than 50K observations. UMAP embedding can take aome time")
          #scale and turn into matrix
          Result_UMAP <- uwot::tumap(DATA_matrix, n_components = 2L, ret_model = TRUE)
          Coords <- uwot::umap_transform(X = Neighbor_patterns_scaled %>% as.matrix(),
                                         model = Result_UMAP)
          DATA_Reduction <-dplyr::bind_cols(DIMENSION_1 = unlist(Coords[,1]), DIMENSION_2 = unlist(Coords[,2]))
        }
      }
    }

    #Generate a specific version of Markers with dimension reduction data for clustering
    if(Cluster_on_Reduced)  Neighbor_patterns_scaled <- DATA_Reduction


    #Define what to do if consensus clustering is required
    if(Strategy == "Consensus_Clustering"){
      #Perform consensus clustering
      Neighborhood_result <- try(ConsensusClusterPlus::ConsensusClusterPlus(t(as.matrix((Neighbor_patterns_scaled))),
                                                                            maxK = Max_N_neighborhoods,
                                                                            reps = Consensus_reps,
                                                                            pItem = Consensus_p_Items,
                                                                            pFeature = 1,
                                                                            title = Consensus_Name,
                                                                            clusterAlg = Consensus_Cluster_Alg,
                                                                            distance = Consensus_Distance,
                                                                            plot = "png",
                                                                            verbose = T)
      )

      if(berryFunctions::is.error(Neighborhood_result)) stop("Consensus Clustering Failed. Your data is probably too large for this method. Please try another strategy")

      #Make the user decide the number of neighborhoods according to results
      N_Neighbors <- menu(choices = as.character(1:Max_N_neighborhoods), title = paste0("Check the results at: ", getwd(), ". Then decide the appropiate number of Neighborhoods"))
      Neighbor_patterns <- Neighbor_patterns %>%dplyr::mutate(Neighborhood_assignment = Neighborhood_result[[as.double(N_Neighbors)]][["consensusClass"]])
    }
    #Define what to do if SOM is required
    if(Strategy == "SOM"){
      print("Executing Self Organizing Map algorithm")
      #Transform data into a scaled matrix and perform Self Organizing Map
      SOM_results <- try(FlowSOM::FlowSOM(as.matrix(Neighbor_patterns_scaled),
                                          scale = F,
                                          colsToUse = 1:ncol(Neighbor_patterns_scaled),
                                          maxMeta = Max_SOM_neighborhoods, #To find optimal meta clusters
                                          silent = F,
                                          seed = 21)
      )
      #Test if SOM returned an error
      if(berryFunctions::is.error(SOM_results)) {
        stop("Data is too large for Self-Organizing Maps. Please try another strategy")
      }
      else{
        #Assign phenotypes to each cell
        Neighbor_patterns <- Neighbor_patterns %>%dplyr::mutate(Neighborhood_assignment = FlowSOM::GetMetaclusters(SOM_results))
      }
    }
    #Then define what to do if graph based clustering is required
    if(Strategy == "Graph_Based"){
      #Generate graphs according to user defined preferences
      if(Graph_type == "Complete"){
        print("Generating the complete graph")
        #Calculate distance matrix and then calculate the graph
        #We define the number and ID of edges of the graph
        Graph_tibble <- as_tibble(expand.grid.unique(1:nrow(Neighbor_patterns_scaled), 1:nrow(Neighbor_patterns_scaled)))
        names(Graph_tibble) <- c("from", "to")
        Graph_tibble <- Graph_tibble %>%dplyr::mutate(ID = stringr::str_c(from, to, sep = "_"))

        #We determine the distance between nodes that will be the features of the edges
        DISTANCE_MATRIX <- as_tibble(as.matrix(dist(Neighbor_patterns_scaled, method = Graph_Distance_method)))
        DISTANCE_MATRIX <- DISTANCE_MATRIX %>%dplyr::mutate(from = as.character(1:nrow(DISTANCE_MATRIX)))
        DISTANCE_MATRIX <- DISTANCE_MATRIX[c(ncol(DISTANCE_MATRIX), 2:(ncol(DISTANCE_MATRIX)-1))] %>% tidyr::pivot_longer(-1, names_to = "to", values_to = "weight") %>%
          dplyr::mutate(ID = stringr::str_c(from, to, sep = "_")) %>% dplyr::select(-from, -to)

        #We bind the edges to their features (distance) and we build the graph
        GRAPH_DF <-dplyr::left_join(Graph_tibble, DISTANCE_MATRIX, by = "ID") %>% dplyr::select(-ID)
        GRAPH_DF <- GRAPH_DF %>%dplyr::mutate(weight = 1/weight)
        Neighborhood_ID <- tibble(Name = as.character(1:nrow(Neighbor_patterns_scaled)))
        Neighborhood_pattern_graph <- igraph::graph_from_data_frame(GRAPH_DF, directed = F, vertices = Neighborhood_ID)
      }
      if(Graph_type == "SNN"){
        print("Generating the SNN graph")
        #Transform data into a nearest neighbor graph
        Neighborhood_pattern_graph <- try(bluster::makeSNNGraph(as.matrix(Neighbor_patterns_scaled),
                                                                k = Nearest_neighbors_for_graph)
        )

        #Test if Graph construction process returned an error
        if(berryFunctions::is.error(Neighborhood_pattern_graph)) {
          stop("Data is too large to build a graph. Please try another strategy")
        }
      }

      print("Performing graph-based clustering")
      #Cluster the graph with louvain or leiden clustering
      if(Graph_Method == "Louvain") {
        Neighbor_patterns <- Neighbor_patterns %>%dplyr::mutate(Neighborhood_assignment = igraph::cluster_louvain(Neighborhood_pattern_graph,
                                                                                                                  weights = NULL,
                                                                                                                  resolution = Graph_Resolution)$membership)
      }

      else if(Graph_Method == "Leiden") {
        Neighbor_patterns <- Neighbor_patterns %>%dplyr::mutate(Neighborhood_assignment = igraph::cluster_leiden(Neighborhood_pattern_graph,
                                                                                                                 objective_function = "modularity",
                                                                                                                 weights = NULL,
                                                                                                                 resolution = Graph_Resolution,
                                                                                                                 beta = 0.01,
                                                                                                                 initial_membership = NULL,
                                                                                                                 n_iterations = 100,
                                                                                                                 vertex_weights = NULL)$membership)
      }

      else if(Graph_Method == "Optimal"){
        Neighbor_patterns <- Neighbor_patterns %>%dplyr::mutate(Neighborhood_assignment = igraph::cluster_optimal(Neighborhood_pattern_graph)$membership)
      }

      else if(Graph_Method == "Greedy"){
        Neighbor_patterns <- Neighbor_patterns %>%dplyr::mutate(Neighborhood_assignment = igraph::cluster_fast_greedy(Neighborhood_pattern_graph)$membership)
      }

      else if(Graph_Method == "WalkTrap"){
        Neighbor_patterns <- Neighbor_patterns %>%dplyr::mutate(Neighborhood_assignment = igraph::cluster_walktrap(Neighborhood_pattern_graph,
                                                                                                                   steps = N_steps,
                                                                                                                   membership = T)$membership)
      }

      else if (Graph_Method == "Spinglass") {
        Neighbor_patterns <- Neighbor_patterns %>%dplyr::mutate(Neighborhood_assignment = igraph::cluster_spinglass(Neighborhood_pattern_graph,
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
        Neighbor_patterns <- Neighbor_patterns %>%dplyr::mutate(Neighborhood_assignment = igraph::cluster_leading_eigen(Neighborhood_pattern_graph)$membership)
      }

      else if(Graph_Method == "Edge_Betweenness"){
        Neighbor_patterns <- Neighbor_patterns %>%dplyr::mutate(Neighborhood_assignment = igraph::cluster_edge_betweenness(Neighborhood_pattern_graph,
                                                                                                                           weights = NULL,
                                                                                                                           directed = FALSE,
                                                                                                                           edge.betweenness = FALSE,
                                                                                                                           merges = FALSE,
                                                                                                                           bridges = FALSE,
                                                                                                                           modularity = FALSE,
                                                                                                                           membership = TRUE)$membership)
      }
    }
    #Define what to do if K means meta clustering is required
    if(Strategy == "K_Means_Meta_clustering"){
      print("Performing initial K-means algorithm")
      if(N_K_centroids >= nrow(Neighbor_patterns_scaled)) stop(paste0("N_K_centroids should be smaller than: ", nrow(Neighbor_patterns_scaled)))

      #First we need to perform K means Clustering
      cl <- try(kmeans(as.matrix(Neighbor_patterns_scaled), #Scale it and turn it into a matrix
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
        DATA_filter_Markers <- Neighbor_patterns_scaled %>%dplyr::mutate(K_means_Cl = cl$cluster)

        #Prepare data for Meta-Clustering
        #Create a tibble with the K means centroids and the format it for Consensus clustering
        K_medoids <- as_tibble(cl$centers) %>%dplyr::mutate(K_means_Cl = 1:nrow(as_tibble(cl$centers)))
        tK_medoids <- K_medoids %>% dplyr::select(-K_means_Cl) %>% as.matrix %>% t

        print("Perorming Consensus Clustering")
        #Perform Consensus clustering with hierarchical clustering
        HC <- try(ConsensusClusterPlus::ConsensusClusterPlus(tK_medoids,
                                                             maxK = Max_N_neighborhoods_Meta,
                                                             reps = Consensus_reps_Meta,
                                                             pItem = Consensus_p_Items_Meta,
                                                             pFeature = 1,
                                                             title = Consensus_Name_Meta,
                                                             distance = "euclidean",
                                                             clusterAlg = "pam",
                                                             plot = "png",
                                                             verbose = T))
        #Test if consensus clustering returned an error
        if(berryFunctions::is.error(HC)) {
          stop("Data is too large for Meta Clustering. Please try another strategy or select a smaller N_K_centroids value")
        }
        else {
          #Make the user decide the number of neighborhoods according to results
          N_Phenotypes<- menu(choices = as.character(1:Max_N_neighborhoods_Meta), title = paste0("Check the results at: ", getwd(), ". Then decide the appropiate number of Neighborhoods"))

          #Bind the final Neighborhood_assignment to the K medoids tibble
          K_medoids <- K_medoids %>%dplyr::mutate(Neighborhood_assignment = HC[[as.double(N_Phenotypes)]][["consensusClass"]])
          K_medoids_for_join <- K_medoids %>% dplyr::select(K_means_Cl, Neighborhood_assignment)

          #Bind The DATA and the K_meoids to obtain the final matrix
          Neighbor_patterns <-dplyr::left_join(DATA_filter_Markers, K_medoids_for_join, by = "K_means_Cl") %>% dplyr::select(-K_means_Cl)
          #Reorder the columns
          Neighbor_patterns <- Neighbor_patterns[c(1:(ncol(Neighbor_patterns)-1), ncol(Neighbor_patterns))]
        }
      }
    }
    #Define what to do if Batch K means is required
    if(Strategy == "Batch_K_means"){
      #First we calculate a metric to decide the number of total phenotypes
      #Specify the params
      params_mbkm <- list(batch_size = Batch_size,
                          init_fraction = 1,
                          early_stop_iter = 10)
      print("Starting Cluster number stimation process")

      if(Batch_size >= nrow(Neighbor_patterns_scaled)) stop(paste0("Batch_size should be smaller than: ", nrow(Neighbor_patterns_scaled)))
      #Run the specified test using each of the number of clusters
      Optimal <- try(ClusterR::Optimal_Clusters_KMeans(Neighbor_patterns_scaled,
                                                       max_clusters = Max_N_neighborhoods_Batch,
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
        N_Phenotypes<- menu(choices = as.character(1:Max_N_neighborhoods_Batch),
                            title = paste0("Look at the plot generated, Then decide the appropiate number of Neighborhoods"))

        print("Performing Batched K means algorithm")
        #Calculate the desired number of clusters with batch k menas
        Batch_k_means <- ClusterR::MiniBatchKmeans(Neighbor_patterns_scaled,
                                                   clusters = as.double(N_Phenotypes),
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

        #Assign the cluster to each observation of MARKER
        pr_mb <- predict(object = Batch_k_means, fuzzy = F, newdata = Neighbor_patterns_scaled)
        pr_mb <- as_tibble(pr_mb)
        names(pr_mb) <- "Neighborhood_assignment"

        #Generate the data phenotypes tibble
        Neighbor_patterns <-dplyr::bind_cols(Neighbor_patterns, pr_mb)
      }
    }
    #Define what to do if GMM is required
    if(Strategy == "GMM"){
      print("Starting Cluster number stimation process")
      #First we calculate a metric to decide the number of total phenotypes
      #Run the specified test using each of the number of clusters
      Optimal <- try(ClusterR::Optimal_Clusters_GMM(Neighbor_patterns_scaled,
                                                    criterion = Quality_metric,
                                                    max_clusters = Max_N_neighborhoods_GMM,
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
        N_Phenotypes<- menu(choices = as.character(1:Max_N_neighborhoods_GMM),
                            title = paste0("Look at the plot generated, Then decide the appropiate number of Phenotypes"))

        print("Calculating Gaussian Mixed Model")
        #Calculate the desired number of clusters with batch k menas
        GMM_model <- ClusterR::GMM(Neighbor_patterns_scaled,
                                   gaussian_comps = as.double(N_Phenotypes),
                                   dist_mode = GMM_Distance,
                                   seed_mode = "random_subset",
                                   km_iter = Max_iterations_km,
                                   em_iter = Max_iterations_em,
                                   verbose = TRUE,
                                   var_floor = 1e-10,
                                   full_covariance_matrices = FALSE
        )

        #Assign the cluster to each observation of MARKER
        pr_mb <- predict(object = GMM_model, fuzzy = F, newdata = Neighbor_patterns_scaled)
        pr_mb <- as_tibble(pr_mb)
        names(pr_mb) <- "Neighborhood_assignment"

        #Generate the data phenotypes tibble
        Neighbor_patterns <-dplyr::bind_cols(Neighbor_patterns, pr_mb)
      }
    }
    #Define what to do if CLARA clustering is required
    if(Strategy == "CLARA_clustering"){
      print("Starting Cluster number stimation process")
      #First we calculate a metric to decide the number of total phenotypes
      Optimal <-  try(ClusterR::Optimal_Clusters_Medoids(Neighbor_patterns_scaled,
                                                         max_clusters = Max_N_neighborhoods_CLARA,
                                                         distance_metric = Distance_CLARA,
                                                         criterion = "silhouette" ,
                                                         clara_samples = Samples_CLARA,
                                                         clara_sample_size = Sample_per_CLARA,
                                                         swap_phase = F,
                                                         threads = N_cores,
                                                         verbose = T,
                                                         plot_clusters = T
      )
      )
      #Test if optimal number of clusters returned an error
      if(berryFunctions::is.error(Optimal)) {
        stop("Could not calculate best cluster number for the data provided. Please try another strategy")
      }
      #Continue if everything OK
      else{
        #Make the user decide the total number of clusters to be used in the final analysis
        N_Phenotypes<- menu(choices = as.character(1:Max_N_neighborhoods_CLARA),
                            title = paste0("Based on the plots generated and you previous choice, decide the appropiate number of final Neighborhoods"))

        print("Performing CLARA (Clustering Large Applications)")
        CLARA_Clustering <- ClusterR::Clara_Medoids(Neighbor_patterns_scaled,
                                                    clusters = as.double(N_Phenotypes),
                                                    samples = Samples_CLARA,
                                                    sample_size = Sample_per_CLARA,
                                                    distance_metric = Distance_CLARA,
                                                    threads = N_cores,
                                                    swap_phase = F,
                                                    fuzzy = FALSE,
                                                    verbose = T,
                                                    seed = 21
        )
        #Assign the cluster to each observation of MARKER
        pr_mb <- predict(object = CLARA_Clustering, fuzzy = F, newdata = Neighbor_patterns_scaled)
        pr_mb <- as_tibble(pr_mb)
        names(pr_mb) <- "Neighborhood_assignment"

        #Generate the data phenotypes tibble
        Neighbor_patterns <-dplyr::bind_cols(Neighbor_patterns, pr_mb)
      }
    }

    print("Generating plots")


    if(Perform_Dimension_reduction){
      #plot dimension reduction according to the number of cells
      if(nrow(DATA_Reduction) <= 100000){
        try(plot(
          DATA_Reduction %>%dplyr::mutate(Neighborhood_assignment = Neighbor_patterns[["Neighborhood_assignment"]]) %>%
            ggplot(aes(x = DIMENSION_1, y = DIMENSION_2, color = as.factor(Neighborhood_assignment))) +
            geom_point(size = 2, alpha = 0.95) +
            cowplot::theme_cowplot() +
            scale_color_manual("Neighborhood_assignment", values = unname(pals::polychrome(length(unique(Neighbor_patterns$Neighborhood_assignment)))))
        )
        )
      }
      if(nrow(DATA_Reduction) > 100000){
        message(">100K observations to generate plots. A random subset containing 10% of the dataset will be selected for Dimension reduction plots")
        try(plot(
          DATA_Reduction %>%dplyr::mutate(Neighborhood_assignment = Neighbor_patterns[["Neighborhood_assignment"]]) %>%
            sample_n(size = 100000) %>%
            ggplot(aes(x = DIMENSION_1, y = DIMENSION_2, color = as.factor(Neighborhood_assignment))) +
            geom_point(size = 2, alpha = 0.95) +
            cowplot::theme_cowplot() +
            scale_color_manual("Neighborhood_assignment", values = unname(pals::polychrome(length(unique(Neighbor_patterns$Neighborhood_assignment)))))
        )
        )
      }
    }

    #Visualize the neighbor composition data for each neighborhood
    plot(Neighbor_patterns %>% tidyr::pivot_longer(cols = -Neighborhood_assignment) %>%
           ggplot(aes(x = as.factor(Neighborhood_assignment), y = value)) +
           geom_violin(aes(color = name, fill = name), alpha=0.3, position=position_dodge(width=0.5)) +
           stat_summary(aes(color = name),
                        fun = median, geom = "crossbar", width = 0.4, linetype = 1, linewidth = 0.5,
                        position = position_dodge(width = 0.5)) +
           cowplot::theme_cowplot()+
           scale_x_discrete("Neighborhood")+
           scale_y_continuous("Neighbors in tha hood")+
           scale_color_discrete("") +
           scale_fill_discrete(""))

    #Visualize the heatmap of mean by neighborhood
    Mean_tibble <- Neighbor_patterns %>% group_by(Neighborhood_assignment) %>% dplyr::summarize_all(mean) %>%dplyr::ungroup() #Obtain mean tibble
    Mean_matrix <- as.matrix(Mean_tibble[-1] %>% scale()) #Scale it and transform it into a  mtrix
    row.names(Mean_matrix) <- Mean_tibble[[1]]

    plot(ComplexHeatmap::Heatmap(Mean_matrix,
                                 name = "Scaled")
    )

    #Generate the final results
    Final_result <-dplyr::bind_cols(DATA_Neighbors, Neighbor_patterns["Neighborhood_assignment"])

    #Print the number of observations by neighborhoods
    print(Final_result %>% dplyr::count(Neighborhood_assignment) %>% dplyr::arrange(desc(n)))
    #Turn into a factor
    Final_result$Neighborhood_assignment <- as.factor(Final_result$Neighborhood_assignment)

    if(Perform_Dimension_reduction) return(list(DATA = Final_result,
                                                Dimension_reduction = DATA_Reduction)
    )
    else return(Final_result)

  }
