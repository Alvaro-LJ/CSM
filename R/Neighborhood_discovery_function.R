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
#' @param Stop_at_preprocessing A logical value indicating if the function should stop after Data pre-processing and return the interim results (see details).
#' @param Pre_processed_data (OPTIONAL) If a pre-processing object is available, it can be provided to skip data pre-processing steps (see details).
#'
#' @param Strategy The clustering strategy. One of the following Consensus_Clustering, SOM, Graph_Based, K_Means_Meta_clustering, Batch_K_means, GMM or CLARA_clustering (see details).
#' @param Force_N_Neighborhoods A logical value indicating if the number of neighborhoods indicated should be used directly to compute the clusters. If TRUE, no estimating process is conducted. Applicable for SOM, Batch_K_means, GMM or CLARA_clustering.
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
#' @param Graph_type If Strategy is Graph_Based: A character value indicating the type of graph to build. Either "Complete" or "SNN" (nearest neighbors).
#' @param Nearest_neighbors_for_graph If strategy is Graph_Based: The number of closest neighbors to calculate the graph.
#' @param Graph_Method If strategy is Graph_Based: One of Louvain, Leiden, Greedy, WalkTrap, Spinglass, Leading_Eigen or Edge_Betweenness.
#' @param Graph_Distance_method If strategy is Graph_Based and Graph_type is complete, a character value indicating the distance method used to build the graph. One of the following: "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski".
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
#' If the dataset is really large, the pre-processing steps (dimension reduction) and clustering can take a while. This can make the tuning of
#' clustering parameters tedious (each try takes a lot of time). To tackle this issue, pre-processing and clustering steps can be split into two independent steps.
#' By setting the Stop_at_preprocessing argument to TRUE, the function will perform only the required pre-processing and will return the result without performing the
#' clustering. Once the user is satisfied with the pre-processing results, the object can be provided to the function to execute the clustering process.
#'
#' Consensus clustering is performed using the ConsensusClusterPlus::ConsensusClusterPlus function.
#'
#' Self Organizing Maps clustering is performed using the FlowSOM::FlowSOM function.
#'
#' For graph based clustering Nearest neighbors graphs (SNNG) are built using bluster::makeSNNGraph and clustered using functions included in the igraph package.
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
#' @examples
#' \dontrun{
#'#Calculate the closest neighbor matrix----------------------------------------
#'DATA_Closest_Neighbors <-
#'Tailored_Closest_neighbor_calculator(
#'    N_cores = 1,
#'    DATA = CSM_Phenotypecell_test,
#'    Strategy = "Distance",
#'    Include_COO_in_neighborhood = TRUE,
#'    Max_dist_allowed = 50,
#'    Cell_Of_Origin = c("TUMOR", "CD8_GZMBneg", "CD8_GZMBpos"),
#'    Target_Cell = c("TUMOR", "CD8_GZMBneg", "CD8_GZMBpos", "OTHER")
#')
#'
#'#Cluster the result to obtain neighborhoods-----------------------------------
#'Neighborhood_discovery_function(
#'    DATA = DATA_Closest_Neighbors$Absolute_count,
#'    Allowed_max_Dist = 51,
#'    Allowed_avg_Dist = 51,
#'    Allowed_median_Dist = 51,
#'
#'    Perform_Dimension_reduction = TRUE,
#'    Dimension_reduction = "UMAP",
#'    Dimension_reduction_prop = 1,
#'    Cluster_on_Reduced = TRUE,
#'
#'    Strategy = "GMM",
#'    Quality_metric = "AIC",
#'    Max_N_neighborhoods_GMM = 5,
#'    Max_iterations_km = 10,
#'    Max_iterations_em = 10,
#'    GMM_Distance = "eucl_dist"
#')
#' }
#'
#' @export

Neighborhood_discovery_function2 <-
  function(DATA,

           Allowed_max_Dist,
           Allowed_avg_Dist,
           Allowed_median_Dist,

           #Dimension reduction
           Perform_Dimension_reduction = FALSE,
           Dimension_reduction = NULL,
           Dimension_reduction_prop = NULL,
           Cluster_on_Reduced = NULL,

           #Perform only Pre_processing
           Stop_at_preprocessing = FALSE,
           Pre_processed_data = NULL,

           #Strategy for clustering
           Strategy,
           Force_N_Neighborhoods = FALSE,

           #Parameters for Consensus Clustering
           Max_N_neighborhoods = NULL,
           Consensus_reps = NULL,
           Consensus_p_Items = NULL,
           Consensus_Cluster_Alg = NULL,
           Consensus_Distance = NULL,
           Consensus_Name = NULL,

           #Parameters for Self-Organizing Maps
           Max_SOM_neighborhoods = NULL, 

           #Parameters for Graph methods
           Graph_type = NULL,
           Graph_Method = NULL,
           Nearest_neighbors_for_graph = NULL,
           Graph_Resolution = NULL,
           Graph_Distance_method = NULL,
           N_steps = NULL,

           #Parameters for K means Meta Clustering
           N_K_centroids = NULL, 
           Max_N_neighborhoods_Meta = NULL, 
           Consensus_reps_Meta = NULL, 
           Consensus_p_Items_Meta = NULL, 
           Consensus_Name_Meta = NULL, 

           #Parameters for Batched K means
           Batch_size = NULL, 
           Max_N_neighborhoods_Batch = NULL,
           N_initiations = NULL, 
           Max_iterations = NULL, 

           #Parameters for Gaussian Mixture Model
           Quality_metric = NULL,
           Max_N_neighborhoods_GMM = NULL, 
           Max_iterations_km = NULL,
           Max_iterations_em = NULL, 
           GMM_Distance = NULL, 

           #Parameters for CLARA clustering
           Samples_CLARA = NULL, 
           Sample_per_CLARA = NULL, 
           Max_N_neighborhoods_CLARA = NULL, 
           Distance_CLARA = NULL, 
           N_cores = NULL 
  ) {
    ##################################GENERAL ARGUMENT CHECK######################################

    #Always check the strategy
    if(!Strategy %in% c("Consensus_Clustering", "SOM", "Graph_Based", "K_Means_Meta_clustering", "Batch_K_means", "GMM", "CLARA_clustering")) stop("Strategy must be one of the following Consensus_Clustering, SOM, Graph_Based, K_Means_Meta_clustering, Batch_K_means, GMM, CLARA_clustering")
    #Check force neighborhoods
    if(!all(length(Force_N_Neighborhoods) == 1, is.logical(Force_N_Neighborhoods))) stop("Force_N_Neighborhoods must be a logical value")
    if(Force_N_Neighborhoods){
      if(!Strategy %in% c("SOM", "Batch_K_means", "GMM", "CLARA_clustering")) message("Force_N_Neighborhoods cannot be used with current strategy, argument will be ignored")
    }
    
    #If NO Pre-processed data provided check if Stop_at_pre_processing is logical and other pre-processing variables
    if(is.null(Pre_processed_data)){
      if(!all(c("min_DIST", "max_DIST", "avg_DIST", "median_DIST") %in% names(DATA))) stop("DATA must be obtained from the Tailored_Closest_neighbor_calculator functions")
      if(!is.logical(Stop_at_preprocessing)) stop("Stop_at_preprocessing should be a logical value")
      if(!all(c(is.numeric(c(Allowed_max_Dist, Allowed_avg_Dist, Allowed_median_Dist)),
                c(Allowed_max_Dist, Allowed_avg_Dist, Allowed_median_Dist) > 0))
      ) stop("Allowed distances must be numeric and > 0")
      if(!is.logical(Perform_Dimension_reduction)) stop("Perform_Dimension_reduction must be a logical value")
      if(Perform_Dimension_reduction){
        if(!Dimension_reduction %in% c("UMAP", "TSNE", "PCA")) stop("Dimension_reduction must be one of the following: UMAP, TSNE, PCA")
        if(!all(is.numeric(Dimension_reduction_prop), Dimension_reduction_prop > 0, Dimension_reduction_prop <= 1)) stop("Dimension_reduction_prop must be a numeric value between 0 and 1")
      }
      if(!is.logical(Cluster_on_Reduced)) stop("Cluster_on_Reduced must be a logical value")
      if(Cluster_on_Reduced){
        if(!Perform_Dimension_reduction) stop("If Clustering needs to be performed on Dimension reduced data please set Perform_Dimension_reduction to TRUE")
      }
    }

    ##################################SUGGESTED PACKAGE AND SPECIFIC ARGUMENT CHECK######################################

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
    #If pre-processed data is not provided, check suggested packages
    if(is.null(Pre_processed_data)){
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

    ##################################DATA PRE-PROCESSING######################################

    #If pre-processed data is not provided, perform data pre-processing
    if(is.null(Pre_processed_data)){
      #Start executing the code
      print("Preparing data for analysis")
      
      #Import closest neighbor data
      DATA_Neighbors <- DATA

      #Filter out cells not whitin a neighborhood
      DATA_Neighbors <- DATA_Neighbors %>% dplyr::filter(max_DIST <= Allowed_max_Dist, avg_DIST <= Allowed_avg_Dist, median_DIST <= Allowed_median_Dist)

      #We select only neighbor columns and the 4 main columns
      Neighbor_patterns <- DATA_Neighbors %>% dplyr::select(-all_of("Phenotype"), #Remove phenotype
                                                            -contains("N_neighbors"), #Remove N_neighbors
                                                            -c((ncol(DATA_Neighbors)-4):(ncol(DATA_Neighbors))) #The last 4 columns containing distances
                                                            )


      #Perform dimension reduction if required
      if(Perform_Dimension_reduction){
        
        #Provide unscaled with the 4 main columns (will be revomed by function)
        DATA_Reduction <- 
          CSM_Dimension_reduction_function(
            Original_data = Neighbor_patterns,
            Dimension_reduction_strategy = Dimension_reduction,
            Dimension_reduction_prop = Dimension_reduction_prop
          )
      }

      #GENERATE THE MARKERS tibble to CLUSTER ACCORDING TO CLUSTER ON REDUCE
      if(Cluster_on_Reduced) MARKERS <- DATA_Reduction[-which(names(DATA_Reduction) == "Cell_no")]
      if(!Cluster_on_Reduced) MARKERS <- Neighbor_patterns %>% dplyr::select(-c(1:4))

      #If Stop_at_preprocessing is true return the interim results
      if(Stop_at_preprocessing){
        #Obtain the arguments
        Pre_processing_argument <- list(
          Perform_Dimension_reduction = Perform_Dimension_reduction,
          Dimension_reduction = Dimension_reduction,
          Dimension_reduction_prop = Dimension_reduction_prop,
          Cluster_on_Reduced = Cluster_on_Reduced
        )

        #No dim reduction required
        if(!Perform_Dimension_reduction){
          #Return the list
          return(list(Pre_processing_argument = Pre_processing_argument,
                      DATA_Neighbors = DATA_Neighbors,
                      Neighbor_patterns = Neighbor_patterns,
                      MARKERS = MARKERS
          ))
        }
        #If dim reduction required
        if(Perform_Dimension_reduction){
          #Print the graph according to results
          if(nrow(DATA_Reduction) <= 100000){
            print("Generating Dimension Reduction Plot")
            plot(
              DATA_Reduction %>% ggplot(aes(x = DIMENSION_1, y = DIMENSION_2)) +
                geom_point(size = 2, alpha = 0.95) +
                cowplot::theme_cowplot()
            )
          }
          if(nrow(DATA_Reduction) > 100000){
            print("More than 100K observations. Generating Dimension Reduction Plot on a random subset of the data")
            plot(
              DATA_Reduction %>% dplyr::slice_sample(prop = 0.1) %>% ggplot(aes(x = DIMENSION_1, y = DIMENSION_2)) +
                geom_point(size = 2, alpha = 0.95) +
                cowplot::theme_cowplot()
            )
          }
          #Return the list
          return(list(Pre_processing_argument = Pre_processing_argument,
                      DATA_Neighbors = DATA_Neighbors,
                      Neighbor_patterns = Neighbor_patterns,
                      MARKERS = MARKERS,
                      DATA_Reduction = DATA_Reduction
          ))
        }
      }
    }

    #If Pre-processed data provided obtain the datasets from the object provided
    if(!is.null(Pre_processed_data)){
      message("Pre_processed_data provided. Pre-processing related arguments will be ignored.")
      if(names(Pre_processed_data)[1] != "Pre_processing_argument") stop("Pre_processing_argument not found in Pre_processed_data object provided")

      Perform_Dimension_reduction <- Pre_processed_data[["Pre_processing_argument"]][["Perform_Dimension_reduction"]]

      #If no dimension reduction required
      if(!Perform_Dimension_reduction){
        DATA_Neighbors <- Pre_processed_data[["DATA_Neighbors"]]
        Neighbor_patterns <- Pre_processed_data[["Neighbor_patterns"]]
        MARKERS <- Pre_processed_data[["MARKERS"]]

      }
      #If dimension reduction is required
      if(Perform_Dimension_reduction){
        DATA_Neighbors <- Pre_processed_data[["DATA_Neighbors"]]
        Neighbor_patterns <- Pre_processed_data[["Neighbor_patterns"]]
        MARKERS <- Pre_processed_data[["MARKERS"]]
        DATA_Reduction <- Pre_processed_data[["DATA_Reduction"]]
      }

    }

    ##################################CLUSTERING######################################

    #Execute the clustering function
    Neighbor_patterns <- 
      CSM_Clustering_function(
        Original_data = Neighbor_patterns,
        MARKERS = MARKERS,
        
        Strategy = Strategy,
        Force_N_Clusters = Force_N_Neighborhoods,
        
        Max_N_clusters_Consensus = Max_N_neighborhoods, 
        Consensus_reps = Consensus_reps, 
        Consensus_p_Items = Consensus_p_Items, 
        Consensus_Cluster_Alg = Consensus_Cluster_Alg,
        Consensus_Distance = Consensus_Distance, 
        Consensus_Name = Consensus_Name, 
        
        Max_SOM_clusters = Max_SOM_neighborhoods, 
        
        Graph_type = Graph_type,
        Graph_Distance_method = Graph_Distance_method,
        Nearest_neighbors_for_graph = Nearest_neighbors_for_graph, 
        Graph_Method = Graph_Method,
        Graph_Resolution = Graph_Resolution, 
        N_steps = N_steps, 
        
        N_K_centroids = N_K_centroids, 
        Max_N_clusters_Meta = Max_N_neighborhoods_Meta, 
        Consensus_reps_Meta = Consensus_reps_Meta, 
        Consensus_p_Items_Meta = Consensus_p_Items_Meta,
        Consensus_Name_Meta = Consensus_Name_Meta, 
        
        Batch_size = Batch_size,
        Max_N_clusters_Batch = Max_N_neighborhoods_Batch, 
        N_initiations = N_initiations, 
        Max_iterations = Max_iterations, 
        
        Quality_metric = Quality_metric, 
        Max_N_clusters_GMM = Max_N_neighborhoods_GMM, 
        Max_iterations_km = Max_iterations_km, 
        Max_iterations_em = Max_iterations_em,
        GMM_Distance = GMM_Distance, 
        
        Samples_CLARA = Samples_CLARA, 
        Sample_per_CLARA = Sample_per_CLARA, 
        Max_N_clusters_CLARA = Max_N_neighborhoods_CLARA, 
        Distance_CLARA = Distance_CLARA,
        N_cores = N_cores 
      )
    
    Neighbor_patterns <- Neighbor_patterns %>% dplyr::rename("Neighborhood_assignment" = "Cluster")

    ##################################RESULT PLOTTING AND FUNCTION EXIT######################################

    print("Generating plots")


    if(Perform_Dimension_reduction){
      #plot dimension reduction according to the number of cells
      if(nrow(DATA_Reduction) <= 100000){
        try(plot(
          DATA_Reduction %>% dplyr::mutate(Neighborhood_assignment = Neighbor_patterns[["Neighborhood_assignment"]]) %>%
            ggplot(aes(x = DIMENSION_1, y = DIMENSION_2, color = as.factor(Neighborhood_assignment))) +
            geom_point(size = 2, alpha = 0.95) +
            cowplot::theme_cowplot() +
            scale_color_manual("Neighborhood_assignment", values = unname(pals::polychrome(length(unique(Neighbor_patterns$Neighborhood_assignment)))))
        )
        )
      }
      if(nrow(DATA_Reduction) > 100000){
        message(">100K observations to generate plots. A random subset containing 10% of the dataset will be selected for Dimension reduction plots.")
        try(plot(
          DATA_Reduction %>% dplyr::mutate(Neighborhood_assignment = Neighbor_patterns[["Neighborhood_assignment"]]) %>%
            dplyr::slice_sample(n = 100000) %>%
            ggplot(aes(x = DIMENSION_1, y = DIMENSION_2, color = as.factor(Neighborhood_assignment))) +
            geom_point(size = 2, alpha = 0.95) +
            cowplot::theme_cowplot() +
            scale_color_manual("Neighborhood_assignment", values = unname(pals::polychrome(length(unique(Neighbor_patterns$Neighborhood_assignment)))))
        )
        )
      }
    }

    #Visualize the neighbor composition data for each neighborhood
    try(
      plot(Neighbor_patterns %>% dplyr::select(-c(1:4)) %>% tidyr::pivot_longer(cols = -Neighborhood_assignment) %>%
             ggplot(aes(x = as.factor(Neighborhood_assignment), y = value)) +
             geom_violin(aes(color = name, fill = name), alpha = 0.3, position=position_dodge(width=0.5)) +
             stat_summary(aes(color = name),
                          fun = median, geom = "crossbar", width = 0.4, linetype = 1, linewidth = 0.5,
                          position = position_dodge(width = 0.5)) +
             cowplot::theme_cowplot()+
             scale_x_discrete("Neighborhood")+
             scale_y_continuous("Neighbors in tha hood")+
             scale_color_discrete("") +
             scale_fill_discrete(""))
    )

    #Visualize the heatmap of mean by neighborhood
    Mean_tibble <- Neighbor_patterns %>% dplyr::select(-c(1:4)) %>% group_by(Neighborhood_assignment) %>% dplyr::summarize_all(mean) %>%dplyr::ungroup() #Obtain mean tibble
    Mean_matrix <- as.matrix(Mean_tibble[-1] %>% scale()) #Scale it and transform it into a  mtrix
    row.names(Mean_matrix) <- Mean_tibble[[1]]

    try(
      plot(ComplexHeatmap::Heatmap(Mean_matrix,
                                   name = "Scaled")
      )
    )


    #Generate the final results
    Final_result <- dplyr::bind_cols(DATA_Neighbors, Neighbor_patterns["Neighborhood_assignment"])

    #Print the number of observations by neighborhoods
    print(Final_result %>% dplyr::count(Neighborhood_assignment) %>% dplyr::arrange(desc(n)))
    #Turn into a factor
    Final_result$Neighborhood_assignment <- as.factor(Final_result$Neighborhood_assignment)

    if(Perform_Dimension_reduction) return(list(DATA = Final_result,
                                                Dimension_reduction = DATA_Reduction)
    )
    else return(Final_result)
  }

