#' Performs phenotype assignment using clustering approaches
#'
#' `Clustering_Phenotyper()` assigns cell phenotypes based on clustering. Phenotype names can be then modidified using [DATA_Phenotype_renamer()]
#'
#' @param DATA A dataframe or tibble containing cell feature data.
#' @param Strategy One of the following Consensus_Clustering, SOM, Graph_Based, K_Means_Meta_clustering, Batch_K_means, GMM or CLARA_clustering (see details).
#'
#' @param Apply_Denoise A logical value. Specify if a denoising filtering is required before clustering (see details).
#' @param Denoising Denoising strategy. One of the following: Quantile, Standard_Deviation, Threshold, Otsu or DimRed_DBscan.
#' @param Percentile A numeric value indicating the percentile for quantile threshold. Cells below percentile for all features will be considered to be noise.
#' @param N_Standard_Deviations A numeric value indicating the number of standard deviations from mean for Standard_Deviation method. Cells below SD for all features will be considered to be noise.
#' @param Selected_threshold A numeric value indicating the threshold to be used for the Threshold method. Cells below the threshold for all features will be considered to be noise.
#' @param Min_cell_no An integer value for the DBscan method. Minimum cell number in distance to consider a cell to be clustered.
#' @param Distance_radius A numeric value for the DBscan method. Distance to be sampled.
#'
#' @param Perform_Dimension_reduction Logical value. Should Dimension Reduction be performed (see details)
#' @param Dimension_reduction Dimension reduction method. One of the following: PCA, TSNE, UMAP
#' @param Dimension_reduction_prop A numeric value between 0 and 1 to indicate the percentage of the cells to be used in dimension computation (applicable for TSNE and UMAP)
#' @param Cluster_on_Reduced A logical value indicating if clustering should be performed on new dimensions
#'
#' @param Max_N_phenotypes If Strategy is Consensus_Clustering: Number of maximum phenotypes that can be identified.
#' @param Consensus_reps If Strategy is Consensus_Clustering: Number of iterations to converge.
#' @param Consensus_p_Items If Strategy is Consensus_Clustering: Percentage of cells that you desire to sample in each iteration.
#' @param Consensus_Cluster_Alg If Strategy is Consensus_Clustering: Clustering algorithm to be used (’hc’ hierarchical (hclust), ’pam’ for paritioning around medoids, ’km’ for k-means).
#' @param Consensus_Distance If Strategy is Consensus_Clustering: Distance metric to be used (pearson(1 - Pearson correlation), spearman(1 - Spearman correlation), euclidean, binary, maximum, canberra, minkowski.
#' @param Consensus_Name If Strategy is Consensus_Clustering: Name of the folder that is going to be created in order to place the resulting graphs.
#'
#' @param Max_SOM_phenotypes If Strategy is SOM: umber of maximum phenotypes that can be identified.
#'
#' @param Nearest_neighbors_for_graph If strategy is Graph_Based: The number of closest neighbors to calculate the graph.
#' @param Graph_Method If strategy is Graph_Based: One of Louvain, Leiden, Greedy, WalkTrap, Spinglass, Leading_Eigen or Edge_Betweenness.
#' @param Graph_Resolution If strategy is Graph_Based: Used for Louvain and Leiden. 1 is default. The smaller the value, the larger the clusters will be.
#' @param N_steps If strategy is Graph_Based: Number of steps given in the WalkTrap algorithm.
#'
#' @param N_K_centroids If strategy is K_Means_Meta_clustering: Number of centroids to perform K means.
#' @param Max_N_phenotypes_Meta If strategy is K_Means_Meta_clustering: Number of maximum phenotypes that can be identified.
#' @param Consensus_reps_Meta If strategy is K_Means_Meta_clustering: Number of iterations to converge.
#' @param Consensus_p_Items_Meta If strategy is K_Means_Meta_clustering: Percentage of cells that you desire to sample in each iteration.
#' @param Consensus_Name_Meta If strategy is K_Means_Meta_clustering: Name of the folder that is going to be created in order to place the resulting graphs.
#'
#' @param Batch_size If strategy is Batch_K_means: Number of cells to be included in each random batch.
#' @param Max_N_phenotypes_Batch If strategy is Batch_K_means: Number of maximum phenotypes that can be identified.
#' @param N_initiations If strategy is Batch_K_means: Number of times the algorithm is going to be tried to find the best clustering result.
#' @param Max_iterations If strategy is Batch_K_means: Max number of iterations in each try.
#'
#' @param Quality_metric If strategy is GMM:T he quality measure used to test the number of clusters ("AIC" or "BIC").
#' @param Max_N_phenotypes_GMM If strategy is GMM: Number of maximum phenotypes that can be identified.
#' @param Max_iterations_km If strategy is GMM: Number of max iterations in the K means clustering performed.
#' @param Max_iterations_em If strategy is GMM: Number of max iterations in the Expectation Maximization algorithm.
#' @param GMM_Distance If strategy is GMM: Distance metric used in the model ("eucl_dist" or "maha_dist").
#'
#' @param Samples_CLARA If strategy is CLARA_clustering: Number of samples the CLARA algorithm is going to use to be calculated.
#' @param Sample_per_CLARA If strategy is CLARA_clustering: Percentage (from 0 to 1) of the total cells that are going to be allocated to each sample.
#' @param Max_N_phenotypes_CLARA If strategy is CLARA_clustering: Number of maximum phenotypes that can be identified.
#' @param Distance_CLARA If strategy is CLARA_clustering: Distance metric used in the model (euclidean, manhattan, chebyshev, canberra, braycurtis, pearson_correlation, simple_matching_coefficient, minkowski, hamming, jaccard_coefficient, Rao_coefficient, mahalanobis, cosine)
#' @param N_cores If strategy is CLARA_clustering: Number of cores to parallelize your computation
#'
#' @details
#' De-noising process does not remove cells from the final output. It rather assigns noise cells to a single phenotype. Otsu thresholding and DBSCAN based denoising are based on EBImage::otsu and dbscan::dbscan functions, respectively.
#'
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
#' @seealso [DATA_Phenotype_renamer()], [ReClustering_function()], [Consensus_phenotype_assigner()], [Concordance_calculator()], [Confusion_matrix_plotter()],
#' [Phenotyping_evaluator_shiny_app_launcher()]
#'
#' @returns Returns a tibble with cell features and a column named 'Phenotype' containing cell labels.
#' If dimension reduction has been performed returns a list with the cell feature dataset as above and a tibble containing dimension reduction coordinates.
#'
#' @export

Clustering_Phenotyper <-
  function(DATA = NULL,
           Strategy = NULL,

           #Denoising parameters
           Apply_Denoise = NULL, #Specify if a denoising filtering is required before clustering
           Denoising = NULL, #Select denoising strategy from: Quantile, Standard_Deviation, Threshold, Otsu or DimRed_DBscan
           Percentile = NULL, #Select the adequate percentile for quantile threshold
           N_Standard_Deviations = NULL, #Select the number of standard deviations from mean for Standard_Deviation method
           Selected_threshold = NULL, #Select the absolute threshold for the Threshold method
           Min_cell_no = NULL, #Parameter for DBscan
           Distance_radius = NULL, #Parameter for DBscan

           #Dimension reduction
           Perform_Dimension_reduction = NULL,
           Dimension_reduction = NULL,
           Dimension_reduction_prop = NULL,
           Cluster_on_Reduced = NULL,

           #Parameters for Consensus Clustering
           Max_N_phenotypes = NULL, #Number of maximum neighborhods that you desire to find
           Consensus_reps = NULL, #Number of iterations of the algorithm to try to converge
           Consensus_p_Items = NULL, #Percentage of the closest neighbor patterns that you desire to sample in each iteration
           Consensus_Cluster_Alg = NULL, #Clustering algorithm to be used (’hc’ hierarchical (hclust), ’pam’ for paritioning around medoids, ’km’ for k-means )
           Consensus_Distance = NULL, #Distance metric to be used (pearson(1 - Pearson correlation), spearman(1 - Spearman correlation), euclidean, binary, maximum, canberra, minkowski
           Consensus_Name = NULL, #Name of the folder that is going to be created in order to place the resulting graphs

           #Parameters for Self-Organizing Maps
           Max_SOM_phenotypes = NULL, #Maximum number of clusters (phenotypes) to try in the algorithm

           #Parameters for Graph-Based approaches
           Nearest_neighbors_for_graph = NULL, #Specify the number of closest neighbors to calculate the graph
           Graph_Method = NULL, #Specify the clustering method
           Graph_Resolution = NULL, #Specify the graph resolution
           N_steps = NULL, #Number of steps given in the WalkTrap algorithm

           #Parameters for K means Meta Clustering
           N_K_centroids = NULL, #Number of centroids to perform K means
           Max_N_phenotypes_Meta = NULL, #Number of maximum clusters (phenotypes) that you desire to find
           Consensus_reps_Meta = NULL, #Number of iterations of the algorithm to try to converge
           Consensus_p_Items_Meta = NULL, #Percentage of cells that you desire to sample in each iteration
           Consensus_Name_Meta = NULL, #Name of the folder that is going to be created in order to place the resulting graphs

           #Parameters for Batched K means
           Batch_size = NULL, #The number of cells to be included in each random batch
           Max_N_phenotypes_Batch = NULL, #Number of maximum clusters (phenotypes) that you desire to find
           N_initiations = NULL, #Number of times the algorithm is going to be tried to find the best clustering result
           Max_iterations = NULL, #Max number of iterations in each try

           #Parameters for Gaussian Mixture Model
           Quality_metric = NULL, #The quality measure used to test the number of clusters ("AIC" or "BIC")
           Max_N_phenotypes_GMM = NULL, #Number of maximum clusters (phenotypes) that you desire to find
           Max_iterations_km = NULL, #Number of max iterations in the K means clustering performed
           Max_iterations_em = NULL, #Number of max iterations in the Expectation Maximization algorithm
           GMM_Distance = NULL, #Distance metric to use in the model ("eucl_dist" or "maha_dist")

           #Parameters for CLARA clustering
           Samples_CLARA = NULL, #Number of samples the CLARA algorithm is going to use to be calculated
           Sample_per_CLARA = NULL, #Percentage (from 0 to 1) of the total cells that are going to be allocated to each sample
           Max_N_phenotypes_CLARA = NULL, #Number of maximum clusters (phenotypes) that you desire to find
           Distance_CLARA = NULL, #euclidean, manhattan, chebyshev, canberra, braycurtis, pearson_correlation,
           #simple_matching_coefficient, minkowski, hamming, jaccard_coefficient, Rao_coefficient, mahalanobis, cosine
           N_cores = NULL #Number of cores to parallelize your computation
  ) {
    on.exit(gc())

    #Check general arguments
    if(!identical(names(DATA)[1:4], c("Cell_no", "X", "Y", "Subject_Names"))) {
      stop("Please generate an appropiate data object using the Data_arrange_function")
    }
    if(!Strategy %in% c("Consensus_Clustering", "SOM", "Graph_Based", "K_Means_Meta_clustering", "Batch_K_means", "GMM", "CLARA_clustering")){
      stop("Strategy must be one of the following: Consensus_Clustering, SOM, Graph_Based, K_Means_Meta_clustering, Batch_K_means, GMM, CLARA_clustering")
    }
    if(!is.logical(Apply_Denoise)) stop("Apply_Denoise must be a logical value")
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
      Argument_checker <- c(Max_N_phenotypes_OK = (Max_N_phenotypes >= 2 & Max_N_phenotypes%%1 == 0),
                            Consensus_reps_OK = (Consensus_reps >= 1 & Consensus_reps%%1 == 0),
                            Consensus_p_Items_OK = (Consensus_p_Items > 0 & Consensus_p_Items <= 1),
                            Consensus_Cluster_Alg_OK = Consensus_Cluster_Alg %in% c("hc", "pam", "km"),
                            Consensus_Distance_OK = Consensus_Distance %in% c("pearson", "spearman", "euclidean", "binary", "maximum", "canberra", "minkowski"),
                            Consensus_Name_OK = is.character(as.character(Consensus_Name))
      )
      Stop_messages <- c(Max_N_phenotypes_OK = "Max_N_Phenotypes must be an integer value > 1",
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
      if(!(Max_SOM_phenotypes > 1 & Max_SOM_phenotypes%%1 == 0)) stop("Max_SOM_phenotypes must be an integer value > 1")
    }
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
      Argument_checker <- c(Nearest_neighbors_for_graph_OK = (Nearest_neighbors_for_graph >= 1 & Nearest_neighbors_for_graph%%1 == 0),
                            Graph_Method_OK = Graph_Method %in% c("Louvain", "Leiden", "Greedy", "WalkTrap", "Spinglass", "Leading_Eigen", "Edge_Betweenness"),
                            Graph_Resolution_OK = all(is.numeric(Graph_Resolution), Graph_Resolution > 0),
                            N_steps_OK = is.null(N_steps) || (N_steps >=1 & N_steps%%1 == 0)

      )
      Stop_messages <- c(Nearest_neighbors_for_graph = "Nearest_neighbors_for_graph must be an integer value > 0",
                         Graph_Method = "Graph_Method must be one of the following: Louvain, Leiden, Greedy, WalkTrap, Spinglass, Leading_Eigen, Edge_Betweenness",
                         Graph_Resolution = "Graph_Resolution must be a numeric value > 0",
                         N_steps = "N_steps must be a integer value > 0"
      )
      #Check arguments and stop if necessary
      if(!all(Argument_checker)){
        stop(cat(Stop_messages[!Argument_checker],
                 fill = sum(!Argument_checker)))
      }
    }
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
      Argument_checker <- c(N_K_centroids_OK = all(nrow(MARKERS) > N_K_centroids, N_K_centroids%%1 == 0, N_K_centroids > 0),
                            Max_N_phenotypes_Meta_OK = (Max_N_phenotypes_Meta >= 2 & Max_N_phenotypes_Meta%%1 == 0),
                            Consensus_reps_Meta_OK = (Consensus_reps_Meta >= 1 & Consensus_reps_Meta%%1 == 0),
                            Consensus_p_Items_Meta_OK = (Consensus_p_Items_Meta > 0 & Consensus_p_Items_Meta <= 1),
                            Consensus_Name_Meta_OK = is.character(as.character(Consensus_Name_Meta))
      )
      Stop_messages <- c(N_K_centroids_OK = "N_K_centroids must be smaller than the number of cells in DATA and a integer value > 0",
                         Max_N_phenotypes_Meta_OK = "Max_N_phenotypes_Meta must be an integer value > 1",
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
    if(Strategy == "Batch_K_means"){
      #Check suggested packages
      if(!requireNamespace("ClusterR", quietly = FALSE)) stop(
        paste0("ClusterR CRAN package is required to execute the function. Please install using the following code: ",
               expression(install.packages("ClusterR")))
      )
      #Check arguments
      Argument_checker <- c(Max_N_phenotypes_Batch_OK = (Max_N_phenotypes_Batch >= 2 & Max_N_phenotypes_Batch%%1 == 0),
                            N_initiations_OK = (N_initiations >= 1 & N_initiations%%1 == 0),
                            Max_iterations_OK = (Max_iterations%%1 == 0)
      )
      Stop_messages <- c(Max_N_phenotypes_Batch_OK = "Max_N_phenotypes_Batch must be an integer value > 1",
                         N_initiations_OK = "N_initiations must be an integer value > 0",
                         Max_iterations_OK = "Max_iterations must be an integer value > 0"
      )
      #Check arguments and stop if necessary
      if(!all(Argument_checker)){
        stop(cat(Stop_messages[!Argument_checker],
                 fill = sum(!Argument_checker)))
      }
    }
    if(Strategy == "GMM"){
      #Check suggested packages
      if(!requireNamespace("ClusterR", quietly = FALSE)) stop(
        paste0("ClusterR CRAN package is required to execute the function. Please install using the following code: ",
               expression(install.packages("ClusterR")))
      )

      #Check arguments
      Argument_checker <- c(Quality_metric_OK = Quality_metric %in% c("AIC", "BIC"),
                            Max_N_phenotypes_GMM_OK = (Max_N_phenotypes_GMM >= 2 & Max_N_phenotypes_GMM%%1 == 0),
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
    if(Strategy == "CLARA_clustering"){
      #Check suggested packages
      if(!requireNamespace("ClusterR", quietly = FALSE)) stop(
        paste0("ClusterR CRAN package is required to execute the function. Please install using the following code: ",
               expression(install.packages("ClusterR")))
      )

      #Check arguments
      Argument_checker <- c(Samples_CLARA_OK = (Samples_CLARA >= 1 & Samples_CLARA%%1 == 0),
                            Sample_per_CLARA_OK = (Sample_per_CLARA > 0 & Sample_per_CLARA <= 1),
                            Max_N_phenotypes_CLARA_OK = (Max_N_phenotypes_CLARA >= 2 & Max_N_phenotypes_CLARA%%1 == 0),
                            Distance_CLARA_OK = Distance_CLARA %in% c("euclidean", "manhattan", "chebyshev", "canberra", "braycurtis",
                                                                      "pearson_correlation", "simple_matching_coefficient", "minkowski",
                                                                      "hamming", "jaccard_coefficient", "Rao_coefficient", "mahalanobis", "cosine"),
                            N_cores_OK = (N_cores >= 1 & N_cores%%1 == 0)
      )
      Stop_messages <- c(Samples_CLARA_OK = "Samples_CLARA must be an integer value > 0",
                         Sample_per_CLARA_OK = "Sample_per_CLARA must be a numeric value between 0 and 1",
                         Max_N_phenotypes_CLARA_OK = "Max_N_phenotypes_CLARA must be an integer value > 1",
                         Distance_CLARA_OK = "Distance_CLARA must be one of the following: euclidean, manhattan, chebyshev, canberra, braycurtis,
                                                                     pearson_correlation, simple_matching_coefficient, minkowski,
                                                                     hamming, jaccard_coefficient, Rao_coefficient, mahalanobis, cosine",
                         N_cores_OK = "N_cores must be an integer value > 0"
      )
      #Check arguments and stop if necessary
      if(!all(Argument_checker)){
        stop(cat(Stop_messages[!Argument_checker],
                 fill = sum(!Argument_checker)))
      }
    }
    if(Apply_Denoise){
      #check denoising argument is correctly stated
      if(!Denoising %in% c("Quantile", "Standard_Deviation", "Threshold", "Otsu", "DimRed_DBscan")) {
        stop("Denoising should be one of Quantile, Standard_Deviation, Threshold, Otsu, DimRed_DBscan")
      }
      if(Denoising == "DimRed_DBscan"){
        #Requires previous dimension reduction
        if(!Perform_Dimension_reduction) stop("DBscan clustering requires Dimension reduction. Please set Perform_Dimension_reduction to TRUE")
      }

      #Check suggested packages
      if(Denoising == "Otsu"){
        if(!requireNamespace("EBImage", quietly = TRUE)) stop(
          paste0("EBImage Bioconductor package is required to execute the function. Please install using the following code: ",
                 expression({
                   if (!require("BiocManager", quietly = TRUE))
                     install.packages("BiocManager")

                   BiocManager::install("EBImage")
                 })
          )
        )
      }
      if(Denoising == "DimRed_DBscan"){
        if(!requireNamespace("dbscan", quietly = FALSE)) stop(
          paste0("dbscan CRAN package is required to execute the function. Please install using the following code: ",
                 expression(install.packages("dbscan")))
        )
      }
    }
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

    #Perform dimension reduction if required
    if(Perform_Dimension_reduction){
      #First PCA
      if(Dimension_reduction == "PCA"){
        if(Dimension_reduction_prop != 1) stop("PCA must be performed using Dimension_reduction_prop = 1")
        print("Generating PCA projections")
        #Scale and turn into matrix
        DATA_matrix <- DATA %>% dplyr::select(-c(1:4)) %>% scale() %>% as.matrix()
        Result_PCA <- svd::propack.svd(DATA_matrix, neig = 2)$u
        DATA_Reduction <- tibble(Cell_no = DATA$Cell_no, DIMENSION_1 = unlist(Result_PCA[,1]), DIMENSION_2 = unlist(Result_PCA[,2]))
      }

      #Second TSNE
      if(Dimension_reduction == "TSNE"){
        if(Dimension_reduction_prop == 1) {
          print("Generating TSNE projections")
          if(nrow(DATA) > 50000) print("Warning! Data set contains more than 50K observations. tSNE embedding can take a long time")
          #scale and turn into matrix
          DATA_matrix <- DATA %>% dplyr::select(-c(1:4)) %>% scale() %>% as.matrix()
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
          DATA_Reduction <-dplyr::bind_cols(DATA["Cell_no"], DIMENSION_1 = unlist(Result_TSNE[,1]), DIMENSION_2 = unlist(Result_TSNE[,2]))
        }

        if(Dimension_reduction_prop != 1) {
          print("Generating TSNE projections")
          DATA_matrix <- DATA %>% dplyr::group_by(Subject_Names) %>% dplyr::sample_frac(size = Dimension_reduction_prop) %>% dplyr::ungroup() %>%
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
                                     new = DATA  %>% dplyr::select(-c(1:4)) %>% scale() %>% as.matrix(),
                                     old = DATA_matrix)
          DATA_Reduction <-dplyr::bind_cols(DATA["Cell_no"], DIMENSION_1 = unlist(Coords[,1]), DIMENSION_2 = unlist(Coords[,2]))
        }
      }

      #Third UMAP
      if(Dimension_reduction == "UMAP"){
        if(Dimension_reduction_prop == 1) {
          print("Generating UMAP projections")
          if(nrow(DATA) > 50000) print("Warning! Data set contains more than 50K observations. UMAP embedding can take some time")
          #scale and turn into matrix
          DATA_matrix <- DATA %>% dplyr::select(-c(1:4)) %>% scale() %>% as.matrix()
          Result_UMAP <- uwot::tumap(DATA_matrix, n_components = 2L)
          DATA_Reduction <-dplyr::bind_cols(DATA["Cell_no"], DIMENSION_1 = unlist(Result_UMAP[,1]), DIMENSION_2 = unlist(Result_UMAP[,2]))
        }

        if(Dimension_reduction_prop != 1) {
          print("Generating UMAP projections")
          DATA_matrix <- DATA %>% dplyr::group_by(Subject_Names) %>% dplyr::sample_frac(size = Dimension_reduction_prop) %>% dplyr::ungroup() %>%
            dplyr::select(-c(1:4)) %>% scale() %>% as.matrix()
          if(nrow(DATA_matrix) > 50000) print("Warning! Data set contains more than 50K observations. UMAP embedding can take aome time")
          #scale and turn into matrix
          Result_UMAP <- uwot::tumap(DATA_matrix, n_components = 2L, ret_model = TRUE)
          Coords <- uwot::umap_transform(X = DATA  %>% dplyr::select(-c(1:4)) %>% scale() %>% as.matrix(),
                                         model = Result_UMAP)
          DATA_Reduction <-dplyr::bind_cols(DATA["Cell_no"], DIMENSION_1 = unlist(Coords[,1]), DIMENSION_2 = unlist(Coords[,2]))
        }
      }
    }

    #If denoising is required apply required function
    if(Apply_Denoise){

      print("Filtering out noisy cells")
      #Identify each cell in the experiment with a unique ID
      DATA <- DATA %>% dplyr::mutate(Unique_ID = 1:nrow(DATA))
      DATA <- DATA[c(ncol(DATA), 1:(ncol(DATA)-1))]

      #Apply desired filters
      if(Denoising == "Quantile") {
        #Check arguments
        if(Percentile < 0.01 | Percentile > 0.99) stop("Percentile must be between 0.01 and 0.99")

        FILTER <-purrr::map_dfc(DATA[-(1:5)], function(x){
          x <= quantile(x, Percentile)
        })
      }

      else if(Denoising == "Standard_Deviation"){
        #Check arguments
        if(!is.numeric(N_Standard_Deviations)) stop("N_Standard_Deviations must be a numeric value")

        FILTER <-purrr::map_dfc(DATA[-(1:5)], function(x){
          x <= (mean(x) - (N_Standard_Deviations*sd(x)))
        })
      }

      else if(Denoising == "Threshold"){
        #Check arguments
        if(!is.numeric(Selected_threshold)) stop("Selected_threshold must be a numeric value")

        FILTER <-purrr::map_dfc(DATA[-(1:5)], function(x){
          x <= Selected_threshold
        })
      }

      else if(Denoising == "Otsu"){
        FILTER <-purrr::map2_df(.x = DATA[-c(1:5)],
                                .y =purrr::map_dbl(DATA[-c(1:5)], function(z){
                                  EBImage::otsu(array(z, dim = c(1, length(z))), range = c(min(z), max(z)), levels = length(unique(z)))
                                }),
                                function(.x, .y) .x <= .y)
      }

      else if(Denoising == "DimRed_DBscan"){
        #Requires previous dimension reduction
        if(!Perform_Dimension_reduction) stop("DBscan clustering requires Dimension reduction. Please set Perform_Dimension_reduction to TRUE")
        #Check other arguments
        if(!all(is.numeric(Min_cell_no), Min_cell_no%%1 == 0, Min_cell_no > 0)) stop("Min_cell_no must be an integer value > 0")
        if(!all(is.numeric(Distance_radius), Distance_radius > 0)) stop("Distance_radius must be a numeric value > 0")

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
              dplyr::sample_n(size = 100000) %>%
              ggplot(aes(x = Dim_1, y = Dim_2, color = Cluster)) + geom_point(size = 1.5) +
              scale_color_manual("", values = c("black", "grey"))+
              cowplot::theme_cowplot()+
              theme(panel.grid = element_blank())
          )
        }

        DB_OK <- menu(choices = c("Proceed", "Abort"), title = "Are the results of the filtering OK?")
        if(DB_OK == 2) stop("Please refine Distance_radius and Min_cell_no parameters and retry")

        #Generate the NOISE column with a logical vector
        DATA <- DATA %>%dplyr::mutate(NOISE = DB_results$cluster) %>%dplyr::mutate(NOISE = case_when(NOISE == 0 ~ TRUE,
                                                                                                     NOISE != 0 ~ FALSE))
      }

      #For no DBscan methods apply a row wise method to determine which cells are noise (a logical vector)
      if(Denoising != "DimRed_DBscan"){
        #Generate the variable to specify if the cell is noise or not
        DATA <- DATA %>%dplyr::mutate(NOISE = unlist(apply(FILTER, MARGIN = 1, function(x) sum(x) == ncol(FILTER))))
      }

      #Generate two tibbles, one with noisy cells and other (DATA) with the actual cells
      NOISE_VECTOR <- DATA[["NOISE"]] #We generate a NOISE_VECTOR if we require to filter noise cells for Dimension reduction methods
      DATA_NOISE <- DATA %>% dplyr::filter(NOISE) %>%dplyr::mutate(Phenotype = 1)
      DATA <- DATA %>% dplyr::filter(!NOISE)
      MARKERS <- DATA %>% dplyr::select(-Unique_ID, -Cell_no, -X, -Y, -Subject_Names, -NOISE)
    }

    #If no denoising required obtain MARKERS directly from DATA
    else{
      DATA <- DATA
      MARKERS <- DATA %>% dplyr::select(-Cell_no, -X, -Y, -Subject_Names)
    }

    #Generate a specific version of Markers with dimension reduction data for Dimension_SNN
    if(Cluster_on_Reduced){
      #Depending on Denoising Obtain directly from DATA_Reduction or filter first
      if(!Apply_Denoise) MARKERS <- DATA_Reduction[-"Cell_no"]
      if(Apply_Denoise) MARKERS <- DATA_Reduction[!NOISE_VECTOR, ] %>% dplyr::select(-Cell_no)
    }

    #Check that MARKERS is a numeric matrix
    if(!is.numeric(as.matrix(MARKERS))) stop("DATA provided must contain marker intensity values")

    #Continue with clustering strategies
    if(Strategy == "Consensus_Clustering"){
      #Perform consensus clustering
      Phenotype_result <- try(ConsensusClusterPlus::ConsensusClusterPlus(t(as.matrix((MARKERS %>% scale()))),
                                                                         maxK = Max_N_phenotypes,
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
      if(berryFunctions::is.error(Phenotype_result)) {
        stop("Data is too large for Consensus Clustering. Please try another strategy")
      }
      else {
        #Make the user decide the number of neighborhoods according to results
        N_Phenotypes<- menu(choices = as.character(1:Max_N_phenotypes), title = paste0("Check the results at: ", getwd(), ". Then decide the appropiate number of Phenotypes"))
        DATA_Phenotypes <- DATA %>%dplyr::mutate(Phenotype = Phenotype_result[[as.double(N_Phenotypes)]][["consensusClass"]])
      }
    }

    else if(Strategy == "SOM"){
      print("Executing Self Organizing Map algorithm")
      #Transform data into a scaled matrix and perform Self Organizing Map
      SOM_results <- try(FlowSOM::FlowSOM(MARKERS %>% scale() %>% as.matrix(),
                                          scale = F,
                                          colsToUse = c(1:ncol(MARKERS)),
                                          maxMeta = Max_SOM_phenotypes,#To find optimal meta clusters
                                          silent = F,
                                          seed = 21)
      )
      #Test if SOM returned an error
      if(berryFunctions::is.error(SOM_results)) {
        stop("Data is too large for Self-Organizing Maps. Please try another strategy")
      }
      else{
        #Assign phenotypes to each cell
        DATA_Phenotypes <- DATA %>%dplyr::mutate(Phenotype = FlowSOM::GetMetaclusters(SOM_results))
      }

    }

    else if(Strategy == "Graph_Based"){
      print("Start Graph building process")

      #Transform data into a nearest neighbor graph
      SNN_graph <- try(bluster::makeSNNGraph(MARKERS %>% scale() %>% as.matrix(),
                                             k = Nearest_neighbors_for_graph)
      )

      #Test if Graph construction process returned an error
      if(berryFunctions::is.error(SNN_graph)) {
        stop("Data is too large to build a graph. Please try another strategy")
      }
      else{
        print("Performing Graph-based clustering")
        #Go for graph clustering
        #Cluster the graph with louvain or leiden clustering
        if(Graph_Method == "Louvain") {
          DATA_Phenotypes <- DATA %>%dplyr::mutate(Phenotype = igraph::cluster_louvain(SNN_graph,
                                                                                       weights = NULL,
                                                                                       resolution = Graph_Resolution)$membership)

        }

        else if(Graph_Method == "Leiden") {
          DATA_Phenotypes <- DATA %>%dplyr::mutate(Phenotype = igraph::cluster_leiden(SNN_graph,
                                                                                      objective_function = "modularity",
                                                                                      weights = NULL,
                                                                                      resolution = Graph_Resolution,
                                                                                      beta = 0.01,
                                                                                      initial_membership = NULL,
                                                                                      n_iterations = 100,
                                                                                      vertex_weights = NULL)$membership)

        }

        else if(Graph_Method == "Greedy"){
          DATA_Phenotypes <- DATA %>%dplyr::mutate(Phenotype = igraph::cluster_fast_greedy(SNN_graph)$membership)
        }

        else if(Graph_Method == "WalkTrap"){
          DATA_Phenotypes <- DATA %>%dplyr::mutate(Phenotype = igraph::cluster_walktrap(SNN_graph,
                                                                                        steps = N_steps,
                                                                                        membership = T)$membership)
        }

        else if (Graph_Method == "Spinglass") {
          DATA_Phenotypes <- DATA %>%dplyr::mutate(Phenotype = igraph::cluster_spinglass(SNN_graph,
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
          DATA_Phenotypes <- DATA %>%dplyr::mutate(Phenotype = igraph::cluster_leading_eigen(SNN_graph,
                                                                                             membership = T)$membership)
        }

        else if(Graph_Method == "Edge_Betweenness"){
          DATA_Phenotypes <- DATA %>%dplyr::mutate(Phenotype = igraph::cluster_edge_betweenness(SNN_graph,
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
        DATA_filter_Markers <- DATA %>%dplyr::mutate(K_means_Cl = cl$cluster)

        #Prepare data for Meta-Clustering
        #Create a tibble with the K means centroids and the format it for Consensus clustering
        K_medoids <- as_tibble(cl$centers) %>%dplyr::mutate(K_means_Cl = 1:nrow(as_tibble(cl$centers)))
        tK_medoids <- K_medoids %>% dplyr::select(-K_means_Cl) %>% as.matrix() %>% t()

        #Perform Consensus clustering with hierarchical clustering
        print("Perorming Consensus Clustering")
        HC <- try(ConsensusClusterPlus::ConsensusClusterPlus(tK_medoids,
                                                             maxK = Max_N_phenotypes_Meta,
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
        if(berryFunctions::is.error(HC)) {
          stop("Data is too large for Meta Clustering. Please try another strategy or select a smaller N_K_centroids value")
        }
        else {
          #Make the user decide the number of neighborhoods according to results
          N_Phenotypes<- menu(choices = as.character(1:Max_N_phenotypes_Meta), title = paste0("Check the results at: ", getwd(), ". Then decide the appropiate number of Phenotypes"))

          #Bind the final Phenotype to the K medoids tibble
          K_medoids <- K_medoids %>%dplyr::mutate(Phenotype = HC[[as.double(N_Phenotypes)]][["consensusClass"]])
          K_medoids_for_join <- K_medoids %>% dplyr::select(K_means_Cl, Phenotype)

          #Bind The DATA and the K_meoids to obtain the final matrix
          DATA_Phenotypes <-dplyr::left_join(DATA_filter_Markers, K_medoids_for_join, by = "K_means_Cl") %>% dplyr::select(-K_means_Cl)
        }
      }

    }

    else if(Strategy == "Batch_K_means"){
      #First we calculate a metric to decide the number of total phenotypes
      #Specify the params
      if(Batch_size > nrow(MARKERS)) stop(paste0("Batch_size cannot be larger than ", nrow(MARKERS)))

      params_mbkm <- list(batch_size = Batch_size,
                          init_fraction = 1,
                          early_stop_iter = 10)
      #Run the specified test using each of the number of clusters
      print("Starting Cluster number stimation process")
      Optimal <- try(ClusterR::Optimal_Clusters_KMeans(MARKERS %>% scale(),
                                                       max_clusters = Max_N_phenotypes_Batch,
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
        N_Phenotypes<- menu(choices = as.character(1:Max_N_phenotypes_Batch),
                            title = paste0("Look at the plot generated, Then decide the appropiate number of Phenotypes"))

        #Calculate the desired number of clusters with batch k menas
        print("Performing Batched K means algorithm")
        Batch_k_means <- ClusterR::MiniBatchKmeans(MARKERS %>% scale(),
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
        pr_mb <- predict(object = Batch_k_means, fuzzy = F, newdata = MARKERS %>% scale())
        pr_mb <- as_tibble(pr_mb)
        names(pr_mb) <- "Phenotype"

        #Generate the data phenotypes tibble
        DATA_Phenotypes <-dplyr::bind_cols(DATA, pr_mb)
      }
    }

    else if(Strategy == "GMM"){
      #First we calculate a metric to decide the number of total phenotypes
      #Run the specified test using each of the number of clusters
      print("Starting Cluster number stimation process")
      Optimal <- try(ClusterR::Optimal_Clusters_GMM(MARKERS %>% scale(),
                                                    criterion = Quality_metric,
                                                    max_clusters = Max_N_phenotypes_GMM,
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
        N_Phenotypes<- menu(choices = as.character(1:Max_N_phenotypes_GMM),
                            title = paste0("Look at the plot generated, Then decide the appropiate number of Phenotypes"))

        #Calculate the desired number of clusters with batch k menas
        print("Calculating Gaussian Mixed Model")
        GMM_model <- ClusterR::GMM(MARKERS %>% scale(),
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
        pr_mb <- predict(object = GMM_model, fuzzy = F, newdata = MARKERS %>% scale())
        pr_mb <- as_tibble(pr_mb)
        names(pr_mb) <- "Phenotype"

        #Generate the data phenotypes tibble
        DATA_Phenotypes <-dplyr::bind_cols(DATA, pr_mb)
      }
    }

    else if(Strategy == "CLARA_clustering"){
      #First we calculate a metric to decide the number of total phenotypes
      print("Starting Cluster number stimation process")
      Optimal <-  try(ClusterR::Optimal_Clusters_Medoids(MARKERS %>% scale(),
                                                         max_clusters = Max_N_phenotypes_CLARA,
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
        N_Phenotypes<- menu(choices = as.character(1:Max_N_phenotypes_CLARA),
                            title = paste0("Based on the plots generated and you previous choice, decide the appropiate number of final Phenotypes"))
        print("Performing CLARA (Clustering Large Applications)")
        CLARA_Clustering <- ClusterR::Clara_Medoids(MARKERS %>% scale(),
                                                    clusters = as.double(N_Phenotypes),
                                                    samples = Samples_CLARA,
                                                    sample_size = Sample_per_CLARA,
                                                    distance_metric = Distance_CLARA,
                                                    threads = N_cores,
                                                    swap_phase = F,
                                                    fuzzy = FALSE,
                                                    verbose = T,
                                                    seed = 21)
        #Assign the cluster to each observation of MARKER
        pr_mb <- predict(object = CLARA_Clustering, fuzzy = F, newdata = MARKERS %>% scale())
        pr_mb <- as_tibble(pr_mb)
        names(pr_mb) <- "Phenotype"

        #Generate the data phenotypes tibble
        DATA_Phenotypes <-dplyr::bind_cols(DATA, pr_mb)
      }
    }

    #If there are noisy and real cells bind both tibbles
    if(Apply_Denoise){
      DATA_Phenotypes <- DATA_Phenotypes %>%dplyr::mutate(Phenotype = as.numeric(as.numeric(Phenotype) + 1))
      DATA_Phenotypes <-dplyr::bind_rows(DATA_NOISE, DATA_Phenotypes) %>% dplyr::arrange(Unique_ID) %>% dplyr::select(-Unique_ID, -NOISE)
      warning("If denoising is applied, Cluster nº1 contains the noisy cells")
    }

    #Turn Phenotype into a factor
    DATA_Phenotypes <- DATA_Phenotypes %>%dplyr::mutate(Phenotype = factor(Phenotype))

    #plot dimension reduction according to the number of cells if required
    if(Perform_Dimension_reduction){
      if(nrow(DATA_Phenotypes) <= 100000){
        try(plot(
          dplyr::left_join(DATA_Phenotypes, DATA_Reduction, by = "Cell_no") %>%
            ggplot(aes(x = DIMENSION_1, y = DIMENSION_2, color = Phenotype)) +
            geom_point(size = 2, alpha = 0.95) +
            cowplot::theme_cowplot() +
            scale_color_manual("Phenotype", values = unname(pals::polychrome(length(unique(DATA_Phenotypes$Phenotype)))))
        )
        )
      }
      if(nrow(DATA_Phenotypes) > 100000){
        message(">100K observations to generate plots. A random subset containing 10% of the dataset will be selected for Dimension reduction plots")
        try(plot(
          dplyr::left_join(DATA_Phenotypes, DATA_Reduction, by = "Cell_no") %>% sample_n(size = 100000) %>%
            ggplot(aes(x = DIMENSION_1, y = DIMENSION_2, color = Phenotype)) +
            geom_point(size = 2, alpha = 0.95) +
            cowplot::theme_cowplot() +
            scale_color_manual("Phenotype", values = unname(pals::polychrome(length(unique(DATA_Phenotypes$Phenotype)))))
        )
        )
      }
    }


    #Skip plot rendering if too many cell features
    if(ncol(DATA_Phenotypes %>% dplyr::select(-c(1:4))) > 100){
      Plot_required <- menu(choices = c("Proceed", "Abort"), title = "More than 100 features. Proceed with plot rendering?")
      if(Plot_required == 2){
        #Print a summary with the results
        print(DATA_Phenotypes %>% dplyr::count(Neighborhood_assignment))

        if(Perform_Dimension_reduction) return(list(DATA = DATA_Phenotypes,
                                                    Dimension_reduction = DATA_Reduction)
        )
        else return(DATA_Phenotypes)
      }
    }

    #Visualize the neighbor composition data for each neighborhood
    plot(DATA_Phenotypes %>% dplyr::select(-c(1:4)) %>% pivot_longer(cols = -Phenotype) %>%
           ggplot(aes(x = as.factor(Phenotype), y = value)) +
           geom_violin(aes(color = name, fill = name), alpha=0.3, position=position_dodge(width=0.5)) +
           stat_summary(aes(color = name),
                        fun = median, geom = "crossbar", width = 0.4, linetype = 1, linewidth = 0.5,
                        position = position_dodge(width = 0.5)) +
           cowplot::theme_cowplot()+
           scale_x_discrete("Phenotype")+
           scale_y_continuous("Marker intensity"))

    #Visualize the heatmap of mean by neighborhood
    Mean_tibble <- DATA_Phenotypes %>% dplyr::select(-c(1:4)) %>% group_by(Phenotype) %>% summarize_all(mean) %>%dplyr::ungroup() #Obtain mean tibble
    Mean_matrix <- as.matrix(Mean_tibble[-1] %>% scale()) #Scale it and transform it into a  mtrix
    row.names(Mean_matrix) <- Mean_tibble[[1]]

    plot(ComplexHeatmap::Heatmap(Mean_matrix,
                                 name = "Scaled")
    )

    #Print the summary
    print(DATA_Phenotypes %>% dplyr::count(Phenotype))
    #Return data and Dimension reduction if generated
    if(Perform_Dimension_reduction) return(list(DATA = DATA_Phenotypes,
                                                Dimension_reduction = DATA_Reduction)
    )
    else return(DATA_Phenotypes)

  }
