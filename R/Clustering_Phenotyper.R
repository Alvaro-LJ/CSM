#' Performs phenotype assignment using clustering approaches
#'
#' `Clustering_Phenotyper()` assigns cell phenotypes based on clustering. Phenotype names can be then modidified using [DATA_Phenotype_renamer()]
#'
#' @param DATA A dataframe or tibble containing cell feature data.
#' @param Strategy The clustering strategy. One of the following Consensus_Clustering, SOM, Graph_Based, K_Means_Meta_clustering, Batch_K_means, GMM or CLARA_clustering (see details).
#' @param Force_N_Phenotypes A logical value indicating if the number of phenotypes indicated should be used directly to compute the phenotypes. If TRUE, no estimating process is conducted. Applicable for SOM, Batch_K_means, GMM or CLARA_clustering.
#'
#' @param Apply_Denoise A logical value. Specify if a denoising filtering is required before clustering (see details).
#' @param Denoising Denoising strategy. A vector containing any of the following: Quantile, Standard_Deviation, Threshold, Otsu or DimRed_DBscan. More than one method can be applied.
#' @param Percentile A numeric value indicating the percentile for quantile denoising. Cells below percentile for all features will be considered to be noise.
#' @param N_Standard_Deviations A numeric value indicating the number of standard deviations from mean for Standard_Deviation method. Cells below SD for all features will be considered to be noise.
#' @param Selected_threshold A numeric value indicating the threshold to be used for the Threshold method. Cells below the threshold for all features will be considered to be noise.
#' @param Min_cell_no An integer value for the DBscan method. Minimum cell number to consider a group of cells to be clustered.
#' @param Distance_radius A numeric value for the DBscan method. Distance to be sampled.
#'
#' @param Perform_Dimension_reduction Logical value. Should Dimension Reduction be performed (see details).
#' @param Dimension_reduction Dimension reduction method. One of the following: PCA, TSNE, UMAP.
#' @param Dimension_reduction_prop A numeric value between 0 and 1 to indicate the percentage of the cells to be used in dimension computation (applicable for TSNE and UMAP).
#' @param Cluster_on_Reduced A logical value indicating if clustering should be performed on new dimensions.
#'
#' @param Stop_at_preprocessing A logical value indicating if the function should stop after Data pre-processing and return the interim results (see details).
#' @param Pre_processed_data (OPTIONAL) If a pre-processing object is available, it can be provided to skip data pre-processing steps (see details).
#'
#' @param Max_N_phenotypes If Strategy is Consensus_Clustering: Number of maximum phenotypes that can be identified.
#' @param Consensus_reps If Strategy is Consensus_Clustering: Number of iterations to converge.
#' @param Consensus_p_Items If Strategy is Consensus_Clustering: Percentage of cells that you desire to sample in each iteration.
#' @param Consensus_Cluster_Alg If Strategy is Consensus_Clustering: Clustering algorithm to be used: hc (hierarchical clustering), pam (paritioning around medoids), km (for k-means).
#' @param Consensus_Distance If Strategy is Consensus_Clustering: Distance metric to be used (pearson(1 - Pearson correlation), spearman(1 - Spearman correlation), euclidean, binary, maximum, canberra, minkowski.
#' @param Consensus_Name If Strategy is Consensus_Clustering: Name of the folder that is going to be created in order to place the resulting graphs.
#'
#' @param Max_SOM_phenotypes If Strategy is SOM: number of maximum phenotypes that can be identified.
#'
#' @param Nearest_neighbors_for_graph If strategy is Graph_Based: The number of closest neighbors to calculate the graph.
#' @param Graph_Method If strategy is Graph_Based: One of Louvain, Leiden, Greedy, WalkTrap, Spinglass, Leading_Eigen or Edge_Betweenness.
#' @param Graph_Resolution If strategy is Graph_Based: Used for Louvain and Leiden. 1 is default. The smaller the value, the larger the clusters will be.
#' @param N_steps If strategy is Graph_Based: Number of steps given in the WalkTrap algorithm.
#'
#' @param N_K_centroids If strategy is K_Means_Meta_clustering: Number of centroids to perform K means.
#' @param Max_N_phenotypes_Meta If strategy is K_Means_Meta_clustering: Number of maximum phenotypes that can be identified.
#' @param Consensus_reps_Meta If strategy is K_Means_Meta_clustering: Number of iterations to converge.
#' @param Consensus_p_Items_Meta If strategy is K_Means_Meta_clustering: Percentage of cells to sample in each iteration.
#' @param Consensus_Name_Meta If strategy is K_Means_Meta_clustering: Name of the folder that is going to be created in order to place the resulting graphs.
#'
#' @param Batch_size If strategy is Batch_K_means: Number of cells to be included in each random batch.
#' @param Max_N_phenotypes_Batch If strategy is Batch_K_means: Number of maximum phenotypes that can be identified.
#' @param N_initiations If strategy is Batch_K_means: Number of times the algorithm is going to be tried to find the best clustering result.
#' @param Max_iterations If strategy is Batch_K_means: Max number of iterations in each try.
#'
#' @param Quality_metric If strategy is GMM: The quality measure used to test the number of clusters ("AIC" or "BIC").
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
#' Dimension reduction can be performed using PCA (svd::propack.svd function), t-SNE (snifter::fitsne function) and UMAP (uwot::tumap function). For t-SNE and UMAP, a model can be fitted using a subset of data and then generalized to all the dataset. This can be more computationally efficient.
#' The first time TSNE is used, snifter will install a CONDA distribution in yout R library using the basilisk package. The path to the library is relevant. Any non-standard characters in you path may lead to errors (non ASCII characters, spaces...).
#'
#' If the dataset is really large, the pre-processing steps (de-noising and dimension reduction) and clustering can take a while. This can make the tuning of
#' clustering parameters tedious (each try takes a lot of time). To tackle this issue, pre-processing and clustering steps can be split into two independent steps.
#' By setting the Stop_at_preprocessing argument to TRUE, the function will perform only the required pre-processing and will return the result without performing the
#' clustering. Once the user is satisfied with the pre-processing results, the object can be provided to the function to execute the clustering process.
#'
#' Consensus clustering is performed using the ConsensusClusterPlus::ConsensusClusterPlus function.
#'
#' Self Organizing Maps clustering is performed using the FlowSOM::FlowSOM function.
#'
#' For graph based clustering, Nearest Neighbors Graphs (SNNG) are built using bluster::makeSNNGraph and clustered using functions included in the igraph package.
#'
#' K_Means_Meta_clustering first summarizes the cell feature matrix observations using K means algorithm and the performs Consensus Clustering. Afterwards results are generalized to all cells.
#'
#' Batch K-means, Gaussian Mixture Models and Clustering Large Applications (CLARA) are all based on the ClusterR package.
#'
#' @seealso [DATA_Phenotype_renamer()], [ReClustering_function()], [Consensus_phenotype_assigner()], [Concordance_calculator()], [Confusion_matrix_plotter()],
#' [Phenotyping_evaluator_shiny_app_launcher()]
#'
#' @returns Returns a tibble with cell features and a column named 'Phenotype' containing cell labels.
#' If dimension reduction has been performed, returns a list with the cell feature dataset as above and a tibble containing dimension reduction coordinates.
#'
#' @examples
#' \dontrun{
#' Clustering_Phenotyper(
#'     DATA = CSM_Arrangedcellfeaturedata_test,
#'     Strategy = "GMM",
#'
#'     Apply_Denoise = TRUE,
#'     Denoising = "DimRed_DBscan",
#'     Min_cell_no = 5,
#'     Distance_radius = 100,
#'
#'     Perform_Dimension_reduction = TRUE,
#'     Dimension_reduction = "UMAP",
#'     Dimension_reduction_prop = 1,
#'     Cluster_on_Reduced = TRUE,
#'
#'     Quality_metric = "AIC",
#'     Max_N_phenotypes_GMM = 5,
#'     Max_iterations_km = 15,
#'     Max_iterations_em = 15,
#'     GMM_Distance = "eucl_dist"
#')
#' }
#'
#'
#' @export

Clustering_Phenotyper <-
  function(DATA,
           Strategy,
           Force_N_Phenotypes = FALSE,

           #Denoising parameters
           Apply_Denoise = FALSE,
           Denoising = NULL,
           Percentile = NULL,
           N_Standard_Deviations = NULL,
           Selected_threshold = NULL,
           Min_cell_no = NULL,
           Distance_radius = NULL,

           #Dimension reduction
           Perform_Dimension_reduction = FALSE,
           Dimension_reduction = NULL,
           Dimension_reduction_prop = NULL,
           Cluster_on_Reduced = FALSE,

           #Perform only Pre_processing
           Stop_at_preprocessing = FALSE,
           Pre_processed_data = NULL,

           #Parameters for Consensus Clustering
           Max_N_phenotypes = NULL,
           Consensus_reps = NULL,
           Consensus_p_Items = NULL,
           Consensus_Cluster_Alg = NULL,
           Consensus_Distance = NULL,
           Consensus_Name = NULL,

           #Parameters for Self-Organizing Maps
           Max_SOM_phenotypes = NULL,

           #Parameters for Graph-Based approaches
           Nearest_neighbors_for_graph = NULL,
           Graph_Method = NULL,
           Graph_Resolution = NULL,
           N_steps = NULL,

           #Parameters for K means Meta Clustering
           N_K_centroids = NULL,
           Max_N_phenotypes_Meta = NULL,
           Consensus_reps_Meta = NULL,
           Consensus_p_Items_Meta = NULL,
           Consensus_Name_Meta = NULL,

           #Parameters for Batched K means
           Batch_size = NULL,
           Max_N_phenotypes_Batch = NULL,
           N_initiations = NULL,
           Max_iterations = NULL,

           #Parameters for Gaussian Mixture Model
           Quality_metric = NULL,
           Max_N_phenotypes_GMM = NULL,
           Max_iterations_km = NULL,
           Max_iterations_em = NULL,
           GMM_Distance = NULL,

           #Parameters for CLARA clustering
           Samples_CLARA = NULL,
           Sample_per_CLARA = NULL,
           Max_N_phenotypes_CLARA = NULL,
           Distance_CLARA = NULL,
           N_cores = NULL
  ) {
    on.exit(gc())

    ##################################GENERAL ARGUMENT CHECK######################################

    #Check general arguments
    if(!Strategy %in% c("Consensus_Clustering", "SOM", "Graph_Based", "K_Means_Meta_clustering", "Batch_K_means", "GMM", "CLARA_clustering")){
      stop("Strategy must be one of the following: Consensus_Clustering, SOM, Graph_Based, K_Means_Meta_clustering, Batch_K_means, GMM, CLARA_clustering")
    }
    if(!all(length(Force_N_Phenotypes) == 1, is.logical(Force_N_Phenotypes))) stop("Force_N_Phenotypes must be a logical value")
    if(Force_N_Phenotypes){
      if(!Strategy %in% c("SOM", "Batch_K_means", "GMM", "CLARA_clustering")) message("Force_N_Phenotypes cannot be used with current strategy, argument will be ignored")
    }

    #If NO Pre-processed data provided check if Stop_at_pre_processing is logical and other pre-processing variables
    if(is.null(Pre_processed_data)){
      if(!identical(names(DATA)[1:4], c("Cell_no", "X", "Y", "Subject_Names"))) {
        stop("Please generate an appropiate data object using the Data_arrange_function")
      }
      if(!is.logical(Stop_at_preprocessing)) stop("Stop_at_preprocessing should be a logical value")
      if(!is.logical(Apply_Denoise)) stop("Apply_Denoise must be a logical value")
      if(!is.logical(Perform_Dimension_reduction)) stop("Perform_Dimension_reduction must be a logical value")
      if(all(Stop_at_preprocessing, !Apply_Denoise, !Perform_Dimension_reduction)) message("Stop_at_processing will be ignored as no data pre-processing steps are required")
      if(Perform_Dimension_reduction){
        if(!Dimension_reduction %in% c("UMAP", "TSNE", "PCA")) stop("Dimension_reduction must be one of the following: UMAP, TSNE, PCA")
        if(!all(is.numeric(Dimension_reduction_prop), Dimension_reduction_prop > 0, Dimension_reduction_prop <= 1)) stop("Dimension_reduction_prop must be a numeric value between 0 and 1")
      }
      if(Apply_Denoise){
        #check denoising argument is correctly stated
        if(!all(Denoising %in% c("Quantile", "Standard_Deviation", "Threshold", "Otsu", "DimRed_DBscan"))) {
          stop("Denoising should be any of the following: Quantile, Standard_Deviation, Threshold, Otsu, DimRed_DBscan")
        }
        if(any(Denoising == "DimRed_DBscan")){
          #Requires previous dimension reduction
          if(!Perform_Dimension_reduction) stop("DBscan clustering requires Dimension reduction. Please set Perform_Dimension_reduction to TRUE")
        }
      }
      if(!is.logical(Cluster_on_Reduced)) stop("Cluster_on_Reduced must be a logical value")
      if(Cluster_on_Reduced){
        if(!Perform_Dimension_reduction) stop("If Clustering needs to be performed on Dimension reduced data please set Perform_Dimension_reduction to TRUE")
      }
    }

    #If Pre-processed data provided obtain the datasets from the object provided
    if(!is.null(Pre_processed_data)){
      message("Pre_processed_data provided. Pre-processing related arguments will be ignored.")
      if(names(Pre_processed_data)[1] != "Pre_processing_argument") stop("Pre_processing_argument not found in Pre_processed_data object provided")

      Apply_Denoise <- Pre_processed_data[["Pre_processing_argument"]][["Apply_Denoise"]]
      Perform_Dimension_reduction <- Pre_processed_data[["Pre_processing_argument"]][["Perform_Dimension_reduction"]]

      if(Apply_Denoise & !Perform_Dimension_reduction){
        DATA <- Pre_processed_data[["DATA"]]
        MARKERS <- Pre_processed_data[["MARKERS"]]
        DATA_NOISE <- Pre_processed_data[["DATA_NOISE"]]

      }
      if(!Apply_Denoise & Perform_Dimension_reduction){
        DATA <- Pre_processed_data[["DATA"]]
        MARKERS <- Pre_processed_data[["MARKERS"]]
        DATA_Reduction <- Pre_processed_data[["DATA_Reduction"]]
      }
      if(Apply_Denoise & Perform_Dimension_reduction){
        DATA <- Pre_processed_data[["DATA"]]
        MARKERS <- Pre_processed_data[["MARKERS"]]
        DATA_NOISE <- Pre_processed_data[["DATA_NOISE"]]
        DATA_Reduction <- Pre_processed_data[["DATA_Reduction"]]
      }
    }

    ##################################SUGGESTED PACKAGE CHECK######################################

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
      Argument_checker <- c(N_K_centroids_OK = all(nrow(DATA) > N_K_centroids, N_K_centroids%%1 == 0, N_K_centroids > 0),
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
    #If pre-processed data is not provided, check suggested packages
    if(is.null(Pre_processed_data)){
      if(Apply_Denoise){
        #Check suggested packages
        if("Otsu" %in% Denoising){
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
        if("DimRed_DBscan" %in% Denoising){
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

    #If pre-processed data is not provided, perform data pre-processing as required by user
    if(is.null(Pre_processed_data)){

      DATA_Reduction <- NULL #Generate a NULL DATA_Reduction just in case it is not generated and evaluated by Denoising function

      #Perform dimension reduction if required
      if(Perform_Dimension_reduction){
        DATA_Reduction <-
          CSM_Dimension_reduction_function(Original_data = DATA,

                                           Dimension_reduction_strategy = Dimension_reduction,
                                           Dimension_reduction_prop = Dimension_reduction_prop
                                           )
      }

      #If denoising is required apply required function
      if(Apply_Denoise){

        print("Filtering out noisy cells")

        Denoising_results <-
          CSM_Denoising_function(
            Original_data = DATA,

            Denoising_strategy = Denoising,

            Percentile = Percentile,
            N_Standard_Deviations = N_Standard_Deviations,
            Selected_threshold = Selected_threshold,

            Perform_Dimension_reduction = Perform_Dimension_reduction,
            DATA_Reduction = DATA_Reduction,
            Min_cell_no = Min_cell_no,
            Distance_radius = Distance_radius
          )

        #Generate two tibbles, one with noisy cells and other (DATA) with the actual cells
        NOISE_VECTOR <- Denoising_results[["NOISE_VECTOR"]]
        DATA_NOISE <- Denoising_results[["DATA_NOISE"]]
        DATA <- Denoising_results[["DATA"]]
        MARKERS <- Denoising_results[["MARKERS"]]
      }

      #If no denoising required obtain MARKERS directly from DATA
      else{
        DATA <- DATA
        MARKERS <- DATA %>% dplyr::select(-Cell_no, -X, -Y, -Subject_Names)
      }

      #Generate a specific version of Markers with dimension reduction data for Dimension_SNN
      if(Cluster_on_Reduced){
        #Depending on Denoising Obtain directly from DATA_Reduction or filter first
        if(!Apply_Denoise) MARKERS <- DATA_Reduction[-which(names(DATA_Reduction) == "Cell_no")]
        if(Apply_Denoise) MARKERS <- DATA_Reduction[!NOISE_VECTOR, ] %>% dplyr::select(-Cell_no)
      }

      #Check that MARKERS is a numeric matrix
      if(!is.numeric(as.matrix(MARKERS))) stop("DATA provided must contain marker intensity values")

      #If Stop_at_preprocessing is true return the interim results
      if(Stop_at_preprocessing){
        #Obtain the arguments
        Pre_processing_argument <- list(
          Apply_Denoise = Apply_Denoise,
          Denoising = Denoising,
          Percentile = Percentile,
          N_Standard_Deviations = N_Standard_Deviations,
          Selected_threshold = Selected_threshold,
          Min_cell_no = Min_cell_no,
          Distance_radius = Distance_radius,
          Perform_Dimension_reduction = Perform_Dimension_reduction,
          Dimension_reduction = Dimension_reduction,
          Dimension_reduction_prop = Dimension_reduction_prop,
          Cluster_on_Reduced = Cluster_on_Reduced
        )

        #If only De-noising needs to be performed
        if(Apply_Denoise & !Perform_Dimension_reduction){
          #Print a message with the results
          message(paste0(nrow(DATA_NOISE), " cells have been identified as noise"))
          #Return the list
          return(list(Pre_processing_argument = Pre_processing_argument,
                      DATA = DATA,
                      MARKERS = MARKERS,
                      DATA_NOISE = DATA_NOISE
          ))
        }
        #If only DIM reduction needs to be performed
        if(!Apply_Denoise & Perform_Dimension_reduction){
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
                      DATA = DATA,
                      MARKERS = MARKERS,
                      DATA_Reduction = DATA_Reduction
          ))
        }
        if(Apply_Denoise & Perform_Dimension_reduction){
          message(paste0(nrow(DATA_NOISE), " cells have been identified as noise"))
          #Print the graph according to results
          if(nrow(DATA_Reduction) <= 100000){
            print("Generating Dimension Reduction Plot")
            plot(
              DATA_Reduction %>% dplyr::mutate(NOISE = NOISE_VECTOR) %>% ggplot(aes(x = DIMENSION_1, y = DIMENSION_2, color = NOISE)) +
                geom_point(size = 2, alpha = 0.95) +
                cowplot::theme_cowplot()
            )
          }
          if(nrow(DATA_Reduction) > 100000){
            print("More than 100K observations. Generating Dimension Reduction Plot on a random subset of the data")
            plot(
              DATA_Reduction %>% dplyr::mutate(NOISE = NOISE_VECTOR) %>% dplyr::slice_sample(prop = 0.1) %>% ggplot(aes(x = DIMENSION_1, y = DIMENSION_2, color = NOISE)) +
                geom_point(size = 2, alpha = 0.95) +
                cowplot::theme_cowplot()
            )
          }
          #Return the list
          return(list(Pre_processing_argument = Pre_processing_argument,
                      DATA = DATA,
                      MARKERS = MARKERS,
                      DATA_NOISE = DATA_NOISE,
                      DATA_Reduction = DATA_Reduction
          ))
        }
      }
    }

    ##################################CLUSTERING######################################

    DATA_Phenotypes <-
      CSM_Clustering_function(
        Original_data = DATA,
        MARKERS = MARKERS,

        Strategy = Strategy,
        Force_N_Clusters = Force_N_Phenotypes,

        Max_N_clusters_Consensus = Max_N_phenotypes,
        Consensus_reps = Consensus_reps,
        Consensus_p_Items = Consensus_p_Items,
        Consensus_Cluster_Alg = Consensus_Cluster_Alg,
        Consensus_Distance = Consensus_Distance,
        Consensus_Name = Consensus_Name,

        Max_SOM_clusters = Max_SOM_phenotypes,

        Graph_type = "SNN", #ONLY SNN supported in for clustering phenotyper
        Nearest_neighbors_for_graph = Nearest_neighbors_for_graph,
        Graph_Method = Graph_Method,
        Graph_Resolution = Graph_Resolution,
        N_steps = N_steps,

        N_K_centroids = N_K_centroids,
        Max_N_clusters_Meta = Max_N_phenotypes_Meta,
        Consensus_reps_Meta = Consensus_reps_Meta,
        Consensus_p_Items_Meta = Consensus_p_Items_Meta,
        Consensus_Name_Meta = Consensus_Name_Meta,

        Batch_size = Batch_size,
        Max_N_clusters_Batch = Max_N_phenotypes_Batch,
        N_initiations = N_initiations,
        Max_iterations = Max_iterations,

        Quality_metric = Quality_metric,
        Max_N_clusters_GMM = Max_N_phenotypes_GMM,
        Max_iterations_km = Max_iterations_km,
        Max_iterations_em = Max_iterations_em,
        GMM_Distance =  GMM_Distance,

        Samples_CLARA = Samples_CLARA,
        Sample_per_CLARA = Sample_per_CLARA,
        Max_N_clusters_CLARA = Max_N_phenotypes_CLARA,
        Distance_CLARA = Distance_CLARA,
        N_cores = N_cores
      )

    #Change the name of the column 'Cluster' for 'Phenotype'
    DATA_Phenotypes <- DATA_Phenotypes %>% dplyr::rename("Phenotype" = "Cluster")

    ##################################RESULT PLOTTING AND FUNCTION EXIT######################################

    #If there are noisy and real cells bind both tibbles
    if(Apply_Denoise){
      #Change the name of the column 'Cluster' for 'Phenotype'
      DATA_NOISE <- DATA_NOISE %>% dplyr::rename("Phenotype" = "Cluster")

      DATA_Phenotypes <- DATA_Phenotypes %>% dplyr::mutate(Phenotype = as.numeric(as.numeric(Phenotype) + 1))
      DATA_Phenotypes <- dplyr::bind_rows(DATA_NOISE, DATA_Phenotypes) %>% dplyr::arrange(Unique_ID) %>% dplyr::select(-Unique_ID)
      warning("If denoising is applied, Cluster number 1 contains the noisy cells")
    }

    #Turn Phenotype into a factor
    DATA_Phenotypes <- DATA_Phenotypes %>% dplyr::mutate(Phenotype = factor(Phenotype))

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
          dplyr::left_join(DATA_Phenotypes, DATA_Reduction, by = "Cell_no") %>% dplyr::slice_sample(n = 100000) %>%
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
    plot(DATA_Phenotypes %>% dplyr::select(-c(1:4)) %>% tidyr::pivot_longer(cols = -Phenotype) %>%
           ggplot(aes(x = as.factor(Phenotype), y = value)) +
           geom_violin(aes(color = name, fill = name), alpha=0.3, position=position_dodge(width=0.5)) +
           stat_summary(aes(color = name),
                        fun = median, geom = "crossbar", width = 0.4, linetype = 1, linewidth = 0.5,
                        position = position_dodge(width = 0.5)) +
           cowplot::theme_cowplot()+
           scale_x_discrete("Phenotype")+
           scale_y_continuous("Marker intensity"))

    #Visualize the heatmap of mean by neighborhood
    Mean_tibble <- DATA_Phenotypes %>% dplyr::select(-c(1:4)) %>% group_by(Phenotype) %>% summarize_all(mean) %>% dplyr::ungroup() #Obtain mean tibble
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

