#' Performs sub-clustering of specific cell populations
#'
#' After phenotyping, specific cell phenotypes can be further subclassified according to previously un-used features. Using this function users may cluster cells in a step wise manner.
#'
#'
#' @param DATA A dataframe or tibble containing cell feature data.
#' @param DATA_aside (OPTIONAL) A dataframe or tibble containing cell feature data that can be used in the re-clustering process.
#' @param Phenotype_variable A character value indicating the name of the column containing cell phenotype labels.
#' @param Phenotype_to_recluster A character value indicating the name of cell phenotype label to be re-clustered.
#' @param Variables_to_recluster A character vector indicating the variables contained in DATA or DATA_aside to be used in the re-clustering process.
#' @param Strategy One of the following: "Arbitrary", "Multilevel", "Consensus_Clustering", "SOM", "Graph_Based", "K_Means_Meta_clustering", "Batch_K_means", "GMM" or "CLARA_clustering" (see details).
#' @param Force_N_Phenotypes A logical value indicating if the number of clusters indicated (using the N_clusters argument) should be used directly to compute the clusters. If TRUE, no estimating process is conducted. Applicable for SOM, Batch_K_means, GMM or CLARA_clustering.
#'
#' @param N_Clusters For Consensus_Clustering, SOM, Graph_Based, K_Means_Meta_clustering, Batch_K_means, GMM or CLARA_clustering, an integer indicating the number of clusters to generate.
#' @param Levels For Multilevel, an integer indicating the number of levels to divide each feature.
#' @param Arbitrary_cutoff For arbitrary, a numeric value or numeric vector indicating the cut-off values to divide each feature.
#'
#' @param Consensus_reps If Strategy is Consensus_Clustering: Number of iterations to converge.
#' @param Consensus_p_Items If Strategy is Consensus_Clustering: Percentage of cells that you desire to sample in each iteration.
#' @param Consensus_Cluster_Alg If Strategy is Consensus_Clustering: Clustering algorithm to be used: hc (hierarchical clustering), pam (paritioning around medoids), km (for k-means).
#' @param Consensus_Distance If Strategy is Consensus_Clustering: Distance metric to be used (pearson(1 - Pearson correlation), spearman(1 - Spearman correlation), euclidean, binary, maximum, canberra, minkowski.
#' @param Consensus_Name If Strategy is Consensus_Clustering: Name of the folder that is going to be created in order to place the resulting graphs.
#'
#' @param Nearest_neighbors_for_graph If strategy is Graph_Based: The number of closest neighbors to calculate the graph.
#' @param Graph_Method If strategy is Graph_Based: One of Louvain, Leiden, Greedy, WalkTrap, Spinglass, Leading_Eigen or Edge_Betweenness.
#' @param Graph_Resolution If strategy is Graph_Based: Used for Louvain and Leiden. 1 is default. The smaller the value, the larger the clusters will be.
#' @param N_steps If strategy is Graph_Based: Number of steps given in the WalkTrap algorithm.
#'
#' @param N_K_centroids If strategy is K_Means_Meta_clustering: Number of centroids to perform K means.
#' @param Consensus_reps_Meta If strategy is K_Means_Meta_clustering: Number of iterations to converge.
#' @param Consensus_p_Items_Meta If strategy is K_Means_Meta_clustering: Percentage of cells to sample in each iteration.
#' @param Consensus_Name_Meta If strategy is K_Means_Meta_clustering: Name of the folder that is going to be created in order to place the resulting graphs.
#'
#' @param Batch_size If strategy is Batch_K_means: Number of cells to be included in each random batch.
#' @param N_initiations If strategy is Batch_K_means: Number of times the algorithm is going to be tried to find the best clustering result.
#' @param Max_iterations If strategy is Batch_K_means: Max number of iterations in each try.
#'
#' @param Quality_metric If strategy is GMM: The quality measure used to test the number of clusters ("AIC" or "BIC").
#' @param Max_iterations_km If strategy is GMM: Number of max iterations in the K means clustering performed.
#' @param Max_iterations_em If strategy is GMM: Number of max iterations in the Expectation Maximization algorithm.
#' @param GMM_Distance If strategy is GMM: Distance metric used in the model ("eucl_dist" or "maha_dist").
#'
#' @param Samples_CLARA If strategy is CLARA_clustering: Number of samples the CLARA algorithm is going to use to be calculated.
#' @param Sample_per_CLARA If strategy is CLARA_clustering: Percentage (from 0 to 1) of the total cells that are going to be allocated to each sample.
#' @param Distance_CLARA If strategy is CLARA_clustering: Distance metric used in the model (euclidean, manhattan, chebyshev, canberra, braycurtis, pearson_correlation, simple_matching_coefficient, minkowski, hamming, jaccard_coefficient, Rao_coefficient, mahalanobis, cosine)
#' @param N_cores If strategy is CLARA_clustering: Number of cores to parallelize your computation
#'
#' @details
#' Multi_level thresholds are calculated using imagerExtra::ThresholdML function.
#'
#' Self Organizing Maps clustering is performed using the FlowSOM::FlowSOM function.
#'
#' Batch K-means and Gaussian Mixture Models are all based on the ClusterR package.
#'
#' @returns A tibble containing cell features and the new column with labels after re-clusterization process.
#'
#' @examples
#' \dontrun{
#' #OPTIONAL: set aside any features that will not be used in the initial phenotyping process--
#' DATA_list <- Data_set_aside(
#'    DATA = CSM_Arrangedcellfeaturedata_test,
#'    Markers_to_set = "GZMB_AVERAGE"
#'    )
#'
#'#Perform initial phenotyping----------------
#'DATA_thresholded <-
#' Thresholding_function(
#'    DATA = DATA_list$DATA,
#'    Strategy = "TriClass_Otsu",
#'    Local_thresholding = FALSE,
#'    number_iterations_TriClass = 20
#' )
#'
#'Phenotype_possibilities <-
#' Marker_combinator_generator(
#'    DATA = DATA_thresholded,
#'    Markers = names(DATA_thresholded)[-c(1:4)]
#')
#'Phenotype_possibilities$Phenotype <- c("TUMOR", "OTHER", "CD8", "CD8")
#'
#'DATA_Phenotypes <-
#' Phenotype_assigner_function(
#'    DATA = DATA_thresholded,
#'    Phenotype_possibilities = Phenotype_possibilities
#'    )
#'
#'#Perform Re-Clustering of a subset of cells according to set aside features--------------
#' ReClustering_function(
#'     DATA = DATA_Phenotypes,
#'     DATA_aside = DATA_list$Aside,
#'     Phenotype_variable = "Phenotype",
#'     Phenotype_to_recluster = "CD8",
#'     Variables_to_recluster = "GZMB_AVERAGE",
#'     Strategy = "Multilevel",
#'     Levels = 3
#'     )
#' }
#'
#' @export

ReClustering_function <-
  function(DATA,
           DATA_aside = NULL,
           Phenotype_variable,
           Phenotype_to_recluster,
           Variables_to_recluster,
           Strategy,

           N_Clusters = NULL,
           Force_N_Clusters = FALSE,

           #Parameters for MULTILEVEL or Arbitrary
           Levels = NULL,
           Arbitrary_cutoff = NULL,

           #Parameters for Consensus Clustering
           Consensus_reps = NULL,
           Consensus_p_Items = NULL,
           Consensus_Cluster_Alg = NULL,
           Consensus_Distance = NULL,
           Consensus_Name = NULL,

           #Parameters for Graph-Based approaches
           Graph_type = NULL,
           Nearest_neighbors_for_graph = NULL,
           Graph_Method = NULL,
           Graph_Resolution = NULL,
           N_steps = NULL,

           #Parameters for K means Meta Clustering
           N_K_centroids = NULL,
           Consensus_reps_Meta = NULL,
           Consensus_p_Items_Meta = NULL,
           Consensus_Name_Meta = NULL,

           #Parameters for Batched K means
           Batch_size = NULL,
           N_initiations = NULL,
           Max_iterations = NULL,

           #Parameters for Gaussian Mixture Model
           Quality_metric = NULL,
           Max_iterations_km = NULL,
           Max_iterations_em = NULL,
           GMM_Distance = NULL,

           #Parameters for CLARA clustering
           Samples_CLARA = NULL,
           Sample_per_CLARA = NULL,
           Distance_CLARA = NULL,
           N_cores = NULL
  ){
    DATA <- DATA
    DATA_aside <- DATA_aside
    if(!all(identical(names(DATA)[1:4], c("Cell_no", "X", "Y", "Subject_Names")),
            ifelse(is.null(DATA_aside), yes = TRUE, no = identical(names(DATA_aside)[1:4], c("Cell_no", "X", "Y", "Subject_Names")))
    )) stop("DATA and DATA_aside must (if not NULL) must be adequately formatted")
    #Check the phenotype variable
    if(!all(length(Phenotype_variable) == 1, Phenotype_variable %in% names(DATA))) stop(paste0(Phenotype_variable, " must be present in DATA"))
    if(!Phenotype_to_recluster %in% unique(DATA[[Phenotype_variable]])) stop(paste0(Phenotype_to_recluster, " not present in ", stringr::str_c(unique(DATA[[Phenotype_variable]]), collapse = ", ")))
    #Check the variables to be re-clusterized
    if(!all(Variables_to_recluster %in% c(names(DATA), names(DATA_aside)))) stop(paste0(stringr::str_c(Variables_to_recluster, collapse = ", "),
                                                                                        " must be present either in DATA or in DATA_aside"))
    #If names in DATA_aside and in DATA are repeated and are required for reclustering then stop
    if(!is.null(DATA_aside)){
      Repeated_names <- names(DATA_aside)[names(DATA_aside) %in% names(DATA)] %in% Variables_to_recluster
      if(any(Repeated_names)){
        Conflictive_names <- names(DATA_aside)[names(DATA_aside) %in% names(DATA)][Repeated_names]
        stop(paste0(stringr::str_c(Conflictive_names, collapse = ", "), " is present in DATA and in DATA_aside. Please check data sources"))
      }
    }

    #Select only the cells that are going to be reclustered
    DATA_selected <- DATA[DATA[[Phenotype_variable]] == Phenotype_to_recluster, ]
    #Select the required variables from data and data aside
    DATA_variables <- DATA_selected %>% dplyr::select(Cell_no, any_of(Variables_to_recluster))
    #Add DATA_aside variables if required
    if(!is.null(DATA_aside)){
      #Select the required variables
      ASIDE_variables <- DATA_aside %>% dplyr::select(Cell_no, any_of(Variables_to_recluster))
      #Check that cells in DATA_aside match the required cells in the DATA
      if(!all(DATA_selected$Cell_no %in% ASIDE_variables$Cell_no)) stop("Cells in DATA and in DATA_aside do not match. Please check before running the re-clustering")
      #Bind to the other variables
      DATA_variables <-dplyr::left_join(DATA_variables, ASIDE_variables, by = "Cell_no")
    }
    #Check that all variables are numeric
    if(!all(purrr::map_lgl(DATA_variables[-1], function(Column) is.numeric(Column)))){
      Problematic_variable <- names(DATA_variables[-1])[!purrr::map_lgl(DATA_variables[-1], function(Column) is.numeric(Column))]
      stop(paste0("The following variables are not numeric and cannot be included in the reclustering process: ", stringr::str_c(Problematic_variable, collapse = ", ")))
    }

    #Check the strategy
    if(!Strategy %in% c("Arbitrary", "Multilevel", "Consensus_Clustering", "SOM", "Graph_Based", "K_Means_Meta_clustering", "Batch_K_means", "GMM", "CLARA_clustering")) stop(paste0("Strategy must be one of the following: ", stringr::str_c(c("Arbitrary", "Multilevel", "Consensus_Clustering", "SOM", "Graph_Based", "K_Means_Meta_clustering", "Batch_K_means", "GMM", "CLARA_clustering"), collapse = ", ")))

    #Check force neighborhoods
    if(!all(length(Force_N_Clusters) == 1, is.logical(Force_N_Clusters))) stop("Force_N_Clusters must be a logical value")
    if(Force_N_Clusters){
      if(!Strategy %in% c("SOM", "Batch_K_means", "GMM", "CLARA_clustering")) message("Force_N_Clusters cannot be used with current strategy, argument will be ignored")
    }

    #Check arguments for Arbitrary
    if(Strategy == "Arbitrary"){
      if(!is.numeric(Arbitrary_cutoff)) stop("Arbitrary_cutoff must be a numeric vector")
      if(length(Variables_to_recluster) > 1) message("The same Arbitrary_cutoff will be applied to all variables used in re-clustering")
    }
    #Check arguments for Multilevel
    if(Strategy == "Multilevel"){
      #Check suggested package
      if(!requireNamespace("imagerExtra", quietly = FALSE)) stop(
        paste0("imagerExtra CRAN package is required to execute the function. Please install using the following code: ",
               expression(install.packages("imagerExtra")))
      )

      #Check argument
      if(!all(is.numeric(Levels), Levels%%1 == 0, Levels > 2)) stop("Levels must be a positive integer value larger than 2")
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
      Argument_checker <- c(Consensus_reps_OK = (Consensus_reps >= 1 & Consensus_reps%%1 == 0),
                            Consensus_p_Items_OK = (Consensus_p_Items > 0 & Consensus_p_Items <= 1),
                            Consensus_Cluster_Alg_OK = Consensus_Cluster_Alg %in% c("hc", "pam", "km"),
                            Consensus_Distance_OK = Consensus_Distance %in% c("pearson", "spearman", "euclidean", "binary", "maximum", "canberra", "minkowski"),
                            Consensus_Name_OK = is.character(as.character(Consensus_Name))
      )
      Stop_messages <- c(Consensus_reps_OK = "Consensus_reps_OK must be an integer value > 0",
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
                            Consensus_reps_Meta_OK = (Consensus_reps_Meta >= 1 & Consensus_reps_Meta%%1 == 0),
                            Consensus_p_Items_Meta_OK = (Consensus_p_Items_Meta > 0 & Consensus_p_Items_Meta <= 1),
                            Consensus_Name_Meta_OK = is.character(as.character(Consensus_Name_Meta))
      )
      Stop_messages <- c(N_K_centroids_OK = "N_K_centroids must be smaller than the number of cells in DATA and a integer value > 0",
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
      Argument_checker <- c(N_initiations_OK = (N_initiations >= 1 & N_initiations%%1 == 0),
                            Max_iterations_OK = (Max_iterations%%1 == 0)
      )
      Stop_messages <- c(N_initiations_OK = "N_initiations must be an integer value > 0",
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
                            Max_iterations_km_OK = (Max_iterations_km >= 1 & Max_iterations_km%%1 == 0),
                            Max_iterations_em_OK = (Max_iterations_em >= 1 & Max_iterations_em%%1 == 0),
                            GMM_Distance_OK = GMM_Distance %in% c("eucl_dist", "maha_dist")
      )
      Stop_messages <- c(Quality_metric_OK = "Quality_metric must be one of the following: AIC, BIC",
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
                            Distance_CLARA_OK = Distance_CLARA %in% c("euclidean", "manhattan", "chebyshev", "canberra", "braycurtis",
                                                                      "pearson_correlation", "simple_matching_coefficient", "minkowski",
                                                                      "hamming", "jaccard_coefficient", "Rao_coefficient", "mahalanobis", "cosine"),
                            N_cores_OK = (N_cores >= 1 & N_cores%%1 == 0)
      )
      Stop_messages <- c(Samples_CLARA_OK = "Samples_CLARA must be an integer value > 0",
                         Sample_per_CLARA_OK = "Sample_per_CLARA must be a numeric value between 0 and 1",
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

    #ALL Strategies will converge in a single solution a DATA frame with a single NEW_Cluster column
    if(Strategy == "Arbitrary"){
      Arbitrary_cutoff <- sort(Arbitrary_cutoff)
      Min_value <- min(as.matrix(DATA_variables[-1]))
      Max_value <- max(as.matrix(DATA_variables[-1]))

      if(any(Arbitrary_cutoff <= Min_value)) stop(paste0("Arbitrary_cutoff value cannot be smaller than: ", Min_value))
      if(any(Arbitrary_cutoff >= Max_value)) stop(paste0("Arbitrary_cutoff value cannot be larger than: ", Max_value))

      Arbitrary_cutoff <- c(Min_value-1, Arbitrary_cutoff, Max_value+1)

      #Decide Label names
      if(length(Arbitrary_cutoff) == 3) Labels <- c("Negative", "Positive")
      if(length(Arbitrary_cutoff) == 4) Labels <- c("Low", "Mid", "High")
      if(length(Arbitrary_cutoff) > 4) Labels <-stringr::str_c("Level_", 0:(length(Arbitrary_cutoff)-2))

      #Cut every variable
      Labels_tibble <- suppressMessages(
        purrr::map_dfc(2:ncol(DATA_variables), function(Column){
          Label_2 <- as.character(cut(DATA_variables[[Column]], breaks = Arbitrary_cutoff, labels = Labels))
          stringr::str_c(names(DATA_variables)[Column], Label_2, sep = "_")
        })
      )
      #Generate the new Cluster variable
      DATA_variables <- DATA_variables %>% dplyr::mutate(NEW_Cluster = purrr::map_chr(1:nrow(Labels_tibble), function(Row) stringr::str_c(Labels_tibble[Row,], collapse = ".")))
    }
    if(Strategy == "Multilevel"){
      #Split every variable
      Labels_tibble <-
        suppressWarnings(
          suppressMessages(
            purrr::map_dfc(2:ncol(DATA_variables), function(Column){
              Label_2 <- as.double(imagerExtra::ThresholdML(imager::cimg(array(DATA_variables[[Column]], dim = c(1, length(DATA_variables[[Column]]), 1, 1))), k = (Levels-1)))
              stringr::str_c(names(DATA_variables)[Column], "Level", Label_2, sep = "_")
            })
          )
        )
      #Generate the new Cluster variable
      DATA_variables <- DATA_variables %>% dplyr::mutate(NEW_Cluster =purrr::map_chr(1:nrow(Labels_tibble), function(Row) stringr::str_c(Labels_tibble[Row,], collapse = ".")))
    }
    if(Strategy %in% c("Consensus_Clustering", "SOM", "Graph_Based", "K_Means_Meta_clustering", "Batch_K_means", "GMM", "CLARA_clustering")){

      DATA_variables <-
        CSM_Clustering_function(
          Original_data = DATA_variables,
          MARKERS = DATA_variables %>% dplyr::select(-Cell_no),

          Strategy = Strategy,
          Force_N_Clusters = Force_N_Clusters,

          Max_N_clusters_Consensus = N_Clusters,
          Consensus_reps = Consensus_reps,
          Consensus_p_Items = Consensus_p_Items ,
          Consensus_Cluster_Alg = Consensus_Cluster_Alg,
          Consensus_Distance = Consensus_Distance,
          Consensus_Name = Consensus_Name,

          Max_SOM_clusters = N_Clusters,

          Graph_type = Graph_type,
          Graph_Distance_method = Graph_Distance_method,
          Nearest_neighbors_for_graph = Nearest_neighbors_for_graph,
          Graph_Method = Graph_Method,
          Graph_Resolution = Graph_Resolution,
          N_steps = N_steps,

          N_K_centroids = N_K_centroids,
          Max_N_clusters_Meta = N_Clusters,
          Consensus_reps_Meta = Consensus_reps_Meta,
          Consensus_p_Items_Meta = Consensus_p_Items_Meta,
          Consensus_Name_Meta = Consensus_Name_Meta,

          Batch_size = Batch_size,
          Max_N_clusters_Batch = N_Clusters,
          N_initiations = N_initiations,
          Max_iterations = Max_iterations,

          Quality_metric = Quality_metric,
          Max_N_clusters_GMM = N_Clusters,
          Max_iterations_km = Max_iterations_km,
          Max_iterations_em = Max_iterations_em,
          GMM_Distance = GMM_Distance,

          Samples_CLARA = Samples_CLARA,
          Sample_per_CLARA = Sample_per_CLARA,
          Max_N_clusters_CLARA = N_Clusters,
          Distance_CLARA = Distance_CLARA,
          N_cores = N_cores
        )

      #Change the name of the column 'Cluster' for 'NEW_Cluster'
      DATA_variables <- DATA_variables %>% dplyr::rename("NEW_Cluster" = "Cluster")

      #Change the values of the NEW_Cluster column
      DATA_variables <- DATA_variables %>% dplyr::mutate(NEW_Cluster = stringr::str_c("Cluster", DATA_variables$NEW_Cluster, sep = "_"))
    }


    #Plot the results
    Mean_tibble <- DATA_variables %>% dplyr::select(-Cell_no) %>% group_by(NEW_Cluster) %>% summarize_all(mean) %>%dplyr::ungroup() #Obtain mean tibble
    Mean_matrix <- as.matrix(Mean_tibble[-1] %>% scale()) #Scale it and transform it into a  mtrix
    row.names(Mean_matrix) <- Mean_tibble[[1]]

    plot(
      ComplexHeatmap::Heatmap(Mean_matrix,
                              name = "Scaled")
    )

    Answer <- menu(choices = c("Proceed", "Abort"), title = "Are the results adequate?")
    if(Answer == 2) stop("Re-clustering has been aborted")

    #If all OK then bind the results with the cell phenotype in origin
    DATA <- dplyr::left_join(DATA, DATA_variables %>% dplyr::select(Cell_no, NEW_Cluster), by = "Cell_no")

    #If NA then NEW cluster does not add any info
    DATA$NEW_Cluster[is.na(DATA$NEW_Cluster)] <- ""
    DATA <- DATA %>% dplyr::mutate(Phenotype = stringr::str_c(Phenotype, NEW_Cluster, sep = "_")) %>% dplyr::select(-NEW_Cluster)
    DATA$Phenotype <- sub("_$", "", DATA$Phenotype)

    #Return the final data
    return(DATA)
  }
