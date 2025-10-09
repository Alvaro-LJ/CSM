#' Summarize images according to their cell-to-cell pairwise spatial interaction patterns
#'
#' The function summarizes the results of the the pairwise cell-to-cell spatial interaction. In addition, the function identifies recurrent patterns of cell-to-cell spatial interactions in images. Before running the function pairwise spatial interactions need to be computed using [Interaction_counter()].
#'
#'
#' @param DATA A list of pairwise cell-to-cell spatial interactions computed using [Interaction_counter()].
#' @param Exclude_non_significant A logical value indicating if non-significant spatial interactions should be removed from the analysis.
#'
#' @param Cluster A logical value indicating if clustering should be performed.
#' @param Max_N_Clusters Number of maximum clusters that can be identified.
#' @param Consensus_reps Number of iterations to converge.
#' @param Consensus_p_Items Percentage of features that you desire to sample in each iteration.
#' @param Consensus_Cluster_Alg Clustering algorithm to be used (’hc’ hierarchical (hclust), ’pam’ for paritioning around medoids, ’km’ for k-means).
#' @param Consensus_Distance Distance metric to be used (pearson(1 - Pearson correlation), spearman(1 - Spearman correlation), euclidean, binary, maximum, canberra, minkowski.
#' @param Consensus_Name Name of the folder that is going to be created in order to place the resulting graphs.
#'
#' @seealso [Interaction_counter()]
#'
#' @returns A list containing a summary of the pairwise cell-to-cell interactions by sample.
#'
#' @export

Interaction_analyzer <-
  function(DATA = NULL,
           Exclude_non_significant = NULL,

           Cluster = TRUE,
           Max_N_Clusters = NULL,
           Consensus_reps = NULL,
           Consensus_p_Items = NULL,
           Consensus_Cluster_Alg = NULL,
           Consensus_Distance = NULL,
           Consensus_Name = NULL){

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
    if(!all(purrr::map_lgl(DATA, function(x) all(c("from_label", "to_label") %in% names(x))))) stop("DATA must be created using the function Interaction_counter")
    if(!is.logical(Exclude_non_significant)) stop("Exclude_non_significant must be a logical value")
    if(!is.logical(Cluster)) stop("Cluster must be a logical value")
    if(Cluster){
      #Check arguments of consensus clustering
      Argument_checker <- c(Max_N_Clusters_OK = (Max_N_Clusters >= 2 & Max_N_Clusters%%1 == 0),
                            Consensus_reps_OK = (Consensus_reps >= 1 & Consensus_reps%%1 == 0),
                            Consensus_p_Items_OK = (Consensus_p_Items > 0 & Consensus_p_Items <= 1),
                            Consensus_Cluster_Alg_OK = Consensus_Cluster_Alg %in% c("hc", "pam", "km"),
                            Consensus_Distance_OK = Consensus_Distance %in% c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski"),
                            Consensus_Name_OK = is.character(as.character(Consensus_Name))
      )
      Stop_messages <- c(Max_N_Clusters_OK = "Max_N_Clusters must be an integer value > 1",
                         Consensus_reps_OK = "Consensus_reps_OK must be an integer value > 0",
                         Consensus_p_Items_OK = "Consensus_p_Items must be a numeric value > 0 and lower than 1",
                         Consensus_Cluster_Alg_OK = "Consensus_Cluster_Alg must be one of the following: hc, pam, km",
                         Consensus_Distance_OK = "Consensus_Distance must be one the following: euclidean, maximum, manhattan, canberra, binary, minkowski",
                         Consensus_Name_OK = "Consensus_Name must ve a character value")
      #Check arguments and stop if necessary
      if(!all(Argument_checker)){
        stop(cat(Stop_messages[!Argument_checker],
                 fill = sum(!Argument_checker)))
      }
    }



    #If significance is not an issue proceed with analysis
    if(!Exclude_non_significant){

      print("Obtaining data")
      Interaction_DF <-purrr::map_dfr(DATA, function(Image){
        Image %>% dplyr::mutate(col_id = stringr::str_c(from_label, to_label, sep = "_")) %>% pivot_wider(id_cols = group_by, names_from = col_id, values_from = ct) %>%
          dplyr::rename("Subject_Names" = "group_by")
      })

      #If no clustering is required return the actual data
      if(!Cluster){
        return(Interaction_DF)
      }

      #If Clustering is required then proceed with clustering
      if(Cluster){
        print("Performing clustering")
        #We will remove the non-clusteable samples using a while loop
        Na_count <- sum(is.na(dist(Interaction_DF[-1], Consensus_Distance))) #Generate the entering condition
        while(Na_count >= 1){
          Na_tibble <- as_tibble(as.matrix(dist(Interaction_DF[-1], "euclidean", diag = T, upper = T)))
          warning(paste0("Some samples had insufficient data. The following samples will be removedd: ",
                         stringr::str_c(Interaction_DF[[1]][
                           which(purrr::map_dbl(Na_tibble, function(Column) sum(is.na(Column))) == max(purrr::map_dbl(Na_tibble, function(Column) sum(is.na(Column)))))],
                           collapse = ", "))) #Print a warning message everytime samples are removed
          Interaction_DF <- Interaction_DF[
            -which(purrr::map_dbl(Na_tibble, function(Column) sum(is.na(Column))) == max(purrr::map_dbl(Na_tibble, function(Column) sum(is.na(Column))))),] #Generate the output tibble
          Na_count <- sum(is.na(dist(Interaction_DF[-1], Consensus_Distance)))#Update the NA counter
        }

        #We need to calculate our own distance matrix to account for potential NA values
        for_dist <- Interaction_DF %>% dplyr::select(-Subject_Names) %>% scale()
        Distance_matrix <- dist(for_dist, Consensus_Distance)
        #Perform consensus clustering
        Clustering_result <- ConsensusClusterPlus::ConsensusClusterPlus(Distance_matrix,
                                                                        maxK = Max_N_Clusters,
                                                                        reps = Consensus_reps,
                                                                        pItem = Consensus_p_Items,
                                                                        pFeature = 1, #For a distance matrix we need to always use 100% of features
                                                                        title = Consensus_Name,
                                                                        clusterAlg = Consensus_Cluster_Alg,
                                                                        plot = "png",
                                                                        verbose = T)

        #Make the user decide the number of neighborhoods according to results
        N_Clusters <- menu(choices = as.character(1:Max_N_Clusters), title = paste0("Check the results at: ", getwd(), ". Then decide the appropiate number of Clusters"))
        Interaction_DF <- Interaction_DF %>%dplyr::mutate(Cluster_assignment = as.character(Clustering_result[[as.double(N_Clusters)]][["consensusClass"]]))

        #Visualize the heatmap of mean by Cluster
        print("Generating Heatmap")
        Mean_tibble <- Interaction_DF %>% dplyr::select(-Subject_Names) %>% dplyr::group_by(Cluster_assignment) %>% dplyr::summarize_all(.funs = function(x) mean(x, na.rm = T)) %>%
          dplyr::ungroup()  #Obtain mean tibble
        Mean_matrix <- as.matrix(Mean_tibble[-1] %>% scale()) #Scale it and transform it into a  mtrix
        row.names(Mean_matrix) <- Mean_tibble[[1]]

        plot(ComplexHeatmap::Heatmap(Mean_matrix,
                                     name = "Scaled")
        )
        #Print the cluster count by sample
        print(Interaction_DF %>% dplyr::count(Cluster_assignment))

        #Return the final data
        return(Interaction_DF)
      }
    }

    #If non significant metrics are to be excluded execute the following
    if(Exclude_non_significant){
      print("Obtaining data")
      Interaction_DF <-purrr::map_dfr(DATA, function(Image){
        Image %>% dplyr::mutate(col_id = stringr::str_c(from_label, to_label, sep = "_")) %>%
          dplyr::filter(sig) %>%
          dplyr::pivot_wider(id_cols = group_by, names_from = col_id, values_from = ct) %>%
          dplyr::rename("Subject_Names" = "group_by")
      })
      #If no clustering is required return the actual data
      if(!Cluster){
        return(Interaction_DF)
      }

      #If Clustering is required then proceed with clustering
      if(Cluster){
        print("Performing clustering")
        #We will remove the non-clusteable samples using a while loop
        Na_count <- sum(is.na(dist(Interaction_DF[-1], Consensus_Distance))) #Generate the entering condition
        while(Na_count >= 1){
          Na_tibble <- as_tibble(as.matrix(dist(Interaction_DF[-1], Consensus_Distance, diag = T, upper = T)))
          warning(paste0("Some samples had insufficient data. The following samples will be removedd: ",
                         stringr::str_c(Interaction_DF[[1]][
                           which(purrr::map_dbl(Na_tibble, function(Column) sum(is.na(Column))) == max(purrr::map_dbl(Na_tibble, function(Column) sum(is.na(Column)))))],
                           collapse = ", "))) #Print a warning message everytime samples are removed
          Interaction_DF <- Interaction_DF[
            -which(purrr::map_dbl(Na_tibble, function(Column) sum(is.na(Column))) == max(purrr::map_dbl(Na_tibble, function(Column) sum(is.na(Column))))),] #Generate the output tibble
          Na_count <- sum(is.na(dist(Interaction_DF[-1], Consensus_Distance)))#Update the NA counter
        }

        #We need to calculate our own distance matrix to account for potential NA values
        for_dist <- Interaction_DF %>% dplyr::select(-Subject_Names) %>% scale()
        Distance_matrix <- dist(for_dist, Consensus_Distance)

        #Perform consensus clustering
        Clustering_result <- ConsensusClusterPlus::ConsensusClusterPlus(Distance_matrix,
                                                                        maxK = Max_N_Clusters,
                                                                        reps = Consensus_reps,
                                                                        pItem = Consensus_p_Items,
                                                                        pFeature = 1, #For a distance matrix we need to always use 100% of features
                                                                        title = Consensus_Name,
                                                                        clusterAlg = Consensus_Cluster_Alg,
                                                                        plot = "png",
                                                                        verbose = T)

        #Make the user decide the number of neighborhoods according to results
        N_Clusters <- menu(choices = as.character(1:Max_N_Clusters), title = paste0("Check the results at: ", getwd(), ". Then decide the appropiate number of Clusters"))
        Interaction_DF <- Interaction_DF %>%dplyr::mutate(Cluster_assignment = as.character(Clustering_result[[as.double(N_Clusters)]][["consensusClass"]]))

        #Visualize the heatmap of mean by Cluster
        Mean_tibble <- Interaction_DF %>% dplyr::select(-Subject_Names) %>% dplyr::group_by(Cluster_assignment) %>% dplyr::summarize_all(.funs = function(x) mean(x, na.rm = T)) %>%
          dplyr::ungroup() #Obtain mean tibble
        Mean_matrix <- as.matrix(Mean_tibble[-1] %>% scale()) #Scale it and transform it into a  mtrix
        row.names(Mean_matrix) <- Mean_tibble[[1]]

        plot(ComplexHeatmap::Heatmap(Mean_matrix,
                                     name = "Scaled")
        )
        #Print the cluster count by sample
        print(Interaction_DF %>% dplyr::count(Cluster_assignment))

        #Return the final data
        return(Interaction_DF)
      }

    }
  }
