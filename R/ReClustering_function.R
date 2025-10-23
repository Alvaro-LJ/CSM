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
#' @param Strategy One of the following: "Arbitrary", "Multilevel", "SOM", "Batch_K_means", "GMM" (see details).
#'
#' @param N_Clusters For SOM, Batched K-means and GMM, an integer indicating the number of clusters to generate.
#' @param Levels For Multilevel, an integer indicating the number of levels to divide each feature.
#' @param Arbitrary_cutoff For arbitrary, a numeric value or numeric vector indicating the cut-off values to divide each feature.
#'
#' @param Batch_size If strategy is Batch_K_means: Number of cells to be included in each random batch.
#' @param N_initiations If strategy is Batch_K_means: Number of times the algorithm is going to be tried to find the best clustering result.
#' @param Max_iterations If strategy is Batch_K_means: Max number of iterations in each try.
#'
#' @param Quality_metric If strategy is GMM:T he quality measure used to test the number of clusters ("AIC" or "BIC").
#' @param Max_iterations_km If strategy is GMM: Number of max iterations in the K means clustering performed.
#' @param Max_iterations_em If strategy is GMM: Number of max iterations in the Expectation Maximization algorithm.
#' @param GMM_Distance If strategy is GMM: Distance metric used in the model ("eucl_dist" or "maha_dist").
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
#' #OPTIONAL, set aside any features that will not be used in the initial phenotyping process----------------
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
#'#Perform Re-Clustering of a subset of cells according to set aside features----------------
#' ReClustering_function(
#'     DATA = DATA_Phenotypes,
#'     DATA_aside = DATA_list$Aside,
#'     Phenotype_variable = "Phenotype",
#'     Phenotype_to_recluster = "CD8",
#'     Variables_to_recluster = "GZMB_AVERAGE",
#'     Strategy = "Multilevel",
#'     Levels = 3
#'     )
#'
#' @export

ReClustering_function <-
  function(DATA = NULL,
           DATA_aside = NULL,
           Phenotype_variable = NULL,
           Phenotype_to_recluster = NULL,
           Variables_to_recluster = NULL,
           Strategy = NULL,

           N_Clusters = NULL,

           Levels = NULL,
           Arbitrary_cutoff = NULL,

           #Parameters for Batched K means
           Batch_size = NULL,
           N_initiations = NULL,
           Max_iterations = NULL,

           #Parameters for Gaussian Mixture Model
           Quality_metric = NULL,
           Max_iterations_km = NULL,
           Max_iterations_em = NULL,
           GMM_Distance = NULL
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
      Problematic_variable <- names(DATA_variables[-1])[!map_lgl(DATA_variables[-1], function(Column) is.numeric(Column))]
      stop(paste0("The following variables are not numeric and cannot be included in the reclustering process: ", stringr::str_c(Problematic_variable, collapse = ", ")))
    }

    #Check the strategy
    if(!Strategy %in% c("Arbitrary", "Multilevel", "SOM", "Batch_K_means", "GMM")) stop(paste0("Strategy must be one of the following: ", stringr::str_c(c("Arbitrary", "Multilevel", "Consensus_Clustering", "SOM", "Batch_K_means", "GMM"), collapse = ", ")))
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
    #Check arguments for SOM
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
      #Check argument
      if(!all(is.numeric(N_Clusters), N_Clusters%%1 == 0, N_Clusters >= 2)) stop("N_clusters must be a positive integer value > 1")
    }
    #Check arguments for Batch_K_means
    if(Strategy == "Batch_K_means"){
      #Check suggested packages
      if(!requireNamespace("ClusterR", quietly = FALSE)) stop(
        paste0("ClusterR CRAN package is required to execute the function. Please install using the following code: ",
               expression(install.packages("ClusterR")))
      )
      #Check arguments
      if(!all(is.numeric(N_Clusters), N_Clusters%%1 == 0, N_Clusters >= 2)) stop("N_clusters must be a positive integer value > 1")
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
    #Check arguments for GMM
    if(Strategy == "GMM"){
      #Check suggested packages
      if(!requireNamespace("ClusterR", quietly = FALSE)) stop(
        paste0("ClusterR CRAN package is required to execute the function. Please install using the following code: ",
               expression(install.packages("ClusterR")))
      )
      #Check arguments
      if(!all(is.numeric(N_Clusters), N_Clusters%%1 == 0, N_Clusters >= 2)) stop("N_clusters must be a positive integer value > 1")
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
      DATA_variables <- DATA_variables %>% dplyr::mutate(NEW_Cluster =purrr::map_chr(1:nrow(Labels_tibble), function(Row) stringr::str_c(Labels_tibble[Row,], collapse = ".")))
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
    if(Strategy == "SOM"){
      #Transform data into a scaled matrix and perform Self Organizing Map
      SOM_results <- try(
        FlowSOM::FlowSOM(DATA_variables %>% dplyr::select(-Cell_no) %>% scale() %>% as.matrix(),
                         scale = F,
                         colsToUse = names(DATA_variables)[-1],
                         nClus = N_Clusters,
                         silent = F,
                         seed = 21)
      )

      if(berryFunctions::is.error(SOM_results)) {
        stop("Number of clusters is too low. Please try increasing the number of nclus")
      }
      else{
        #Assign phenotypes to each cell
        DATA_variables <- DATA_variables %>% dplyr::mutate(NEW_Cluster =stringr::str_c("Cluster", FlowSOM::GetMetaclusters(SOM_results), sep = "_"))
      }
    }
    if(Strategy == "Batch_K_means"){
      if(nrow(DATA_variables) < Batch_size) stop(paste0("Batch_size must be smaller than ", nrow(DATA_variables)))
      Batch_k_means <- ClusterR::MiniBatchKmeans(DATA_variables %>% dplyr::select(-Cell_no) %>% scale(),
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

      #Assign the cluster to each observation of MARKER
      pr_mb <- suppressWarnings(predict(object = Batch_k_means, fuzzy = F, newdata = DATA_variables %>% dplyr::select(-Cell_no) %>% scale()))
      pr_mb <- as_tibble(pr_mb)
      names(pr_mb) <- "NEW_Cluster"
      pr_mb$NEW_Cluster <-stringr::str_c("Cluster", pr_mb$NEW_Cluster, sep = "_")

      #Generate the data phenotypes tibble
      DATA_variables <-dplyr::bind_cols(DATA_variables, pr_mb)
    }
    if(Strategy == "GMM"){
      GMM_model <- ClusterR::GMM(DATA_variables %>% dplyr::select(-Cell_no) %>% scale(),
                                 gaussian_comps = as.double(N_Clusters),
                                 dist_mode = GMM_Distance,
                                 seed_mode = "random_subset",
                                 km_iter = Max_iterations_km,
                                 em_iter = Max_iterations_em,
                                 verbose = TRUE,
                                 var_floor = 1e-10,
                                 full_covariance_matrices = FALSE
      )

      pr_mb <- suppressWarnings(predict(object = GMM_model, fuzzy = F, newdata = DATA_variables %>% dplyr::select(-Cell_no) %>% scale()))
      pr_mb <- as_tibble(pr_mb)
      names(pr_mb) <- "NEW_Cluster"
      pr_mb$NEW_Cluster <-stringr::str_c("Cluster", pr_mb$NEW_Cluster, sep = "_")

      #Generate the data phenotypes tibble
      DATA_variables <-dplyr::bind_cols(DATA_variables, pr_mb)
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
    DATA <-dplyr::left_join(DATA, DATA_variables %>% dplyr::select(Cell_no, NEW_Cluster), by = "Cell_no")

    #If NA then NEW cluster does not add any info
    DATA$NEW_Cluster[is.na(DATA$NEW_Cluster)] <- ""
    DATA <- DATA %>% dplyr::mutate(Phenotype = stringr::str_c(Phenotype, NEW_Cluster, sep = "_")) %>% dplyr::select(-NEW_Cluster)
    DATA$Phenotype <- sub("_$", "", DATA$Phenotype)

    #Return the final data
    return(DATA)
  }
