#' Identifies cell neighborhoods using the pipeline present in the SPIAT package
#'
#' The function identifies cell neighborhoods using the SPIAT::identify_neighborhoods function.
#'
#' @param N_cores Integer. Number of cores to parallelize your computation.
#' @param DATA_SPIAT A list containing SPIAT objects generated using [SPIAT_object_generator()].
#' @param Strategy A character value indicating the method to find the clusters. One of "hierarchical" or "dbscan".
#' @param Cell_types_of_interest A character vector indicating the cell phenotype labels included in the analysis.
#' @param Radius A numeric value indicating the distance radius used to calculate neighborhoods.
#' @param Min_neighborhood_size A integer value indicating the minimum number of cells to constitute a neighborhood.
#' @param No_Phenotype_name (OPTIONAL) A character value indicating the cell phenotype label corresponding to non-phenotyped cells.
#'
#' @returns A tibble containing a summary of the analysis result by image.
#'
#' @seealso [SPIAT_object_generator()].
#' @export

SPIAT_neighborhood_identifier <-
  function(N_cores = NULL,
           DATA_SPIAT = NULL,
           Strategy = NULL,
           Cell_types_of_interest = NULL,
           Radius = NULL,
           Min_neighborhood_size = NULL,
           K_phenograph = NULL, #Still under development
           No_Phenotype_name = NULL) {

    #Check suggested packages
    if(!requireNamespace("SPIAT", quietly = TRUE)) stop(
      paste0("SPIAT Bioconductor package is required to execute the function. Please install using the following code: ",
             expression({
               if (!require("BiocManager", quietly = TRUE))
                 install.packages("BiocManager")

               BiocManager::install("SPIAT")
             })
      )
    )

    #Check arguments
    if(!all(N_cores >= 1 & N_cores%%1 == 0)) stop("N_cores must be an integer value > 0")

    if(!exists(DATA_SPIAT, envir = .GlobalEnv)) stop("DATA_SPIAT must be an existing object")
    DATA_SPIAT <- get(DATA_SPIAT, envir = .GlobalEnv)
    if(!all(purrr::map_lgl(DATA_SPIAT, function(Image) class(Image) == "SpatialExperiment"))) stop("DATA_SPIAT must be created using the SPIAT_object_generator function")

    if(!Strategy %in% c("hierarchical", "dbscan")) stop("Strategy must be one of the following: hierarchical, dbscan")
    if(!all(Cell_types_of_interest %in% unique(unlist(purrr::map(DATA_SPIAT, function(x) x@colData@listData$Phenotype))))) {
      Absent_cell_types <- Cell_types_of_interest[!Cell_types_of_interest %in% unique(unlist(purrr::map(DATA_SPIAT, function(x) x@colData@listData$Phenotype)))]
      stop(paste0(stringr::str_c(Absent_cell_types, collapse = ", "), " not found in DATA_SPIAT"))
    }
    if(!all(is.numeric(Radius), Radius > 0)){
      stop("Radius must be a numeric value > 0")
    }
    if(!all(is.numeric(Min_neighborhood_size), Min_neighborhood_size%%1 == 0, Min_neighborhood_size > 0)) stop("Min_neighborhood_size must be a integer value > 0")
    if(!any(is.null(No_Phenotype_name), is.character(No_Phenotype_name))) stop("No_Phenotype_name must be Null or a character value")


    #save exit function if parallelization fails
    on.exit({
      future::plan("future::sequential")
      gc()
    })
    #Proceed with analysis
    future::plan("future::multisession", workers = N_cores)
    options(future.globals.maxSize = Inf, future.rng.onMisuse = "ignore")
    furrr::furrr_options(scheduling = Inf)

    RESULTS <-
      suppressMessages(
        furrr::future_map(DATA_SPIAT, function(Image){
          library(SPIAT)
          #First calculate the neighborhood
          SPIAT_neighborhoods <-  identify_neighborhoods(
            Image,
            method = Strategy, #Method to be implemented
            cell_types_of_interest = Cell_types_of_interest, #Cell types of interest to find the neighbors
            radius = Radius, #Distance threshold to consider two cells interacting (required for hierarchical and DBSCAN)
            min_neighborhood_size = Min_neighborhood_size, #Minimum size of neighborhoods (required for hierarchical and DBSCAN)
            k = K_phenograph, #A parameter required by phenographr
            feature_colname = "Phenotype",
            no_pheno = No_Phenotype_name #Name for the cells without a phenotype
          )

          #We analyze neighborhood composition
          SPIAT_neighborhoods_Composition <- composition_of_neighborhoods(SPIAT_neighborhoods,
                                                                          feature_colname = "Phenotype")

          #We plot neighborhood composition
          plot_composition_heatmap(SPIAT_neighborhoods_Composition,
                                   feature_colname = "Phenotype",
                                   log_values = T)

          #We calculate the Average_Nearest_Neighbour_Index ANN_index by each cluster
          Cluster_Interaction_analysis <-purrr::map_dfr(str_subset(unique(SPIAT_neighborhoods$Neighborhood), "Cluster"), function(Cluster) {
            Results <- average_nearest_neighbor_index(SPIAT_neighborhoods,
                                                      reference_celltypes = Cluster,
                                                      feature_colname = "Neighborhood",
                                                      p_val = 0.05 #select threshold p value
            )
            return(c(Cluster = Cluster,
                     ANN_index = Results$ANN_index,
                     Pattern = Results$pattern,
                     p_val = Results$`p-value`))
          })

          #We list our results
          return(list(Neighborhoods = SPIAT_neighborhoods,
                      Neighborhoods_Composition = SPIAT_neighborhoods_Composition,
                      Neighborhoods_Composition_plot = plot_composition_heatmap(SPIAT_neighborhoods_Composition,
                                                                                feature_colname = "Phenotype",
                                                                                log_values = T),
                      Average_Nearest_Neighbour_Index <- Cluster_Interaction_analysis
          ))
        },
        .progress = TRUE)
      )
    future::plan("future::sequential")
    gc()

    return(RESULTS)
  }
