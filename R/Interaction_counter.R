#' Analyzes the spatial interaction pattern between all cell types of all the images in the dataset
#'
#' The function calculates the cell to cell spatial interaction pattern between cell type pairs. It is based on the imcRtools::countInteractions function.
#' Spatial interactions are quantified on graph representation of datasets and can be computed in different ways. Afterwards, recurrent spatial interaction patterns
#' can be identified in the dataset by running the [Interaction_analyzer()] function.
#'
#'
#' @param DATA A dataframe or tibble containing a column named 'Phenotype' containing cell phenotype labels.
#' @param Phenotypes_included A character vector indicating the phenotype labels that will be included in the analysis.
#' @param N_cores Integer. Number of cores to parallelize your computation.
#'
#' @param Graph_type A character value indicating how to build the cell graph. One of the following: "expansion", "knn" or "delaunay".
#' @param K_number If Graph_type is knn, an integer value indicating the number of K nearest neighbors to build the graph.
#' @param Dist_threshold If Graph_type is expansion, a numeric value indicating the distance threshold.
#'
#' @param Method A character value indicating the approach to spatial interaction quantification. One of the following: "classic", "histocat" or "patch" (see details).
#' @param patch_size If Method is "patch", a integer value indicating the size of the patch.
#'
#' @param Perform_significance_testing A logical value indicating if significance testing should be performed.
#' @param N_iterations If significance testing needs to be performed, the number of iterations to calculate p values.
#' @param p_threshold If significance testing needs to be performed, the p value threshold for significance.
#'
#' @details
#' For additional information please visit imcRtools Vignette: https://www.bioconductor.org/packages/release/bioc/vignettes/imcRtools/inst/doc/imcRtools.html#5_Spatial_analysis.
#'
#'
#' @seealso [Interaction_analyzer()]
#'
#' @returns A list containing pairwise cell spatial interactions by sample.
#'
#' @examples
#' \dontrun{
#' Interaction_counter(
#'     DATA = CSM_Phenotypecell_test,
#'     Phenotypes_included = unique(CSM_Phenotypecell_test$Phenotype),
#'     N_cores = 1,
#'
#'     Graph_type = "knn",
#'     K_number = 20
#'
#'     Method = "classic",
#'     patch_size = 125,
#'
#'     Perform_significance_testing = TRUE,
#'     N_iterations = 10,
#'     p_threshold = 0.05
#')
#' }
#'
#' @export

Interaction_counter <-
  function(DATA = NULL,
           Phenotypes_included = NULL,
           N_cores = NULL,

           Graph_type = NULL,
           K_number = NULL,
           Dist_threshold = NULL,

           Method = NULL,
           patch_size = NULL,

           Perform_significance_testing = NULL,
           N_iterations = NULL,
           p_threshold = NULL){

    #Check suggested packages
    {
      if(!requireNamespace("imcRtools", quietly = TRUE)) stop(
        paste0("imcRtools Bioconductor package is required to execute the function. Please install using the following code: ",
               expression({
                 if (!require("BiocManager", quietly = TRUE))
                   install.packages("BiocManager")

                 BiocManager::install("imcRtools")
               })
        )
      )
      if(!requireNamespace("SpatialExperiment", quietly = TRUE)) stop(
        paste0("SpatialExperiment Bioconductor package is required to execute the function. Please install using the following code: ",
               expression({
                 if (!require("BiocManager", quietly = TRUE))
                   install.packages("BiocManager")

                 BiocManager::install("SpatialExperiment")
               })
        )
      )
    }

    #Check arguments
    if(!all(N_cores >= 1 & N_cores%%1 == 0)) stop("N_cores must be an integer value > 0")
    #Import data
    DATA <- DATA

    if(!all(Phenotypes_included %in% unique(DATA$Phenotype))){
      stop(paste0("Phenotypes provided must be one of: ", stringr::str_c(unique(DATA$Phenotype), collapse = ", ")))
    }
    if(!is.logical(Perform_significance_testing)) stop("Perform_significance_testing must be a logical value")
    if(!Graph_type %in% c("expansion", "knn", "delaunay")) stop("Graph_type must be one of the following: expansion, knn, delaunay")
    if(!all(is.numeric(K_number), K_number%%1 == 0, K_number >= 1)) stop("K_numer must be an integer value > 0")
    if(Graph_type == "expansion"){
      if(!all(is.numeric(Dist_threshold), Dist_threshold > 0)) stop("Dist_threshold must be a numeric value > 0")
    }
    if(!Method %in% c("classic", "histocat", "patch")) stop("Method must be one of the following: classic, histocat, patch")
    if(!all(is.numeric(patch_size), patch_size > 0)) stop("patch_size must be a numeric value > 0")
    if(Perform_significance_testing){
      if(!all(is.numeric(N_iterations), N_iterations%%1 == 0, N_iterations > 0)) stop("N_iterations must be a integer value > 0")
      if(!all(is.numeric(p_threshold), p_threshold <= 1, p_threshold > 0)) stop("p_threshold must be a numeric value between 0 and 1")
    }

    DATA <- DATA %>% dplyr::filter(Phenotype %in% Phenotypes_included)

    #save exit function if parallelization fails
    on.exit({
      future::plan("future::sequential")
      gc()
    })

    #Generate the cluster
    future::plan("future::multisession", workers = N_cores)
    options(future.globals.maxSize = Inf, future.rng.onMisuse = "ignore")
    furrr::furrr_options(scheduling = Inf)

    Final_list <-
      suppressMessages(
        furrr::future_map(unique(DATA$Subject_Names), function(Image){
          #Get the Image data
          Interim <- DATA %>% dplyr::filter(Subject_Names == Image)

          #First we create a SpatialExperiment object before we can proceed
          Spatial_object <- SpatialExperiment::SpatialExperiment(sample_id = Image,
                                                                 colData = Interim,
                                                                 spatialCoordsNames = names(Interim)[2:3]
          )

          #Now we build the graph to the calculate interactions
          Graph <- imcRtools::buildSpatialGraph(Spatial_object,
                                                img_id = "Subject_Names",
                                                type = Graph_type,
                                                directed = TRUE,
                                                threshold = Dist_threshold,
                                                name = "spatialcontext_graph",
                                                k = K_number,
                                                coords = names(Interim)[2:3]
          )

          #Now we calculate the output

          if(!Perform_significance_testing){
            out <- imcRtools::countInteractions(
              Graph,
              group_by = "Subject_Names",
              label = "Phenotype",
              colPairName = "spatialcontext_graph",
              method =  Method,
              patch_size = patch_size
            )
            out <- as_tibble(out)
            return(out)
          }

          if(Perform_significance_testing){
            out <-
              imcRtools::testInteractions(Graph,
                                          group_by = "Subject_Names",
                                          label = "Phenotype",
                                          method = Method,
                                          colPairName = "spatialcontext_graph",
                                          BPPARAM = BiocParallel::SerialParam(RNGseed = 123),
                                          iter = N_iterations,
                                          p_threshold = p_threshold)

            out <- as_tibble(out)
            return(out)
          }

        },
        .progress = TRUE)
      )

    future::plan("future::sequential")
    gc()

    names(Final_list) <- unique(DATA$Subject_Names)
    return(Final_list)
  }
