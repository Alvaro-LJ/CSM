#' Analyzes the spatial heterogeneity following the pipeline in the SPIAT package
#'
#' The function divides the images into tiles, calculates entropy within tile and then counts the number of tiles above a user defined threshold.
#' Results may be similar to those obtained using [Tiled_image_heterogeneity_analyzer()].
#'
#'
#' @param N_cores Integer. Number of cores to parallelize your computation.
#' @param DATA_SPIAT A list containing SPIAT objects generated using [SPIAT_object_generator()].
#' @param DATA_Phenotypes A dataframe or tibble containing cell feature data with a column named 'Phenotype' containing cell phenotype labels.
#' @param Tile_size A numeric value indicating the size of the tile.
#' @param Phenotypes_included A character vector indicating the phenotype labels that will be included in the analysis.
#' @param Autocorrelation_metric Either "globalmoran" or "GearyC".
#' @param Entropy_threshold A numeric value indicating the the arbitrary threshold to be used.
#'
#' @returns A tibble containing a summary of the analysis result by image.
#'
#' @seealso [SPIAT_object_generator()], [Tiled_image_heterogeneity_analyzer()].
#'
#'
#' @examples
#' #Generate SPIAT object list----------------------------
#' DATA_SPIAT <-
#' SPIAT_object_generator(
#'     DATA_Intensities = CSM_Arrangedcellfeaturedata_test,
#'     DATA_Phenotypes = CSM_Phenotypecell_test
#' )
#'
#' #Analyze heterogeneity--------------------------------
#' SPIAT_Heterogeneity_Analyzer(
#'     N_cores = 2,
#'     DATA_SPIAT = DATA_SPIAT,
#'     DATA_Phenotypes = CSM_Phenotypecell_test,
#'     Tile_size = 125,
#'     Phenotypes_included = c("TUMOR", "CD8_GZMBneg", "CD8_GZMBpos", "OTHER"),
#'     Entropy_threshold = 0.5,
#'     Autocorrelation_metric = "globalmoran"
#' )
#'
#'
#' @export

SPIAT_Heterogeneity_Analyzer <-
  function(N_cores = NULL,
           DATA_SPIAT = NULL,
           DATA_Phenotypes = NULL,
           Tile_size = NULL,
           Phenotypes_included = NULL,
           Autocorrelation_metric = NULL,
           Entropy_threshold = NULL) {

    #Check suggested packages
    {
      if(!requireNamespace("SPIAT", quietly = TRUE)) stop(
        paste0("SPIAT Bioconductor package is required to execute the function. Please install using the following code: ",
               expression({
                 if (!require("BiocManager", quietly = TRUE))
                   install.packages("BiocManager")

                 BiocManager::install("SPIAT")
               })
        )
      )
      if(!requireNamespace("elsa", quietly = FALSE)) stop(
        paste0("elsa CRAN package is required to execute the function. Please install using the following code: ",
               expression(install.packages("elsa")))
      )
    }


    #Check arguments by generating a argument check vector and message vector
    Argument_checker <- c(N_cores_OK = (N_cores >= 1 & N_cores%%1 == 0),
                          DATA_SPIAT_OK = all(purrr::map_lgl(DATA_SPIAT, function(Image) class(Image) == "SpatialExperiment")),
                          DATA_Phenotypes_OK = "Phenotype" %in% names(DATA_Phenotypes),
                          Tile_size_OK = all(is.numeric(Tile_size), Tile_size > 0),
                          Phenotypes_included_OK = all(Phenotypes_included %in% DATA_Phenotypes$Phenotype),
                          Autocorrelation_metric_OK = Autocorrelation_metric %in% c("globalmoran", "GearyC"),
                          Entropy_threshold_OK = is.numeric(Entropy_threshold)
    )

    Stop_messages <- c(N_cores_OK = "N_cores must be an integer value > 0",
                       DATA_SPIAT_OK = "DATA_SPIAT must be generated with the SPIAT_object_generator function",
                       DATA_Phenotypes_OK = "DATA_Phenotypes must contain a Phenotype object",
                       Tile_size_OK = "Tile_size must be a numeric value > 0",
                       Phenotypes_included_OK =stringr::str_c("Phenotypes must be any of the following: ", stringr::str_c(unique(DATA_Phenotypes$Phenotype), collapse = ", ")),
                       Autocorrelation_metric_OK = "Autocorrelation_metric must be one of the following: globalmoran, GearyC",
                       Entropy_threshold_OK = "Entropy_threshold must be a numeric value")

    #Check arguments and stop if necessary
    if(!all(Argument_checker)){
      stop(cat(Stop_messages[!Argument_checker],
               fill = sum(!Argument_checker)))
    }

    #Import list containing SPIAT OBJECTS
    DATA_SPIAT <- DATA_SPIAT

    #Generate list containing Data_phenotypes objects
    DATA_Phenotypes <-purrr::map(names(DATA_SPIAT), function(Image) DATA_Phenotypes %>% dplyr::filter(Subject_Names == Image))

    #save exit function if parallelization fails
    on.exit({
      future::plan("future::sequential")
      gc()
    })

    #prepare the cluster
    future::plan("future::multisession", workers = N_cores)
    options(future.globals.maxSize = Inf, future.rng.onMisuse = "ignore")
    furrr::furrr_options(scheduling = Inf)

    RESULTS <-
      furrr::future_map(seq_along(1:length(DATA_SPIAT)), function(Index){
        suppressPackageStartupMessages(library(SPIAT))
        #Define functions inside the cores
        Own_grid_metrics <-
          function(spe_object, FUN, n_split, ...) {
            split <- SPIAT::image_splitter(spe_object, n_split)
            list.metric <- list()
            for (i in seq_len(length(split))) {
              spe <- split[[i]]
              if (is.na(spe)) {
                spe <- NULL
              }
              if (methods::is(spe, "SpatialExperiment")) {
                metric <- quiet_basic(FUN(spe, ...))
                if (length(metric) == 0) {
                  metric <- 0
                }
                list.metric[[i]] <- metric
              }
              else {
                list.metric[[i]] <- 0
              }
            }
            x <- raster::raster(ncol = n_split, nrow = n_split, xmn = 0,
                                ymn = 0, xmx = max(SpatialExperiment::spatialCoords(spe_object)[,
                                                                                                "Cell.X.Position"]), ymx = max(SpatialExperiment::spatialCoords(spe_object)[,
                                                                                                                                                                            "Cell.X.Position"]))
            raster::values(x) <- unlist(list.metric)
            y <- raster::flip(x, direction = "y")
            return(y)
          }
        quiet_basic <- function(x) {
          sink(tempfile())
          on.exit(sink())
          invisible(force(x))
        }

        #obtain individual SPIAT object
        SPIAT_Object <- DATA_SPIAT[[Index]]

        #Obtain Individual Phenotype object
        Phenotype_Object <- DATA_Phenotypes[[Index]]

        #Calculate Image area (approximation)
        Width <- max(Phenotype_Object$X) - min(Phenotype_Object$X)
        Height <- max(Phenotype_Object$Y) - min(Phenotype_Object$Y)
        Area <- Width * Height

        #Calculate the number of constant size tiles needed to cover the total area
        N_tiles <- ceiling(Area / (Tile_size^2))

        #Perform grid metrics
        grid_object <- Own_grid_metrics(SPIAT_Object, FUN = calculate_entropy, n_split = N_tiles,
                                        cell_types_of_interest = Phenotypes_included,
                                        feature_colname = "Phenotype")

        #Calculate percentage of tiles above a defined entropy threshold
        Percent <- SPIAT::calculate_percentage_of_grids(grid_object, threshold = Entropy_threshold, above = TRUE)
        Spatial_autocorrelation <- SPIAT::calculate_spatial_autocorrelation(grid_object, metric = Autocorrelation_metric)

        c(Subject_Names = names(DATA_SPIAT)[Index],
          PER_above_threshold = Percent,
          N_tiles = N_tiles*N_tiles,
          Spatial_autocorrelation = Spatial_autocorrelation)
      },
      .progress = TRUE)
    future::plan("future::sequential")
    gc()

    return(purrr::map_df(RESULTS,dplyr::bind_rows))
  }
