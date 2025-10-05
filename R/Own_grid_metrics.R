#' Interal function to RUN SPIAT safely and without bugs
#'
#' Intended for internal use only
#'
#' @details
#' Used in [SPIAT_Heterogeneity_Analyzer()]
#'
#'
#' @keywords Internal

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
