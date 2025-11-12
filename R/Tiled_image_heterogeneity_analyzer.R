#' Generates a summary of the heterogeneity by tile.
#'
#' @param Tiled_heterogeneity_DATA A list containing tiled images obtained using [Tiled_image_heterogeneity_calculator()].
#' @param Strategy A character value indicating the analytical approach. Either 'Quantify_by_Threshold' or 'Overall_Summary' (see details).
#' @param Metric A character value indicating the heterogeneity metric to be plotted.
#' @param Threshold If Strategy is Quantify_by_Threshold, the arbitrary threshold to be used.
#'
#' @details
#' Quantify_by_Threshold counts the number of tiles above the arbitrary threshold and spatial autocorrelation (Moran's I) of the positive tiles.
#'
#' Overall_summary quantifies the central tendency and dispersion metrics of the tile heterogeneity metrics and their spatial autocorrelation.
#'
#' @returns A tibble containing the summary of the analysis
#'
#' @seealso [Tiled_image_heterogeneity_calculator()], [Tiled_image_heterogeneity_graph_maker()]
#'
#' @examples
#' \dontrun{
#' #Divide cells into tiles---------
#' Tiled_Images <-
#' Image_tiling_processing_function(
#'    N_cores = 2,
#'    DATA = CSM_Phenotypecell_test,
#'    Tile_width = 125,
#'    Tile_height = 125,
#'    Variables_to_keep = "Phenotype"
#')
#'
#' #Calculate heterogeneity by tile----
#' Tiled_image_heterogeneity_calculator(
#'     Tiled_images = Tiled_Images,
#'     Minimum_cell_no_per_tile = 3,
#'     Phenotypes_included = c("TUMOR", "CD8_GZMBneg", "CD8_GZMBpos", "OTHER")
#')
#'
#' #Generate the summary
#' Tiled_image_heterogeneity_analyzer(
#'     Tiled_heterogeneity_DATA = Heterogeneity_by_tile,
#'     Strategy = "Quantify_by_Threshold",
#'     Metric = "Shannon",
#'     Threshold = 0.5
#')
#' }
#'
#' @export


Tiled_image_heterogeneity_analyzer <-
  function(Tiled_heterogeneity_DATA,
           Strategy,
           Metric,
           Threshold = NULL) {

    #Check required packages
    if(!requireNamespace("ape", quietly = FALSE)) stop(
      paste0("ape CRAN package is required to execute the function. Please install using the following code: ",
             expression(install.packages("ape")))
    )

    #Check arguments
    if (!(Strategy %in% c("Quantify_by_Threshold", "Overall_Summary"))) {
      stop("Strategy should either Quantify_by_Threshold or Overall_Summary")
    }
    else if(!all(purrr::map_lgl(seq_along(1:length(Tiled_heterogeneity_DATA)), function(Image) Metric %in% names(Tiled_heterogeneity_DATA[[Image]])))) {
      stop("Metric should be contained in the Tiled_heterogeneity_DATA object")
    }


    else if(Strategy == "Quantify_by_Threshold") {
      if(!is.numeric(Threshold)) stop("Threshold must be a numeric value")

      purrr::map_dfr(seq_along(1:length(Tiled_heterogeneity_DATA)), function(Image) {

        #Import tiled image data and select requested heterogeneity metric
        Image_filtered <- Tiled_heterogeneity_DATA[[Image]] %>% dplyr::select(1:7, all_of(Metric))
        names(Image_filtered)[8] <- "value"

        #Prepare data to perform Morans I
        For_Morans <- Image_filtered %>% dplyr::mutate(value = case_when(value >= Threshold ~ 1,
                                                                        TRUE ~ 0))
        if(length(unique(For_Morans$value)) > 1) {

          #Prepare data for Morans I calculation
          DIST_Morans <- as.matrix(dist(cbind(For_Morans$tile_X_centroid, For_Morans$tile_Y_centroid), method = "euclidean"))#We calculate distance matrix
          INVERSE_DIST <- 1/DIST_Morans #We calculate inverse of distance
          diag(INVERSE_DIST) <- 0
          results <- ape::Moran.I(For_Morans$value, INVERSE_DIST, na.rm = T)#Calculate MORAN's I

          #Generate result tibble
          Results_tibble <- tibble(Subject_Names = names(Tiled_heterogeneity_DATA)[Image],
                                   Tiles_above_Threshold = sum(Image_filtered$value >= Threshold),
                                   Total_tiles = nrow(Image_filtered),
                                   Percentage_above_Threshold = Tiles_above_Threshold / Total_tiles,
                                   MoransI_obseved = results$observed,
                                   MoransI_expected = results$expected,
                                   MoransI_sd = results$sd,
                                   MoransI_pval = results$p.value,
                                   Threshold = Threshold,
                                   Metric = Metric)
        }

        else if(length(unique(For_Morans$value)) == 1){

          #No need to calculate Morans I
          Results_tibble <- tibble(Subject_Names = names(Tiled_heterogeneity_DATA)[Image],
                                   Tiles_above_Threshold = sum(Image_filtered$value >= Threshold),
                                   Total_tiles = nrow(Image_filtered),
                                   Percentage_above_Threshold = Tiles_above_Threshold / Total_tiles,
                                   MoransI_obseved = NA,
                                   MoransI_expected = NA,
                                   MoransI_sd = NA,
                                   MoransI_pval = NA,
                                   Threshold = Threshold,
                                   Metric = Metric)
        }
        return(Results_tibble)

      },
      .progress = list(clear = F,
                       name = "Calculating heterogeneity summary acording to threshold",
                       show_after = 2,
                       type = "iterator"))

    }

    else if(Strategy == "Overall_Summary") {
      purrr::map_dfr(seq_along(1:length(Tiled_heterogeneity_DATA)), function(Image) {
        Image_filtered <- Tiled_heterogeneity_DATA[[Image]] %>% dplyr::select(1:7, all_of(Metric))

        names(Image_filtered)[8] <- "value"

        if(length(unique(Image_filtered$value)) > 1) {

          #Prepare data for Morans I calculation
          DIST_Morans <- as.matrix(dist(cbind(Image_filtered$tile_X_centroid, Image_filtered$tile_Y_centroid), method = "euclidean"))#We calculate distance matrix
          INVERSE_DIST <- 1/DIST_Morans #We calculate inverse of distance
          diag(INVERSE_DIST) <- 0
          results <- ape::Moran.I(Image_filtered$value, INVERSE_DIST, na.rm = T)#Calculate MORAN's I

          Results_tibble <- tibble(Subject_Names = names(Tiled_heterogeneity_DATA)[Image],
                                   Min_value = min(Image_filtered$value, na.rm = TRUE),
                                   p25 = quantile(Image_filtered$value, 0.25, na.rm = TRUE),
                                   Average = mean(Image_filtered$value, na.rm = TRUE),
                                   p50 = quantile(Image_filtered$value, 0.50, na.rm = TRUE),
                                   p75 = quantile(Image_filtered$value, 0.75, na.rm = TRUE),
                                   Max_value = max(Image_filtered$value, na.rm = TRUE),
                                   sd_value = sd(Image_filtered$value, na.rm = TRUE),
                                   Total_tiles = nrow(Image_filtered),
                                   Evaluable_tiles = nrow(Image_filtered) - sum(is.na(Image_filtered$value)),
                                   MoransI_obseved = results$observed,
                                   MoransI_expected = results$expected,
                                   MoransI_sd = results$sd,
                                   MoransI_pval = results$p.value,
                                   Metric = Metric
          )
        }

        else if(length(unique(Image_filtered$value)) == 1) {

          #no need to calculate Morans I
          Results_tibble <- tibble(Subject_Names = names(Tiled_heterogeneity_DATA)[Image],
                                   Min_value = min(Image_filtered$value, na.rm = TRUE),
                                   p25 = quantile(Image_filtered$value, 0.25, na.rm = TRUE),
                                   Average = mean(Image_filtered$value, na.rm = TRUE),
                                   p50 = quantile(Image_filtered$value, 0.50, na.rm = TRUE),
                                   p75 = quantile(Image_filtered$value, 0.75, na.rm = TRUE),
                                   Max_value = max(Image_filtered$value, na.rm = TRUE),
                                   sd_value = sd(Image_filtered$value, na.rm = TRUE),
                                   Total_tiles = nrow(Image_filtered),
                                   Evaluable_tiles = nrow(Image_filtered) - sum(is.na(Image_filtered$value)),
                                   MoransI_obseved = NA,
                                   MoransI_expected = NA,
                                   MoransI_sd = NA,
                                   MoransI_pval = NA,
                                   Metric = Metric
          )
        }
        return(Results_tibble)

      },
      .progress = list(clear = F,
                       name = "Calculating heterogeneity overall summary",
                       show_after = 2,
                       type = "iterator"))
    }
  }
