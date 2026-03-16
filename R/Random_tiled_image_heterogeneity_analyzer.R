#' Generates a summary of the heterogeneity by tile using tiled images generated with random parameters
#'
#'
#' @param Random_tiled_heterogeneity_DATA A list containing tiled images obtained using [Random_tiled_image_heterogeneity_calculator()].
#' @param Strategy A character value indicating the analytical approach. Either 'Quantify_by_Threshold' or 'Overall_Summary' (see details).
#' @param Metric A character value indicating the heterogeneity metric to be plotted.
#' @param Threshold If Strategy is Quantify_by_Threshold, the arbitrary threshold to be used.
#'
#' @details
#' The function is aimed to estimate impact of the Modifiable areal unit problem (MAUP). For a given set of tiled images generated using random parameters, the function will compute by sample heterogeneity metrics in a similar
#' way as it's homologue function [Tiled_image_heterogeneity_analyzer()].
#'
#' The first two random iterations will always contain tiling strategies in the extremes of the user defined tile size variation factor.
#'
#' @returns A list containing the parameters used to generate the tiled images and the actual results.
#'
#' @seealso [Random_image_tiling_processing_function()], [Random_tiled_image_heterogeneity_calculator()], [Random_texture_features_calculator()]
#'
#' @examples
#' \dontrun{
#' #Generate tiled image data using 10 randomly generated parameter set---------
#' Random_tile_DATA <-
#'   Random_image_tiling_processing_function(
#'     N_cores =  1,
#'     DATA = CSM::CSM_Phenotypecell_test,
#'     Tile_width = 100,
#'     Tile_height = 100,
#'     Variables_to_keep = "Phenotype",
#'
#'     N_iterations = 10,
#'     Max_tile_size_variation = 1.1,
#'     Force_squared_tiles = TRUE
#'   )
#'
#' #Compute the heterogeneity by tile using this random tile dataset------------
#' Tiled_heterogeneity_DATA <-
#'   Random_tiled_image_heterogeneity_calculator(
#'     N_cores = 1,
#'     Random_tiled_images = Random_tile_DATA,
#'     Minimum_cell_no_per_tile = 1,
#'     Phenotypes_included = unique(CSM_Phenotypecell_test$Phenotype)
#'   )
#'
#' #Analyze the data-----------------------------------------------------------
#' Random_tiled_image_heterogeneity_analyzer(
#'     Random_tiled_heterogeneity_DATA = Tiled_heterogeneity_DATA,
#'     Strategy = "Overall_Summary",
#'     Metric = "Shannon"
#'   )
#'}
#'
#'
#' @export

Random_tiled_image_heterogeneity_analyzer <-
  function(
    Random_tiled_heterogeneity_DATA,
    Strategy,
    Metric,
    Threshold = NULL){

    ####Check required packages####
    if(!requireNamespace("ape", quietly = FALSE)) stop(
      paste0("ape CRAN package is required to execute the function. Please install using the following code: ",
             expression(install.packages("ape")))
    )

    ####Check arguments####
    if(!(Strategy %in% c("Quantify_by_Threshold", "Overall_Summary"))) {
      stop("Strategy should either Quantify_by_Threshold or Overall_Summary")
    }
    if(Strategy == "Quantify_by_Threshold"){
      if(!all(length(Threshold) == 1, is.numeric(Threshold))) stop("Threshold must be a single numeric value")
    }

    if(!identical(names(Random_tiled_heterogeneity_DATA), c("Iteration_Parameters", "Heterogeneity_by_tile"))) stop("Random_tiled_heterogeneity_DATA should be generated using the Random_tiled_image_heterogeneity_calculator")

    Available_metrics <- unique(unlist(purrr::map(Random_tiled_heterogeneity_DATA[[2]], function(Image){
      purrr::map(Image, ~names(.)[-c(1:7)])
    })))
    if(!Metric %in% Available_metrics) stop(paste0("Metric must be one of the following: \n", stringr::str_c(Available_metrics, collapse = "\n")))


    #####Compute the metrics
    #Iterate over the available images. A single tibble per image will be created containing the image iterations, then these will be summarized
    Iteration_results <-
      purrr::map(Random_tiled_heterogeneity_DATA[[2]], function(Image){
        #Get the results
        return(
          Tiled_image_heterogeneity_analyzer(
            Tiled_heterogeneity_DATA = Image,
            Strategy = Strategy,
            Metric = Metric,
            Threshold = Threshold
          )
        )
    })

    #Generate a summary of the results according to the selected strategy (generated tables are different)
    #If Quantify_by_Threshold
    if(Strategy == "Quantify_by_Threshold"){
      #Generate the average results
      Average_tibble <- purrr::map_dfr(Iteration_results, function(Image){
        purrr::map_dfc(Image[2:6], function(column) mean(column, na.rm = TRUE))
      })
      names(Average_tibble) <- stringr::str_c("Average_", names(Average_tibble), sep = "")
      Average_tibble <- dplyr::bind_cols(tibble(Subject_Names = names(Iteration_results),
                                                N_iterations = nrow(Iteration_results[[1]])),
                                         Average_tibble)
      Average_tibble$MoransI_pval_below_05 <- purrr::map_dbl(Iteration_results, function(Image) sum(Image[[8]] < 0.05, na.rm = TRUE)/length(Image[[8]]))
      Average_tibble$Threshold <- Threshold
      Average_tibble$Metric <- Metric

      #Generate the median results
      Median_tibble <- purrr::map_dfr(Iteration_results, function(Image){
        purrr::map_dfc(Image[2:6], function(column) quantile(column, 0.5, na.rm = TRUE))
      })
      names(Median_tibble) <- stringr::str_c("Median_", names(Median_tibble), sep = "")
      Median_tibble <- dplyr::bind_cols(tibble(Subject_Names = names(Iteration_results),
                                                N_iterations = nrow(Iteration_results[[1]])),
                                        Median_tibble)
      Median_tibble$MoransI_pval_below_0.05 <- Average_tibble$MoransI_pval_below_05
      Median_tibble$Threshold <- Threshold
      Median_tibble$Metric <- Metric

      #Return all the results in a list
      return(
        list(
          Iteration_Parameters = Random_tiled_heterogeneity_DATA[["Iteration_Parameters"]],
          Raw_results = Iteration_results,
          Average_by_sample = Average_tibble,
          Median_by_sample = Median_tibble
        )
      )
    }

    #If Overall summary
    if(Strategy == "Overall_Summary"){
      #Generate the average results
      Average_tibble <- purrr::map_dfr(Iteration_results, function(Image){
        purrr::map_dfc(Image[2:13], function(column) mean(column, na.rm = TRUE))
      })
      names(Average_tibble) <- stringr::str_c("Average_", names(Average_tibble), sep = "")
      Average_tibble <- dplyr::bind_cols(tibble(Subject_Names = names(Iteration_results),
                                                N_iterations = nrow(Iteration_results[[1]])),
                                         Average_tibble)
      Average_tibble$MoransI_pval_below_05 <- purrr::map_dbl(Iteration_results, function(Image) sum(Image[[14]] < 0.05, na.rm = TRUE)/length(Image[[8]]))
      Average_tibble$Metric <- Metric

      #Generate the median results
      Median_tibble <- purrr::map_dfr(Iteration_results, function(Image){
        purrr::map_dfc(Image[2:13], function(column) quantile(column, 0.5, na.rm = TRUE))
      })
      names(Median_tibble) <- stringr::str_c("Median_", names(Median_tibble), sep = "")
      Median_tibble <- dplyr::bind_cols(tibble(Subject_Names = names(Iteration_results),
                                               N_iterations = nrow(Iteration_results[[1]])),
                                        Median_tibble)
      Median_tibble$MoransI_pval_below_0.05 <- Average_tibble$MoransI_pval_below_05
      Median_tibble$Metric <- Metric

      #Return all the results in a list
      return(
        list(
          Iteration_Parameters = Random_tiled_heterogeneity_DATA[["Iteration_Parameters"]],
          Raw_results = Iteration_results,
          Average_by_sample = Average_tibble,
          Median_by_sample = Median_tibble
        )
      )
    }

  }

