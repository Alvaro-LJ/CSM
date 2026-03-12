#' Calculates heterogeneity by tile using tiled images generated with random parameters
#'
#'
#' @param N_cores Integer. Number of cores to parallelize your computation.
#' @param Random_tiled_images A list containing tiled images obtained using [Random_tiled_image_heterogeneity_calculator()].
#' @param Minimum_cell_no_per_tile An integer indicating the minimum number of cells that a tile must contain. Tiles below the limit will not be included in the analysis.
#' @param Phenotypes_included A character vector indicating the phenotype labels that will be included in the analysis.
#'
#' @details
#' The function is aimed to estimate impact of the Modifiable areal unit problem (MAUP). For a given set of tiled images generated using random parameters, the function will compute by-tile heterogeneity metrics in a similar
#' way as it's homologue function [Tiled_image_heterogeneity_calculator()].
#'
#' @returns A list containing the parameters used to generate the tiled images and the actual results.
#'
#' @seealso [Random_image_tiling_processing_function()], [Random_tiled_image_heterogeneity_analyzer()], [Random_texture_features_calculator()]
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


Random_tiled_image_heterogeneity_calculator <-
  function(N_cores = 1,
           Random_tiled_images,
           Minimum_cell_no_per_tile = 1,
           Phenotypes_included){

    ####Check suggested packages####
    {
      if(!requireNamespace("picante", quietly = FALSE)) stop(
        paste0("picante CRAN package is required to execute the function. Please install using the following code: ",
               expression(install.packages("picante")))
      )
      if(!requireNamespace("vegan", quietly = FALSE)) stop(
        paste0("vegan CRAN package is required to execute the function. Please install using the following code: ",
               expression(install.packages("vegan")))
      )
      if(!requireNamespace("DescTools", quietly = FALSE)) stop(
        paste0("DescTools CRAN package is required to execute the function. Please install using the following code: ",
               expression(install.packages("DescTools")))
      )
      if(!requireNamespace("philentropy", quietly = FALSE)) stop(
        paste0("philentropy CRAN package is required to execute the function. Please install using the following code: ",
               expression(install.packages("philentropy")))
      )
    }

    ####Check arguments####
    #Check that Random_tiled_images are adequately formatted
    if(!identical(names(Random_tiled_images), c("Iteration_Parameters", "Results"))) stop("Random_tiled_images should be generated using the Random_image_tiling_processing_function")

    #Check that the phenotypes included are present in the data
    Unique_phenotypes_in_dataset <- unique(unlist(purrr::map(Random_tiled_images[[2]], function(Image) Image[[2]][["Phenotype"]])))
    if(!all(Phenotypes_included %in% Unique_phenotypes_in_dataset)) stop(paste0("Phenotypes included must be any of: ", stringr::str_c(Unique_phenotypes_in_dataset, collapse = ", ")))

    #Check that Min cells are integer values
    if(!all(Minimum_cell_no_per_tile%%1 == 0, Minimum_cell_no_per_tile > 0)) stop("Minimum_cell_no_per_tile must be an integer value > 0")

    #Check cores
    if(!all(N_cores%%1 == 0, N_cores > 0, length(N_cores) == 1)) stop("N_cores must be an integer value > 0")

    ####Iterate over the elements in the Results list (EVERY IMAGE)####
    #Set the on exit statement
    on.exit({
      future::plan("future::sequential")
      gc()
    })

    #Set the cores
    future::plan("future::multisession", workers = N_cores)
    options(future.globals.maxSize = Inf, future.rng.onMisuse = "ignore")
    furrr::furrr_options(scheduling = Inf)

    Heterogeneity_results <-
      furrr::future_map(Random_tiled_images[["Results"]], function(Image_results){

        #For every image we will build a pseudo tiled version list two-element list items
        Image_iteration_list <-
          purrr::map(1:length(Image_results[["Iteration_tile_list"]]), function(Iteration_index){
            #Get the tile info and the cell_info
            Tile_info <- Image_results[["Iteration_tile_list"]][[Iteration_index]]
            Tile_info <- Tile_info %>% dplyr::rename("tile_id" = all_of(stringr::str_c("tile_id_iteration_", Iteration_index, sep = "")))
            Tile_info <- Tile_info %>% dplyr::select(tile_X_centroid, tile_Y_centroid, tile_id, tile_xmin, tile_xmax, tile_ymin, tile_ymax)

            Cell_info <- Image_results[["Cell_information"]] %>% dplyr::select(Cell_no, Subject_Names, Phenotype, dplyr::all_of(stringr::str_c("tile_id_iteration_", Iteration_index, sep = "")))
            Cell_info <- Cell_info %>% dplyr::rename("tile_id" = all_of(stringr::str_c("tile_id_iteration_", Iteration_index, sep = "")))
            Cell_info <- Cell_info %>% dplyr::select(Cell_no, Phenotype, tile_id)

            return(list(Tile_info, Cell_info))
          })
        #Add names
        names(Image_iteration_list) <- stringr::str_c("Iteration", 1:length(Image_iteration_list), sep = "_")

        #We will compute the phenotypes included for every image just in case any image does not contain all the required phenotypes and
        #Then the function returns an error
        Image_phenotypes <- unique(unlist(map(Image_iteration_list, ~.[[2]]$Phenotype)))
        Image_phenotypes <- Image_phenotypes[Image_phenotypes %in% Phenotypes_included]


        #Now we will run the heterogeneity analysis by tile as usual
        Heterogeneity_list <-
          Tiled_image_heterogeneity_calculator(
            Tiled_images = Image_iteration_list,
            Minimum_cell_no_per_tile = Minimum_cell_no_per_tile,
            Phenotypes_included = Image_phenotypes
          )
        Heterogeneity_list <- purrr::map(Heterogeneity_list, function(results) results %>%
                                           dplyr::select(tile_X_centroid, tile_Y_centroid, tile_id, tile_xmin, tile_xmax, tile_ymin, tile_ymax,
                                                         Shannon, Simpson, Inverse_simpson, Renyi_Scale_Inf, Rao_Dkk, Gini, Kullback_Leibler, Jensen_Shannon))
        names(Heterogeneity_list) <- stringr::str_c("Iteration", 1:length(Heterogeneity_list), sep = "_")
        return(Heterogeneity_list)
      },
      .progress = TRUE)

    #Stop the multicore
    future::plan("future::sequential")
    gc()

    #Return the final results
    return(
      list(Iteration_Parameters = Random_tiled_images[["Iteration_Parameters"]],
           Heterogeneity_by_tile = Heterogeneity_results)
    )
  }


