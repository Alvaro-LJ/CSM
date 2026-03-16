#' Calculates the average texture features of a single cell type using tiled images generated with random parameters
#'
#' The function is aimed to estimate impact of the Modifiable areal unit problem (MAUP). The function calculates the average texture features in a similar
#' way as it's homologue function [Texture_features_calculator()].
#'
#' @param N_cores Integer. Number of cores to parallelize your computation.
#' @param Random_tiled_images A list containing tiled images obtained using [Random_image_tiling_processing_function()].
#' @param Phenotype_included A character value indicating which cell phenotype to analyze.
#'
#' @returns A list containing the parameters used to generate the tiled images and the actual results.
#'
#' @seealso [Random_image_tiling_processing_function()], [Random_tiled_image_heterogeneity_calculator()], [Random_tiled_image_heterogeneity_analyzer()],
#'
#' @examples
#' \dontrun{
#'
#' #Generate tiled image data using 10 randomly generated parameter set
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
#' #Compute texture features in the dataset
#' Random_texture_features_calculator(
#'           N_cores = 1,
#'           Random_tiled_images = Random_tile_DATA,
#'           Phenotype_included = "CD8_GZMBneg"
#'   )
#' }
#'
#' @export

Random_texture_features_calculator <-
  function(N_cores = 1,
           Random_tiled_images,
           Phenotype_included){

    ####Check suggested packages####
    {
      if(!requireNamespace("raster", quietly = FALSE)) stop(
        paste0("raster CRAN package is required to execute the function. Please install using the following code: ",
               expression(install.packages("raster")))
      )
      if(!requireNamespace("glcm", quietly = FALSE)) stop(
        paste0("glcm CRAN package is required to execute the function. Please install using the following code: ",
               expression(install.packages("glcm")))
      )
      if(!requireNamespace("tabularaster", quietly = FALSE)) stop(
        paste0("tabularaster CRAN package is required to execute the function. Please install using the following code: ",
               expression(install.packages("tabularaster")))
      )
    }

    ####Check arguments####
    #Check that Random_tiled_images are adequately formatted
    if(!identical(names(Random_tiled_images), c("Iteration_Parameters", "Results"))) stop("Random_tiled_images should be generated using the Random_image_tiling_processing_function")

    #Check that the phenotypes included are present in the data
    Unique_phenotypes_in_dataset <- unique(unlist(purrr::map(Random_tiled_images[[2]], function(Image) Image[[2]][["Phenotype"]])))
    if(!all(length(Phenotype_included) == 1, Phenotype_included %in% Unique_phenotypes_in_dataset)) stop(paste0("Phenotype included must be one of: ", stringr::str_c(Unique_phenotypes_in_dataset, collapse = ", ")))

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

    Texture_results <-
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

        #Now we will run the heterogeneity analysis by tile as usual
        Texture_tibble <-
          Texture_features_calculator(Tiled_images = Image_iteration_list,
                                      Variable_name = "Phenotype",
                                      Phenotype_included = Phenotype_included
                                      )
        return(Texture_tibble)
      },
      .progress = TRUE)

    #Stop the multicore
    future::plan("future::sequential")
    gc()

    #Generate the summaries by sample
    #Generate the average results
    Average_tibble <- purrr::map_dfr(Texture_results, function(Image){
      purrr::map_dfc(Image[2:9], function(column) mean(column, na.rm = TRUE))
    })
    names(Average_tibble) <- stringr::str_c("Average_", names(Average_tibble), sep = "")
    Average_tibble <- dplyr::bind_cols(tibble(Subject_Names = names(Texture_results),
                                              N_iterations = nrow(Texture_results[[1]])),
                                       Average_tibble)

    #Generate the median results
    Median_tibble <- purrr::map_dfr(Texture_results, function(Image){
      purrr::map_dfc(Image[2:9], function(column) quantile(column, 0.5, na.rm = TRUE))
    })
    names(Median_tibble) <- stringr::str_c("Median_", names(Median_tibble), sep = "")
    Median_tibble <- dplyr::bind_cols(tibble(Subject_Names = names(Texture_results),
                                             N_iterations = nrow(Texture_results[[1]])),
                                      Median_tibble)

    #Return all the results in a list
    return(
      list(
        Iteration_Parameters = Random_tiled_images[["Iteration_Parameters"]],
        Raw_results = Texture_results,
        Average_by_sample = Average_tibble,
        Median_by_sample = Median_tibble
      )
    )
  }

