#' Performs image tiling
#'
#' The function assigns cells to unique tiles. Afterwards, tile information can be used to run subsequent analyses.
#'
#' @param N_cores Integer. Number of cores to parallelize your computation.
#' @param DATA A dataframe or tibble containing cell feature data.
#' @param Tile_width A numeric value indicating the width of the tiles.
#' @param Tile_height A numeric value indicating the height of the tiles.
#' @param Variables_to_keep A character vector indicating the names of the columns to be kept when tiling the image.
#'
#' @returns A list containing tile information for every image in DATA.
#'
#' @details
#' There is a one unit overlap between tile edges to account for potential decimal numbers in cell coordinates. Cells are assigned
#' to the first tile where they are located. In most scenarios, this 1 unit overlap wont impact the final results.
#'
#'
#' @seealso [Image_length_calculator()], [Image_tiling_processing_function()]
#'
#' @examples
#' \dontrun{
#' Image_tiling_processing_function(
#'    N_cores = 2,
#'    DATA = CSM_Phenotypecell_test,
#'    Tile_width = 125,
#'    Tile_height = 125,
#'    Variables_to_keep = "Phenotype"
#')
#' }
#'
#' @export


Image_tiling_processing_function <-
  function(N_cores = 1,
           DATA,
           Tile_width,
           Tile_height,
           Variables_to_keep) {

    #Load Data phenotypes and tile generator function to the environment created in the function
    DATA <- DATA
    if(!(is.numeric(Tile_width) & is.numeric(Tile_height))) stop("Both Tile_width and Tile_height must be numeric values")
    if(!(Tile_width%%1 == 0 & Tile_height%%1 ==0)) message("It is highly recommended that Tile_width and Tile_height are integer values")
    if(!(Tile_width > 0 & Tile_height > 0)) stop("Tile_width and Tile_height must be > 0")
    if(N_cores%%1 != 0) stop("N_cores must be an integer value")
    if(!all(Variables_to_keep %in% names(DATA))) stop(paste0(stringr::str_c(Variables_to_keep[!Variables_to_keep %in% names(DATA)], collapse = ", "),
                                                             " not found in DATA variables"))

    else {
      if (!all(c("X", "Y", "Subject_Names", Variables_to_keep) %in% names(DATA))) {
        stop("DATA supplied must be follow the structure specified in step 0 and contain Variables_to_keep")
      } else {

        #save exit function if parallelization fails
        on.exit({
          future::plan("future::sequential")
          gc()
        })

        future::plan("future::multisession", workers = N_cores)
        options(future.globals.maxSize = Inf, future.rng.onMisuse = "ignore")
        furrr::furrr_options(scheduling = Inf)
        Tiled_Images <-
          furrr::future_map(unique(DATA$Subject_Names), function(x)
            Tile_generator_function(
              Image_name = x,
              DATA = DATA,
              Tile_width = Tile_width,
              Tile_height = Tile_height,
              Variables_to_keep = Variables_to_keep
            ),
            .progress = TRUE)

        future::plan("future::sequential")
        gc()

        names(Tiled_Images) <- unique(DATA$Subject_Names)
        return(Tiled_Images)
      }
    }
  }
