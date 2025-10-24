#' Calculates heterogeneity by tile
#'
#' The function calculates various heterogeneity metrics based on cell composition distribution of every tile in an image.
#'
#' @param Tiled_images A list containing tiled images obtained using [Image_tiling_processing_function()].
#' @param Minimum_cell_no_per_tile An integer indicating the minimum number of cells that a tile must contain. Tiles below the limit will not be included in the analysis.
#' @param Phenotypes_included A character vector indicating the phenotype labels that will be included in the analysis.
#'
#' @details
#' Shannon, Simpson,Inverse Simpson and Renyi entropy indexes are calculated using the vegan package.
#'
#' Rao index is calculated using the picante::raoD function.
#'
#' Gini index is calculated using the DescTools::Gini function.
#'
#' Kullback_Leibler and Jensen_Shannon indexes are calculated using the philentropy package. They compare the observed cell composition within each tile against the global cell composition of the sample.
#'
#' @returns A list containing image information with by-tile heterogeneity metrics.
#'
#' @examples
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
#' @export

Tiled_image_heterogeneity_calculator <-
  function(Tiled_images = NULL,
           Minimum_cell_no_per_tile = NULL,
           Phenotypes_included = NULL) {

    #Check suggested packages
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

    #Check arguments
    #Check that the phenotypes included are present in the data
    if(!all(Phenotypes_included %in% unique(unlist(purrr::map(Tiled_images, function(df) df[[2]]$Phenotype))))) {
      stop(paste0("Phenotypes included must be any of: ", stringr::str_c(unique(unlist(purrr::map(Tiled_images, function(df) df[[2]]$Phenotype))), collapse = ", ")))
    }
    #Check that Min cells are integer values
    if(!all(Minimum_cell_no_per_tile%%1 == 0, Minimum_cell_no_per_tile > 0)) stop("Minimum_cell_no_per_tile must be an integer value > 0")

    #Else proceed with analysis
    Results <-
      purrr::map(1:length(Tiled_Images), function(x) {

        Image <- Tiled_Images[[x]]
        Interim <- Image[[2]] %>% dplyr::filter(Phenotype %in% Phenotypes_included)

        #Generate the global cell count expected proportions (required for KL and JSD)
        Cell_counts_by_image <- Interim %>% count(Phenotype)
        Expected_proportions <- Cell_counts_by_image[["n"]] / sum(Cell_counts_by_image[["n"]])
        names(Expected_proportions) <-stringr::str_c("PER_", Cell_counts_by_image[["Phenotype"]], sep = "")
        Expected_proportions <- Expected_proportions[sort(names(Expected_proportions))]

        #Filter out tiles with less than minimum cells per tile
        Filtered_tiles <-
          Interim  %>%
          group_by(tile_id) %>% dplyr::count() %>%dplyr::ungroup() %>% dplyr::filter(n >= Minimum_cell_no_per_tile)

        #If not enough tiles are present in the image print a warning message
        if(nrow(Filtered_tiles) == 0) {
          warning(paste0("No valid tiles present in: ", names(Tiled_Images)[x], ", hence it will be removed from analysis. Please consider lowering Minimum_cell_no_per_tile threshold or increasing the cell types included in analysis"))
          return(NULL)
        }

        else{

          #Build cell count by tile matrix
          Interim2 <-
            Interim %>% dplyr::filter(tile_id %in% Filtered_tiles[[1]]) %>% group_by(tile_id, Phenotype) %>% dplyr::count() %>%
            tidyr::pivot_wider(id_cols = tile_id,
                        names_from = Phenotype,
                        values_from = n)
          Interim2[is.na(Interim2)] <- 0

          #Calculate the percentage of the total cells in the tile that belong to each phenotype and build a tibble with the info
          Interim_per <-purrr::map_dfc(Interim2[-1], function(Column) Column / rowSums(Interim2[-1]))
          names(Interim_per) <-stringr::str_c("PER_", names(Interim2)[-1], sep = "")
          Interim_per <- Interim_per[sort(names(Interim_per))]

          #Add missing columns if any cell type has gone missing after tile removal process (is required for KL and JSD analyses that are based on proportion)
          if(!all(names(Expected_proportions) %in% names(Interim_per))){
            Interim_missing <- suppressWarnings(as_tibble(matrix(0, nrow = nrow(Interim_per), ncol = sum(!names(Expected_proportions) %in% names(Interim_per)))))
            names(Interim_missing) <- names(Expected_proportions)[!names(Expected_proportions) %in% names(Interim_per)]
            Interim_per <-dplyr::bind_cols(Interim_per, Interim_missing)
            Interim_per <- Interim_per[sort(names(Interim_per))]
          }

          #Calculate the total cells per tile
          Total_cells_tibble <- tibble(n_cells = rowSums(Interim2[-1]))

          #Calculate the metrics and generate the result tibble
          Results <-dplyr::bind_cols(
            Interim2[1],
            tibble(Shannon = apply(Interim2[-1], MARGIN = 1, function(row) vegan::diversity(row, index = "shannon")),
                   Simpson = apply(Interim2[-1], MARGIN = 1, function(row) vegan::diversity(row, index = "simpson")),
                   Inverse_simpson = apply(Interim2[-1], MARGIN = 1, function(row) vegan::diversity(row, index = "invsimpson"))
            ),
            tibble(renyi = as.double(vegan::renyi(Interim2[-1], scales = Inf))),
            tibble(rao_Dkk = picante::raoD(Interim2[-1])$Dkk),
            tibble(Gini = unlist(apply(Interim2[-1], MARGIN = 1, function(Sample) DescTools::Gini(Sample, conf.level = NA)))),
            tibble(Kullback_Leibler = unlist(apply(Interim_per, MARGIN = 1, function(Row){
              Observed <- Row
              Expected <- Expected_proportions
              suppressMessages(philentropy::KL(rbind(Observed, Expected), unit = "log2"))
            }))),
            tibble(Jensen_Shannon = unlist(apply(Interim_per, MARGIN = 1, function(Row){
              Observed <- Row
              Expected <- Expected_proportions
              suppressMessages(philentropy::JSD(rbind(Observed, Expected),  test.na = FALSE, unit = "log2"))
            }))),
            Interim2[-1],
            Total_cells_tibble,
            Interim_per
          )

          #Change names
          names(Results)[2:9] <-
            c("Shannon",
              "Simpson",
              "Inverse_simpson",
              "Renyi_Scale_Inf",
              "Rao_Dkk",
              "Gini",
              "Kullback_Leibler",
              "Jensen_Shannon")

          #Bind results to the tile info matrix and eliminate rows with NA
          Image_results <-dplyr::left_join(Image[[1]], Results, by = "tile_id")
          #Safely remove images with less than one evaluable metric
          Image_results <- Image_results[apply(Image_results, MARGIN = 1, function(row) sum(is.na(row)) < 9),]
          #reorder the output and return
          Image_results <- Image_results %>% dplyr::select(1:7,
                                                           all_of(c("Shannon",
                                                                    "Simpson",
                                                                    "Inverse_simpson",
                                                                    "Renyi_Scale_Inf",
                                                                    "Rao_Dkk",
                                                                    "Gini",
                                                                    "Kullback_Leibler",
                                                                    "Jensen_Shannon")),
                                                           n_cells,
                                                           names(Image_results)[names(Image_results) %in% Phenotypes_included],
                                                           contains("PER_"))
          return(Image_results)

        }
      },
      .progress = list(clear = F,
                       name = "Calculating heterogeneity by tile",
                       show_after = 2,
                       type = "iterator"))

    names(Results) <- names(Tiled_Images)

    #Return results except for NULL values
    return(Results[!purrr::map_lgl(Results, is.null)])
  }
