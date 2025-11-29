#' Summarizes cell phenotypes for every image in a dataset
#'
#' The function generates cell counts, percentages and optionally, cell densities.
#'
#' @param DATA A dataframe or tibble containing cell feature data and a column named 'Phenotype' containing cell labels.
#' @param Calculate_Density A logical value indicating if density should be calculated.
#' @param DATA_Area If Calculate_Density is TRUE, a dataframe or tibble containing image names and tissue area information. It can be calculated using [Image_size_calculator()].
#'
#' @returns Returns a tibble with cell counts, percentages and optionally cell densities per image.
#'
#' @examples
#' \dontrun{
#' #Calculate tissue surface----------------
#' DATA_AREA <-
#' Image_size_calculator(
#'    DATA = CSM_Phenotypecell_test,
#'    Strategy = "Tiling",
#'    Tile_accuracy = 100
#'  )
#'
#'#Calculate the cell counts by image-------
#' Phenotype_quantifier(
#'     DATA = CSM_Phenotypecell_test,
#'     Calculate_Density = TRUE,
#'     DATA_Area = DATA_AREA
#')
#' }
#'
#' @export

Phenotype_quantifier <-
  function(DATA,
           Calculate_Density = FALSE,
           DATA_Area = NULL) {

    #Check arguments
    if(!is.logical(Calculate_Density)){
      stop("Calculate_Density must be a logical value")
    }

    #Obtain the number of cells by image according to the phenotypes
    Results <- DATA %>% group_by(Subject_Names, Phenotype) %>% dplyr::count() %>% dplyr::ungroup() %>%
      tidyr::pivot_wider(names_from = Phenotype, values_from = n)
    Results[is.na(Results)] <- 0

    #Calculate the number of total cells
    N_cells <- apply(Results[-1], MARGIN = 1, sum)

    #Calculate proportions
    Prop_tibble <- as_tibble(Results[-1] / N_cells)
    names(Prop_tibble) <-stringr::str_c("PROP_", names(Prop_tibble))
    Results$N_cells <- N_cells

    if(Calculate_Density){
      #Check arguments
      if(names(DATA_Area)[ncol(DATA_Area)] != "Area"){
        stop("The last column of the DATA_Area must be named Area")
      }

      #Select cell counts
      For_density <- Results %>% dplyr::select(-N_cells, -contains("PROP_"))

      #Join the cell counts with the area tibbles
      For_density <-dplyr::left_join(For_density, DATA_Area, by = "Subject_Names")

      #Calculate densities
      Density_results <- as_tibble(For_density[c(2:(ncol(For_density)-1))] / For_density[[ncol(For_density)]])
      names(Density_results) <-stringr::str_c("Density_", names(Density_results))

      #Bind the results and return the tibble
      return(dplyr::bind_cols(Results, Prop_tibble, Density_results, For_density[ncol(For_density)]))
    }

    else{
      #Return the bind results
      return(dplyr::bind_cols(Results, Prop_tibble))
    }
  }
