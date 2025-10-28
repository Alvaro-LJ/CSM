#' Summarizes cell neighborhoods for every image in the dataset
#'
#' The function generates neighborhood counts, percentages and optionally neighborhood densities.
#'
#' @param DATA A dataframe or tibble containing cell feature data and a column named 'Neighborhood_assignment'.
#' @param Calculate_Density A logical value indicating if density should be calculated.
#' @param DATA_Area If Calculate_Density is TRUE, a dataframe or tibble containing image names and tissue area information. It can be calculated using [Image_size_calculator()].
#'
#' @returns Returns a tibble with neighborhood counts, percentages and optionally densities per image.
#'
#' @examples
#' \dontrun{
#' Neighborhood_Quantifier(
#'     DATA = CSM_Neighborhoods_test,
#'     Calculate_Density = FALSE
#')
#' }

#'
#' @export

Neighborhood_Quantifier <-
  function(DATA,
           Calculate_Density = FALSE,
           DATA_Area = NULL
  ) {

    #Check arguments
    if(!is.logical(Calculate_Density)){
      stop("Calculate_Density must be a logical value")
    }
    #Check that DATA has a neighborhood_assignment variable
    if(!"Neighborhood_assignment" %in% names(DATA)) stop("DATA must be generated using Neighborhood_discovery_function or UTAG_Neighborhood_identifier")

    #Import Neighborhood data
    DATA_neighborhoods <- DATA

    #Generate a tibble where the number of neighborhoods by sample are specified
    RESULTS_TIBBLE <-purrr::map_dfr(unique(DATA_neighborhoods$Subject_Names), function(Image) {
      #Select individual images
      Interim <- DATA_neighborhoods %>% dplyr::filter(Subject_Names == Image)

      #Count the number of neighborhoods
      Results <- Interim %>% dplyr::count(Neighborhood_assignment) %>% tidyr::pivot_wider(names_from = Neighborhood_assignment, values_from = n)
      #Quantify the proportion
      Results_prop <- Results/nrow(Interim)
      names(Results_prop) <- stringr::str_c("PROP_", names(Results))
      Results$N_cells <- nrow(Interim)
      Results$Subject_Names <- Image

      #arrange the tibble and exit loop
      Results <- Results[c(ncol(Results), 1:(ncol(Results)-1))]
      dplyr::bind_cols(Results, Results_prop)

    }, .progress = list(clear = F,
                        name = "Calculating neighborhood counts",
                        show_after = 1,
                        type = "iterator"))

    #Substitute NA by 0
    RESULTS_TIBBLE[is.na(RESULTS_TIBBLE)] <- 0

    #Graph without clinical data
    plot(RESULTS_TIBBLE %>% dplyr::select(Subject_Names, contains("PROP_")) %>% tidyr::pivot_longer(-Subject_Names) %>%
           ggplot(aes(x = Subject_Names, y = value)) + facet_wrap(~name, "free_y", ncol = 1, nrow = length(unique(DATA_neighborhoods$Neighborhood_assignment))) +
           geom_col(width = 0.5, color = "black") +
           cowplot::theme_cowplot() +
           scale_y_continuous("% of neighborhood in sample")+
           scale_x_discrete("")+
           theme(axis.text.x = element_text(angle = -90,
                                            vjust = 0.5))
    )

    #If density computing is required, execute the following code
    if(Calculate_Density){
      #Check arguments
      if(names(DATA_Area)[ncol(DATA_Area)] != "Area"){
        stop("The last column of the DATA_Area must be named Area")
      }

      print("Calculating densities")
      #Select cell counts
      For_density <- RESULTS_TIBBLE %>% dplyr::select(-N_cells, -contains("PROP_"))

      #Join the cell counts with the area tibbles
      For_density <-dplyr::left_join(For_density, DATA_Area, by = "Subject_Names")

      #Calculate densities
      Density_results <- as_tibble(For_density[c(2:(ncol(For_density)-1))] / For_density[[ncol(For_density)]])
      names(Density_results) <-stringr::str_c("Density_", names(Density_results))

      #Graph without clinical data
      plot(bind_cols(For_density[1], Density_results) %>% dplyr::select(Subject_Names, contains("Density_")) %>% tidyr::pivot_longer(-Subject_Names) %>%
             ggplot(aes(x = Subject_Names, y = value)) + facet_wrap(~name, "free_y", ncol = 1, nrow = length(unique(DATA_neighborhoods$Neighborhood_assignment))) +
             geom_col(width = 0.5, color = "black") +
             cowplot::theme_cowplot() +
             scale_y_continuous("Density of neighborhood in sample")+
             scale_x_discrete("")+
             theme(axis.text.x = element_text(angle = -90,
                                              vjust = 0.5))
      )
      #Bind the results and return the tibble
      return(bind_cols(RESULTS_TIBBLE, Density_results, For_density[ncol(For_density)]))
    }

    #If not required exit the function returning the common results
    else{
      #Exit the function with the final results
      return(RESULTS_TIBBLE)
    }
  }
