#' Summarizes thresholded data by image
#'
#' The function generates positive cell counts, percentages and optionally cell densities according to thresholded feature data.
#'
#' @param DATA A dataframe or tibble containing cell feature data that has been thresholded.
#' @param Calculate_Density A logical value indicating if density should be calculated.
#' @param DATA_Area If Calculate_Density is TRUE, a dataframe or tibble containing image names and tissue area information.
#' @returns Returns a tibble with cell counts, pertentages and optionally cell densities per image.
#'
#'
#' @examples
#' \dontrun{
#'#Calculate tissue area
#'DATA_AREA <- Image_size_calculator(
#'   DATA = CSM_Arrangedcellfeaturedata_test,
#'   Strategy = "Tiling",
#'   Image_to_plot = NULL,
#'   Tile_accuracy = 60
#')
#'
#'#Threshold data
#'DATA_thresholded <- Thresholding_function(
#'   DATA = CSM_Arrangedcellfeaturedata_test,
#'   Strategy = "EBI_Otsu",
#'   Local_thresholding = FALSE,
#'   Method_autothreshold = "Otsu",
#'   number_iterations_TriClass = 20,
#'   Percentile = 0.5,
#'   Defined_threshold = 0.1,
#'   Levels = 3
#' )
#'
#'#Generate the summary tibble
#'Thresholding_exploration_function(
#'   DATA = DATA_thresholded,
#'   Calculate_Density = TRUE,
#'   DATA_Area = DATA_AREA
#' )
#' }
#'
#' @export

Thresholding_exploration_function <-
  function(DATA,
           Calculate_Density = FALSE ,
           DATA_Area = NULL) {

    #Check arguments
    if(!is.logical(Calculate_Density)){
      stop("Calculate_Density must be a logical value")
    }

    #Import the data
    DATA <- DATA

    #Perform the quantification of variables across cases but taking into account if data is dichotomic or polychotomic
    Results <-
      purrr::map_dfr(unique(DATA$Subject_Names), function(Image){
        #Select variables
        DATA_vars <- DATA[-c(1:4)]

        #Select dichotomic variables
        DATA_Dichom <- dplyr::bind_cols(DATA[c(1:4)], DATA_vars[purrr::map_lgl(DATA_vars, function(Var) length(unique(Var)) <= 2)])
        #Select polychotomic variables
        DATA_Poly <- dplyr::bind_cols(DATA[c(1:4)], DATA_vars[purrr::map_lgl(DATA_vars, function(Var) length(unique(Var)) > 2)])

        #Start with dichotomic variables if these are present in the data
        if(length(names(DATA_Dichom)) > 4){
          #Start with dichomotic variables
          Results_Dichom <- DATA_Dichom %>% dplyr::filter(Subject_Names == Image) %>% tidyr::pivot_longer(-c(1:4)) %>%
            group_by(name) %>% dplyr::count(value) %>% dplyr::ungroup() %>% dplyr::filter(value) %>%
            tidyr::pivot_wider(names_from = name, values_from = n) %>% dplyr::select(-value)

          #If no cells are positive prime the tibble with the first column
          if(nrow(Results_Dichom) == 0){
            Results_Dichom<- tibble(value = 0)
            names(Results_Dichom) <- names(DATA_Dichom)[5]
          }
        }

        #Continue with polychotomic variables if these are present
        if(length(names(DATA_Poly)) > 4){
          Results_Poly <- DATA_Poly %>% dplyr::filter(Subject_Names == Image) %>% tidyr::pivot_longer(-c(1:4)) %>%
            group_by(name) %>% dplyr::count(value) %>% dplyr::ungroup() %>% dplyr::filter(value != 0) %>%
            dplyr::mutate(Variable = stringr::str_c(name, "Level", value, sep = "_")) %>% dplyr::select(-name,-value) %>%
            tidyr::pivot_wider(names_from = Variable, values_from = n)
          #If no cells are positive prime the tibble with the first column
          if (nrow(Results_Poly) == 0) {
            Results_Poly <- tibble(value = 0)
            names(Results_Poly) <- stringr::str_c(names(DATA_Poly)[5], "Level", "1", sep = "_")
          }
        }

        #If both polychotomic and dichotomic variables are present in th data return both results
        if(all(exists("Results_Dichom"), exists("Results_Poly"))){
          return(bind_cols(Results_Dichom, Results_Poly))
        }
        #If just dichotomic variables then return Dichom
        else if(exists("Results_Dichom")) {
          return(Results_Dichom)
        }
        #If just polychotomic variables then return Poly
        else if(exists("Results_Poly")){
          return(Results_Poly)
        }
      },
      .progress = list(clear = F,
                       name = "Counting positive cells by marker",
                       show_after = 2,
                       type = "iterator"))

    #Substitute NA values for 0
    Results[is.na(Results)] <- 0
    #Arrange columns in alphabetic order
    Results <- Results[sort(names(Results))]

    #Bind the subject names with the results
    Results <-dplyr::bind_cols(unique(DATA$Subject_Names), Results)
    names(Results)[1] <- "Subject_Names"

    #Calculate the total number of cells
    N_cells <-purrr::map_dbl(unique(DATA$Subject_Names), function(Image) {
      nrow(DATA %>% dplyr::filter(Subject_Names == Image))
    })

    #Calculate the proportion of cells that are positive
    print("Calculating cell percentages")
    Prop_tibble <- as_tibble(Results[-1] / N_cells)
    names(Prop_tibble) <-stringr::str_c("PROP_", names(Prop_tibble))

    Results$N_cells <- N_cells

    #If required calculate the density and return the results
    if(Calculate_Density){
      #Check arguments
      if(names(DATA_Area)[ncol(DATA_Area)] != "Area"){
        stop("The last column of the DATA_Area must be named Area")
      }

      print("Calculating cell densities")
      #Select cell counts
      For_density <- Results %>% dplyr::select(-N_cells, -contains("PROP_"))

      #Join the cell counts with the area tibbles
      For_density <-dplyr::left_join(For_density, DATA_Area, by = "Subject_Names")

      #Calculate densities
      Density_results <- as_tibble(For_density[c(2:(ncol(For_density)-1))] / For_density[[ncol(For_density)]])
      names(Density_results) <-stringr::str_c("Density_", names(Density_results))

      #Bind the results and return the tibble
      return(bind_cols(Results, Prop_tibble, Density_results, For_density[ncol(For_density)]))
    }

    #If not require return only the proportions
    else{
      #Bind the results and return the tibble
      return(bind_cols(Results, Prop_tibble))
    }
  }
