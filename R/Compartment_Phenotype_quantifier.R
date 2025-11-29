#' Quantifies cell proportions and cell densities by compartment
#'
#' The function calculates cell proportions and optionally cell densities in tumor and stromal compartments.
#' Compartment areas must be first computed using [Advanced_Tumor_Stroma_identifier()] function.
#'
#' @param DATA A list containing cell compartment location and compartment area generated with [Advanced_Tumor_Stroma_identifier()].
#' @param Calculate_Density A logical value indicating if density should be calculated. Currently computing cell density in border is not supported.
#' @param DATA_Area If density needs to be calculated, a tibble or dataframe containing total tissue size information. It may be generated using [Image_size_calculator()].
#'
#' @returns Returns tibble with information of cell percentages and (optionally) densities by image and compartment.
#'
#' @seealso [Compartment_Phenotype_quantifier()], [Image_size_calculator()].
#'
#' @examples
#' \dontrun{
#' #Calculate tumor, stroma and border compartments----------------------------
#' DATA_structured <-
#' Advanced_Tumor_Stroma_identifier(
#'     DATA_Phenotypes = CSM_Phenotypecell_test,
#'     Index_phenotype = "TUMOR",
#'
#'     Filtering_Method = "DBSCAN",
#'     Min_cell_no = 10,
#'     Distance_radius = 100,
#'
#'     Hull_ratio = 0.05,
#'     Calculate_border = TRUE,
#'     Dist_to_border = 10
#' )
#'
#' #Calculate cell percentages in each compartment-----------------------------
#'Compartment_Phenotype_quantifier(
#'       DATA = DATA_structured,
#'       Calculate_Density = FALSE
#')
#' }
#'
#' @export

Compartment_Phenotype_quantifier <-
  function(DATA,
           Calculate_Density = FALSE,
           DATA_Area = NULL){
    DATA <- DATA
    #Check that DATA has been created with the advanced tumor stroma identifier
    if(!identical(names(DATA), c("DATA_Phenotypes", "DATA_Compartment_Area"))){
      stop("DATA must have been created with the Advanced_Tumor_Stroma_identifier")
    }

    if(!is.logical(Calculate_Density)){
      stop("Calculate_Density must be a logical value")
    }
    if(Calculate_Density){
      if(names(DATA_Area)[ncol(DATA_Area)] != "Area"){
        stop("The last column of the DATA_Area must be named Area")
      }
    }


    #Import our DATA
    DATA_Phenotypes <- DATA[[1]]
    DATA_Tumor_Area <- DATA[[2]]
    names(DATA_Tumor_Area)[2] <- "Tumor_Area"


    #Generate a list of the cell types divided by compartment
    Cells_by_compartment <-purrr::map(unique(DATA_Phenotypes$Compartment), function(Unique_Compartment){
      DATA_Phenotypes %>% dplyr::filter(Compartment == Unique_Compartment)
    })
    names(Cells_by_compartment) <- unique(DATA_Phenotypes$Compartment)

    #Calculate the cell counts and percentages by compartment
    Cell_prop_by_compartment <-purrr::map(Cells_by_compartment, function(Unique_Compartment){
      #Obtain the number of cells by image according to the phenotypes
      Results <- Unique_Compartment %>% group_by(Subject_Names, Phenotype) %>% dplyr::count() %>%dplyr::ungroup() %>%
        tidyr::pivot_wider(names_from = Phenotype, values_from = n)
      Results[is.na(Results)] <- 0

      #Calculate the number of total cells
      N_cells <- apply(Results[-1], MARGIN = 1, sum)

      #Calculate proportions
      Prop_tibble <- as_tibble(Results[-1] / N_cells)
      names(Prop_tibble) <-stringr::str_c("PROP_", names(Prop_tibble))
      Results$N_cells <- N_cells

      #Arrange columns by alphabetic order
      Results <- Results[c("Subject_Names", sort(names(Results)[-c(1, ncol(Results))]), "N_cells")]
      Prop_tibble <- Prop_tibble[sort(names(Prop_tibble))]

      return(bind_cols(Results, Prop_tibble))
    })

    #If densities need to be calculated execute the following code
    if(Calculate_Density){
      #Check arguments
      if(names(DATA_Area)[ncol(DATA_Area)] != "Area"){
        stop("The last column of the DATA_Area must be named Area")
      }

      #Obtain the total area and tumor area to calculate Tumor / Stroma areas
      DATA_Overall_Area <- DATA_AREA
      DATA_Areas <-dplyr::left_join(DATA_Overall_Area, DATA_Tumor_Area, by = "Subject_Names") %>%dplyr::mutate(Stroma_Area = Area - Tumor_Area)


      #Start with tumor densities
      TUMOR <-dplyr::left_join(Cell_prop_by_compartment[["Tumor"]], DATA_Areas[c("Subject_Names", "Tumor_Area")], by = "Subject_Names")

      #Select cell counts
      For_density_TUMOR <- TUMOR %>% dplyr::select(-N_cells, -contains("PROP_"))

      #Calculate tumor densities
      Density_results_TUMOR <- as_tibble(For_density_TUMOR[c(2:(ncol(For_density_TUMOR)-1))] / For_density_TUMOR[[ncol(For_density_TUMOR)]])
      names(Density_results_TUMOR) <-stringr::str_c("Density_", names(Density_results_TUMOR))
      TUMOR <-dplyr::bind_cols(TUMOR, Density_results_TUMOR)


      #Continue with stromal densities
      STROMA <-dplyr::left_join(Cell_prop_by_compartment[["Stroma"]], DATA_Areas[c("Subject_Names", "Stroma_Area")], by = "Subject_Names")

      #Select cell counts
      For_density_STROMA <- STROMA %>% dplyr::select(-N_cells, -contains("PROP_"))

      #Calculate stromal densities
      Density_results_STROMA <- as_tibble(For_density_STROMA[c(2:(ncol(For_density_STROMA)-1))] / For_density_STROMA[[ncol(For_density_STROMA)]])
      names(Density_results_STROMA) <-stringr::str_c("Density_", names(Density_results_STROMA))
      STROMA <-dplyr::bind_cols(STROMA, Density_results_STROMA)

      Cell_prop_by_compartment[["Tumor"]] <- TUMOR
      Cell_prop_by_compartment[["Stroma"]] <- STROMA

      #If the input data has stromal compartment, then print a warning
      if(length(Cell_prop_by_compartment) > 2) warning("Border densities are not stimated in the current CSM version")
    }

    #Return the final result
    return(Cell_prop_by_compartment)
  }
