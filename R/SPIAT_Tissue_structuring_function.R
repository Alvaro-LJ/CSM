#' Generates tissue structures according to a cell phenotype label using the SPIAT package.
#'
#' The function generates tissue structure according to a cell phenotype label (usually tumor cells). Then all the cells are asigned a label according to their distance to the tissue structures.
#'
#' @param N_cores Integer. Number of cores to parallelize your computation.
#' @param DATA_SPIAT A list containing SPIAT objects generated using [SPIAT_object_generator()].
#' @param DATA_Phenotypes A dataframe or tibble containing cell feature data with a column named 'Phenotype' containing cell phenotype labels.
#' @param Cell_type_to_define_cluster A character value indicating which cell phenotype label should be used  to build tissue structures.
#' @param Minimum_number_cells_cluster A integer value indicating the minimum number of cells required to be considered a structure.
#' @param Cell_types_of_interest Cell types to be related to these structures.
#' @param Simplify_result A logical value indicating if the final result should be simplified.
#'
#' @returns If Simplify_result is TRUE, the function returns a tibble with simplify cell location labels.
#' If Simplify_result is FALSE, the function returns a list containing cell location labels for every image.
#'
#' @seealso [SPIAT_object_generator()].
#'
#' @examples
#' #Generate SPIAT object list----------------------------
#' DATA_SPIAT <-
#' SPIAT_object_generator(
#'     DATA_Intensities = CSM_Arrangedcellfeaturedata_test,
#'     DATA_Phenotypes = CSM_Phenotypecell_test
#' )
#'
#' #Calculate tissue structures-------------
#' SPIAT_Tissue_structuring_function(
#'     N_cores = 2,
#'     DATA_SPIAT = DATA_SPIAT,
#'     DATA_Phenotypes = CSM_Phenotypecell_test,
#'     Cell_type_to_define_cluster = "TUMOR",
#'     Minimum_number_cells_cluster = 5,
#'     Cell_types_of_interest = c("TUMOR", "CD8_GZMBneg", "CD8_GZMBpos", "OTHER"),
#'     Layers_margin = 5,
#'     Simplify_result = T
#')
#'
#' @export

SPIAT_Tissue_structuring_function <-
  function(N_cores = NULL,
           DATA_SPIAT = NULL,
           DATA_Phenotypes = NULL,
           Cell_type_to_define_cluster = NULL,
           Minimum_number_cells_cluster = NULL,
           Cell_types_of_interest = NULL,
           Layers_margin = NULL,
           Simplify_result = NULL
  ){

    {
      #Check suggested packages
      if(!requireNamespace("SPIAT", quietly = TRUE)) stop(
        paste0("SPIAT Bioconductor package is required to execute the function. Please install using the following code: ",
               expression({
                 if (!require("BiocManager", quietly = TRUE))
                   install.packages("BiocManager")

                 BiocManager::install("SPIAT")
               })
        )
      )
      if(!requireNamespace("alphahull", quietly = FALSE)) stop(
        paste0("alphahull CRAN package is required to execute the function. Please install using the following code: ",
               expression(install.packages("alphahull")))
      )
    }



    #Check arguments
    if(!all(N_cores >= 1 & N_cores%%1 == 0)) stop("N_cores must be an integer value > 0")

    DATA_SPIAT <- DATA_SPIAT
    if(!all(purrr::map_lgl(DATA_SPIAT, function(Image) class(Image) == "SpatialExperiment"))) stop("DATA_SPIAT must be created using the SPIAT_object_generator function")

    DATA_Phenotypes <- DATA_Phenotypes
    if(!identical(names(DATA_Phenotypes)[c(1:4)], c("Cell_no", "X", "Y", "Subject_Names"))) stop("DATA_Phenotypes must have an adequate format")
    if(!"Phenotype" %in% names(DATA_Phenotypes)) stop("DATA_Phenotypes must have an adequate format")

    if(!all(Cell_type_to_define_cluster %in% unique(unlist(purrr::map(DATA_SPIAT, function(x) x@colData@listData$Phenotype))))){
      stop(paste0(Cell_type_to_define_cluster, " not found in DATA_SPIAT"))
    }
    if(!all(Cell_types_of_interest %in% unique(unlist(purrr::map(DATA_SPIAT, function(x) x@colData@listData$Phenotype))))) {
      Absent_cell_types <- Cell_types_of_interest[!Cell_types_of_interest %in% unique(unlist(purrr::map(DATA_SPIAT, function(x) x@colData@listData$Phenotype)))]
      stop(paste0(stringr::str_c(Absent_cell_types, collapse = ", "), " not found in DATA_SPIAT"))
    }
    if(!all(is.numeric(Minimum_number_cells_cluster), Minimum_number_cells_cluster > 0, Minimum_number_cells_cluster%%1 == 0)){
      stop("Minimum_number_cell_cluster must be a integer value > 0")
    }
    if(!all(is.numeric(Layers_margin), Layers_margin%%1 == 0)) stop("Layers_margin must be a integer value")
    if(!is.logical(Simplify_result)) stop("Simplify_result must be a logical value")

    #Proceed with execution
    Phenotypes_list <-purrr::map(names(DATA_SPIAT), function(Image) DATA_Phenotypes %>% dplyr::filter(Subject_Names == Image))

    #save exit function if parallelization fails
    on.exit({
      future::plan("future::sequential")
      gc()
    })

    future::plan("future::multisession", workers = N_cores)
    options(future.globals.maxSize = Inf, future.rng.onMisuse = "ignore")
    furrr::furrr_options(scheduling = Inf)
    RESULTS <-
      suppressMessages(
        furrr::future_map(seq_along(1:length(DATA_SPIAT)), function(Image) {
          library(SPIAT)


          pdf(NULL) #Avoids anoying PDF generation


          Formatted_image <- DATA_SPIAT[[Image]] #Get formatted image

          #Calculate the ratio of target cells that are in a border
          Border_Ratio <- R_BC(Formatted_image,
                               cell_type_of_interest = Cell_type_to_define_cluster, #Cell type that defines the tumor
                               feature_colname = "Phenotype")

          #Calculate the clusters
          formatted_border <- identify_bordering_cells(Formatted_image,
                                                       reference_cell = Cell_type_to_define_cluster,
                                                       feature_colname = "Phenotype",
                                                       ahull_alpha = NULL, #Controls the size of cell clusters
                                                       n_to_exclude = Minimum_number_cells_cluster, #Controls the minimum amount of cells to be a cluster
                                                       plot_final_border = FALSE)

          #Calculate the number of cluster islands found in the image
          N_tumor_clusters <- attr(formatted_border, "n_of_clusters")

          #Now we calculate distances to this borders defined previously
          formatted_distance <- calculate_distance_to_margin(formatted_border)

          #Assign a structure to the tissue
          formatted_structure <- define_structure(formatted_distance,
                                                  cell_types_of_interest = Cell_types_of_interest, #Select our cell types of interest
                                                  feature_colname = "Phenotype",
                                                  n_margin_layers = Layers_margin #Layers of cells that define the internal and external border of the tumor margin
          )
          categories <- unique(formatted_structure$Structure)


          #Calculate proportions of cells in each structure
          immune_proportions <- calculate_proportions_of_cells_in_structure(
            spe_object = formatted_structure,
            cell_types_of_interest = Cell_types_of_interest,
            feature_colname = "Phenotype")

          #Calculate distances
          immune_distances <- calculate_summary_distances_of_cells_to_borders(
            spe_object = formatted_structure,
            cell_types_of_interest = Cell_types_of_interest,
            feature_colname = "Phenotype")

          dev.off()

          #Assign location to cell types
          return(
            list(Cells = Phenotypes_list[[Image]] %>% dplyr::mutate(Location = formatted_structure$Structure),
                 Border_Ratio = Border_Ratio,
                 N_tumor_clusters = N_tumor_clusters,
                 Immune_proportions = immune_proportions,
                 Immune_distances = immune_distances)
          )
        },
        .progress = TRUE)
      )
    future::plan("future::sequential")
    gc()

    #Change names of the results
    names(RESULTS) <- names(DATA_SPIAT)

    #Simplify the results if desired
    if(Simplify_result) {
      return(purrr::map_dfr(RESULTS, function(x) x[[1]]) %>% dplyr::select(1:4, Phenotype, Location) %>%
               dplyr::mutate(Location = case_when(Location == "Border" ~ "Border",
                                                  Location == "Infiltrated.CoI" ~ "Core",
                                                  Location == "Inside" ~ "Core",
                                                  Location == "Stromal.CoI" ~ "Stroma",
                                                  Location == "Outside" ~ "Stroma",
                                                  Location == "Internal.margin.CoI" ~ "Internal_Border",
                                                  Location == "Internal.margin" ~ "Internal_Border",
                                                  Location == "External.margin.CoI" ~ "External_Border",
                                                  Location == "External.margin" ~ "External_Border")))
    }

    #If not return the complete results
    else {return(RESULTS)}
  }
