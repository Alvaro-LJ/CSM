#' Calculates tumor and stromal compartment based on Local Indicators of Spatial Association
#'
#' The function calculates tissue compartments (usually tumor and stroma) based on Local Indicators of Spatial Association.
#' The function is based on the lisaClust::lisaClustfunction.
#'
#' @param DATA_Phenotypes A dataframe or tibble containing a column named 'Phenotype' containing cell phenotype labels.
#' @param Index_phenotype A character value indicating the cell phenotype to be used to calculate the tissue compartment.
#' @param Image_preview A character value indicating the name of the image to be used in the preview.
#' @param Association_Dist_min A numeric value indicating the minimum distance to test cell-cell spatial interaction.
#' @param Association_Dist_max A numeric value indicating the maximum distance to test cell-cell spatial interaction.
#' @param Type_of_LISA_function A character value indicating the type of LISA function used: Either 'L' or 'K'.
#' @param Window_type A character value indicating the interaction window shape: ’square’, ’convex’ or ’concave’.
#' @param N_cores Integer. Number of cores to parallelize your computation.
#'
#' @returns Returns a tibble with cell features and a column named 'Compartment' containing cell location.
#'
#' @examples
#' \dontrun{
#' LISA_Tumor_Stroma_identifier(
#'    DATA_Phenotypes = CSM_Phenotypecell_test,
#'    Index_phenotype = "TUMOR",
#'    Association_Dist_min = 5,
#'    Association_Dist_max = 50,
#'    Type_of_LISA_function = "L",
#'    Window_type = "convex",
#'    N_cores = 1
#')
#' }
#'
#' @export

LISA_Tumor_Stroma_identifier <-
  function(DATA_Phenotypes = NULL,
           Index_phenotype = NULL,
           Image_preview = NULL,
           Association_Dist_min = NULL,
           Association_Dist_max = NULL,
           Type_of_LISA_function = NULL,
           Window_type = NULL,
           N_cores = 1){

    #Check suggested packages
    {
      if(!requireNamespace("SpatialExperiment", quietly = TRUE)) stop(
        paste0("SpatialExperiment Bioconductor package is required to execute the function. Please install using the following code: ",
               expression({
                 if (!require("BiocManager", quietly = TRUE))
                   install.packages("BiocManager")

                 BiocManager::install("SpatialExperiment")
               })
        )
      )
      if(!requireNamespace("lisaClust", quietly = TRUE)) stop(
        paste0("lisaClust Bioconductor package is required to execute the function. Please install using the following code: ",
               expression({
                 if (!require("BiocManager", quietly = TRUE))
                   install.packages("BiocManager")

                 BiocManager::install("lisaClust")
               })
        )
      )
    }
    #Import all required Data
    DATA_Phenotypes <- DATA_Phenotypes

    #Check more arguments
    if(!Index_phenotype %in% unique(DATA_Phenotypes$Phenotype)) stop(paste0("Index phenotype must be one of the following: ", stringr::str_c(unique(DATA_Phenotypes$Phenotype), collapse = ", ")))
    if(is.null(Image_preview)) Image_preview <- sample(unique(DATA_Phenotypes$Subject_Names), size = 1)
    if(!Image_preview %in% unique(DATA_Phenotypes$Subject_Names)) stop(paste0(Image_preview, " not found in Subject_Names"))

    if(!all(is.numeric(Association_Dist_min), Association_Dist_min > 0)) stop("Association_Dist_min must be a numeric value > 0")
    if(!all(is.numeric(Association_Dist_max), Association_Dist_max > 0)) stop("Association_Dist_max must be a numeric value > 0")
    if(!Association_Dist_max > Association_Dist_min) stop("Association_Dist_max must be greater than Association_Dist_min")
    if(!Type_of_LISA_function %in% c("L", "K")) stop("Type_of_LISA_function must be one of the following: L, K")
    if(!Window_type %in% c("square", "convex", "concave")) stop("Window_type must be one of the following: square, convex, concave")
    if(!all(N_cores >= 1 & N_cores%%1 == 0)) stop("N_cores must be an integer value > 0")

    DATA_Phenotypes <- DATA_Phenotypes %>% dplyr::rename("x" = "X", "y" = "Y")

    #Select Image to be previewed
    For_test <- DATA_Phenotypes %>% dplyr::filter(Subject_Names == Image_preview)

    #Change all non-tumor cells to OTHER phenotype
    For_test <- For_test %>%dplyr::mutate(Cell_ref = case_when(Phenotype == Index_phenotype ~ "TUMOR",
                                                               TRUE ~ "OTHER"))


    #Turn our data into a SpatialExperiment object
    Spatial_object_Test <- SpatialExperiment::SpatialExperiment(sample_id = Image_preview,
                                                                colData = For_test,
                                                                spatialCoordsNames = c("x", "y")
    )

    #Perform LISA clustering
    Lisa_test <- lisaClust::lisaClust(Spatial_object_Test,
                                      k = 3,
                                      Rs = seq(from = Association_Dist_min, to = Association_Dist_max, by = floor((Association_Dist_max - Association_Dist_min)/5)),
                                      imageID = "Subject_Names",
                                      cellType = "Cell_ref",
                                      spatialCoords = c("x", "y"),
                                      window = Window_type,
                                      lisaFunc = Type_of_LISA_function)

    #Add results to original data
    For_test <- For_test %>% dplyr::mutate(Region = SummarizedExperiment::colData(Lisa_test)$region)

    rm(Lisa_test)
    gc()

    #Generate the tumor proportion results and add labels
    Proportion_results <- as_tibble(prop.table(table(For_test$Region, For_test$Cell_ref), margin = 1), .name_repair = "universal")
    names(Proportion_results) <- c("Region", "Cell", "Proportion")
    Proportion_results <- Proportion_results %>% tidyr::pivot_wider(names_from = Cell, values_from = Proportion)

    #Generate name results
    Tumor_region <- Proportion_results[[which(Proportion_results[["TUMOR"]] == max(Proportion_results[["TUMOR"]])), 1]]
    Stroma_region <- Proportion_results[[which(Proportion_results[["TUMOR"]] == min(Proportion_results[["TUMOR"]])), 1]]
    Interface_region <- Proportion_results[["Region"]][which((Proportion_results[["Region"]] != Tumor_region & Proportion_results[["Region"]] != Stroma_region))]

    For_test <- For_test %>%dplyr::mutate(Region = case_when(Region == Tumor_region ~ "Tumor",
                                                             Region == Stroma_region ~ "Stroma",
                                                             Region == Interface_region ~ "Border"))

    PLOT1 <- For_test %>% ggplot(aes(x = x, y = y, color = Cell_ref )) +
      geom_point(size = 1.5) +
      theme_minimal() +
      scale_color_manual(values = c(alpha("grey", 0.05), alpha("black", 0.8)))+
      ggtitle("Tumor / Non tumor cells") +
      theme(panel.grid = element_blank(),
            axis.text = element_blank(),
            axis.title = element_blank(),
            legend.position = "bottom",
            plot.title = element_text(hjust = 0.5))

    PLOT2 <-  For_test %>% ggplot(aes(x = x, y = y, color = Region)) +
      geom_point(size = 1.5) +
      theme_minimal() +
      scale_color_manual(values = c(alpha("grey", 0.05), alpha("blue", 0.5), alpha("red", 0.5)))+
      ggtitle("Tissue segmented cells") +
      theme(panel.grid = element_blank(),
            axis.text = element_blank(),
            axis.title = element_blank(),
            legend.position = "bottom",
            plot.title = element_text(hjust = 0.5))
    plot(patchwork::wrap_plots(PLOT1, PLOT2, nrow = 1))

    #Ask user if they want to proceed
    User_answer <- menu(choices = c("Proceed", "Abort"), title = "Evaluate images and decide if analysis should proceed")
    if(User_answer == 2){
      stop("Algorithm has been aborted. Refine parameters to improve performance (Association_Dist, LISA_function, Window_type)")
    }
    #If they want to continue
    else{
      #save exit function if parallelization fails
      on.exit({
        future::plan("future::sequential")
        gc()
      })

      future::plan("future::multisession", workers = N_cores)
      options(future.globals.maxSize = Inf, future.rng.onMisuse = "ignore")
      furrr::furrr_options(scheduling = Inf)

      Final_list <-
        suppressMessages(furrr::future_map(unique(DATA_Phenotypes$Subject_Names), function(Image){
          #Select Image to be previewed
          For_test <- DATA_Phenotypes %>% dplyr::filter(Subject_Names == Image)

          #Change all non-tumor cells to OTHER phenotype
          For_test <- For_test %>%dplyr::mutate(Cell_ref = case_when(Phenotype == Index_phenotype ~ "TUMOR",
                                                                     TRUE ~ "OTHER"))
          #Turn our data into a SpatialExperiment object
          Spatial_object_Test <- SpatialExperiment::SpatialExperiment(sample_id = Image,
                                                                      colData = For_test,
                                                                      spatialCoordsNames = names(For_test)[2:3]
          )

          #Perform LISA clustering
          Lisa_test <- lisaClust::lisaClust(Spatial_object_Test,
                                            k = 3,
                                            Rs = seq(from = Association_Dist_min, to = Association_Dist_max, by = floor((Association_Dist_max - Association_Dist_min)/5)),
                                            imageID = "Subject_Names",
                                            cellType = "Cell_ref",
                                            spatialCoords = c("x", "y"),
                                            window = Window_type,
                                            lisaFunc = Type_of_LISA_function)

          #Add results to original data
          For_test <- For_test %>%dplyr::mutate(Compartment = SummarizedExperiment::colData(Lisa_test)$region)

          rm(Lisa_test)
          gc()

          #Generate the tumor proportion results and add labels
          Proportion_results <- as_tibble(prop.table(table(For_test$Compartment, For_test$Cell_ref), margin = 1), .name_repair = "universal")
          names(Proportion_results) <- c("Region", "Cell", "Proportion")
          Proportion_results <- Proportion_results %>% tidyr::pivot_wider(names_from = Cell, values_from = Proportion)

          #Generate name results
          Tumor_region <- Proportion_results[[which(Proportion_results[["TUMOR"]] == max(Proportion_results[["TUMOR"]])), 1]]
          Stroma_region <- Proportion_results[[which(Proportion_results[["TUMOR"]] == min(Proportion_results[["TUMOR"]])), 1]]
          Interface_region <- Proportion_results[["Region"]][which((Proportion_results[["Region"]] != Tumor_region & Proportion_results[["Region"]] != Stroma_region))]

          For_test <- For_test %>%dplyr::mutate(Compartment = dplyr::case_when(Compartment == Tumor_region ~ "Tumor",
                                                                               Compartment == Stroma_region ~ "Stroma",
                                                                               Compartment == Interface_region ~ "Border"))
          For_test <- For_test %>% dplyr::select(-Cell_ref)
          return(For_test)
        },
        .progress = TRUE))


      future::plan("future::sequential")
      gc()

      FINAL_result <-purrr::map_dfr(Final_list,dplyr::bind_rows)
      return(FINAL_result %>% dplyr::rename("X" = "x", "Y" = "y"))
    }
  }
