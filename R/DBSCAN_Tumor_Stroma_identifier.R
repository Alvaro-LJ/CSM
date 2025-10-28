#' Calculates tumor and stromal compartment based on density based clustering
#'
#' The function calculates tissue compartments (usually tumor and stroma) based on density based clustering. Cells are considered to be in the tissue compartment if they are clustered in proximity to index cell types.
#' The function is based on the dbscan::dbscan function.
#'
#' @param DATA_Phenotypes A dataframe or tibble containing a column named 'Phenotype' containing cell phenotype labels.
#' @param Index_phenotype A character value indicating the cell phenotype to be used to calculate the tissue compartment.
#' @param Image_preview A character value indicating the name of the image to be used in the preview.
#' @param Min_cells A integer value indicating the minimum number of cells required to be a cluster.
#' @param Distance_radius A numeric value indicating the distance to be sampled.
#' @param N_cores Integer. Number of cores to parallelize your computation.
#'
#' @returns Returns a tibble with cell features and a column named 'Compartment' containing cell location.
#'
#'@examples
#'\dontrun{
#'DBSCAN_Tumor_Stroma_identifier(
#'     DATA_Phenotypes = CSM_Phenotypecell_test,
#'     Index_phenotype = "TUMOR",
#'     Min_cells = 3,
#'     Distance_radius = 50,
#'     N_cores = 1
#')
#'}
#'
#' @export

DBSCAN_Tumor_Stroma_identifier <-
  function(DATA_Phenotypes,
           Index_phenotype,
           Image_preview = NULL,
           Min_cells,
           Distance_radius,
           N_cores = 1
  ){

    #Check suggested packages
    if(!requireNamespace("dbscan", quietly = FALSE)) stop(
      paste0("dbscanCRAN package is required to execute the function. Please install using the following code: ",
             expression(install.packages("dbscan")))
    )

    #Import all required Data
    DATA_Phenotypes <- DATA_Phenotypes

    #Check more arguments
    if(!Index_phenotype %in% unique(DATA_Phenotypes$Phenotype)) stop(paste0("Index phenotype must be one of the following: ", stringr::str_c(unique(DATA_Phenotypes$Phenotype), collapse = ", ")))
    if(!all(is.numeric(Min_cells), Min_cells%%1 == 0, Min_cells > 0)) stop("Min_cells must be an integer value > 0")
    if(is.null(Image_preview)) Image_preview <- sample(unique(DATA_Phenotypes$Subject_Names), size = 1)
    if(!Image_preview %in% unique(DATA_Phenotypes$Subject_Names)) stop(paste0(Image_preview, " not found in Subject_Names"))
    if(!all(is.numeric(Distance_radius), Distance_radius > 0)) stop("Distance_radius must be an integer value > 0")
    if(!all(N_cores >= 1 & N_cores%%1 == 0)) stop("N_cores must be an integer value > 0")



    print("Performing initial test")
    #First generate the DBSCAN model for the Image preview
    #Select Image to be previewed
    For_test <- DATA_Phenotypes %>% dplyr::filter(Subject_Names == Image_preview)
    #Obtain individual data sets containing Index cell and other cells
    For_test_Index <- For_test %>% dplyr::filter(Phenotype == Index_phenotype)
    For_test_no_Index <- For_test %>% dplyr::filter(Phenotype != Index_phenotype)

    #Calculate the DBSCAN model with the required parameters
    DB_results <- dbscan::dbscan(For_test_Index[c("X", "Y")], eps = Distance_radius, minPts = Min_cells, borderPoints = FALSE)

    #Obtain position in the tumor and the stroma for both index cells and other cells
    For_test_Index <- For_test_Index %>% dplyr::mutate(Cluster = as.character(DB_results$cluster)) %>%
      dplyr::mutate(Cluster = dplyr::case_when(Cluster == "0" ~ "Stroma",
                                               T ~ "Tumor"))
    For_test_no_Index <- For_test_no_Index %>%
      dplyr::mutate(Cluster = predict(DB_results, newdata = For_test_no_Index[c("X", "Y")], data = For_test_Index[c("X", "Y")])) %>%
      dplyr::mutate(Cluster = dplyr::case_when(Cluster == "0" ~ "Stroma",
                                               T ~ "Tumor"))

    #Obtain final result
    Test_final <- dplyr::bind_rows(For_test_Index, For_test_no_Index) %>% dplyr::mutate(Cell_arrange = as.numeric(sub(".*CELL_([^_]+)__.*", "\\1", c(For_test_Index$Cell_no, For_test_no_Index$Cell_no)))) %>%
      dplyr::arrange(Cell_arrange) %>% dplyr::select(-Cell_arrange)

    #Obtain three plots
    PLOT1 <- For_test %>%dplyr::mutate(Phenotype = case_when(Phenotype == Index_phenotype ~ "Tumor Cell",
                                                             TRUE ~ "OTHER")) %>%
      ggplot(aes(x = X, y = Y, color = Phenotype)) +
      geom_point(size = 1) +
      theme_minimal() +
      scale_color_manual("", values = c(alpha("grey", 0.1), alpha("black", 0.5))) +
      ggtitle("Cells in analysis") +
      theme(panel.grid = element_blank(),
            axis.text = element_blank(),
            axis.title = element_blank(),
            legend.position = "bottom",
            plot.title = element_text(hjust = 0.5))

    PLOT2 <- For_test_Index %>% dplyr::filter(Cluster != "Stroma") %>%
      ggplot(aes(x = X, y = Y)) +
      geom_point() +
      theme_minimal() +
      ggtitle("Tumor cells that are clustered") +
      theme(panel.grid = element_blank(),
            axis.text = element_blank(),
            axis.title = element_blank(),
            legend.position = "bottom",
            plot.title = element_text(hjust = 0.5))

    PLOT3 <- Test_final %>%
      ggplot(aes(x = X, y = Y, color = Cluster)) +
      geom_point(size = 1) +
      theme_minimal() +
      scale_color_manual("", values = c(alpha("grey", 0.1), alpha("black", 0.5))) +
      ggtitle("Final Result") +
      theme(panel.grid = element_blank(),
            axis.text = element_blank(),
            axis.title = element_blank(),
            legend.position = "bottom",
            plot.title = element_text(hjust = 0.5))

    #plot final result
    plot(patchwork::wrap_plots(PLOT1, PLOT2, nrow = 1))

    #Ask user if they want to proceed
    User_answer <- menu(choices = c("Proceed", "Abort"), title = "Evaluate images and decide if analysis should proceed")
    if(User_answer == 2){
      stop("Algorithm has been aborted. Refine Min_cells and Distance radius analysis")
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
        furrr::future_map(unique(DATA_Phenotypes$Subject_Names), function(Image){
          #Select Image to be previewed
          For_test <- DATA_Phenotypes %>% dplyr::filter(Subject_Names == Image)
          #Obtain individual data sets containing Index cell and other cells
          For_test_Index <- For_test %>% dplyr::filter(Phenotype == Index_phenotype)
          For_test_no_Index <- For_test %>% dplyr::filter(Phenotype != Index_phenotype)

          #Calculate the DBSCAN model with the required parameters
          DB_results <- dbscan::dbscan(For_test_Index[c("X", "Y")], eps = Distance_radius, minPts = Min_cells, borderPoints = FALSE)

          #Obtain position in the tumor and the stroma for both index cells and other cells
          For_test_Index <- For_test_Index %>%dplyr::mutate(Compartment = as.character(DB_results$cluster)) %>%
            dplyr::mutate(Compartment = dplyr::case_when(Compartment == "0" ~ "Stroma",
                                                         T ~ "Tumor"))
          For_test_no_Index <- For_test_no_Index %>%
            dplyr::mutate(Compartment = predict(DB_results, newdata = For_test_no_Index[c("X", "Y")], data = For_test_Index[c("X", "Y")])) %>%
            dplyr::mutate(Compartment = dplyr::case_when(Compartment == "0" ~ "Stroma",
                                                         T ~ "Tumor"))

          #Obtain final result
          Test_final <- dplyr::bind_rows(For_test_Index, For_test_no_Index) %>% dplyr::mutate(Cell_arrange = as.numeric(sub(".*CELL_([^_]+)__.*", "\\1", c(For_test_Index$Cell_no, For_test_no_Index$Cell_no)))) %>%
            dplyr::arrange(Cell_arrange) %>% dplyr::select(-Cell_arrange)
        },
        .progress = TRUE)
      future::plan("future::sequential")
      gc()

      return(
        purrr::map_dfr(Final_list,dplyr::bind_rows)
      )
    }
  }
