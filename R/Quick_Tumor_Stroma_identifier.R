#' Calculates tumor and stromal compartment based on tiles
#'
#' The function calculates tissue compartments (usually tumor and stroma) based on tile occupation by tumor cells.
#'
#' @param DATA A dataframe or tibble containing a column named 'Phenotype' containing cell phenotype labels.
#' @param Index_phenotype A character value indicating the cell phenotype to be used to calculate the tissue compartment.
#' @param Accuracy A numeric value indicating the size of the tile. Smaller sizes give more accurate result.
#' @param Min_cell_no A integer value indicating the minimum number of cells within a tile to consider the tile positive for a compartment.
#' @param Image_preview A character value indicating the name of the image to be used in the preview.
#' @param N_cores Integer. Number of cores to parallelize your computation.
#'
#'
#' @returns Returns a tibble with cell features and a column named 'Compartment' containing cell location.
#'
#' @export

Quick_Tumor_Stroma_identifier <-
  function(DATA_Phenotypes = NULL,
           Index_phenotype = NULL,
           Accuracy = NULL,
           Min_cell_no = NULL,
           Image_preview = NULL,
           N_cores = NULL){
    #Check arguments
    if(!all(is.character(DATA_Phenotypes), exists(DATA_Phenotypes, envir = .GlobalEnv))) stop("DATA_Phenotypes must be the name of an existing object")
    #Import all required Data from the environment
    DATA_Phenotypes <- get(DATA_Phenotypes, envir = globalenv())

    #Check more arguments
    if(!Index_phenotype %in% unique(DATA_Phenotypes$Phenotype)) stop(paste0("Index phenotype must be one of the following: ", stringr::str_c(unique(DATA_Phenotypes$Phenotype), collapse = ", ")))
    if(!all(is.numeric(Accuracy), Accuracy > 0)) stop("Accuracy must be a numeric value > 0")
    if(!all(is.numeric(Min_cell_no), Min_cell_no%%1 == 0, Min_cell_no > 0)) stop("Min_cell_no must be an integer value > 0")
    if(!Image_preview %in% unique(DATA_Phenotypes$Subject_Names)) stop(paste0(Image_preview, " not found in Subject_Names"))
    if(!all(N_cores >= 1 & N_cores%%1 == 0)) stop("N_cores must be an integer value > 0")



    print("Performing initial test")
    For_plot <- DATA_Phenotypes %>% dplyr::filter(Subject_Names == Image_preview)

    #Generate tiled image
    Test_plot_source <- For_plot %>% dplyr::filter(Phenotype == Index_phenotype) %>%
      ggplot(aes(x = X, y = Y)) + geom_bin2d(binwidth = Accuracy)

    plot(
      as_tibble(layer_data(Test_plot_source)) %>% dplyr::filter(value >= Min_cell_no) %>%
        ggplot() +
        geom_rect(aes(xmin = xmin, ymin = ymin, xmax = xmax, ymax = ymax), fill = "grey", color = "black") +
        theme_minimal() +
        geom_point(aes(x = X, y = Y), color = "red", alpha = 0.2, data = For_plot %>% dplyr::filter(Phenotype == Index_phenotype)) +
        theme(panel.grid = element_blank())
    )

    #Ask the user if the algorihtm should proceed
    answer <- menu(c("Proceed", "Abort"), title = "The plot depics the suggested area occupied by tumor cells. Should the algorithm proceed with the analysis?")
    #If user decides to stop then abort function and return stop message
    if(answer == 2) {
      stop("The function has been stopped. Please tune the parameters for a better result")
    }

    #If not proceed with computation
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
          #Select individual images
          Base <- DATA_Phenotypes %>% dplyr::filter(Subject_Names == Image)

          #Generate tiled image
          Plot_source <- Base %>% dplyr::filter(Phenotype == Index_phenotype) %>%
            ggplot(aes(x = X, y = Y)) + geom_bin2d(binwidth = Accuracy)

          #Filter layer information according to pre-specified requirements
          Interim <- as_tibble(layer_data(Plot_source)) %>% dplyr::filter(value >= Min_cell_no) %>% dplyr::select(xbin, ybin, x, y, xmin, xmax, ymin, ymax)

          #Assign cells to a specific compartment
          Final <- Base %>% dplyr::mutate(Compartment =purrr::map2_lgl(.x = Base$X, .y = Base$Y, function(.x, .y) {
            x_position <- (.x >= Interim$xmin) & (.x < Interim$xmax)
            y_position <- (.y >= Interim$ymin) & (.y < Interim$ymax)
            Final_pos <- any(x_position & y_position)
          })) %>% dplyr::mutate(Compartment = case_when(Compartment ~ "Tumor",
                                                       !Compartment ~ "Stroma"))
          return(Final)
        },
        .progress = TRUE)
      future::plan("future::sequential")
      gc()

      return(
        purrr::map_dfr(Final_list,dplyr::bind_rows)
      )
    }
  }
