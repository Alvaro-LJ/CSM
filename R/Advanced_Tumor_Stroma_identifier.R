#' Calculates tumor and stromal compartment based on concave hull
#'
#' The function calculates tissue compartments (usually tumor and stroma) based on concave hull.
#' The border is calculated based on a user defined distance to the tissue structure. Before the concave hull is calculated a filtering method must be applied to remove non-spatially clustered cells.
#' The concave hull is calculated using th sf R package.
#' DBscan filtering is based on the dbscan::dbscan function.
#'
#' The function results can be used to run the [Compartment_Phenotype_quantifier()] function.
#'
#' @param DATA_Phenotypes A dataframe or tibble containing a column named 'Phenotype' containing cell phenotype labels.
#' @param Index_phenotype A character value indicating the cell phenotype to be used to calculate the tissue compartment.
#' @param Filtering_Method A character value indicating the filtering method type. Either "Tiling" or "DBSCAN".
#'
#' @param Accuracy If Filtering_Method is 'Tiling': A numeric value indicating the size of the tile. Smaller sizes give more accurate result.
#' @param Min_cells If Filtering_Method is 'DBSCAN': A integer value indicating the minimum number of cells required to be a cluster.
#' @param Distance_radius If Filtering_Method is 'DBSCAN': A numeric value indicating the distance to be sampled.
#'
#' @param Image_preview A character value indicating the name of the image to be used in the preview.
#' @param N_cores Integer. Number of cores to parallelize your computation.
#'
#' @returns Returns a list with a tibble with cell features including a column named 'Compartment' containing cell location and a tibble with tumor compartment area.
#'
#' @seealso [Compartment_Phenotype_quantifier()]
#'
#' @examples
#' \dontrun{
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
#' }
#'
#' @import dplyr
#' @export

Advanced_Tumor_Stroma_identifier <-
  function(DATA_Phenotypes = NULL,
           Index_phenotype = NULL,
           Filtering_Method = NULL,

           Accuracy = NULL,
           Min_cell_no = NULL,
           Distance_radius = NULL,



           Hull_ratio = NULL,
           Calculate_border = NULL,
           Dist_to_border = NULL,

           Image_preview = NULL,
           N_cores = 1){
    #check suggested packages
    {
      if(!requireNamespace("sf", quietly = FALSE)) stop(
        paste0("sf CRAN package is required to execute the function. Please install using the following code: ",
               expression(install.packages("sf")))
      )
      if(Filtering_Method == "DBSCAN"){
        if(!requireNamespace("dbscan", quietly = FALSE)) stop(
          paste0("dbscan CRAN package is required to execute the function. Please install using the following code: ",
                 expression(install.packages("dbscan")))
        )
      }
    }
    #Import all required Data
    DATA_Phenotypes <- DATA_Phenotypes

    #Check more arguments
    if(!Index_phenotype %in% unique(DATA_Phenotypes$Phenotype)) stop(paste0("Index phenotype must be one of the following: ", stringr::str_c(unique(DATA_Phenotypes$Phenotype), collapse = ", ")))
    if(!Filtering_Method %in% c("Tiling", "DBSCAN")) stop("Filtering_Method must be one of the following: Tiling, DBSCAN")
    if(Filtering_Method == "Tiling"){
      if(!all(is.numeric(Accuracy), Accuracy > 0)) stop("Accuracy must be a numeric value > 0")
      if(!all(is.numeric(Min_cell_no), Min_cell_no%%1 == 0, Min_cell_no > 0)) stop("Min_cell_no must be an integer value > 0")
    }
    if(Filtering_Method == "DBSCAN"){
      if(!all(is.numeric(Min_cell_no), Min_cell_no%%1 == 0, Min_cell_no > 0)) stop("Min_cell_no must be an integer value > 0")
      if(!all(is.numeric(Distance_radius), Distance_radius > 0)) stop("Distance_radius must be a numeric value > 0")
    }

    if(!all(is.numeric(Hull_ratio), Hull_ratio >= 0, Hull_ratio <= 1)) stop("Hull_ratio must be a numeric value between 0 and 1")
    if(!is.logical(Calculate_border)) stop("Calculate_border must be a logical value")
    if(Calculate_border) if(!all(is.numeric(Dist_to_border), Dist_to_border > 0)) stop("Dist_to_border must be a numeric value > 0")
    if(is.null(Image_preview)) Image_preview <- sample(unique(DATA_Phenotypes$Subject_Names), size = 1)
    if(!Image_preview %in% unique(DATA_Phenotypes$Subject_Names)) stop(paste0(Image_preview, " not found in Subject_Names"))
    if(!all(N_cores >= 1 & N_cores%%1 == 0)) stop("N_cores must be an integer value > 0")


    print("Performing initial test")
    #First work on the tiling filtering method
    if(Filtering_Method == "Tiling"){
      #Select the image to be previewer
      For_plot <- DATA_Phenotypes %>% dplyr::filter(Subject_Names == Image_preview)

      #Generate tiled image
      Test_plot_source <- For_plot %>% dplyr::filter(Phenotype == Index_phenotype) %>%
        ggplot(aes(x = X, y = Y)) + geom_bin2d(binwidth = Accuracy)

      #Plot adequate tiles and overlay the tumor cells
      Tile_plot <- as_tibble(layer_data(Test_plot_source)) %>% dplyr::filter(value >= Min_cell_no) %>%
        ggplot() +
        geom_rect(aes(xmin = xmin, ymin = ymin, xmax = xmax, ymax = ymax), fill = "grey", color = "black") +
        theme_minimal() +
        scale_x_continuous("") +
        scale_y_continuous("") +
        geom_point(aes(x = X, y = Y), color = "red", alpha = 0.2, data = For_plot %>% dplyr::filter(Phenotype == Index_phenotype)) +
        theme(panel.grid = element_blank(),
              axis.text = element_blank())

      #Filter layer information according to pre-specified requirements
      Interim <- as_tibble(layer_data(Test_plot_source)) %>% dplyr::filter(value >= Min_cell_no) %>% dplyr::select(xbin, ybin, x, y, xmin, xmax, ymin, ymax)

      #Select the tumor cells
      Base <- For_plot %>% dplyr::filter(Phenotype == Index_phenotype)

      #Filter tumor cells according to their belonging or not to and adequate tile
      Final_tumor_cells <- Base %>% dplyr::mutate(Compartment = purrr::map2_lgl(.x = Base$X, .y = Base$Y, function(.x, .y) {
        x_position <- (.x >= Interim$xmin) & (.x < Interim$xmax)
        y_position <- (.y >= Interim$ymin) & (.y < Interim$ymax)
        Final_pos <- any(x_position & y_position)
      })) %>% dplyr::filter(Compartment)

      #Now we are going to build a polygon with these points
      #Prepare the require objects (sf object with tumor cells, sf object with all cells, polygon object and line object)
      Tumor_cells_sf <- sf::st_as_sf(Final_tumor_cells, coords = c("X", "Y"))
      All_cells_sf <- sf::st_as_sf(For_plot, coords = c("X", "Y"))
      Final_tumor_cells_polygon <- sf::st_cast((Tumor_cells_sf %>% summarise() %>% sf::st_concave_hull(ratio = Hull_ratio) %>% summarise), "POLYGON")
      Final_tumor_cells_line <- sf::st_cast((Tumor_cells_sf %>% summarise() %>% sf::st_concave_hull(ratio = Hull_ratio) %>% summarise), "LINESTRING")

      #Define compartments
      For_plot$Compartment <- "Stroma"
      For_plot[For_plot$Cell_no %in% as_tibble(sf::st_filter(All_cells_sf, Final_tumor_cells_polygon))[[1]],
               "Compartment"] <- "Tumor"
      #If border compartment needs to be calculated perform appropriate changes in the data
      if(Calculate_border){
        For_plot[unlist(sf::st_is_within_distance(All_cells_sf, Final_tumor_cells_line, sparse = F, dist = Dist_to_border)),
                 "Compartment"] <- "Border"
      }

      #Generate the final plot of the cell assignment
      Final_plot <- For_plot %>% ggplot(aes(x = X, y = Y, color = Compartment)) + geom_point() +
        scale_color_manual(values = c("red", "blue", "green")) +
        theme_minimal() +
        scale_x_continuous("") +
        scale_y_continuous("") +
        theme(panel.grid = element_blank(),
              axis.text = element_blank(),
              legend.position = "bottom")

      #Generate the polygon layout plot
      Layout_plot <- Final_tumor_cells_polygon %>% ggplot() + geom_sf(linewidth = 1.5, fill = "white", color = "black") +
        theme_minimal() +
        scale_x_continuous("") +
        scale_y_continuous("") +
        theme(panel.grid = element_blank(),
              axis.text = element_blank())

      #Plot all the results and then ask the user
      plot(patchwork::wrap_plots(Tile_plot, Final_plot, Layout_plot, nrow = 2))


      #Ask the user if the algorihtm should proceed
      answer <- menu(c("Proceed", "Abort"), title = "The plot depics the suggested area occupied by tumor cells. Should the algorithm proceed with the analysis?")
      #If user decides to stop then abort function and return stop message
      if(answer == 2) {
        stop("The function has been stopped. Please tune the parameters for a better result")
      }

      #If everything is OK proceed with computation
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
            #Select the image to be previewer
            For_plot <- DATA_Phenotypes %>% dplyr::filter(Subject_Names == Image)

            #Generate tiled image
            Test_plot_source <- For_plot %>% dplyr::filter(Phenotype == Index_phenotype) %>%
              ggplot(aes(x = X, y = Y)) + geom_bin2d(binwidth = Accuracy)

            #Filter layer information according to pre-specified requirements
            Interim <- as_tibble(layer_data(Test_plot_source)) %>% dplyr::filter(value >= Min_cell_no) %>% dplyr::select(xbin, ybin, x, y, xmin, xmax, ymin, ymax)

            #Select the tumor cells
            Base <- For_plot %>% dplyr::filter(Phenotype == Index_phenotype)

            #Filter tumor cells according to their belonging or not to and adequate tile
            Final_tumor_cells <- Base %>%dplyr::mutate(Compartment =purrr::map2_lgl(.x = Base$X, .y = Base$Y, function(.x, .y) {
              x_position <- (.x >= Interim$xmin) & (.x < Interim$xmax)
              y_position <- (.y >= Interim$ymin) & (.y < Interim$ymax)
              Final_pos <- any(x_position & y_position)
            })) %>% dplyr::filter(Compartment)

            #Now we are going to build a polygon with these points
            #Prepare the require objects (sf object with tumor cells, sf object with all cells, polygon object and line object)
            Tumor_cells_sf <- sf::st_as_sf(Final_tumor_cells, coords = c("X", "Y"))
            All_cells_sf <- sf::st_as_sf(For_plot, coords = c("X", "Y"))
            Final_tumor_cells_polygon <- sf::st_cast((Tumor_cells_sf %>% summarise() %>% sf::st_concave_hull(ratio = Hull_ratio) %>% summarise), "POLYGON")
            Final_tumor_cells_line <- sf::st_cast((Tumor_cells_sf %>% summarise() %>% sf::st_concave_hull(ratio = Hull_ratio) %>% summarise), "LINESTRING")

            #Define compartments
            For_plot$Compartment <- "Stroma"

            For_plot[For_plot$Cell_no %in% as_tibble(sf::st_filter(All_cells_sf, Final_tumor_cells_polygon))[[1]],
                     "Compartment"] <- "Tumor"

            #If border compartment needs to be calculated perform appropriate changes in the data
            if(Calculate_border){
              For_plot[unlist(sf::st_is_within_distance(All_cells_sf, Final_tumor_cells_line, sparse = F, dist = Dist_to_border)),
                       "Compartment"] <- "Border"
            }

            return(list(DATA = For_plot,
                        Area = sf::st_area(Final_tumor_cells_polygon)))


          },
          .progress = TRUE)
        future::plan("future::sequential")
        gc()

        #Name the list accorrdingly
        names(Final_list) <- unique(DATA_Phenotypes$Subject_Names)
        #Collapse data into a single DF
        Final_DATA <-purrr::map_dfr(Final_list, function(Image) Image[[1]])
        #Collapse areas into a single DF
        Final_Areas <- tibble(Subject_Names = unique(DATA_Phenotypes$Subject_Names),
                              Area =purrr::map_dbl(Final_list, function(Image) Image[[2]]))

        #Return all the results
        return(list(DATA_Phenotypes = Final_DATA,
                    DATA_Compartment_Area = Final_Areas))
      }
    }

    #If DBSCAN is the filtering method of choice, proceed accordingly
    if(Filtering_Method == "DBSCAN"){
      #Select Image to be previewed
      For_test <- DATA_Phenotypes %>% dplyr::filter(Subject_Names == Image_preview)
      #Obtain individual data sets containing Index cell and other cells
      For_test_Index <- For_test %>% dplyr::filter(Phenotype == Index_phenotype)
      For_test_no_Index <- For_test %>% dplyr::filter(Phenotype != Index_phenotype)

      #Calculate the DBSCAN model with the required parameters
      DB_results <- dbscan::dbscan(For_test_Index[c("X", "Y")], eps = Distance_radius, minPts = Min_cell_no, borderPoints = FALSE)

      #Obtain position in the tumor and the stroma for both index cells and other cells
      For_test_Index <- For_test_Index %>% dplyr::mutate(Compartment = as.character(DB_results$cluster)) %>%
        dplyr::mutate(Compartment = case_when(Compartment == "0" ~ "Stroma",
                                              T ~ "Tumor"))
      #Obtain the final tumor cells
      Final_tumor_cells <- For_test_Index %>% dplyr::filter(Compartment == "Tumor")

      #Now we are going to build a polygon with these points
      #Prepare the require objects (sf object with tumor cells, sf object with all cells, polygon object and line object)
      Tumor_cells_sf <- sf::st_as_sf(Final_tumor_cells, coords = c("X", "Y"))
      All_cells_sf <- sf::st_as_sf(For_test, coords = c("X", "Y"))
      Final_tumor_cells_polygon <- sf::st_cast((Tumor_cells_sf %>% summarise() %>% sf::st_concave_hull(ratio = Hull_ratio) %>% summarise), "POLYGON")
      Final_tumor_cells_line <- sf::st_cast((Tumor_cells_sf %>% summarise() %>% sf::st_concave_hull(ratio = Hull_ratio) %>% summarise), "LINESTRING")

      #Define compartments
      For_test$Compartment <- "Stroma"
      For_test[For_test$Cell_no %in% as_tibble(sf::st_filter(All_cells_sf, Final_tumor_cells_polygon))[[1]],
               "Compartment"] <- "Tumor"
      #If border compartment needs to be calculated perform appropriate changes in the data
      if(Calculate_border){
        For_test[unlist(sf::st_is_within_distance(All_cells_sf, Final_tumor_cells_line, sparse = F, dist = Dist_to_border)),
                 "Compartment"] <- "Border"
      }

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

      PLOT2 <- For_test_Index %>% dplyr::filter(Compartment != "Stroma") %>%
        ggplot(aes(x = X, y = Y)) +
        geom_point() +
        theme_minimal() +
        ggtitle("Tumor cells that are clustered") +
        theme(panel.grid = element_blank(),
              axis.text = element_blank(),
              axis.title = element_blank(),
              legend.position = "bottom",
              plot.title = element_text(hjust = 0.5))

      #Generate the final plot of the cell assignment
      Final_plot <- For_test %>% ggplot(aes(x = X, y = Y, color = Compartment)) + geom_point() +
        scale_color_manual(values = c("red", "blue", "green")) +
        theme_minimal() +
        scale_x_continuous("") +
        scale_y_continuous("") +
        theme(panel.grid = element_blank(),
              axis.text = element_blank(),
              legend.position = "bottom")

      #Generate the polygon layout plot
      Layout_plot <- Final_tumor_cells_polygon %>% ggplot() + geom_sf(linewidth = 1.5, fill = "white", color = "black") +
        theme_minimal() +
        scale_x_continuous("") +
        scale_y_continuous("") +
        theme(panel.grid = element_blank(),
              axis.text = element_blank())

      plot(patchwork::wrap_plots(PLOT1, PLOT2, Final_plot, Layout_plot, ncol = 2, nrow = 2, widths = 1))

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

            #Calculate the DBSCAN model with the required parameters
            DB_results <- dbscan::dbscan(For_test_Index[c("X", "Y")], eps = Distance_radius, minPts = Min_cell_no, borderPoints = FALSE)

            #Obtain position in the tumor and the stroma for both index cells and other cells
            For_test_Index <- For_test_Index %>%dplyr::mutate(Compartment = as.character(DB_results$cluster)) %>%
              dplyr::mutate(Compartment = case_when(Compartment == "0" ~ "Stroma",
                                                    T ~ "Tumor"))
            #Obtain the final tumor cells
            Final_tumor_cells <- For_test_Index %>% dplyr::filter(Compartment == "Tumor")

            #Now we are going to build a polygon with these points
            #Prepare the require objects (sf object with tumor cells, sf object with all cells, polygon object and line object)
            Tumor_cells_sf <- sf::st_as_sf(Final_tumor_cells, coords = c("X", "Y"))
            All_cells_sf <- sf::st_as_sf(For_test, coords = c("X", "Y"))
            Final_tumor_cells_polygon <- sf::st_cast((Tumor_cells_sf %>% summarise() %>% sf::st_concave_hull(ratio = Hull_ratio) %>% summarise), "POLYGON")
            Final_tumor_cells_line <- sf::st_cast((Tumor_cells_sf %>% summarise() %>% sf::st_concave_hull(ratio = Hull_ratio) %>% summarise), "LINESTRING")

            #Define compartments
            For_test$Compartment <- "Stroma"
            For_test[For_test$Cell_no %in% as_tibble(sf::st_filter(All_cells_sf, Final_tumor_cells_polygon))[[1]],
                     "Compartment"] <- "Tumor"
            #If border compartment needs to be calculated perform appropriate changes in the data
            if(Calculate_border){
              For_test[unlist(sf::st_is_within_distance(All_cells_sf, Final_tumor_cells_line, sparse = F, dist = Dist_to_border)),
                       "Compartment"] <- "Border"
            }

            #Return the final results
            return(list(DATA = For_test,
                        Area = sf::st_area(Final_tumor_cells_polygon)))
          },
          .progress = TRUE)
        future::plan("future::sequential")
        gc()

        #Name accordingly
        names(Final_list) <- unique(DATA_Phenotypes$Subject_Names)

        #Collapse data into a single DF
        Final_DATA <-purrr::map_dfr(Final_list, function(Image) Image[[1]])

        #Collaps areas into a single DF
        Final_Areas <- tibble(Subject_Names = unique(DATA_Phenotypes$Subject_Names),
                              Area =purrr::map_dbl(Final_list, function(Image) Image[[2]]))

        #Return the final results
        return(list(DATA_Phenotypes = Final_DATA,
                    DATA_Compartment_Area = Final_Areas))

      }
    }
  }
