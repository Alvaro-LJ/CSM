#' Executes cell phenotype label prediction
#'
#' Given a cell feature dataset and a cell phenotype label prediction model, the function assigns cell phenotype labels to all cells in the dataset.
#' The model must have been fitted using [Image_based_phenotyper_App_launcher()].
#'
#' @param DATA A dataframe or tibble containing cell feature data.
#' @param Model_parameters The model obtained using [Image_based_phenotyper_App_launcher()].
#' @param N_cores Integer. Number of cores to parallelize your computation.
#'
#' @seealso [Image_based_phenotyper_App_launcher()]
#'
#' @export

Model_cell_phenotyper <-
  function(DATA = NULL,
           Model_parameters = NULL,
           N_cores = NULL){

    #Check suggested packages
    {
      if(!requireNamespace("rtree", quietly = FALSE)) stop(
        paste0("rtree GitHub package is required to execute the function. Please install using the following code: ",
               expression(remotes::install_github("akoyabio/rtree")))
      )
      if(!requireNamespace("magick", quietly = FALSE)) stop(
        paste0("magick CRAN package is required to execute the function. Please install using the following code: ",
               expression(install.packages("magick")))
      )
      if(!requireNamespace("tidymodels", quietly = FALSE)) stop(
        paste0("tidymodels CRAN package is required to execute the function. Please install using the following code: ",
               expression(install.packages("tidymodels")))
      )
      if(!requireNamespace("randomForest", quietly = FALSE)) stop(
        paste0("randomForest CRAN package is required to execute the function. Please install using the following code: ",
               expression(install.packages("randomForest")))
      )
      if(!requireNamespace("xgboost", quietly = FALSE)) stop(
        paste0("xgboost CRAN package is required to execute the function. Please install using the following code: ",
               expression(install.packages("xgboost")))
      )
      if(!requireNamespace("brulee", quietly = FALSE)) stop(
        paste0("brulee CRAN package is required to execute the function. Please install using the following code: ",
               expression(install.packages("brulee")))
      )
    }

    #On exit always return to sequential computing
    on.exit({
      future::plan("future::sequential")
      gc()
    })

    print("Argument check...")
    #Get the data
    DATA <- DATA
    if(!identical(names(DATA)[1:4],  c("Cell_no", "X", "Y", "Subject_Names"))) stop("DATA provided should have an adecuate format")

    #Check the cores
    if(!all(N_cores >= 1 & N_cores%%1 == 0)) stop("N_cores must be an integer value > 0")

    #If the DATA does not contain the required features return an error
    Model_parameters_list <- Model_parameters
    if(!all(Model_parameters_list$Model_parameters$Model_features %in% names(DATA))){
      Absent_features <- Model_parameters_list$Model_parameters$Model_features[!Model_parameters_list$Model_parameters$Model_features %in% names(DATA)]
      stop(paste0("The following features are not present in data: ", stringr::str_c(Absent_features, collapse = ", ")))
    }


    #Proceed if features are OK (no need to check arguments because they have been pre-processed by the app)
    #If NO SPATIAL CONTEXT REQUIRED then proceed
    if(!Model_parameters_list$Model_parameters$Spatial_context){
      print("Generating predictions")
      #Generate a list of images with the required features
      Features_DATA_list <-purrr::map(unique(DATA$Subject_Names), function(Image){
        DATA %>% dplyr::filter(Subject_Names == Image) %>% dplyr::select(all_of(Model_parameters_list$Model_parameters$Model_features))
      })
      names(Features_DATA_list) <- unique(DATA$Subject_Names)

      #Get the Model
      Model <- Model_parameters_list[["Model"]]

      #Predict for every Image
      future::plan("future::multisession", workers = N_cores)
      options(future.globals.maxSize = Inf, future.rng.onMisuse = "ignore")
      furrr::furrr_options(scheduling = Inf)
      PREDICTIONS <-
        furrr::future_map_dfr(Features_DATA_list, function(Image_features){
          #Load tidymodels to solve method dispatch issues
          suppressMessages(
            library(tidymodels)
          )
          #Generate the predictions
          Predictions <- predict(Model, new_data = Image_features, type = "prob")
          #Obtain the label with the highest probability
          Col_index <- max.col(Predictions, ties.method = "random")
          Predictions <- tibble(Label = colnames(Predictions)[Col_index],
                                Probability =purrr::map2_dbl(.x = 1:nrow(Predictions), .y = Col_index, function(.x, .y) Predictions[[.x, .y]]))
          Predictions$Label <- substr(Predictions$Label, start = 7, stop = nchar(Predictions$Label))
          names(Predictions)[1] <- "Phenotype"

          #If the label with the highest probability does not reach the threshold change to unnasigned
          Predictions$Phenotype[Predictions$Probability < Model_parameters_list$Model_parameters$Model_threshold] <- "Unassigned"
          #Return only the predicted phenotype
          return(Predictions %>% dplyr::select(Phenotype))
        },
        .progress = TRUE)
      future::plan("future::sequential")
      gc()
    }

    #If SPATIAL CONTEXT REQUIRED FIRST CALCULATE THE NEIGHBORHOOD DATA
    if(Model_parameters_list$Model_parameters$Spatial_context){

      print("Obtaining cell neighbors related features")
      DATA_Neighbors <- DATA %>% dplyr::select(1:4, all_of(Model_parameters_list$Model_parameters$Model_features))
      DATA_Neighbors <- UTAG_message_passing_Image_based_phenotyper(DATA = DATA_Neighbors,
                                                                    COO_to_visit = NULL,
                                                                    Neighbor_strategy = Model_parameters_list$Model_parameters$Neighbor_strategy,
                                                                    Message_strategy = Model_parameters_list$Model_parameters$Message_strategy,
                                                                    N_neighbors = Model_parameters_list$Model_parameters$N_neighbors,
                                                                    Max_dist_allowed = Model_parameters_list$Model_parameters$Max_dist_allowed,
                                                                    Weighting_Strategy = Model_parameters_list$Model_parameters$Weighting_Strategy,
                                                                    N_cores = N_cores)

      DATA_Neighbors <- DATA_Neighbors %>% dplyr::select(-X, -Y, -Subject_Names, -mean_DIST, -max_DIST, -N_neighbors)
      names(DATA_Neighbors)[-1] <-stringr::str_c("Neighbor_", names(DATA_Neighbors)[-1], sep = "")

      Features_DATA_list <- DATA %>% dplyr::select(Cell_no, Subject_Names, all_of(Model_parameters_list$Model_parameters$Model_features))

      Features_DATA_list <-dplyr::left_join(Features_DATA_list,  DATA_Neighbors, by = "Cell_no")


      print("Generating predictions")
      #Generate a list of images with the required features
      Features_DATA_list <-purrr::map(unique(DATA$Subject_Names), function(Image){
        Features_DATA_list %>% dplyr::filter(Subject_Names == Image) %>% dplyr::select(-Cell_no, -Subject_Names)
      })
      names(Features_DATA_list) <- unique(DATA$Subject_Names)

      #Get the Model
      Model <- Model_parameters_list[["Model"]]

      #Predict for every Image
      future::plan("future::multisession", workers = N_cores)
      options(future.globals.maxSize = Inf, future.rng.onMisuse = "ignore")
      furrr::furrr_options(scheduling = Inf)
      PREDICTIONS <-
        furrr::future_map_dfr(Features_DATA_list, function(Image_features){
          #Load tidymodels to solve method dispatch issues
          suppressMessages(
            library(tidymodels)
          )
          #Generate the predictions
          Predictions <- predict(Model, new_data = Image_features, type = "prob")
          #Obtain the label with the highest probability
          Col_index <- max.col(Predictions, ties.method = "random")
          Predictions <- tibble(Label = colnames(Predictions)[Col_index],
                                Probability =purrr::map2_dbl(.x = 1:nrow(Predictions), .y = Col_index, function(.x, .y) Predictions[[.x, .y]]))
          Predictions$Label <- substr(Predictions$Label, start = 7, stop = nchar(Predictions$Label))
          names(Predictions)[1] <- "Phenotype"

          #If the label with the highest probability does not reach the threshold change to unnasigned
          Predictions$Phenotype[Predictions$Probability < Model_parameters_list$Model_parameters$Model_threshold] <- "Unassigned"
          #Return only the predicted phenotype
          return(Predictions %>% dplyr::select(Phenotype))
        },
        .progress = TRUE)
      future::plan("future::sequential")
      gc()
    }

    #Bind the predicted phenotype with the original data
    FINAL_DATA <- dplyr::bind_cols(DATA, PREDICTIONS)
    return(FINAL_DATA)
  }
