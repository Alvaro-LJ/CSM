#' Performs image tiling according to randomly iterating parameters
#'
#'
#' @param N_cores Integer. Number of cores to parallelize your computation.
#' @param DATA A dataframe or tibble containing cell feature data.
#' @param Tile_width A numeric value indicating the reference width of the tiles.
#' @param Tile_height A numeric value indicating the reference height of the tiles.
#' @param Variables_to_keep A character vector indicating the names of the columns to be kept when tiling the image.
#'
#' @param N_iterations A integer value indicating the number of random parameter combinations that will generate tiled images.
#' @param Max_tile_size_variation A numeric value >= 1, indicating the max ratio of tile size that will be used in the random parameters.
#' @param Force_squared_tiles A logical value indicating if tiles width and height should be forced to be equal in every random parameter iteration.
#'
#' @details
#' The function is aimed to estimate impact of the Modifiable areal unit problem (MAUP) (please see the 'Calculating tissue heterogeneity' Vignette for more details).
#' The function will generate tiled information for every image randomly iterating certain image parameters (tile dimension, cell coordinate angle rotating and position flipping).
#' These tiled images can be used to compute the influence of the MAU on spatial heterogeneity analysis (using the [Random_tiled_image_heterogeneity_calculator()] and
#' [Random_texture_features_calculator()] functions).
#'
#' The first two random iterations will always contain tiling strategies in the extremes of the user defined tile size variation factor.
#'
#' @returns A list containing the parameters used to generate the tiled images and the actual results.
#'
#' @seealso [Random_tiled_image_heterogeneity_calculator()], [Random_tiled_image_heterogeneity_analyzer()], [Random_texture_features_calculator()]
#'
#' @examples
#' \dontrun{
#' #Generate tiled image data using 10 randomly generated parameter set---------
#' Random_tile_DATA <-
#'   Random_image_tiling_processing_function(
#'     N_cores =  1,
#'     DATA = CSM::CSM_Phenotypecell_test,
#'     Tile_width = 100,
#'     Tile_height = 100,
#'     Variables_to_keep = "Phenotype",
#'
#'     N_iterations = 10,
#'     Max_tile_size_variation = 1.1,
#'     Force_squared_tiles = TRUE
#'   )
#'
#' #Compute the heterogeneity by tile using this random tile dataset------------
#' Tiled_heterogeneity_DATA <-
#'   Random_tiled_image_heterogeneity_calculator(
#'     N_cores = 1,
#'     Random_tiled_images = Random_tile_DATA,
#'     Minimum_cell_no_per_tile = 1,
#'     Phenotypes_included = unique(CSM_Phenotypecell_test$Phenotype)
#'   )
#'
#' #Analyze the data-----------------------------------------------------------
#' Random_tiled_image_heterogeneity_analyzer(
#'     Random_tiled_heterogeneity_DATA = Tiled_heterogeneity_DATA,
#'     Strategy = "Overall_Summary",
#'     Metric = "Shannon"
#'   )
#'}
#'
#'
#' @export

Random_image_tiling_processing_function <-
  function(N_cores = 1,
           DATA,
           Tile_width,
           Tile_height,
           Variables_to_keep,

           N_iterations = 10,
           Max_tile_size_variation = 1.5,
           Force_squared_tiles = FALSE){

    ####ARGUMENT CHECK####
    if(!(is.numeric(Tile_width) & is.numeric(Tile_height))) stop("Both Tile_width and Tile_height must be numeric values")
    if(!(Tile_width%%1 == 0 & Tile_height%%1 ==0)) message("It is highly recommended that Tile_width and Tile_height are integer values")
    if(!(Tile_width > 0 & Tile_height > 0)) stop("Tile_width and Tile_height must be > 0")
    if(!all(N_cores%%1 == 0, N_cores > 0, length(N_cores) == 1)) stop("N_cores must be an integer value > 0")
    if(!all(Variables_to_keep %in% names(DATA))) stop(paste0(stringr::str_c(Variables_to_keep[!Variables_to_keep %in% names(DATA)], collapse = ", "),
                                                             " not found in DATA variables"))
    if(!all(c("X", "Y", "Subject_Names", Variables_to_keep) %in% names(DATA))) stop("DATA supplied must be follow the structure specified in step 0 and contain Variables_to_keep")

    if(!all(length(N_iterations) == 1, N_iterations%%1 == 0, N_iterations >= 3)) stop("N_iterations must be a single integer value > 2")
    if(!all(length(Max_tile_size_variation) == 1, is.numeric(Max_tile_size_variation), Max_tile_size_variation >= 1)) stop("Max_tile_size_variation must be a single numeric value >= 1")
    if(!all(length(Force_squared_tiles) == 1, is.logical(Force_squared_tiles))) stop("Force_squared_tiles must be a single logical value")

    #if Squared tiles need to be forced then check that Tiled_width and _height are equal
    if(Force_squared_tiles){
      if(!identical(Tile_width, Tile_height)) stop("If Force_squared_tiles is TRUE, Tiled_width and Tile_height must be identical")
    }

    #Filter the data to save memory use
    DATA <- DATA %>% dplyr::select(Cell_no, X, Y, Subject_Names, dplyr::all_of(Variables_to_keep))

    ####Generate the random variation look-up table####
    #The first two rows of the look-up table will always contain the extremes of tile size (to guarantee extreme size representation)

    Random_iteration_lookup_table <- tibble(Iteration = stringr::str_c("Iteration", 1:N_iterations, sep = "_"))

    #If no force squared tiles then proceed as usual
    if(!Force_squared_tiles){
      #Tile width
      Random_iteration_lookup_table$Tile_width <-
        c(Tile_width * Max_tile_size_variation, #First element is always the max
          Tile_width / Max_tile_size_variation, #Second element is always the min

          purrr::map_dbl(seq_along(1:(N_iterations-2)), function(Iter){
            Multi_or_divide <- sample(c("Multiply", "Divide"), size = 1) #Randomly select if multiply or divide

            if(Multi_or_divide == "Multiply") return(Tile_width * runif(n = 1, min = 1, max = Max_tile_size_variation)) #Return multiplication between tile width and a random number between 1 and max size
            if(Multi_or_divide == "Divide") return(Tile_width / runif(n = 1, min = 1, max = Max_tile_size_variation)) #Return division between tile width and a random number between 1 and max size

          })
      )
      Random_iteration_lookup_table$Tile_width <- ceiling(Random_iteration_lookup_table$Tile_width)

      #Tile height
      Random_iteration_lookup_table$Tile_height <-
        c(Tile_height * Max_tile_size_variation, #First element is always the max
          Tile_height / Max_tile_size_variation, #Second element is always the min

          purrr::map_dbl(seq_along(1:(N_iterations-2)), function(Iter){
            Multi_or_divide <- sample(c("Multiply", "Divide"), size = 1) #Randomly select if multiply or divide

            if(Multi_or_divide == "Multiply") return(Tile_height * runif(n = 1, min = 1, max = Max_tile_size_variation)) #Return multiplication between tile width and a random number between 1 and max size
            if(Multi_or_divide == "Divide") return(Tile_height / runif(n = 1, min = 1, max = Max_tile_size_variation)) #Return division between tile width and a random number between 1 and max size

          })
        )
      Random_iteration_lookup_table$Tile_height <- ceiling(Random_iteration_lookup_table$Tile_height)
    }
    if(Force_squared_tiles){
      #Tile width
      Random_iteration_lookup_table$Tile_width <-
        c(Tile_width * Max_tile_size_variation, #First element is always the max
          Tile_width / Max_tile_size_variation, #Second element is always the min

          purrr::map_dbl(seq_along(1:(N_iterations-2)), function(Iter){
            Multi_or_divide <- sample(c("Multiply", "Divide"), size = 1) #Randomly select if multiply or divide

            if(Multi_or_divide == "Multiply") return(Tile_width * runif(n = 1, min = 1, max = Max_tile_size_variation)) #Return multiplication between tile width and a random number between 1 and max size
            if(Multi_or_divide == "Divide") return(Tile_width / runif(n = 1, min = 1, max = Max_tile_size_variation)) #Return division between tile width and a random number between 1 and max size

          })
        )
      Random_iteration_lookup_table$Tile_width <- ceiling(Random_iteration_lookup_table$Tile_width)

      #Tile height is the same as width
      Random_iteration_lookup_table$Tile_height <- Random_iteration_lookup_table$Tile_width
    }

    Random_iteration_lookup_table$X_flip <- sample(c(TRUE, FALSE), size = N_iterations, replace = TRUE)
    Random_iteration_lookup_table$Y_flip <- sample(c(TRUE, FALSE), size = N_iterations, replace = TRUE)
    Random_iteration_lookup_table$Rotation <- floor(runif(n = N_iterations, min = 0, max = 360))

    #Print the result
    print(Random_iteration_lookup_table, n = 100)

    #Ask the user if everything is OK
    Test_OK <- menu(choices = c("Proceed", "Abort"), title = "Check the random iteration plan generated. Should the analysis proceed?")
    if(Test_OK == 2) stop("Analysis aborted")
    if(N_iterations > 20) cat("\n More than 20 iterations will be generated, this can be computationally demanding...")


    ####Compute the random tiling iterations for every sample####
    Result_list <-
      purrr::map(seq_along(1:length(unique(DATA$Subject_Names))), function(Image_index){

        #Obtain the image name
        Image_name <- unique(DATA$Subject_Names)[Image_index]

        #Print a message
        cat(paste0("\n", "Generating tiles for image ", Image_index, "/", length(unique(DATA$Subject_Names)), " - ", Image_name))

        #Proceed with analysis
        Interim <- DATA %>% dplyr::filter(Subject_Names == Image_name)

        #Iterate over random iterations
        #Set the on exit statement
        on.exit({
          future::plan("future::sequential")
          gc()
        })

        #Set the cores
        future::plan("future::multisession", workers = N_cores)
        options(future.globals.maxSize = Inf, future.rng.onMisuse = "ignore")
        furrr::furrr_options(scheduling = Inf)

        #Compute the results list
        Iteration_results <-
          furrr::future_map(seq_along(1:N_iterations), function(Iteration_index){
          #Get the parameters
          Iteration_N <- Random_iteration_lookup_table$Iteration[Iteration_index]
          Tile_width_current_iter <- Random_iteration_lookup_table$Tile_width[Iteration_index]
          Tile_height_current_iter <- Random_iteration_lookup_table$Tile_height[Iteration_index]
          Perform_X_flip <- Random_iteration_lookup_table$X_flip[Iteration_index]
          Perform_Y_flip <- Random_iteration_lookup_table$Y_flip[Iteration_index]
          Rotation_angle <- Random_iteration_lookup_table$Rotation[Iteration_index]

          #Prepare the modified data
          Modified_data <- Interim

          #If required perform X flip
          if(Perform_X_flip){
            Modified_data$X <- -Modified_data$X
            Modified_data$X <- Modified_data$X + (-min(Modified_data$X)) + min(Interim$X)
          }
          #If required perform Y flip
          if(Perform_X_flip){
            Modified_data <- Interim
            Modified_data$Y <- -Modified_data$Y
            Modified_data$Y <- Modified_data$Y + (-min(Modified_data$Y)) + min(Interim$Y)
          }

          #Perform angle rotation
          New_x <- (Modified_data$X * round(cos(Rotation_angle * pi / 180), digits = 4)) - (Modified_data$Y * round(sin(Rotation_angle * pi / 180), digits = 4))
          New_y <- (Modified_data$X * round(sin(Rotation_angle * pi / 180), digits = 4)) + (Modified_data$Y * round(cos(Rotation_angle * pi / 180), digits = 4))

          Modified_data$X <- New_x + (-min(New_x)) + min(Modified_data$X)
          Modified_data$Y <- New_y + (-min(New_y)) + min(Modified_data$Y)

          #Compute the tiling results
          Tiles <- Image_tiling_processing_function(
            N_cores = 1,
            DATA = Modified_data,
            Tile_width = Tile_width_current_iter,
            Tile_height = Tile_height_current_iter,
            Variables_to_keep = Variables_to_keep
          )

          #Get the results and change names
          Tile_info <- Tiles[[1]][[1]]
          names(Tile_info)[3] <- stringr::str_c("tile_id", "iteration", Iteration_index, sep = "_")

          Cell_info <- Tiles[[1]][[2]] %>% dplyr::select(Cell_no, Subject_Names, tile_id, dplyr::all_of(Variables_to_keep))
          names(Cell_info)[3] <- stringr::str_c("tile_id", "iteration", Iteration_index, sep = "_")

          return(list(Tile_info = Tile_info,
                      Cell_info = Cell_info))

          },
          .progress = TRUE)

        #Stop the multicore
        future::plan("future::sequential")
        gc()

        #Add iteration_name
        names(Iteration_results) <- stringr::str_c("Iteration", 1:N_iterations, sep = "_")

        return(Iteration_results)
      })
    #Add Subject_Names
    names(Result_list) <- unique(DATA$Subject_Names)

    ####Build the final list in an adequate format####
    #Build the final list
    Result_list_formatted <-
      purrr::map(Result_list, function(Image){
        #Get tile information (a list)
        Iteration_tiles <- purrr::map(Image, ~.[[1]])

        #Generate a single tibble with ALL the tile IDs
        Cell_info_tibble <- Image[[1]][[2]][c(1,2,4)]
        Tile_IDs <- purrr::map_dfc(Image, ~.[[2]][3])
        Cell_info_tibble <- dplyr::bind_cols(Cell_info_tibble, Tile_IDs)

        Final_list <- list(Iteration_tile_list = Iteration_tiles,
                           Cell_information = Cell_info_tibble)
        return(Final_list)
      })

    #Return the final result
    return(list(Iteration_Parameters = Random_iteration_lookup_table,
                Results = Result_list_formatted)
    )
  }


