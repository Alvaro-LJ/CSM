#' Calculates concordance between different methods.
#'
#' The function calculates concordance between two or more methods that generate cell labels.
#'
#' @param ... Two or more dataframes or tibbles containing a column with cell labels. A name must be provided for every tibble.
#' @param Variable A character value of the column name containing the labels. The name must be shared across datasets provided.
#' @param Strategy One of the following: 'Rand' or 'FM' (Fowlkes–Mallows) (see details).
#'
#' @details
#' Rand index is calculated using catsim::rand_index function.
#'
#' Fowlkes-Mallows index is calculated using dendextend::FM_index function.
#'
#' @returns Returns a tibble containing concordance values by sample. Generates summary plots.
#'
#' @export

Concordance_calculator <-
  function(...,
           Variable = NULL,
           Strategy = NULL){
    #Check suggested packages
    {
      if(Strategy == "Rand"){
        if(!requireNamespace("catsim", quietly = FALSE)) stop(
          paste0("catsim CRAN package is required to execute the function. Please install using the following code: ",
                 expression(install.packages("catsim")))
        )
      }
      if(Strategy == "FM"){
        if(!requireNamespace("dendextend", quietly = FALSE)) stop(
          paste0("dendextend CRAN package is required to execute the function. Please install using the following code: ",
                 expression(install.packages("dendextend")))
        )
      }

    }

    on.exit(gc())
    Tibble_list <- list(...)

    print("Argument check...")
    #Check the variable argument
    if(!is.character(Variable)) stop("Variable must be a character value")
    if(!Strategy %in% c("Rand", "FM")) stop("Strategy must be one of the following: Rand, FM")

    #Check the user has provided names to the DATA provided
    if(any(names(Tibble_list) == "")) stop("User must provide a name for each DATA provided")

    #Check that names not are repeated
    if(length(unique(names(Tibble_list))) != length(Tibble_list)) stop("DATA provided must have each a different name")

    #Check that names do not contain the  '_' character. If they do modify with '.'
    names(Tibble_list) <- str_replace_all(names(Tibble_list), pattern = "_", replacement = ".")

    #Check the length of the list
    if(length(Tibble_list) < 2) stop("At least 2 DATA sources must be provided to calculate concordance")

    #Check that the first four columns in the data are the sampe
    if(!all(purrr::map_lgl(Tibble_list, function(tibble){
      identical(names(tibble)[c(1:4)], c("Cell_no", "X", "Y", "Subject_Names"))
    }))
    ) stop("DATA provided must have been formatted appropiately. Names of DATA provided do not match required standards")

    #Check that data have the same number of rows
    if(length(unique(purrr::map_dbl(Tibble_list, function(tibble) nrow(tibble)))) != 1) stop("All DATA provided must have the same number of rows")

    #Check that Variable is present in all data
    Variable_in_tibbles <-purrr::map_lgl(Tibble_list, function(tibble) as.character(Variable) %in% names(tibble))
    if(!all(Variable_in_tibbles)){
      stop(paste0("The following DATA provided does not contain the ", Variable, " column: ",
                  stringr::str_c(names(Tibble_list)[!Variable_in_tibbles], collapse = ", ")))
    }

    #check the presence of NA values in the Variable column
    Tibbles_with_NA <-purrr::map_lgl(Tibble_list, function(tibble) sum(is.na(tibble[[Variable]])) >= 1)
    if(any(Tibbles_with_NA)){
      stop(paste0(stringr::str_c(names(Tibble_list)[Tibbles_with_NA], collapse = ", "), " contains NA in ", Variable, " column. Please remove NA values"))
    }

    #Check that the same subject names are the same in all the datasets
    Subject_names_tibble <-purrr::map_df(Tibble_list, function(tibble){
      tibble[["Subject_Names"]]
    })
    Subject_names_tibble <- Subject_names_tibble %>% distinct
    if(
      !all(purrr::map_lgl(seq_along(1:nrow(Subject_names_tibble)), function(Row){
        Row <- unlist(Subject_names_tibble [Row,])
        length(unique(Row)) == 1
      }))
    ) stop("Subject_Names must be the same in all DATA provided")

    #Check that all Cells_ID are matched
    Cell_no_tibble <-purrr::map_df(Tibble_list, function(tibble){
      tibble[["Cell_no"]]
    })
    Cell_no_tibble <- Cell_no_tibble %>% distinct
    if(
      !all(purrr::map_lgl(seq_along(1:nrow(Cell_no_tibble)), function(Row){
        Row <- unlist(Cell_no_tibble[Row,])
        length(unique(Row)) == 1
      }))
    ) stop("Rows in each DATA provided must have the information of the same cell")
    #Remove the cell_no_tibble
    rm(Cell_no_tibble)

    #select the variable and turn it into a numeric value
    Tibble_list <-purrr::map(Tibble_list, function(tibble){
      Interim_tibble <- tibble %>% dplyr::select(1:4, all_of(Variable))
      Interim_tibble[5] <- as.numeric(as.factor(Interim_tibble[[5]]))
      Interim_tibble
    })

    #Start performing the computation
    By_subject_results <-purrr::map(unique(Tibble_list[[1]][["Subject_Names"]]),
                                    function(Subject){
                                      #Generate a tibble that has a column for every method
                                      purrr::map_dfc(seq_along(1:length(Tibble_list)), function(Method){
                                        Interim <- Tibble_list[[Method]] %>% dplyr::filter(Subject_Names == Subject) %>% dplyr::select(all_of(Variable))
                                        names(Interim) <- names(Tibble_list)[Method]
                                        return(Interim)
                                      })
                                    }, .progress = list(clear = F,
                                                        name = paste0("Obtaining ", Variable, " information for each Subject Name"),
                                                        show_after = 1,
                                                        type = "iterator"))

    names(By_subject_results) <- unique(Tibble_list[[1]][["Subject_Names"]])
    By_subject_results

    #Calculate the rand index for every 2 by 2 method
    Concordance_index_result <-purrr::map(By_subject_results, function(Subject){
      #generate a comparison plan (it will indicate thepurrr::map function what 2 by 2 comparisons need to be performed)
      Comparison_plan <- expand.grid.unique(names(Subject), names(Subject))
      Comparison_plan <- tibble(Method_A = Comparison_plan[,1], Method_B = Comparison_plan[,2])

      #Execute the comparison plan according to the strategy selected
      if(Strategy == "Rand"){
        Result <- Comparison_plan %>%dplyr::mutate(Concordance_result = pmap_dbl(Comparison_plan, function(Method_A, Method_B){
          catsim::rand_index(Subject[[Method_A]], Subject[[Method_B]])
        }))
      }
      if(Strategy == "FM"){
        Result <- Comparison_plan %>%dplyr::mutate(Concordance_result = pmap_dbl(Comparison_plan, function(Method_A, Method_B){
          dendextend::FM_index(Subject[[Method_A]], Subject[[Method_B]])[[1]]
        }))
      }
      return(Result)
    }, .progress = list(clear = F,
                        name = paste0("Calculating ", as.character(Strategy), " index for each Subject"),
                        show_after = 1,
                        type = "iterator"))

    #Generate a tibble with the results by patient
    By_Patient_results <-purrr::map_dfr(seq_along(1:length(Concordance_index_result)),
                                        function(Index){
                                          Interim_tibble <- Concordance_index_result[[Index]] %>%dplyr::mutate(Comparison = stringr::str_c(Method_A, Method_B, sep = "_")) %>%
                                            dplyr::select(Concordance_result, Comparison) %>%
                                            pivot_wider(names_from = Comparison, values_from = Concordance_result) %>%
                                            dplyr::mutate(Subject_Names = names(Concordance_index_result)[Index])
                                          Interim_tibble[c(ncol(Interim_tibble), 1:(ncol(Interim_tibble)-1))]
                                        })
    #Calculate the median and p25 and q75 for each comparsion and arrange it all in a single list
    q50_rand <-purrr::map_dbl(By_Patient_results[-1], function(Comparison) quantile(Comparison, 0.50))
    q25_rand <-purrr::map_dbl(By_Patient_results[-1], function(Comparison) quantile(Comparison, 0.25))
    q75_rand <-purrr::map_dbl(By_Patient_results[-1], function(Comparison) quantile(Comparison, 0.75))
    Result_list <- list(Quantile_25 = q25_rand,
                        Quantile_50 = q50_rand,
                        Quantile_75 = q75_rand)

    #Use this list to generate heatmaps of comparsion results
    Heatmap_plots <-  purrr::map(seq_along(1:length(Result_list)), function(Index){
      #Get the result list and change the names of the tibble, then add a column specifying the comparison
      Interim <- as_tibble(Result_list[[Index]])
      names(Interim) <- "Value"
      Interim$Comparison <- names(By_Patient_results)[-1]

      #Split the comparison back into the two methods being compared and bind it to the interim tibble
      Comparisons <- str_split(Interim$Comparison, pattern = "_", simplify = TRUE)
      colnames(Comparisons) <- c("Method_A", "Method_B")
      Interim <-dplyr::bind_cols(Interim, as_tibble(Comparisons))

      #Start generating the tibble for the final ggplot (it will require to columns with all possible comparisons, including duplicated ones)
      All_comparison_tibble <- as_tibble(expand.grid(unique(c(Interim$Method_A, Interim$Method_B)), unique(c(Interim$Method_A, Interim$Method_B))))
      names(All_comparison_tibble) <- c("Method_A", "Method_B")
      All_comparison_tibble$Comparison <- stringr::str_c(All_comparison_tibble$Method_A, All_comparison_tibble$Method_B, sep = "_")
      All_comparison_tibble <- All_comparison_tibble

      #Generate another tibble with the same information but reverse labels
      Interim2 <- tibble(Value = Interim$Value,
                         Comparison = stringr::str_c(Interim$Method_B, Interim$Method_A, sep = "_"),
                         Method_A = Interim$Method_A,
                         Method_B = Interim$Method_B)


      Interim <-dplyr::bind_rows(Interim, Interim2) %>% dplyr::select(-Method_A, -Method_B)
      Tibble_plot <-dplyr::left_join(All_comparison_tibble, Interim, by = "Comparison")

      Tibble_plot %>%
        ggplot(aes(x = Method_A, y = Method_B, fill = Value)) +
        geom_tile() +
        scale_x_discrete("") +
        scale_y_discrete("", limits = rev) +
        scale_fill_gradient2(limits = c(0,1), low = "#f23a3a", mid = "white", high = "#3af256", midpoint = 0.5) +
        cowplot::theme_cowplot() +
        geom_text(aes(label = as.character(round(Value, 2))), size = 7) +
        guides(fill = "none") +
        ggtitle(names(Result_list)[Index]) +
        theme(plot.title = element_text(size = 12, hjust = 0.5),
              axis.text.x = element_text(size = 9, face = "bold", color = "black"),
              axis.text.y = element_text(size = 9, face = "bold", color = "black"))
    })
    #Print the plots
    suppressWarnings(plot(cowplot::plot_grid(plotlist = Heatmap_plots, nrow = 1, ncol = 3)))

    #Generate a heatmap with the results by patient
    plot(By_Patient_results %>% pivot_longer(-1) %>%
           ggplot(aes(x = Subject_Names, y = name, fill = value)) +
           geom_tile() +
           scale_x_discrete("") +
           scale_y_discrete("", limits = rev) +
           scale_fill_gradient2(as.character(Strategy), limits = c(0,1), low = "#f23a3a", mid = "white", high = "#3af256", midpoint = 0.5) +
           cowplot::theme_cowplot() +
           theme(axis.text.x = element_text(size = 8, color = "black", angle = -85, hjust = 0, vjust = 0.5),
                 axis.text.y = element_text(size = 10, face = "bold", color = "black")))

    #Return the Rand by patient
    return(By_Patient_results)
  }
