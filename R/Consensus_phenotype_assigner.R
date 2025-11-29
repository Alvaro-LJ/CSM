#' Assigns phenotype labels by reaching a consensus between various methods
#'
#' The function allows various cell phenotyping methods to be combined into a final consensus label.
#'
#' @param ... Two or more dataframes or tibbles containing a column with cell labels. A name must be provided for every tibble.
#' @param Win_strategy One of the following: 'Majority', 'Arbitrary', 'Difference' (see details).
#' @param Arbitrary_threshold If Win_strategy is Arbitrary or Difference, a numeric threshold between 0 and 1.
#' @param Weights (OPTIONAL) A numeric vector of length equal to the number of datasets provided.
#' @param No_consensus_value A character value indicating which label to use if no consensus is reached (default is NA_character_).
#' @param Tie_break_method If Win_strategy is Majority, the method to break potential ties. Either 'Random' or 'No_consensus' (see details).
#' @param N_cores Integer. Number of cores to parallelize your computation.
#'
#' @details
#' Majority voting assigns the label with the highest votes. If there is a tie between two or more labels, the tie can be resolved by randomly selecting a label, or by using the No_consensus_value.
#'
#' Arbitrary voting assigns the label with the highest votes as long as it has received enough votes. Arbitrary threshold must be between 0.5 and 1.
#'
#' Difference voting assigns the label with the highest votes as long as the difference between the selected label and the second most voted option is above an arbitrary threshold.
#'
#' @returns Returns a tibble containing consensus cell labels.
#'
#' @examples
#' \dontrun{
#' #Generate datasets with random cell phenotype labels----------------------
#' CSM_Phenotypecell_test_1 <-
#'  CSM_Phenotypecell_test %>%
#'      mutate(Phenotype = sample(unique(CSM_Phenotypecell_test$Phenotype),
#'                                size = nrow(CSM_Phenotypecell_test),
#'                                replace = TRUE))
#' CSM_Phenotypecell_test_2 <-
#'  CSM_Phenotypecell_test %>%
#'      mutate(Phenotype = sample(unique(CSM_Phenotypecell_test$Phenotype),
#'                                size = nrow(CSM_Phenotypecell_test),
#'                                replace = TRUE))
#' CSM_Phenotypecell_test_3 <-
#'  CSM_Phenotypecell_test %>%
#'      mutate(Phenotype = sample(unique(CSM_Phenotypecell_test$Phenotype),
#'                                size = nrow(CSM_Phenotypecell_test),
#'                                replace = TRUE))
#'
#' #Find consensus phenotype labels between the three results--------------
#'DATA_Phenotypes_consensus <-
#'Consensus_phenotype_assigner(
#'  Random_1 = CSM_Phenotypecell_test_1,
#'  Random_2 = CSM_Phenotypecell_test_2,
#'  Random_3 = CSM_Phenotypecell_test_3,
#'
#'  Win_strategy = "Majority",
#'  Weights = c(2,1,1),
#'  No_consensus_value = NA_character_,
#'  Tie_break_method = "Random",
#'  N_cores = 2
#'  )
#' }
#'
#' @export

Consensus_phenotype_assigner <-
  function(...,
           Win_strategy,
           Arbitrary_threshold = NULL,
           Weights = NULL,
           No_consensus_value = NA_character_,
           Tie_break_method = NULL,
           N_cores = 1) {
    on.exit({
      future::plan("future::sequential")
      gc()
    })

    print("Checking provided function arguments")
    #Build a list with the methods selected
    Phenotype_tibble_list <- list(...)
    #Check that Phenotype_tibble_list is larger or equal to 2
    if (length(Phenotype_tibble_list) < 2)
      stop("DATA provided must include at least 2 methods")
    #Argument quality check
    if (!Win_strategy %in% c("Majority", "Arbitrary", "Difference"))
      stop("Win_strategy must be one of the following: Majority, Absolute_winner, Arbitrary")
    if (Win_strategy == "Arbitrary") {
      if (!all(
        is.numeric(Arbitrary_threshold),
        Arbitrary_threshold > 0.5,
        Arbitrary_threshold < 1
      ))
        stop("Arbitrary_threshold must be a numeric value larger than 0.5 and lower than 1")
    }
    if (Win_strategy == "Difference") {
      if (!all(
        is.numeric(Arbitrary_threshold),
        Arbitrary_threshold > 0,
        Arbitrary_threshold < 1
      ))
        stop("Arbitrary_threshold must be a numeric value larger than 0 and lower than 1")
    }
    if (Win_strategy == "Majority") {
      if(!Tie_break_method %in% c("Random", "No_consensus")) stop("Tie_break_method must be one of the following: Random, No_consensus")
    }
    if (!is.null(Weights)) {
      if (length(Weights) != length(Phenotype_tibble_list))
        stop(paste0("Weights must be of length ", length(Phenotype_tibble_list)))
      if (!all(is.numeric(Weights),
               dplyr::if_else(Weights > 0, true = TRUE, false = FALSE)))
        stop("Weights must be a numeric vector of values larger than 0")
      #Modify weights to sum 1
      Weights <- Weights / sum(Weights)
    }
    if (!all(!is.null(N_cores), N_cores >= 1 &
             N_cores %% 1 == 0))
      stop("N_cores must be an integer value > 0")
    #List quality check
    #Check that list has names
    if (any(names(Phenotype_tibble_list) == ""))
      stop("User must provide a name for each DATA provided")
    #Check that list contains adequate names
    Conflictive_tibbles <- !purrr::map_lgl(Phenotype_tibble_list, function(Tibble)
      identical(names(Tibble)[1:4], c("Cell_no", "X", "Y", "Subject_Names")))
    if (any(Conflictive_tibbles))
      stop(paste0(
        "The following DATA are not adequately formatted: ",
        stringr::str_c(names(Phenotype_tibble_list)[Conflictive_tibbles], collapse = ", ")
      ))
    Conflictive_tibbles <- !purrr::map_lgl(Phenotype_tibble_list, function(Tibble)
      "Phenotype" %in% names(Tibble))
    if (any(Conflictive_tibbles))
      stop(
        paste0(
          "The following DATA do not contain a column called Phenotype: ",
          stringr::str_c(names(Phenotype_tibble_list)[Conflictive_tibbles], collapse = ", ")
        )
      )
    #Check intersect between cells in each list
    All_cells <- unique(unlist(purrr::map(Phenotype_tibble_list, ~ .[["Cell_no"]])))
    Intersect_cells <- purrr::reduce(purrr::map(Phenotype_tibble_list, ~.[["Cell_no"]]), function(x, y) intersect(x, y))
    if(!all(All_cells %in% Intersect_cells)){
      Phenotype_tibble_list <-purrr::map(Phenotype_tibble_list, function(Tibble){
        return(Tibble %>% dplyr::filter(Cell_no %in% Intersect_cells))
      })
    }
    #Check intersect between phenotype labels within each list
    All_labels <- unique(unlist(purrr::map(Phenotype_tibble_list, ~.[["Phenotype"]])))
    Phenotypes_in_DATA <-
      t(
        purrr::map_dfc(Phenotype_tibble_list, function(Tibble){
          All_labels %in% unique(Tibble[["Phenotype"]])
        })
      )
    colnames(Phenotypes_in_DATA) <- All_labels
    Phenotypes_in_DATA <- as_tibble(Phenotypes_in_DATA)
    Phenotypes_in_DATA$DATA <- names(Phenotype_tibble_list)
    Phenotypes_in_DATA <- Phenotypes_in_DATA[c(ncol(Phenotypes_in_DATA), 1:(ncol(Phenotypes_in_DATA)-1))]
    Conflictive_rows <- apply(Phenotypes_in_DATA[-1], MARGIN = 1, function(Row) sum(Row) < (ncol(Phenotypes_in_DATA)-1))

    #print summary and ask the user if the process should continue
    if(!all(All_cells %in% Intersect_cells)){
      message(paste0("Some cells are not present in all DATA provided. Only cells present across all DATA provided will be used. ",
                     length(All_cells) - length(Intersect_cells),
                     " cells will be removed."
      )
      )

    }
    if(any(Conflictive_rows)){
      message("The following DATA provided do not harbor all the potential cell types")
      print(Phenotypes_in_DATA[Conflictive_rows, ])
    }
    if(!is.null(Weights)){
      names(Weights) <- names(Phenotype_tibble_list)
      print("The following weights will be used")
      print(Weights)
    }
    if(Win_strategy == "Majority"){
      if(Tie_break_method == "Random") message(paste0("A Random label between winners will be selected if no consensus is reached"))
      else message(paste0("The following value will be used if no consensus is reached: ", No_consensus_value))
    }


    #ASK THE USER
    answer <- menu(c("Proceed", "Abort"), title = "Should the analysis proceed")
    #If user decides to stop then abort function and return stop message
    if(answer == 2) stop("The function has been stopped.")

    #PROCEED WITH ELECTIONS
    print("Celebrating phenotype elections...")
    #Generate a tibble where each cell is a row and columns are each of the methods votes
    Vote_tibble <- suppressMessages(purrr::map_dfc(Phenotype_tibble_list, function(Tibble) Tibble[match(Intersect_cells, Tibble$Cell_no),] %>% dplyr::select(Phenotype)))
    names(Vote_tibble) <- names(Phenotype_tibble_list)
    Vote_tibble <-dplyr::bind_cols(tibble(Cell_no = Intersect_cells),
                                   Vote_tibble)
    #Turn your vote tibble into a list where each method asssigns a vote to each possibility
    Vote_tibble <-purrr::map(2:ncol(Vote_tibble), function(Column){
      #Generate the matrix of votes (1 if voted 0 if no voted)
      Tibble <-dplyr::bind_cols(Vote_tibble[1], Vote_tibble[Column])
      names(Tibble)[2] <- "Method"
      Tibble <- Tibble %>%dplyr::mutate(Vote = 1) %>% tidyr::pivot_wider(names_from = Method, values_from = Vote, id_cols = Cell_no)
      Tibble[is.na(Tibble)] <- 0

      #Add columns not present in the method selected
      Missing_values <- All_labels[which(!All_labels %in% names(Tibble))]
      if(length(Missing_values) >= 1){
        Cero_tibble <- matrix(0, nrow = nrow(Tibble), ncol = length(Missing_values))
        colnames(Cero_tibble) <- Missing_values
        Cero_tibble <- as_tibble(Cero_tibble)
        Tibble <-dplyr::bind_cols(Tibble, Cero_tibble)
      }

      #Reorder the final vote tibble so that all tibbles have the same order
      Tibble[c("Cell_no", All_labels)]
    })
    names(Vote_tibble) <- names(Phenotype_tibble_list)

    #Apply weights if required
    if(!is.null(Weights)){
      Vote_tibble <-purrr::map(names(Weights), function(Index){
        Weight <- Weights[Index]
        Tibble <- Vote_tibble[[Index]]
        Tibble[-1] <- Tibble[-1]*Weight
        return(Tibble)
      })
      names(Vote_tibble) <- names(Phenotype_tibble_list)
    }
    #If no weights then divide by the number of voters
    if(is.null(Weights)){
      Vote_tibble <-purrr::map(Vote_tibble, function(Tibble){
        Tibble <- Tibble
        Tibble[-1] <- Tibble[-1]/length(Vote_tibble)
        return(Tibble)
      })
    }

    #Generate a summary of votes (Used to graph the results)
    Vote_count_summary <-purrr::map_dfr(Vote_tibble, function(Tibble){
      purrr::map_dbl(Tibble[-1], sum)
    })
    Vote_count_summary$Method <- names(Vote_tibble)

    #Add the results for all methods
    Vote_tibble <-
      dplyr::bind_cols(
        Vote_tibble[[1]]["Cell_no"],
        purrr::reduce(purrr::map(Vote_tibble, ~.[-1] %>% as.matrix()),
                      function(x, y) as_tibble(x + y))
      )

    print("Counting votes...")
    #Apply Win_strategy to decide final label
    future::plan("future::multisession", workers = N_cores)
    options(future.globals.maxSize = Inf, future.rng.onMisuse = "ignore")
    furrr::furrr_options(scheduling = Inf)

    #MAJORITY
    if(Win_strategy == "Majority"){
      Vote_tibble <-
        furrr::future_map_dfr(1:nrow(Vote_tibble), function(Row){
          #Calculate the results for each cell
          Cell <- Vote_tibble[Row, ]
          Cell <- Cell %>% tidyr::pivot_longer(-1) %>% dplyr::arrange(desc(value))

          #If there is an absolute winner return a tibble with the results
          if(Cell$value[1] > Cell$value[[2]]){
            return(tibble(Cell_no = Cell$Cell_no[1],
                          Result = "Winner",
                          N_winners = 1,
                          Voting_result = Cell$value[1],
                          Phenotype = Cell$name[1]
            )
            )
          }

          #If there is a Tie return according to user preferences
          if(Cell$value[1] == Cell$value[[2]]){
            #Filter rows with winning values
            Winner_rows <- which(Cell$value == Cell$value[1])
            Cell <- Cell[Winner_rows, ]

            #obtain results
            Cell_no <- Cell$Cell_no[1]
            Result <- "Tie"
            N_winners <- nrow(Cell)
            Voting_result <- Cell$value[1]
            Phenotype <- if(Tie_break_method == "Random") sample(Cell$name, size = 1) else No_consensus_value
            return(tibble(Cell_no = Cell_no,
                          Result = Result,
                          N_winners = N_winners,
                          Voting_result = Voting_result,
                          Phenotype = Phenotype
            )
            )
          }

        }, .progress = TRUE)
    }
    #Difference
    if(Win_strategy == "Difference"){
      Vote_tibble <-
        furrr::future_map_dfr(1:nrow(Vote_tibble), function(Row){
          #Calculate the results for each cell
          Cell <- Vote_tibble[Row, ]
          Cell <- Cell %>% tidyr::pivot_longer(-1) %>% dplyr::arrange(desc(value))
          #If there is an winner above threshold with respect to the second highest rival
          if((Cell$value[1] - Cell$value[2]) >= Arbitrary_threshold){
            return(tibble(Cell_no = Cell$Cell_no[1],
                          Result = "Winner",
                          N_winners = 1,
                          Voting_result = Cell$value[1],
                          Phenotype = Cell$name[1]
            )
            )
          }
          #If no winner above threshold return NA and no winners will be asigned
          else{
            return(tibble(Cell_no = Cell$Cell_no[1],
                          Result = "No_winner",
                          N_winners = 0,
                          Voting_result = 0,
                          Phenotype = No_consensus_value
            )
            )
          }
        }, .progress = TRUE)
    }
    #ARBITRARY
    if(Win_strategy == "Arbitrary"){
      Vote_tibble <-
        furrr::future_map_dfr(1:nrow(Vote_tibble), function(Row){
          #Calculate the results for each cell
          Cell <- Vote_tibble[Row, ]
          Cell <- Cell %>% tidyr::pivot_longer(-1) %>% dplyr::arrange(desc(value))
          #If there is an winner above threshold return the result
          if(Cell$value[1] >= Arbitrary_threshold){
            return(tibble(Cell_no = Cell$Cell_no[1],
                          Result = "Winner",
                          N_winners = 1,
                          Voting_result = Cell$value[1],
                          Phenotype = Cell$name[1]
            )
            )
          }
          #If no winner above threshold return NA and no winners will be assigned
          else{
            return(tibble(Cell_no = Cell$Cell_no[1],
                          Result = "No_winner",
                          N_winners = 0,
                          Voting_result = 0,
                          Phenotype = No_consensus_value
            )
            )
          }
        }, .progress = TRUE)
    }
    future::plan("future::sequential")
    gc()

    #Generate a summary and graphs
    print(Vote_tibble %>% dplyr::count(Result))
    print(Vote_tibble %>% dplyr::filter(Voting_result > 0) %>% dplyr::summarise(Min_result = min(Voting_result),
                                                                                Max_result = max(Voting_result),
                                                                                Median_result = quantile(Voting_result, 0.5),
                                                                                Average_result = mean(Voting_result)))
    print(Vote_tibble %>% dplyr::count(Phenotype) %>% dplyr::arrange(desc(n)))
    plot(
      Vote_count_summary %>% tidyr::pivot_longer(cols = -Method) %>%
        ggplot(aes(x = name, y = value, color = Method, group = Method)) +
        geom_line(linewidth = 2) +
        geom_point(size = 3) +
        scale_y_continuous("Votes") +
        scale_x_discrete("")+
        theme_bw() +
        theme(legend.position = "bottom",
              legend.direction = "horizontal",
              axis.text = element_text(color = "black"))
    )

    #Bind the 1st element of the Phenotype tibble list to the consensus phenotype result
    Consensus_result <-
      dplyr::left_join(Phenotype_tibble_list[[1]] %>% dplyr::filter(Cell_no %in% Intersect_cells) %>% dplyr::select(-Phenotype),
                       Vote_tibble %>% dplyr::select(Cell_no, Phenotype),
                       by = "Cell_no")

    return(Consensus_result)
  }
