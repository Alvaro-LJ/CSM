#' Analyzes the expression of features within a cell category
#'
#' The expression level of various features that have been set aside can be analyzed within cell categories (either a cell phenotype or neighborhood label).
#' These features can be set aside with [Data_set_aside()] function.
#'
#' @param DATA A dataframe or tibble containing cell feature data.
#' @param Target_Variable A character indicating the name of the column containing cell categories (for example cell phenotypes).
#' @param Targets_Included A character vector indicating which cell categories should be analyzed (for example a specific cell phenotype to be included).
#' @param DATA_Aside A dataframe or tibble containing the functional features to be analyzed.
#' @param Threshold_functional_Markers A logical value indicating if features should be thresholded before analysis.
#' @param Levels If thresholding is required, a integer value indicating the number of desired levels to be calculated for multilevel thresholding (see details).
#'
#' @details
#' Multilevel thresholding is calculated using imagerExtra::ThresholdML function
#'
#' @returns A list containing cell feature information and a summary of feature expression in target cells. Also generates summary plots.
#'
#' @export

Cell_functional_assessment <-
  function(DATA = NULL,
           Target_Variable = NULL,
           Targets_Included = NULL,
           DATA_Aside = NULL,
           Threshold_functional_Markers = NULL,
           Levels = NULL){

    #Check suggested packages
    if(!requireNamespace("imagerExtra", quietly = FALSE)) stop(
      paste0("imagerExtra CRAN package is required to execute the function. Please install using the following code: ",
             expression(install.packages("imagerExtra")))
    )


    #Check variables
    if(!all(Target_Variable %in% names(DATA))) {
      Missing_arguments <- Target_Variable[!Target_Variable %in% names(DATA)]
      stop(paste0(stringr::str_c(Missing_arguments, collapse = ", "), " not found in DATA"))
    }
    if(!is.logical(Threshold_functional_Markers)) {
      stop("Threshold_functional_Markers must be a logical value")
    }

    #Import phenotype DATA and generate unique ID
    DATA <- DATA %>%dplyr::mutate(Unique_ID = stringr::str_c(Cell_no, Subject_Names, sep = "_"))
    DATA <- DATA %>% dplyr::select(1:4, Unique_ID, all_of(Target_Variable))

    #Check that target variable contains the desired targets before proceeding
    if(!any(Targets_Included %in% DATA[[6]])) {
      stop(paste0(stringr::str_c(Targets_Included, collapse = ", "), " not found in ", Target_Variable))
    }

    #Filter desired cells to be further analyzed
    DATA <- DATA[DATA[[6]] %in% Targets_Included, ]

    #Now select the Data aside and join the results to the original data
    DATA_Aside <- DATA_Aside %>%dplyr::mutate(Unique_ID = stringr::str_c(Cell_no, Subject_Names, sep = "_"))
    DATA_Aside <- DATA_Aside %>% dplyr::select(-c(1:4))
    DATA_Aside <- DATA_Aside[c(ncol(DATA_Aside), 1:(ncol(DATA_Aside)-1))]

    if(!Threshold_functional_Markers){
      #Join both tibbles
      DATA_Joined <-dplyr::left_join(DATA, DATA_Aside, by = "Unique_ID") %>% dplyr::select(-Unique_ID)
      #Change name to allow function generalizaion
      names(DATA_Joined)[5] <- "Variable"

      #Plot the functional marker by each Target variable
      plot(
        DATA_Joined %>% tidyr::pivot_longer(-c(1:5)) %>%
          ggplot(aes(x = Variable, y = log10(value+0.00001), fill = Variable)) + facet_wrap(~name, "free", nrow = 1, ncol = (ncol(DATA_Joined) - 5)) +
          geom_boxplot() +
          cowplot::theme_cowplot() +
          guides(fill = "none")+
          scale_x_discrete("") +
          scale_y_continuous("Log 10 functional Marker expression across target cells")
      )

      #Generate a summary tibble with the functional markers across the Target variables
      By_sample_Results <-
        DATA_Joined %>% dplyr::select(-Cell_no, -X, -Y) %>% group_by(Subject_Names, Variable) %>% summarize_at(vars(-group_cols()), list(min = ~min(.x, na.rm = T),
                                                                                                                                         p25 = ~quantile(.x, 0.25, na.rm = T),
                                                                                                                                         Average = ~mean(.x, na.rm = T),
                                                                                                                                         p50 = ~quantile(.x, 0.5, na.rm = T),
                                                                                                                                         p75 = ~quantile(.x, 0.75, na.rm = T),
                                                                                                                                         max = ~max(.x, na.rm = T))) %>%
        ungroup() %>% tidyr::pivot_longer(-c(1:2)) %>% dplyr::mutate(Final_Variable = stringr::str_c(Variable, name, sep = "_")) %>%
        dplyr::select(-Variable, -name) %>%
        tidyr::pivot_wider(names_from = Final_Variable, values_from = value)

      #Return the target Variable name to its original value
      names(DATA_Joined)[5] <- Target_Variable

      return(list(DATA_functional_markers = DATA_Joined,
                  By_sample_results = By_sample_Results))
    }

    else if(Threshold_functional_Markers){
      if(Levels %%1 != 0){
        stop("Levels must be and integer value")
      }
      print("Thresholding Markers")
      DATA_Aside <-dplyr::bind_cols(DATA_Aside[1],
                                    purrr::map_df(DATA_Aside[-1],
                                                  function(z){
                                                    if(length(unique(z))>Levels){ #requires at least n Levels to be calculated
                                                      as.double(imagerExtra::ThresholdML(imager::cimg(array(z, dim = c(1, length(z), 1, 1))), k = (Levels-1))) #K to specify the amount of cut-off points
                                                    }else(NA)
                                                  })
      )
      #Join both tibbles
      DATA_Joined <-dplyr::left_join(DATA, DATA_Aside, by = "Unique_ID") %>% dplyr::select(-Unique_ID)
      #Change name to allow function generalizaion
      names(DATA_Joined)[5] <- "Variable"

      #Plot the results
      plot(
        DATA_Joined %>% tidyr::pivot_longer(-c(1:5)) %>%
          ggplot(aes(x = Variable, y = 1)) + facet_wrap(~name, "free", nrow = 1, ncol = (ncol(DATA_Joined) - 5)) +
          geom_col(aes(fill = as.factor(value)), position = position_fill(reverse = T),  linewidth = 1) +
          cowplot::theme_cowplot() +
          scale_x_discrete("") +
          scale_y_continuous("%") +
          scale_fill_discrete("Level")
      )

      #Prepare by sample results

      #First generate a tibble with the Subject_Names the Target variable and the functional Markers
      For_Counts <- DATA_Joined %>% dplyr::select(-Cell_no, -X, -Y)

      #Prepare the overall target cell counts to calculate percentages
      Cell_count_tibble <- For_Counts %>% group_by(Subject_Names) %>% dplyr::count(Variable) %>% ungroup() %>%
        tidyr::pivot_wider(names_from = Variable, values_from = n)
      Cell_count_tibble[is.na(Cell_count_tibble)] <- 0
      #Make sure cell count tibble is arranged by subject_names
      Cell_count_tibble <- Cell_count_tibble %>% dplyr::arrange(Subject_Names)


      #Now count the functional marker for each subject name and target cell
      By_sample_Results <-
        #Iterate by functional marker
        purrr::map_dfc(3:ncol(For_Counts), function(var){
          #Generate a simple tibble with 3 columns, Subject_Names, Target cell and functional marker
          Interim <-dplyr::bind_cols(For_Counts[1:2], For_Counts[var])

          #Change the name of the functional marker to make it generizable
          names(Interim)[[3]] <- "Functional_Marker"
          #Prepare the cell count data
          Interim <- Interim %>% group_by(Subject_Names, Variable) %>% dplyr::count(Functional_Marker) %>% ungroup() %>%
            dplyr::mutate(Final_Variable = stringr::str_c(Variable, names(For_Counts)[var], Functional_Marker, sep = "_")) %>%
            dplyr::select(-Variable,-Functional_Marker) %>%
            tidyr::pivot_wider(names_from = Final_Variable, values_from = n)
          Interim[is.na(Interim)] <- 0

          #Make sure both tibbles are ordered in the same way with respect to the Subject Names column
          Interim <- Interim %>% dplyr::arrange(Subject_Names)

          #Generate the proportion tibble with the cell counts
          Prop_tibble <-purrr::map_dfc(names(Cell_count_tibble)[-1], function(Target) {
            Absolute_cell_counts <- Interim %>% dplyr::select(contains(Target))
            For_percentage <- Cell_count_tibble %>% dplyr::select(contains(Target))

            Prop_tibble <- as_tibble(Absolute_cell_counts/For_percentage[[1]])
            names(Prop_tibble) <- stringr::str_c("PROP_", names(Prop_tibble), sep ="")
            Prop_tibble
          })

          dplyr::bind_cols(Interim[-1], Prop_tibble)
        })

      #Generate the final tibble with an adequate order
      By_sample_Results <- dplyr::bind_cols(
        Cell_count_tibble[1],
        purrr::map_dfc(names(Cell_count_tibble)[-1], function(Target){
          dplyr::bind_cols(Cell_count_tibble %>% dplyr::select(contains(Target)),
                           By_sample_Results %>% dplyr::select(contains(Target)))
        })

      )

      #Return the target Variable name to its original value
      names(DATA_Joined)[5] <- Target_Variable

      return(list(DATA_functional_markers = DATA_Joined,
                  By_sample_results = By_sample_Results))

    }
  }
