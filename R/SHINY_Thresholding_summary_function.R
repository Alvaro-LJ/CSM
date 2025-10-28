#' Generates a summary of threshold results designed to work with dedicated shiny app
#'
#' Intended for internal use only

#' @param DATA Dataframe or tibble containing cell features.
#' @param DATA_Thresholded Dataframe or tibble containing thresholded cell features.
#' @param LOCAL Logical value indicating if thresholds are calculated locally.
#' @param CASE A character indicating which image needs to be summarized.
#'
#' @details
#' Used in [Thresholding_tester_app()]
#'
#'
#' @returns A list containing the histogram, the cut-off values and the cells above and below threshold for the selected case
#' @keywords Internal

SHINY_Thresholding_summary_function <-
  function(DATA = NULL,
           DATA_Thresholded = NULL,
           LOCAL = NULL,
           CASE = NULL) {

    #First GLOBAL thresholds (will not vary from case to case)
    if(!as.logical(LOCAL)){

      #First GLOBAL 2 level thresholding
      if(length(unique(DATA_Thresholded$Marker)) < 3){

        #Find the threshold
        Thresholds <- min(DATA$Marker[DATA_Thresholded[[5]]])
        #Generate tibble
        Threshold_results <- tibble(Threshold_level =stringr::str_c("Threshold_Level_", 1:length(Thresholds)),
                                    value = Thresholds)

        #Draw the histogram
        Histogram_plot <- DATA %>% dplyr::filter(Subject_Names == CASE) %>%
          ggplot(aes(x = Marker)) +
          geom_histogram(bins = 1000) +
          cowplot::theme_cowplot() +
          geom_vline(aes(xintercept = value), data = Threshold_results, color = "red", linewidth = 1.2) +
          scale_x_continuous("Marker intensity") +
          scale_y_discrete("Number of cells")

        #Calculate positive cells and total number of cells in the sample and place data in a tibble
        DATA_case <- DATA_Thresholded %>% dplyr::filter(Subject_Names == CASE)

        Sample_summary <- tibble(Cell = c("Positive Cells", "Total Cells"), number = c(sum(DATA_case[[5]]), nrow(DATA_case)))

        #Return the list
        return(list(Histogram = Histogram_plot,
                    Threshold_tibble = Threshold_results,
                    Summary = Sample_summary
        )
        )
      }

      #Then GLOBAL multi-level thresholding
      else if(length(unique(DATA_Thresholded$Marker)) >= 3){
        #Find the threshold
        Thresholds <-purrr::map_dbl(unique(DATA_Thresholded$Marker), function(Level) {
          min(DATA$Marker[DATA_Thresholded$Marker == Level])
        })
        #Generate tibble
        Threshold_results <- tibble(Threshold_level =stringr::str_c("Threshold_Level_", 1:length(Thresholds)),
                                    value = Thresholds)
        #Draw the histogram
        Histogram_plot <- DATA %>% dplyr::filter(Subject_Names == CASE) %>%
          ggplot(aes(x = Marker)) +
          geom_histogram(bins = 1000) +
          cowplot::theme_cowplot() +
          geom_vline(aes(xintercept = value), data = Threshold_results, color = "red", linewidth = 1.2) +
          scale_x_continuous("Marker intensity") +
          scale_y_discrete("Number of cells")

        #Calculate positive cells and total number of cells in the sample and place data in a tibble
        DATA_case <- DATA_Thresholded %>% dplyr::filter(Subject_Names == CASE)

        Sample_summary <- DATA_case %>% dplyr::count(Marker)
        names(Sample_summary) <- c("Cell", "number")
        Sample_summary[[1]] <-stringr::str_c("Cell_", as.character(Sample_summary[[1]]))
        Nrow_Sample_summary <- nrow(Sample_summary)
        Sample_summary[Nrow_Sample_summary + 1, 1] <- "Total Cells"
        Sample_summary[Nrow_Sample_summary + 1, 2] <- nrow(DATA_case)

        #Return the list
        return(list(Histogram = Histogram_plot,
                    Threshold_tibble = Threshold_results,
                    Summary = Sample_summary)
        )
      }
    }

    #Then LOCAL thresholds (will vary from case to case)
    else if(as.logical(LOCAL)){
      #Filter cases by Case_id
      DATA <- DATA %>% dplyr::filter(Subject_Names == CASE)
      DATA_Thresholded <- DATA_Thresholded %>% dplyr::filter(Subject_Names == CASE)

      #First LOCAL 2 level thresholding
      if(length(unique(DATA_Thresholded$Marker)) < 3){

        #Find the threshold
        Thresholds <- min(DATA$Marker[DATA_Thresholded[[5]]])
        #Generate tibble
        Threshold_results <- tibble(Threshold_level =stringr::str_c("Threshold_Level_", 1:length(Thresholds)),
                                    value = Thresholds)

        #Draw the histogram
        Histogram_plot <- DATA %>% dplyr::filter(Subject_Names == CASE) %>%
          ggplot(aes(x = Marker)) +
          geom_histogram(bins = 1000) +
          cowplot::theme_cowplot() +
          geom_vline(aes(xintercept = value), data = Threshold_results, color = "red", linewidth = 1.2) +
          scale_x_continuous("Marker intensity") +
          scale_y_discrete("Number of cells")

        #Calculate positive cells and total number of cells in the sample and place data in a tibble
        DATA_case <- DATA_Thresholded

        Sample_summary <- tibble(Cell = c("Positive Cells", "Total Cells"), number = c(sum(DATA_case[[5]]), nrow(DATA_case)))

        #Return the list
        return(list(Histogram = Histogram_plot,
                    Threshold_tibble = Threshold_results,
                    Summary = Sample_summary)
        )
      }

      #Then LOCAL multi-level thresholding
      else if(length(unique(DATA_Thresholded$Marker)) >= 3){


        #Return the list
        return(list(Histogram = ggplot() + geom_text(aes(x = 1, y = 1, label = "LOCAL MULTI-LEVEL THRESHOLDING\nIS\nSTRONGLY DISCOURAGED"),
                                                     color = "white",
                                                     size = 8) +
                      cowplot::theme_cowplot()+
                      theme(axis.line = element_blank(),
                            axis.ticks = element_blank(),
                            axis.text = element_blank(),
                            axis.title = element_blank(),
                            panel.background = element_rect(fill = "black")),
                    Threshold_tibble = tibble(Warning = "Local Multi-level thresholding is not a good idea"),
                    Summary = tibble(Warning = "Local Multi-level thresholding is not a good idea"))
        )
      }
    }
  }
