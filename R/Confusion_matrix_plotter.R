#' Generates a confusion matrix between two methods used to label cells.
#'
#' The function generates a confusion matrix between two cell labeling methods (cell phenotype or neighborhood labels for example).
#'
#' @param ... Two dataframes or tibbles containing a column with cell labels.
#' @param Variable A character value of the column name containing the labels. The name must be shared across the two datasets provided.
#'
#' @returns Generates a confusion matrix plot
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
#'
#' #Generate the confusion matrix-------------------------------------------
#'
#'Confusion_matrix_plotter(
#'    Random_1 = CSM_Phenotypecell_test_1,
#'    Random_2 = CSM_Phenotypecell_test_2,
#'    Variable = "Phenotype"
#'    )
#' }
#'
#' @export

Confusion_matrix_plotter <-
  function(...,
           Variable){

    DATA_list <- list(...)

    if(length(DATA_list) != 2) stop("Only two datasets must be provided")
    #Check DATA_1 is adequately formatted
    if(!identical(names(DATA_list[[1]])[1:4],  c("Cell_no", "X", "Y", "Subject_Names"))) stop(paste0(names(DATA_list)[1], " should have an adecuate format"))
    #Check DATA_2 is adequately formatted
    if(!identical(names(DATA_list[[2]])[1:4],  c("Cell_no", "X", "Y", "Subject_Names"))) stop(paste0(names(DATA_list)[2], " should have an adecuate format"))

    #Check if "Phenotype" is present in both dataset
    if(!Variable %in% names(DATA_list[[1]])) stop(paste0("A column named", Variable,  " is not present in ", names(DATA_list)[1]))
    if(!Variable %in% names(DATA_list[[2]])) stop(paste0("A column named", Variable,  " is not present in ", names(DATA_list)[1]))

    #Then proceed calculating the intersect
    Common_cells <- intersect(DATA_list[[1]]$Cell_no, DATA_list[[2]]$Cell_no)

    #If any method is devoid of cells then print a warning
    Methods_without_cells <-purrr::map_lgl(DATA_list, ~all(Common_cells %in% .$Cell_no))
    if(any(!Methods_without_cells)) warning("Only Cells present in both methods will be used")

    #Generate the final tibble with the interesect
    Tibble <- tibble(Cell_no = Common_cells)
    Tibble <-dplyr::left_join(Tibble, DATA_list[[1]] %>% dplyr::select(Cell_no, all_of(Variable)), by = "Cell_no")
    names(Tibble)[2] <- names(DATA_list)[1]
    Tibble <-dplyr::left_join(Tibble, DATA_list[[2]] %>% dplyr::select(Cell_no, all_of(Variable)), by = "Cell_no")
    names(Tibble)[3] <- names(DATA_list)[2]

    #Generate the count dataframe
    Tibble <-
      suppressWarnings(
        Tibble %>%
          group_by_(names(Tibble)[2], names(Tibble)[3]) %>% count() %>%
          tidyr::pivot_wider(names_from = names(Tibble)[2], values_from = n)
      )
    Tibble[is.na(Tibble)] <- 0
    Tibble <- Tibble %>% tidyr::pivot_longer(-1, names_to = names(DATA_list)[1])

    #Define the order
    Order <- sort(unique(Tibble[[names(Tibble)[2]]]))
    #Change the plot
    names(Tibble)[c(1,2)] <- c("Method_2", "Method_1")

    #Generate the final plot
    Final_plot <-
      Tibble %>%
      ggplot(aes(x = Method_1, y = Method_2)) +
      geom_tile(aes(fill = value), data = Tibble %>% dplyr::filter(Method_1 != Method_2)) +
      geom_tile(fill = "#aef2a0", data = Tibble %>% dplyr::filter(Method_1 == Method_2)) +
      geom_text(aes(label = value), size = 8) +
      scale_x_discrete(names(DATA_list)[1], limits = Order)+
      scale_y_discrete(names(DATA_list)[2], limits = rev(Order)) +
      scale_fill_viridis_c("N", option = "C") +
      theme_minimal() +
      theme(panel.grid = element_blank(),
            axis.text = element_text(size = 12, color = "black"),
            axis.title = element_text(size = 12, face = "bold"),
            legend.title = element_text(size = 12, hjust = 0.5))

    #Plot
    plot(
      Final_plot
    )

    #Return
    Final_plot
  }
