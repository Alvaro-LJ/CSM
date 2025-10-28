#' Calculates distance between a Cell of Origin and two target cell populations
#'
#' Intended for internal use only
#'
#' @param DATA A dataframe or tibble containing cell phenotype labels from a single image.
#' @param cell_A A character value indicating the cell phenotype label of the Cell of Origin.
#' @param cell_B A character value indicating the cell phenotype label of the Target cell 1.
#' @param cell_C A character value indicating the cell phenotype label of the Target cell 2.
#'
#' @details
#' Used in [Trio_Distance_matrix_generator()]
#'
#'
#' @returns A list containing two tibbles. Tibble A_to_B contains the distance matrix between COO and Target cell 1. Tibble A_to_C contains the distance matrix between COO and Target cell 2
#' @keywords Internal

Advanced_Distance_function_A_B_C <-
  function(DATA,
           cell_A,
           cell_B,
           cell_C){
    Tibble_A <- DATA %>% dplyr::filter(Phenotype == cell_A)
    Tibble_B <- DATA %>% dplyr::filter(Phenotype == cell_B)
    Tibble_C <- DATA %>% dplyr::filter(Phenotype == cell_C)

    #First tibble A to B cell types
    Tibble_interim <- as_tibble(expand.grid(Tibble_A$Cell_no, Tibble_B$Cell_no))
    names(Tibble_interim) <- c("Cell_no", "NOTHING")
    Tibble_interim <- dplyr::left_join(Tibble_interim, Tibble_A, by = "Cell_no") %>% dplyr::select(1:4)
    names(Tibble_interim) <- c("Cell_A", "Cell_no", "X_Cell_A", "Y_Cell_A")
    Tibble_interim <- dplyr::left_join(Tibble_interim, Tibble_B, by = "Cell_no") %>% dplyr::select(1:6)
    names(Tibble_interim) <- c("Cell_A_no", "Cell_B_no", "X_Cell_A", "Y_Cell_A", "X_Cell_B", "Y_Cell_B")

    Final_tibble_A_B <- Tibble_interim %>%dplyr::mutate(DIST = sqrt((X_Cell_A - X_Cell_B)^2 + (Y_Cell_A - Y_Cell_B)^2)) %>%
      dplyr::select(Cell_A_no, Cell_B_no, DIST) %>% tidyr::pivot_wider(id_cols = Cell_A_no, names_from = Cell_B_no, values_from = DIST) %>%
      dplyr::rename("Cell_Of_Origin_no" = "Cell_A_no")

    #Second tibble A to C cell types
    Tibble_interim <- as_tibble(expand.grid(Tibble_A$Cell_no, Tibble_C$Cell_no))
    names(Tibble_interim) <- c("Cell_no", "NOTHING")
    Tibble_interim <-dplyr::left_join(Tibble_interim, Tibble_A, by = "Cell_no") %>% dplyr::select(1:4)
    names(Tibble_interim) <- c("Cell_A", "Cell_no", "X_Cell_A", "Y_Cell_A")
    Tibble_interim <-dplyr::left_join(Tibble_interim, Tibble_C, by = "Cell_no") %>% dplyr::select(1:6)
    names(Tibble_interim) <- c("Cell_A_no", "Cell_C_no", "X_Cell_A", "Y_Cell_A", "X_Cell_C", "Y_Cell_C")
    Final_tibble_A_C <- Tibble_interim %>% dplyr::mutate(DIST = sqrt((X_Cell_A - X_Cell_C)^2 + (Y_Cell_A - Y_Cell_C)^2)) %>%
      dplyr::select(Cell_A_no, Cell_C_no, DIST) %>% tidyr::pivot_wider(id_cols = Cell_A_no, names_from = Cell_C_no, values_from = DIST) %>%
      dplyr::rename("Cell_Of_Origin_no" = "Cell_A_no")

    #Final result
    return(
      list(A_to_B <- Final_tibble_A_B,
           A_to_C <- Final_tibble_A_C)
    )
  }
