#' Calculates distance between two cell types for a given image selecting a random cell of origin
#'
#' Intended for internal use only

#' @param DATA A dataframe or tibble containing cell phenotype labels from a single image.
#' @param cell_A A character value indicating the cell phenotype label of the Cell of Origin.
#' @param cell_B A character value indicating the cell phenotype label of the Target cell.
#'
#' @details
#' Used in [Random_Distance_matrix_generator()]
#'
#'
#' @returns A tibble containing the distance matrix for a single Cell of Origin and all target cells.
#' @keywords Internal

Advanced_Distance_function_A_B_single <-
  function(DATA,
           cell_A,
           cell_B){
    Tibble_A <- DATA %>% dplyr::filter(Phenotype == cell_A) %>% dplyr::slice_sample(n = 1)
    Tibble_B <- DATA %>% dplyr::filter(Phenotype == cell_B)
    Tibble_interim <- as_tibble(expand.grid(Tibble_A$Cell_no, Tibble_B$Cell_no))
    names(Tibble_interim) <- c("Cell_no", "NOTHING")
    Tibble_interim <- dplyr::left_join(Tibble_interim, Tibble_A, by = "Cell_no") %>% dplyr::select(1:4)
    names(Tibble_interim) <- c("Cell_A", "Cell_no", "X_Cell_A", "Y_Cell_A")
    Tibble_interim <- dplyr::left_join(Tibble_interim, Tibble_B, by = "Cell_no") %>% dplyr::select(1:6)
    names(Tibble_interim) <- c("Cell_A_no", "Cell_B_no", "X_Cell_A", "Y_Cell_A", "X_Cell_B", "Y_Cell_B")
    return(
      Tibble_interim %>% dplyr::mutate(DIST = sqrt((X_Cell_A - X_Cell_B)^2 + (Y_Cell_A - Y_Cell_B)^2)) %>%
        dplyr::select(Cell_A_no, Cell_B_no, DIST) %>% tidyr::pivot_wider(id_cols = Cell_A_no, names_from = Cell_B_no, values_from = DIST)
    )
  }
