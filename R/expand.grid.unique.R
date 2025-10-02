#' Generates a tibble with unique combinations between two vectors.
#'
#' Intended for internal use only. In contrast to expand.grid only unique combinations are preserved.
#'
#' @param x A vector
#' @param y A vector
#'
#' @returns A dataframe containing unique combinations
#'
#' @keywords Internal

expand.grid.unique <-
  function(x, y, include.equals=FALSE)
  {
    x <- unique(x)
    y <- unique(y)
    g <- function(i)
    {
      z <- setdiff(y, x[seq_len(i-include.equals)])

      if(length(z)) cbind(x[i], z, deparse.level=0)
    }
    do.call(rbind, lapply(seq_along(x), g))
  }
