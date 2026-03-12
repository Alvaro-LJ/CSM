#' Prints text in color based in the cat function
#'
#' Intended for internal use only
#'
#' @param text The text being printed.
#' @param color One of the following: "black", "red", "green", "yellow", "blue", "magenta", "cyan", "white".
#'
#' @details
#' Used in [Barrier_effect_calculator()]
#' 
#' This code is based on Joseph Crispell's answer to the stack overflow forum question 
#' https://stackoverflow.com/questions/10802806/is-there-a-way-to-output-text-to-the-r-console-in-color
#'
#'
#' @returns the text with the formatted color.
#' 
#' @keywords Internal

Colored_print <- function(text, color = "green") {
  
  #Note ANSI color codes
  colour_codes <- list(
    "black" = 30,
    "red" = 31,
    "green" = 32,
    "yellow" = 33,
    "blue" = 34,
    "magenta" = 35,
    "cyan" = 36,
    "white" = 37
  )
  # Create ANSI version of colour
  ansi_colour <- paste0("\033[", colour_codes[[color]], "m")
  
  return(cat(ansi_colour, text, "\033[0m\n"))
  }

