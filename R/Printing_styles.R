#' Prints text in color based in the cat function and print centered text
#'
#' Intended for internal use only
#'
#' @param text The text being printed.
#' @param color One of the following: "black", "red", "green", "yellow", "blue", "magenta", "cyan", "white".
#'
#' @details
#' Used in [Barrier_effect_calculator()], [Cell_in_edge_remover()], [Cell_to_pixel_distance_calculator()], [Distance_matrix_generator()], [Gcross_calculator()], [Random_Distance_matrix_generator()],
#' [Structure_morphometry_calculator()], [Tiled_Image_Clustering_function()], [Trio_Distance_matrix_generator()], [Trio_Random_Distance_matrix_generator()], [Cell_aligner_function()]
#'
#' This code for Colored_print is based on Joseph Crispell's answer to the stack overflow forum question
#' https://stackoverflow.com/questions/10802806/is-there-a-way-to-output-text-to-the-r-console-in-color
#'
#'
#' @returns the text with the formatted color and or horizontally centerede in the console.
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


#' @rdname Colored_print
Centered_print <- function(text) {
  #Find the width of the console
  width <- getOption("width")

  #Divide message into different lines
  lines <- strsplit(text, "\n")[[1]]

  #Print the lines
  for (line in lines) {
    padding <- max(0, floor((width - nchar(line)) / 2))
    cat(strrep(" ", padding), line, "\n", sep = "")
  }
}

