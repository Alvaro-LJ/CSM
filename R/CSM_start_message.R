#' This regulates the package start-up message
#'
#'
#' @param libname Used to obtain library name
#' @param pkgname Used to obtain package name
#'
#' @keywords Internal
.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    paste(
      stringr::str_c("Comprehensive Spatial Methods - CSM (version ", utils::packageDescription("CSM")$Version, ")", " loaded"),
      "",
      "Please see package help with help(package = 'CSM')",
      "",
      "To see examples of use please see CSM vignettes with browseVignettes(package = 'CSM')",
      "",
      "For further information please visit our github repository: https://github.com/Alvaro-LJ/CSM",
      "",
      "Bug and suggestions can be reported through our github repository: https://github.com/Alvaro-LJ/CSM/issues",
      "",
      "Thanks for using CSM!",
      sep = "\n"
    )
  )
}



