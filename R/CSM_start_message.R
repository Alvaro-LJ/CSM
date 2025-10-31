#' This regulates the package start-up message
#'
#'
#' @param libname Used to obtain library name
#' @param pkgname Used to obtain package name
#'
#' @keywords Internal
.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    cat(
      stringr::str_c("Comprehensive Spatial Methods - CSM (version ", utils::packageDescription("CSM")$Version, ")", " loaded"),
      "",
      "Please see package help with help(package = 'CSM')",
      "",
      "For further information please visit our github repositories: https://github.com/Alvaro-LJ",
      "",
      "Bug and suggestions can be reported through our github profile",
      "",
      "Thanks for using CSM!",
      fill = 4
    )
  )
}



