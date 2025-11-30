#' Obtain CPU information and available RAM
#'
#' The function returns machine characteristics like CPU model and number of cores. It will return the available RAM. It will also return the OS and R version.
#'
#' @return A list containing information regarding CPU, RAM, R version being used and Operating system.
#'
#' @examples
#' \dontrun{
#' Describe_my_computing_unit()
#' }
#'
#' @export

Describe_my_computing_unit <-
  function() {

    #Check that benchmarkme is installed
    if(!requireNamespace("benchmarkme", quietly = FALSE)) stop(
      paste0("benchmarkme CRAN package is required to execute the function. Please install using the following code: ",
             expression(install.packages("benchmarkme")))
    )

    #Proceed
    if(stringr::str_detect(benchmarkme::get_platform_info()$OS.type, "windows")){
      list(
        CPU =  try(paste0("CPU: ", benchmarkme::get_cpu()$model_name, ". Number of available threads (cores): ", benchmarkme::get_cpu()$no_of_cores)),
        RAM =  try(paste0("Total available RAM memory: approximately ",
                      shell('powershell -Command "& {\\"{0:N2} GB\\" -f ((Get-CimInstance -ClassName Win32_ComputerSystem).TotalPhysicalMemory / 1GB)}"', intern = TRUE))),
        R_VERSION =  try(paste0(benchmarkme::get_r_version()$version.string, " - ",
                            benchmarkme::get_r_version()$nickname, ". Version type: ", benchmarkme::get_r_version()$arch)),
        OS = try(paste0("You are running R on ", benchmarkme::get_platform_info()$OS.type))
      )
    }

    else (
      benchmarkme::get_sys_details(byte_compiler = F, linear_algebra = F,  locale = F, installed_packages = F, machine = F)
    )}



