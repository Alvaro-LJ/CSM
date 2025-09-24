#' Obtain CPU information and available RAM
#' @return A list containing information regarding CPU, RAM, R version being used and Operating system

Describe_my_computing_unit <-
  function() {
    if(stringr::str_detect(benchmarkme::get_platform_info()$OS.type, "windows")){
      list(
        CPU =  paste0("CPU: ", benchmarkme::get_cpu()$model_name, ". Number of available threads (cores): ", benchmarkme::get_cpu()$no_of_cores),
        RAM =  paste0("Total available RAM memory: approximately ",
                      round(
                        as.double(substr(shell('wmic path Win32_ComputerSystem get TotalPhysicalMemory/value', translate = F, intern = TRUE)[[3]],
                                         start = 21,
                                         stop = nchar(shell('wmic path Win32_ComputerSystem get TotalPhysicalMemory/value', translate = F, intern = TRUE)[[3]])-2))/104857600,
                        digits = 2),
                      " Gb"),
        R_VERSION =  paste0(benchmarkme::get_r_version()$version.string, " - ",
                            benchmarkme::get_r_version()$nickname, ". Version type: ", benchmarkme::get_r_version()$arch),
        OS = paste0("You are running R on ", benchmarkme::get_platform_info()$OS.type)
      )
    }

    else (
      benchmarkme::get_sys_details(byte_compiler = F, linear_algebra = F,  locale = F, installed_packages = F, machine = F)
    )}
