#' Analyzes the spatial interaction between cell types according to the SPIAT package
#'
#' The function calculates cell interaction pattern based on entropy gradients. For more information see related publication: https://doi.org/10.1038/s41467-023-37822-0.
#'
#' @param DATA_SPIAT A list containing SPIAT objects generated using [SPIAT_object_generator()].
#' @param Gradient_start A numeric value indicating where should the sampling process start.
#' @param Gradient_stop A numeric value indicating where should the sampling process end.
#' @param Gradient_sampling A numeric value indicating the distance to be sampled.
#' @param Phenotypes_included A character vector indicating the phenotype labels that will be included in the analysis.
#'
#' @returns A tibble containing a summary of the analysis result by image.
#'
#' @seealso [SPIAT_object_generator()].
#' @export

SPIAT_entropy_gradient_generator <-
  function(DATA_SPIAT = NULL,
           Gradient_start = NULL,
           Gradient_stop = NULL,
           Gradient_sampling = NULL,
           Phenotypes_included = NULL
  ) {

    #Check suggested packages
    if(!requireNamespace("SPIAT", quietly = TRUE)) stop(
      paste0("SPIAT Bioconductor package is required to execute the function. Please install using the following code: ",
             expression({
               if (!require("BiocManager", quietly = TRUE))
                 install.packages("BiocManager")

               BiocManager::install("SPIAT")
             })
      )
    )

    #Check arguments by generating a argument check vector and message vector
    Argument_checker <- c(DATA_SPIAT_OK = all(purrr::map_lgl(DATA_SPIAT, function(Image) class(Image) == "SpatialExperiment")),
                          Gradient_start_OK = all(is.numeric(Gradient_start), Gradient_start >= 0, Gradient_start < Gradient_stop),
                          Gradient_stop = all(is.numeric(Gradient_stop), Gradient_stop >= 0),
                          Gradient_sampling = all(is.numeric(Gradient_sampling), (Gradient_stop - Gradient_start) >= 0) ,
                          Phenotypes_included_OK = all(Phenotypes_included %in% unique(unlist(purrr::map(DATA_SPIAT, function(x) x@colData@listData$Phenotype))))
    )
    Stop_messages <- c(DATA_SPIAT_OK = "DATA_SPIAT must be generated with the SPIAT_object_generator function",
                       Gradient_start_OK = "Gradient_start must be a positive value smaller than Gradient_stop",
                       Gradient_stop = "Gradient_stop must be a positive value greater than Gradient_start",
                       Gradient_sampling = "Gradient sampling must be a positive value smaller than the differentce between Gradient_stop and Gradient_start",
                       Phenotypes_included_OK =stringr::str_c("Phenotypes_included must be any of the following: ",
                                                              stringr::str_c(unique(unlist(purrr::map(DATA_SPIAT, function(x) x@colData@listData$Phenotype))), collapse = ", "),
                                                              collapse = "")
    )
    #Check arguments and stop if necessary
    if(!all(Argument_checker)){
      stop(cat(Stop_messages[!Argument_checker],
               fill = sum(!Argument_checker)))
    }

    #Get SPIAT DATA
    DATA_SPIAT <- DATA_SPIAT

    #Select the gradients for entropy computation
    gradient_pos <- seq(from = Gradient_start, to = Gradient_stop, by = Gradient_sampling)
    names(gradient_pos) <- stringr::str_c("Radius_", gradient_pos, sep = "")
    #Calculate the results
    RESULTS <-
      purrr::map_dfc(gradient_pos, function(pos) {
        purrr::map_dbl(DATA_SPIAT, function(Image){
          Result <- try(mean(SPIAT::calculate_entropy(Image, cell_types_of_interest = Phenotypes_included,
                                                      radius = pos)[[13]]),#Select entropy cells of origin (first) and target (second)
                        silent = T)
          if(berryFunctions::is.error(Result)){
            return(NA)
          }
          else{
            return(Result)
          }
        })
      }, .progress = list(clear = F,
                          name = "Calculating entropy gradient",
                          show_after = 1,
                          type = "iterator"))

    RESULTS$Subject_Names <- names(DATA_SPIAT)
    RESULTS <- RESULTS[c(ncol(RESULTS), 1:(ncol(RESULTS)-1))]

    plot(RESULTS %>% pivot_longer(-1) %>%
           dplyr::mutate(name = factor(name, levels =stringr::str_c("Radius_", gradient_pos))) %>%
           ggplot(aes(x = name, y = value, color = Subject_Names, group = Subject_Names)) + geom_line() +
           cowplot::theme_cowplot() +
           scale_x_discrete("Radius", labels = as.character(gradient_pos)) +
           scale_y_continuous("Entropy")
    )
    return(RESULTS)
  }
