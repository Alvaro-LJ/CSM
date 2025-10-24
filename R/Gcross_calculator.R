#' Analyzes the spatial interaction between two cell types across images by using the Gcross function
#'
#' The function calculates the probability of spatial encounter between two cell types as distance sampling increases using the Gcross functions.
#' Steep rise in probability as distance increases means high chances of spatial interaction. The function is calculated using the spatstat R metapackage
#' (specifically the spatstat.explore::Gcross function). The Area Under the Curve (AUC) of the resulting function can be used as a metric of spatial association.
#'
#' @param DATA A dataframe or tibble containing a column named 'Phenotype' containing cell phenotype labels.
#' @param Cell_Of_Origin A character value indicating the cell phenotype label of the Cell of Origin.
#' @param Target_Cell A character value indicating the cell phenotype label of the Target cell.
#' @param Stop_at A numeric value indicating the maximum distance to be analyzed.
#' @param Sampling_frequency A numeric value indicating the distance sampling frequency.
#' @param Use_Clinical A logical value indicating if clinical data should be used in the plot generation process.
#' @param DATA_Clinical A dataframe or tibble containing image metadata.
#' @param Clinical_var A character value indicating the column name used in the clinical data analysis.
#'
#' @returns A list containing two outputs: Raw_DATA contains the expected interaction probability for every distance point. Simplified_DATA cotains the AUC data by image.
#'
#' @examples
#'Gcross_calculator(
#'    DATA = CSM_Phenotypecell_test,
#'    Cell_Of_Origin = "CD8_GZMBneg",
#'    Target_Cell = "TUMOR",
#'    Stop_at = 200,
#'    Sampling_frequency = 10,
#'    Use_Clinical = FALSE
#')
#'
#' @export

Gcross_calculator <-
  function(DATA = NULL,
           Cell_Of_Origin = NULL,
           Target_Cell = NULL,
           Stop_at = NULL, #Final distance where cumulative interaction analysis will stop
           Sampling_frequency = NULL, #Sampling distance interval
           Use_Clinical = NULL,
           DATA_Clinical = NULL,
           Clinical_var = NULL
  ) {

    #Check suggested packages
    {
      if(!requireNamespace("spatstat", quietly = FALSE)) stop(
        paste0("spatstat CRAN package is required to execute the function. Please install using the following code: ",
               expression(install.packages("spatstat")))
      )
    }

    #Argument check
    if(is.null(Cell_Of_Origin) || !is.character(Cell_Of_Origin)) {
      stop("Cell_Of_Origin must be a non-null character value.")
    }
    if(is.null(Target_Cell) || !is.character(Target_Cell)) {
      stop("Target_Cell must be a non-null character value.")
    }
    if(!all(c(Cell_Of_Origin, Target_Cell) %in% unique(DATA$Phenotype))){
      stop(paste0("Cell_Of_Origin and Target_Cell must be one of the following: ", stringr::str_c(unique(DATA$Phenotype), collapse = ", ")))
    }
    if(is.null(Stop_at) || !is.numeric(Stop_at) || Stop_at <= 0 || Stop_at < Sampling_frequency) {
      stop("Stop_at must be a positive numeric value larger than Sampling frequency.")
    }
    if(is.null(Sampling_frequency) || !is.numeric(Sampling_frequency) || Sampling_frequency <= 0) {
      stop("Sampling_frequency must be a positive numeric value.")
    }
    if(is.null(Use_Clinical) || !is.logical(Use_Clinical)) {
      stop("Use_Clinical must be a logical (TRUE/FALSE) value.")
    }
    if(Use_Clinical) {
      if(is.null(Clinical_var) || !is.character(Clinical_var)) {
        stop("Clinical_var must be a non-null character value when Use_Clinical is TRUE.")
      }
      # Check if Clinical_var exists in DATA_Clinical
      if(!all(Clinical_var %in% colnames(DATA_Clinical))) {
        stop("Clinical_var must be a valid column in DATA_Clinical.")
      }
    }

    #First obtain our data
    DATA_Phenotypes <- DATA

    #Calculate if any sample needs to be removed from the analysis
    Samples_to_remove <- purrr::map_lgl(unique(DATA_Phenotypes$Subject_Names), function(Image){
      nrow(DATA_Phenotypes %>% dplyr::filter(Subject_Names == Image, Phenotype %in% c(Cell_Of_Origin, Target_Cell)) %>% dplyr::count(Phenotype)) < 2
    })


    #If required, remove the troublesome samples and print a warning
    if(sum(Samples_to_remove) > 0){
      DATA_Phenotypes <- DATA_Phenotypes %>% dplyr::filter(!(Subject_Names %in% unique(DATA_Phenotypes$Subject_Names)[Samples_to_remove]))
      warning(paste0("Samples without COO or target cells will be removed from the analysis. ",
                     "The following samples will be removed: ", stringr::str_c(unique(DATA$Subject_Names)[Samples_to_remove], collapse = ", ")
      ))
    }

    #If not require continue with the analysis
    DATA_Phenotypes <- DATA_Phenotypes %>% dplyr::filter(Phenotype %in% c(Cell_Of_Origin, Target_Cell))


    #Generate our basic result tibble containing the Km calculation for each data
    RESULTS <-purrr::map_dfr(unique(DATA_Phenotypes$Subject_Names), function(x) {
      Interim <- DATA_Phenotypes %>% dplyr::filter(Subject_Names == x)  #select cell types

      #Prepare point pattern
      point_pattern <-
        spatstat.geom::ppp(x = Interim$X, y = Interim$Y,
                           n = nrow(Interim),
                           window = spatstat.geom::owin(xrange = c(min(Interim$X), max(Interim$X)),
                                                        yrange = c(min(Interim$Y), max(Interim$Y))),
                           marks = factor(Interim$Phenotype, levels = c(Cell_Of_Origin, Target_Cell))
        )

      #obtain the results
      RESULTS <- as_tibble(spatstat.explore::Gcross(point_pattern,
                                                    i = Cell_Of_Origin, #Select cell of origin
                                                    j = Target_Cell, #Select target cell
                                                    r = seq(from = 0, to = Stop_at, by = Sampling_frequency), #select radius to calculate the Gcross function
                                                    correction=c("km")
      )
      )

      RESULTS <- RESULTS %>%dplyr::mutate(Subject_Names = x) %>% dplyr::select(Subject_Names, r, theo, km)
    }, .progress = list(clear = F,
                        name = "Calculating G-cross function",
                        show_after = 1,
                        type = "iterator"))

    #If no clinical information provided then proceed with exit from the function
    if(!Use_Clinical) {
      try(plot(RESULTS %>% ggplot(aes(x = r, y = km, group = Subject_Names, color = Subject_Names)) + geom_line(linewidth = 1.2)+
                 cowplot::theme_cowplot() +
                 scale_x_continuous("Radius size") +
                 scale_y_continuous("Probability of encountering a Target Cell")+
                 scale_color_manual(values = unname(pals::polychrome(length(unique(DATA_Phenotypes$Subject_Names)))))
      ))

      AUC_results <-purrr::map_dfr(unique(RESULTS$Subject_Names), function(Image) {
        Interim <- RESULTS %>% dplyr::filter(Subject_Names == Image)

        AUC <-
          integrate(
            approxfun(Interim$r, Interim$km),
            lower = 0,
            upper = max(Interim$r),
            subdivisions = 2000
          )$value

        AUC_tibble <- tibble(
          Subject_Names = Image,
          AUC = AUC,
          COO = Cell_Of_Origin,
          Target_Cell = Target_Cell,
          Radius = Stop_at
        )
      })

      return(list(Raw_DATA = RESULTS,
                  Simplified_DATA = AUC_results))
    }

    #If clinical information provided then we can use the clinical information to  generate a graph
    else if(Use_Clinical) {
      print("Calculating association with clinical data")
      DATA_Clinical <- DATA_Clinical %>% dplyr::select(Subject_Names, all_of(Clinical_var))

      RESULTS <-dplyr::left_join(RESULTS, DATA_Clinical, by = "Subject_Names")
      names(RESULTS)[5] <- "Clin_Var"

      try(
        plot(RESULTS %>% ggplot(aes(x = r, y = km, group = Subject_Names, color = Clin_Var)) + geom_line(linewidth = 1.2)+
               cowplot::theme_cowplot() +
               scale_x_continuous("Radius size") +
               scale_y_continuous("Probability of encountering a Target Cell")+
               scale_color_viridis_d(Clinical_var))
      )

      names(RESULTS)[5] <- Clinical_var


      AUC_results <-purrr::map_dfr(unique(RESULTS$Subject_Names), function(Image) {
        Interim <- RESULTS %>% dplyr::filter(Subject_Names == Image)

        AUC <-
          integrate(
            approxfun(Interim$r, Interim$km),
            lower = 0,
            upper = max(Interim$r),
            subdivisions = 2000
          )$value

        AUC_tibble <- tibble(
          Subject_Names = Image,
          AUC = AUC,
          COO = Cell_Of_Origin,
          Target_Cell = Target_Cell,
          Radius = Stop_at,
          Clinical_var = unique(Interim[Clinical_var])[[1]]
        )
      })
      names(AUC_results)[6] <- Clinical_var
      return(list(Raw_DATA = RESULTS,
                  Simplified_DATA = AUC_results))

    }
  }
