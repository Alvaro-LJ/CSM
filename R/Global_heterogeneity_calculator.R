#' Calculates cell composition heterogeneity.
#'
#' The function calculates various heterogeneity metrics based on cell composition distribution. It does not take into account spatial information.
#'
#' @param DATA A dataframe or tibble containing a column named 'Phenotype' containing cell phenotype labels.
#' @param Phenotypes_included A character vector indicating the phenotype labels that will be included in the analysis.
#'
#' @details
#' Shannon, Simpson,Inverse Simpson and Renyi entropy indexes are calculated using the vegan package.
#'
#' Rao index is calculated using the picante::raoD function.
#'
#' Gini index is calculated using the DescTools::Gini function.
#'
#' Kullback_Leibler and Jensen_Shannon indexes are calculated using the philentropy package. They are computed in two different ways:
#' First '_Sample' indexes compare the observed cell composition against an homogeneous distribution where all cells have the same prevalence in the sample.
#' Second '_Experiment' indexes compare the observed cell composition against the cell distribution observed in all the cells from DATA.
#'
#'
#' @returns Generates a tibble containing heterogeneity metrics by image.
#'
#' @examples
#' \dontrun{
#' Global_heterogeneity_calculator(
#'     DATA = CSM_Phenotypecell_test,
#'     Phenotypes_included = unique(CSM_Phenotypecell_test$Phenotype)
#' )
#' }
#'
#' @export

Global_heterogeneity_calculator <-
  function(DATA,
           Phenotypes_included){
    #Check suggested packages
    {
      if(!requireNamespace("picante", quietly = FALSE)) stop(
        paste0("picante CRAN package is required to execute the function. Please install using the following code: ",
               expression(install.packages("picante")))
      )
      if(!requireNamespace("vegan", quietly = FALSE)) stop(
        paste0("vegan CRAN package is required to execute the function. Please install using the following code: ",
               expression(install.packages("vegan")))
      )
      if(!requireNamespace("DescTools", quietly = FALSE)) stop(
        paste0("DescTools CRAN package is required to execute the function. Please install using the following code: ",
               expression(install.packages("DescTools")))
      )
      if(!requireNamespace("philentropy", quietly = FALSE)) stop(
        paste0("philentropy CRAN package is required to execute the function. Please install using the following code: ",
               expression(install.packages("philentropy")))
      )
    }

    #Check if the phenotypes included are present in DATA
    if(!all(Phenotypes_included %in% unique(DATA$Phenotype))) {
      stop(paste0("Phenotypes included must be any of: ", stringr::str_c(unique(DATA$Phenotype), collapse = ", ")))
    }
    #If everything is OK perform the analysis
    else{

      #Generate the cell count data
      Interim <- DATA %>% dplyr::filter(Phenotype %in% Phenotypes_included) %>%
        dplyr::group_by(Subject_Names, Phenotype) %>% dplyr::count() %>% dplyr::ungroup() %>%
        tidyr::pivot_wider(id_cols = Subject_Names, names_from = Phenotype, values_from = n)
      #Substitute NA for 0
      Interim[is.na(Interim)] <- 0

      #Arrange the Interim matrix according to column name
      Interim <-dplyr::bind_cols(Interim[1], Interim[sort(names(Interim[-1]))])

      #Calculate experiment wise cell proportion distribution (for KL and JS)
      Global_counts <- DATA %>% dplyr::filter(Phenotype %in% Phenotypes_included) %>%
        dplyr::count(Phenotype) %>% dplyr::mutate(n = n/sum(n)) %>% tidyr::pivot_wider(names_from = Phenotype, values_from = n)
      #Sort Global_counts according to column name
      Global_counts <- Global_counts[sort(names(Global_counts))]
      Global_counts <- unlist(Global_counts)

      #Calculate the different metrics
      Results <-
        dplyr::bind_cols(Interim[1],
                         tibble(Shannon = apply(Interim[-1], MARGIN = 1, function(row) vegan::diversity(row, index = "shannon")),
                                Simpson = apply(Interim[-1], MARGIN = 1, function(row) vegan::diversity(row, index = "simpson")),
                                Inverse_simpson = apply(Interim[-1], MARGIN = 1, function(row) vegan::diversity(row, index = "invsimpson"))
                         ),
                         tibble(renyi = as.double(vegan::renyi(Interim[-1], scales = Inf))),
                         tibble(Rao_Dkk = picante::raoD(Interim[-1])$Dkk),
                         tibble(Gini = unlist(apply(Interim[-1], MARGIN = 1, function(Sample) DescTools::Gini(Sample, conf.level = NA)))),
                         tibble(Kullback_Leibler_Sample = unlist(apply(Interim[-1], MARGIN = 1, function(Row){
                           Observed <- Row/sum(Row)
                           Expected <- rep(1/length(Row), times = length(Row))
                           suppressMessages(philentropy::KL(rbind(Observed, Expected), unit = "log2"))
                         }))),
                         tibble(Jensen_Shannon_Sample = unlist(apply(Interim[-1], MARGIN = 1, function(Row){
                           Observed <- Row/sum(Row)
                           Expected <- rep(1/length(Row), times = length(Row))
                           suppressMessages(philentropy::JSD(rbind(Observed, Expected),  test.na = FALSE, unit = "log2"))
                         }))),
                         tibble(Kullback_Leibler_Experiment = unlist(apply(Interim[-1], MARGIN = 1, function(Row){
                           Observed <- Row/sum(Row)
                           Expected <- Global_counts
                           suppressMessages(philentropy::KL(rbind(Observed, Expected), unit = "log2"))
                         }))),
                         tibble(Jensen_Shannon_Experiment = unlist(apply(Interim[-1], MARGIN = 1, function(Row){
                           Observed <- Row/sum(Row)
                           Expected <- Global_counts
                           suppressMessages(philentropy::JSD(rbind(Observed, Expected),  test.na = FALSE, unit = "log2"))
                         })))
        )
      names(Results)[-1] <- c("Shannon", "Simpson", "Inverse_simpson", "Renyi_Scale_Inf", "Rao_Dkk", "Gini", "KL_Sample", "JS_Sample", "KL_Experiment", "JS_Experiment")
      #Return the result
      return(Results)
    }
  }
