#' Calculates association between image feature data and image metadata
#'
#' `Clinical_Data_analyzer()` will calculate the association between image feature data (a numeric variable) and image metadata.
#' Image feature data can be obtained in various ways using CSM, for example running [Phenotype_quantifier()], [Global_heterogeneity_calculator()], [Distance_analyzer()] among others.
#'
#'
#' @param DATA A dataframe or tibble containing image feature data.
#' @param DATA_var A character value indicating which feature from DATA will be analyzed.
#' @param DATA_Clinical A dataframe or tibble containing image metadata, formatted using [Clinical_Data_arrange_function()].
#' @param Clinical_var A character value indicating the name of the column to be analyzed from image metadata.
#' @param Perform_time_to_event A logical value indicating if time to event analysis should be performed (see details).
#' @param Time_variable If time to event is required this argument is used to provide the column name specifying the time variable (Numeric).
#' @param Event_variable If time to event is required this argument is used to provide the column name specifying the event variable (0 or 1).
#' @returns Returns a summary of association and generates summary plots.
#'
#' @details If image metadata variable is a character, the function will calculate t-tests or ANOVA depending on the amount of categories.
#' If image metadata is numeric, Pearson correlation will be calculated. For time-to event analysis, cut-off points for the image feature will be calculated using the survminer::surv_cutpoint function.
#'
#' @examples
#' \dontrun{
#' #Calculate the area of TMA samples-------------------------------------------
#' DATA_AREA <-
#'  Image_size_calculator(
#'     DATA = CSM_PhenotypeTMA_test,
#'     Strategy = "Concave_hull",
#'     Hull_ratio = 0.4
#')
#'
#' #Calculate the cell densities by samples-------------------------------------
#' Cells_by_sample <-
#'  Phenotype_quantifier(
#'     CSM_PhenotypeTMA_test,
#'     Calculate_Density = TRUE,
#'     DATA_Area = DATA_AREA
#' )
#'
#' #Obtain clinical data--------------------------------------------------------
#' DATA_CLINICAL <-
#' Clinical_Data_arrange_function(
#'      DATA = CSM_ClinicalTMA_test,
#'      Subject_Names = "Sample",
#'      Outcomes_to_keep = c("AGE", "MMRP_status", "DEATH", "OS_m")
#' )
#'
#' #Analyze sample features-----------------------------------------------------
#'Clinical_Data_analyzer(
#'    DATA = Cells_by_sample,
#'    DATA_var = c("Density_CD8", "Density_MACROPHAGE"),
#'    DATA_Clinical = DATA_CLINICAL,
#'    Clinical_var = "MMRP_status",
#'    Perform_time_to_event = FALSE
#')
#'
#'Clinical_Data_analyzer(
#'    DATA = Cells_by_sample,
#'    DATA_var = c("Density_CD8", "Density_MACROPHAGE"),
#'    DATA_Clinical = DATA_CLINICAL,
#'    Clinical_var = "AGE",
#'    Perform_time_to_event = FALSE
#')
#'
#'Clinical_Data_analyzer(
#'    DATA = Cells_by_sample,
#'    DATA_var = c("Density_CD8", "Density_MACROPHAGE"),
#'    DATA_Clinical = DATA_CLINICAL,
#'    Perform_time_to_event = TRUE,
#'    Time_variable = "OS_m",
#'    Event_variable = "DEATH"
#')
#' }
#'
#' @export

Clinical_Data_analyzer <-
  function(DATA,
           DATA_var,
           DATA_Clinical,
           Clinical_var,
           Perform_time_to_event = FALSE,
           Time_variable = NULL,
           Event_variable = NULL) {

    #Check suggested packages
    {if(Perform_time_to_event){
      if(!requireNamespace("survival", quietly = FALSE)) stop(
        paste0("survival CRAN package is required to execute the function. Please install using the following code: ",
               expression(install.packages("survival")))
      )
      if(!requireNamespace("survminer", quietly = FALSE)) stop(
        paste0("survminer CRAN package is required to execute the function. Please install using the following code: ",
               expression(install.packages("survminer")))
      )
    }
    }

    #Check if Subject names is in both Data sources
    if(!all("Subject_Names" %in% names(DATA), "Subject_Names" %in% names(DATA_CLINICAL))) {
      stop("Subject_Names variable should be present in DATA and DATA_Clinical")
    }
    #Check if at least one Subject Name is present in both Data sources
    if(!any(DATA_Clinical$Subject_Names %in% DATA$Subject_Names)) {
      stop("No match between samples from DATA and DATA_Clinical. Check that at least a single subject is present in both data sources")
    }
    #Check DATA_variables included in analysis
    if(!all(DATA_var %in% names(DATA))) {
      Missing_arguments <- DATA_var[!DATA_var %in% names(DATA)]
      stop(paste0(stringr::str_c(Missing_arguments, collapse = ", "), " not found in DATA"))
    }
    if(!is.logical(Perform_time_to_event)){
      stop("Perform_time_to_event should be a logical value")
    }

    #Proceed with regular analysis (NO TIME TO EVENT)
    else {
      if(!Perform_time_to_event) {
        #Check Clinical Variables included in Data Clinical
        if(!all(Clinical_var %in% names(DATA_Clinical))) {
          Missing_arguments <- Clinical_var[!Clinical_var %in% names(DATA_Clinical)]
          stop(paste0(stringr::str_c(Missing_arguments, collapse = ", "), " not found in DATA_Clinical"))
        }
        else {
          #Import Data to be correlated with clinical findings
          DATA <- DATA
          DATA <- DATA %>% dplyr::select(Subject_Names, all_of(DATA_var))

          #Import and arrange clinical Data
          DATA_Clinical <- DATA_Clinical
          DATA_Clinical <- DATA_Clinical %>% dplyr::select(Subject_Names, all_of(Clinical_var))
          names(DATA_Clinical)[2] <- "Clin_var"

          #Join the clinical and the actual data
          JOINED <- dplyr::left_join(DATA_Clinical, DATA, by = "Subject_Names")

          #If clinical data is a character or a factor
          if(any(is.character(JOINED$Clin_var), is.factor(JOINED$Clin_var))) {
            #Rearrange Subject_Names according to clinical variable
            JOINED <- JOINED %>% dplyr::arrange(Clin_var)
            JOINED <- JOINED %>% dplyr::mutate(Subject_Names = factor(Subject_Names, levels = JOINED$Subject_Names))

            #Plot samples by clinical variable
            plot(JOINED %>% tidyr::pivot_longer(-c(1:2)) %>%
                   ggplot(aes(x = Subject_Names, fill = Clin_var, y = value)) +
                   facet_wrap(~name, "free_y", ncol = 1, nrow = ncol(JOINED)-2) + geom_col(width = 0.5, color = "black") +
                   cowplot::theme_cowplot() +
                   scale_x_discrete("") +
                   scale_y_continuous("value")+
                   scale_fill_viridis_d(Clinical_var) +
                   theme(axis.text.x = element_text(angle = -90, hjust = 0.5))
            )

            #Plot summary
            plot(JOINED %>% tidyr::pivot_longer(-c(1:2)) %>%
                   ggplot(aes(x = Clin_var, fill = Clin_var, y = value)) +
                   facet_wrap(~name, "free_y", ncol = 1, nrow = ncol(JOINED)-2) + geom_bar(width = 0.5, color = "black", stat = "summary", fun = "mean") +
                   cowplot::theme_cowplot() +
                   scale_x_discrete(Clinical_var) +
                   scale_y_continuous("Mean value") +
                   scale_fill_viridis_d(Clinical_var) +
                   ggpubr::stat_compare_means(label.x.npc = "center", label.y.npc = 'bottom')
            )

            #Prepare the summary for everyone of the DATA vars
            RESULTS <- purrr::map(3:ncol(JOINED), function(Variable){
              JOINED <-dplyr::bind_cols(JOINED[1:2], JOINED[Variable])
              #Calculate mean
              MEAN_tibble <- JOINED %>% group_by(Clin_var) %>% summarize_at(vars(-group_cols(), -Subject_Names), .funs = "mean", na.rm = TRUE) %>% ungroup()
              names(MEAN_tibble)[2] <- stringr::str_c("Average_", names(MEAN_tibble)[2])
              #Calculate sd
              SD_tibble <- JOINED %>% group_by(Clin_var) %>% summarize_at(vars(-group_cols(), -Subject_Names), .funs = "sd", na.rm = TRUE) %>% ungroup()
              names(SD_tibble)[2] <- stringr::str_c("SD_", names(SD_tibble)[2])
              #calculate p25
              p25_tibble <- JOINED %>% group_by(Clin_var) %>% summarize_at(vars(-group_cols(), -Subject_Names), .funs = function(x) quantile(x, 0.25, na.rm = TRUE)) %>% ungroup()
              names(p25_tibble)[2] <- stringr::str_c("p25_", names(p25_tibble)[2])
              #Calculate p50
              Median_tibble <- JOINED %>% group_by(Clin_var) %>% summarize_at(vars(-group_cols(), -Subject_Names), .funs = function(x) quantile(x, 0.5, na.rm = TRUE)) %>% ungroup()
              names(Median_tibble)[2] <- stringr::str_c("Median_", names(Median_tibble)[2])
              #Calculate p75
              p75_tibble <- JOINED %>% group_by(Clin_var) %>% summarize_at(vars(-group_cols(), -Subject_Names), .funs = function(x) quantile(x, 0.75, na.rm = TRUE)) %>% ungroup()
              names(p75_tibble)[2] <- stringr::str_c("p75_", names(p75_tibble)[2])
              #Bind the resulting tibbles
              RESULTS_tibble <-dplyr::bind_cols(MEAN_tibble, SD_tibble[-1], p25_tibble[-1], Median_tibble[-1], p75_tibble[-1])

              #Generate the t.test or the ANOVA test according to the number of levels of the clinical variable
              if(length(unique(JOINED$Clin_var)) == 2) {
                if(berryFunctions::is.error(stats::t.test(JOINED[[3]]~JOINED$Clin_var))) {t_test_results <- NA}
                else{t_test_results <- stats::t.test(JOINED[[3]]~JOINED$Clin_var)$p.value}
                return(list(Summary = RESULTS_tibble,
                            t_test_result = t_test_results))
              }
              else {
                if(berryFunctions::is.error(stats::oneway.test(JOINED[[3]]~JOINED$Clin_var))) {ANOVA_test_results <- NA}
                else{ANOVA_test_results <- stats::oneway.test(JOINED[[3]]~JOINED$Clin_var)$p.value}
                return(list(Summary = RESULTS_tibble,
                            ANOVA_result = ANOVA_test_results))
              }

            })

            #Change the name of the RESULTS list
            names(RESULTS) <- names(JOINED)[3:ncol(JOINED)]
            return(RESULTS)

          }

          #If its numeric then execute the following code
          if(is.numeric(JOINED$Clin_var)) {
            #Generate the pearson correlation coefficient and the associated pvalue
            RESULTS <-
              purrr::map_dfr(JOINED[-(1:2)], function(Variable) {
                Correlation_Result <- stats::cor.test(Variable, JOINED$Clin_var, method = "pearson", alternative = "two.sided")
                c(Pearson_cor = Correlation_Result$estimate,
                  Pearson_pval = Correlation_Result$p.value)
              })
            #Add the names of the pearson
            RESULTS$name <- names(JOINED)[-(1:2)]

            #Change the names and the order of the variables
            names(RESULTS) <- c("Pearson_cor", "Pearson_pval", "name")
            RESULTS <- RESULTS[c(3, 1, 2)]

            #Plot the individual variables and outcomes
            plot(
              JOINED %>% tidyr::pivot_longer(-c(1:2)) %>%
                ggplot(aes(x = Clin_var, y = value)) +
                facet_wrap(~name, "free_y", ncol = 1, nrow = ncol(JOINED)-2) +
                geom_point(size = 1.2) + geom_smooth(method = "lm", color = "red", se = F) +
                cowplot::theme_cowplot() +
                scale_x_continuous(Clinical_var) +
                scale_y_continuous("Variable")
            )

            #Plot the summary
            plot(
              RESULTS %>% ggplot(aes(x = forcats::fct_reorder(name, Pearson_cor), y = Pearson_cor, fill = Pearson_pval)) +
                geom_col(width = 0.5, color = "black") +
                cowplot::theme_cowplot() +
                scale_x_discrete("") +
                scale_y_continuous("Pearson Correlation") +
                geom_hline(yintercept = 0, color ="black", linewidth = 0.9) +
                scale_fill_gradient2(high = "white", low = "red", mid = "white", midpoint = 0.05, limits = c(0, 1))+
                theme(axis.text.x = element_text(angle = -90, vjust = 0.5))
            )

            #Return the results
            return(RESULTS)
          }
        }
      }

      #Proceed with time to event analysis if required
      else if(Perform_time_to_event) {
        #Check time to event and event variables
        if(!all(c(Time_variable, Event_variable) %in% names(DATA_Clinical))) {
          Missing_arguments <- c(Time_variable, Event_variable)[!c(Time_variable, Event_variable) %in% names(DATA_Clinical)]
          stop(paste0(stringr::str_c(Missing_arguments, collapse = ", "), " not found in DATA_Clinical"))
        }

        #Import Data to be correlated with clinical findings
        DATA <- DATA
        DATA <- DATA %>% dplyr::select(Subject_Names, all_of(DATA_var))

        #Import and arrange clinical Data
        DATA_Clinical <- DATA_Clinical
        DATA_Clinical <- DATA_Clinical %>% dplyr::select(Subject_Names, dplyr::all_of(c(Time_variable, Event_variable)))
        names(DATA_Clinical)[2:3] <- c("Time_variable", "Event_variable")

        #Join the clinical and the actual data
        JOINED <-dplyr::left_join(DATA_Clinical, DATA, by = "Subject_Names")

        res.cut <- survminer::surv_cutpoint(JOINED,
                                            time = "Time_variable",
                                            event = "Event_variable",
                                            variables = names(JOINED)[-c(1:3)],
                                            minprop = 0.1,
                                            progressbar=TRUE
        )
        JOINED_Cat <- as_tibble(survminer::surv_categorize(res.cut))

        Survival_models <-
          purrr::map(3:ncol(JOINED_Cat), function(Var) {
            Interim <-dplyr::bind_cols(JOINED_Cat[1:2], JOINED_Cat[Var])
            names(Interim)[3] <- "Target"
            Model <- survival::survfit(survival::Surv(as.numeric(Time_variable), as.numeric(Event_variable)) ~ Target, data = Interim)
            Plot <- survminer::ggsurvplot(Model, data = Interim,
                                          pval = T, risk.table = F, conf.int = FALSE,
                                          linewidth = 1,
                                          palette = "npg",
                                          title = names(JOINED_Cat)[Var],
                                          xlab = Time_variable,
                                          ylab = Event_variable)$plot
            CoxPHModel <- survival::coxph(survival::Surv(as.numeric(Time_variable), as.numeric(Event_variable)) ~ Target, data = Interim)

            return(list(Model = CoxPHModel,
                        Plot = Plot))
          })
        #Plot the KM plots
        plot(cowplot::plot_grid(plotlist =purrr::map(Survival_models, ~.[["Plot"]]), ncol = ceiling( (ncol(JOINED_Cat)-3) / 2)))

        #Return a list with the model
        Model_list <- purrr::map(Survival_models, ~summary(.[["Model"]]))
        names(Model_list) <- names(JOINED_Cat[3:ncol(JOINED_Cat)]) #Add names
        return(Model_list)
      }
    }
  }
