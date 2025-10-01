#' Calculates cut-off points of thresholded data and generates summary plots
#'
#' Given the original and thresholded dataset, the function calculates cut-off values and generates summary graphs.
#'
#' @param DATA A dataframe or tibble containing cell feature data.
#' @param DATA_thresholded A dataframe or tibble containing cell feature data that has been thresholded.
#'
#' @returns Returns a tibble with cut-off values and generates summary plots.
#'
#' @export

Thresholding_summary_function <-
  function(DATA = NULL, DATA_thresholded = NULL){

    #Generate the two basic information sources for this analysis
    Thresholded <- DATA_thresholded[-c(1:4)]
    Markers <- DATA[-c(1:4)]

    #Generate two tibbles for dichotomic thresholds
    Thresholded_Dichom <- Thresholded[map_lgl(Thresholded, function(Var) length(unique(Var)) <= 2)]
    Markers_Dichom <- Markers[names(Thresholded_Dichom)]
    #If dichotomic variables are present calculate the thresholds
    if(ncol(Thresholded_Dichom) >= 1){
      #Find thresholds for dichotomic thresholds
      Thresholds_Dichom <-purrr::map2_dbl(.x = Markers_Dichom, .y = Thresholded_Dichom, function(.x, .y){
        min(.x[.y])
      },
      .progress = list(clear = F,
                       name = "Retrieving Dichotomic marker thresholds",
                       show_after = 2,
                       type = "iterator"))
    }

    #Generate two tibbles for Polychotomic thresholds
    Thresholded_Poly <- Thresholded[map_lgl(Thresholded, function(Var) length(unique(Var)) > 2)]
    Markers_Poly <- Markers[names(Thresholded_Poly)]
    #If polychotomic variables are present, calculate the thresholds
    if(ncol(Thresholded_Poly) >= 1){
      #Find thresholds for Polychotomic thresholds an turn it into a tibble
      Thresholds_Poly <-purrr::map2(.x = Markers_Poly, .y = Thresholded_Poly, function(.x, .y){
        purrr::map_dbl(unique(.y), function(Level) {
          min(.x[.y == Level])
        })[-1]
      },
      .progress = list(clear = F,
                       name = "Retrieving Polychotomic marker thresholds",
                       show_after = 2,
                       type = "iterator"))

      #The calculate the maximum amount of levels in the sample
      Max_level_size <- max(purrr::map_dbl(Thresholds_Poly, function(Individual) length(Individual)))
      #Generate a tibble with the levels
      Thresholds_Poly <-purrr::map_dfc(Thresholds_Poly, function(Individual) {
        c(Individual, rep(NA, times = Max_level_size - length(Individual)))
      })
    }

    #Generate the final tibble if only dichotomic variables are present
    if(ncol(Thresholded_Poly) < 1){
      #Generate the result tibble
      Thresholds_results <- matrix(Thresholds_Dichom, byrow = T, ncol = length(Thresholds_Dichom))
      colnames(Thresholds_results) <- names(Thresholds_Dichom)
      Thresholds_results <- as_tibble(Thresholds_results)
      Thresholds_results$Threshold_level <-stringr::str_c("Threshold_Level_", 1:nrow(Thresholds_results))
      Thresholds_results <- Thresholds_results[c(ncol(Thresholds_results), 1:(ncol(Thresholds_results)-1))]
      For_histogram <- Thresholds_results %>% pivot_longer(-1)
    }
    #Generate the final tibble if only polychotomic variables are present
    else if (ncol(Thresholded_Dichom) < 1){
      #Generate the result tibble
      Thresholds_results <- Thresholds_Poly
      Thresholds_results$Threshold_level <-stringr::str_c("Threshold_Level_", 1:nrow(Thresholds_results))
      Thresholds_results <- Thresholds_results[c(ncol(Thresholds_results), 1:(ncol(Thresholds_results)-1))]
      For_histogram <- Thresholds_results %>% pivot_longer(-1)
    }
    #Generate the final tibble if both type of variables are present
    else if (ncol(Thresholded_Poly) >= 1 & ncol(Thresholded_Dichom) >= 1){
      Thresholds_results <- matrix(Thresholds_Dichom, byrow = T, ncol = length(Thresholds_Dichom))
      colnames(Thresholds_results) <- names(Thresholds_Dichom)
      Thresholds_results <- as_tibble(Thresholds_results)
      Thresholds_results[2:Max_level_size, ] <- NA
      Thresholds_results <-dplyr::bind_cols(Thresholds_results, Thresholds_Poly)
      Thresholds_results$Threshold_level <-stringr::str_c("Threshold_Level_", 1:nrow(Thresholds_results))
      Thresholds_results <- Thresholds_results[c(ncol(Thresholds_results), 1:(ncol(Thresholds_results)-1))]
      For_histogram <- Thresholds_results %>% pivot_longer(-1)
    }

    #If not many markers present in the analysis, then print it in the plot window
    if(ncol(Markers) <= 10){
      plot(Markers %>% pivot_longer(1:ncol(Markers)) %>%
             ggplot(aes(x = value)) + facet_wrap(~name, "free", ncol = 1, nrow = ncol(Markers)) +
             geom_histogram(bins = 1000) +
             cowplot::theme_cowplot() +
             geom_vline(aes(xintercept = value), data = For_histogram, color = "red", linewidth = 1.2) +
             scale_x_continuous("Marker intensity") +
             scale_y_discrete("Number of cells")
      )
    }
    #Else directly generate a graph in the wd
    if(ncol(Markers) > 10) {
      warning(paste0("If the amount of Markers included is more than 10, the resulting graph will be exported to ", getwd()))
      png("Thresholds_for_Markers.png", width = 5000, height = 10000)
      plot(Markers %>% pivot_longer(1:ncol(Markers)) %>%
             ggplot(aes(x = value)) + facet_wrap(~name, "free", ncol = 1, nrow = ncol(Markers)) +
             geom_histogram(bins = 1000) +
             cowplot::theme_cowplot() +
             geom_vline(aes(xintercept = value), data = For_histogram, color = "red", linewidth = 1.2) +
             scale_x_continuous("Marker intensity") +
             scale_y_discrete("Number of cells"))
      dev.off()
    }

    return(Thresholds_results)

  }
