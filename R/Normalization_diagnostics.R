#' Compares normalized and non-normalized data
#'
#' The function will analyze normalization performance by comparing the Coefficient of Variation of Otsu thresholds across samples. In theory, if normalization is successful,
#' the coefficient of variation of the Otsu threshold should be reduced. In addition, the function will also generate smoothened histograms comparing feature expression distribution before and after normalization.
#'
#' @param Original_DATA A dataframe or tibble containing cell feature data before normalization.
#' @param Normalized_DATA A dataframe or tibble containing normalized feature expression data.
#'
#' @returns Returns summary plots.
#'
#' @seealso [Normalization_function()], [Normalization_function_parallel()]
#'
#' @examples
#' \dontrun{
#' #Run normalization----------------------
#' mxnorm_Parameters <-
#' list(
#'   slide_id = "Subject_Names",
#'   image_id = "Subject_Names",
#'   marker_cols = names(CSM_Arrangedcellfeaturedata_test[-c(1:4)]),
#'   transform = "None",
#'   method = "ComBat",
#'   method_override = NULL,
#'   method_override_name = NULL
#')
#'
#' Normalized_data <-  Normalization_function(
#'   CSM_Arrangedcellfeaturedata_test,
#'   Strategy = "mxnorm",
#'   Parameters = mxnorm_Parameters
#' )
#'
#' #Check normalization results---------------
#'Normalization_diagnostics(
#' Original_DATA = CSM_Arrangedcellfeaturedata_test,
#' Normalized_DATA = Normalized_data
#' )
#' }
#'
#'
#' @export

Normalization_diagnostics <-
  function(Original_DATA,
           Normalized_DATA){
    #Check suggested packages
    if(!requireNamespace("EBImage", quietly = TRUE)) stop(
      paste0("EBImage Bioconductor package is required to execute the function. Please install using the following code: ",
             expression({
               if (!require("BiocManager", quietly = TRUE))
                 install.packages("BiocManager")

               BiocManager::install("EBImage")
             })
      )
    )

    on.exit(gc())

    Original_DATA <- Original_DATA
    Normalized_DATA <- Normalized_DATA

    #Check arguments
    if(!identical(names(Original_DATA)[1:4], c("Cell_no", "X", "Y", "Subject_Names"))){
      stop("Original_DATA not correctly specified, please format appropiatetly (see Step 0)")
    }
    #Check arguments
    if(!identical(names(Normalized_DATA)[1:4], c("Cell_no", "X", "Y", "Subject_Names"))){
      stop("Normalized_DATA not correctly specified, please format appropiatetly (see Step 0)")
    }
    #Then check that both names are the same and in the same order
    if(!identical(names(Original_DATA), names(Normalized_DATA))) stop("Names in both datasets don't match. Please review")
    #Then check that both Subject_Names columns are identical
    if(!identical(Original_DATA$Subject_Names, Normalized_DATA$Subject_Names)) stop("Subject_Names columns are not identical in provided datasets")


    #Start by calculating the otsu index for each variable for original and normalized
    Thresholds_original <-dplyr::bind_cols(
      tibble(Subject_Names = unique(Original_DATA$Subject_Names), Data = "Original"),
      purrr::map_dfr(unique(Original_DATA$Subject_Names), function(Image){
        Interim <- Original_DATA %>% dplyr::filter(Subject_Names == Image)
        purrr::map_dbl(Interim[-c(1:4)], function(z){
          EBImage::otsu(array(z, dim = c(1, length(z))), range = c(min(z), max(z)), levels = length(unique(z)))
        })
      }, .progress = list(clear = F,
                          name = "Calculating Otsu thresholds by sample in original data",
                          show_after = 2,
                          type = "iterator"))
    )


    Thresholds_normalized <-dplyr::bind_cols(
      tibble(Subject_Names = unique(Normalized_DATA$Subject_Names), Data = "Normalized"),
      purrr::map_dfr(unique(Normalized_DATA$Subject_Names), function(Image){
        Interim <- Normalized_DATA %>% dplyr::filter(Subject_Names == Image)
        purrr::map_dbl(Interim[-c(1:4)], function(z){
          EBImage::otsu(array(z, dim = c(1, length(z))), range = c(min(z), max(z)), levels = length(unique(z)))
        })
      }, .progress = list(clear = F,
                          name = "Calculating Otsu thresholds by sample in normalied data",
                          show_after = 2,
                          type = "iterator"))
    )

    #Generate the plot of the Otsu thresholds by sample
    plot(
      dplyr::bind_rows(Thresholds_original, Thresholds_normalized) %>% tidyr::pivot_longer(cols = -c(1:2)) %>%
        ggplot(aes(x = factor(Data, levels= c("Original", "Normalized")), y = value, color = Data, fill = Data)) + facet_wrap(~name, scales = "free_y") +
        geom_boxplot(outlier.color=NA, alpha = 0.2, coef = 0) +
        geom_point(position=position_jitter(width = 0.1), alpha = 0.5) +
        guides(color = "none", fill = "none")+
        scale_x_discrete("") +
        scale_y_continuous("Otsu threshold by sample") +
        cowplot::theme_cowplot() +
        theme(axis.text.x = element_text(size = 10, color = "black"),
              strip.background =element_rect(fill="white"),
              strip.text = element_text(size = 12, color = "black", face = "bold"))
    )

    #Generate a tibble with the CV and plot them
    CV_tibble <-purrr::map_dfr(list(Thresholds_original,
                                    Thresholds_normalized),
                               function(Tibble){
                                 dplyr::bind_cols(unique(Tibble[2]),
                                                  purrr::map_df(Tibble[-c(1:2)], function(variable){
                                                    round((sd(variable) / mean(variable))*100, 2)
                                                  }))
                               })

    plot(
      CV_tibble %>%
        tidyr::pivot_longer(-Data) %>%
        ggplot(aes(x = name, y = value, fill = factor(Data, levels = c("Original", "Normalized")))) +
        geom_col(width = 0.5, color = "black", linewidth = 1.2, position = position_dodge(width = 0.6)) +
        scale_x_discrete("") +
        scale_y_continuous("CV of Otsu thresholds") +
        scale_fill_discrete("") +
        guides(fill = guide_legend(title = "")) +
        cowplot::theme_cowplot() +
        theme(axis.text = element_text(size = 12, color = "black"),
              legend.text = element_text(size = 12, color = "black"),
              legend.position = "bottom")
    )

    #Generate scaled versions of Original and normalized data in order to plot the histograms
    Histo_Original <-dplyr::bind_cols(tibble(Data = "Original"),
                                      Original_DATA[-c(1:4)] %>% scale())
    Histo_Normalized <-dplyr::bind_cols(tibble(Data = "Normalized"),
                                        Normalized_DATA[-c(1:4)] %>% scale())

    suppressWarnings(
      plot(
      ggplot() +
        geom_density(aes(y = value, x = -after_stat(density), fill = "Original"), color = "black", linewidth = 0.8, data = Histo_Original %>% tidyr::pivot_longer(-1)) +
        geom_density(aes(y = value, x = after_stat(density), fill = "Normalized"), color = "black", linewidth = 0.8, data = Histo_Normalized %>% tidyr::pivot_longer(-1))+
        facet_wrap(~name, scales = "free") +
        geom_vline(xintercept = 0, linewidth = 0.8, color = "black") +
        scale_x_continuous("", labels = NULL) +
        scale_y_continuous("",
                           limits = c(-3, 3), #Limit the Y scale to avoid over-informating outsiders
                           labels = NULL)+
        cowplot::theme_cowplot() +
        scale_fill_discrete("", breaks = c("Original", "Normalized"))+
        guides(fill = guide_legend(title = "")) +
        theme(legend.position = "bottom",
              strip.background =element_rect(fill="white"),
              strip.text = element_text(size = 12, color = "black", face = "bold"))
    )
    )
  }
