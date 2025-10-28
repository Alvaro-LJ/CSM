#' Thresholds cell features
#'
#' `Thresholding_function()` thresholds cell features following several approaches. The user can test this approaches using the [Thresholding_tester_app()].
#'  Thresholded features can be further processed using the [Marker_combinator_generator()] and [Phenotype_assigner_function()] to obtain cell phenotypes.
#'
#' @param DATA A dataframe or tibble containing cell feature data.
#' @param Strategy A character indicating the thresholding approach. It mus be one of the following: "EBI_Otsu", "Kmeans", "Kmeans_Otsu", "Autothreshold", "TriClass_Otsu", "Mean", "Quantile", "Arbitrary", "Multi_level" (see details).
#' @param Local_thresholding A logical value indicating if local (per image) threshold should be calculated (see details).
#' @param Method_autothreshold If strategy is "Autothreshold", indicate the type of threshold to be used. It must be one of "IJDefault", "Huang", "Huang2", "Intermodes", "IsoData", "Li", "MaxEntropy", "Mean",
#' "MinErrorI", "Minimum", "Moments", "Otsu", "RenyiEntropy", "Shanbhag", "Triangle", "Yen".
#' @param number_iterations_TriClass If strategy is "TriClass_Otsu", an integer value indicating the maximum number of iterations.
#' @param Percentile If strategy is "Quantile", a numeric value indicating the quantile to be used as cut-off value.
#' @param Defined_threshold If strategy is "Arbitrary", a numeric value indicating the cut-off value.
#' @param Levels If strategy is "Multi_level", an integer indicating the number of levels.
#'
#' @details
#' Local thresholding is generally discouraged. It may capture image-specific artifacts or patterns and not true feature positivity.
#'
#' EBI_Otsu calculates Otsu threshold using EBImage::otsu function. Kmeans calculates threshold using the mmand::threshold function. Kmeans_Otsu calculates the threhold using imager::threshold function.
#' Autothreshold calculates thresholds using the autothresholdr::auto_thresh_mask. TriClass_Otsu calculates thresholds using imagerExtra::ThresholdTriclass function.
#' Multi_level thresholds are calculated using imagerExtra::ThresholdML function.
#'
#' @seealso [Thresholding_tester_app()], [Thresholding_function_tailored()], [Thresholding_exploration_function()], [Thresholding_summary_function()], [Thresholding_graphicator_function()].
#'
#' @returns Returns a tibble with thresholded cell features. All features are thresholded using the same approach. For binary thresholding features are converted to logical vectors, for multi-level are converted to numeric vectors.
#'
#' @examples
#'\dontrun{
#' Thresholding_function(
#'   DATA = CSM_Arrangedcellfeaturedata_test,
#'   Strategy = "EBI_Otsu",
#'   Local_thresholding = TRUE,
#'   Method_autothreshold = "Otsu",
#'   number_iterations_TriClass = 20,
#'   Percentile = 0.5,
#'   Defined_threshold = 0.1,
#'   Levels = 3
#' )
#'}
#'
#' @export

Thresholding_function <-
  function(DATA,
           Strategy,
           Local_thresholding = FALSE,
           Method_autothreshold = "Otsu",
           number_iterations_TriClass = 20,
           Percentile = 0.5,
           Defined_threshold = 0.1,
           Levels = 3) {

    #Check suggested packages
    {
      if(Strategy == "EBI_Otsu"){
        if(!requireNamespace("EBImage", quietly = TRUE)) stop(
          paste0("EBImage Bioconductor package is required to execute the function. Please install using the following code: ",
                 expression({
                   if (!require("BiocManager", quietly = TRUE))
                     install.packages("BiocManager")

                   BiocManager::install("EBImage")
                 })
          )
        )
      }
      if(Strategy == "Kmeans"){
        if(!requireNamespace("mmand", quietly = FALSE)) stop(
          paste0("mmand CRAN package is required to execute the function. Please install using the following code: ",
                 expression(install.packages("mmand")))
        )
      }
      if(Strategy == "Kmeans_Otsu"){
        if(!requireNamespace("imager", quietly = FALSE)) stop(
          paste0("imager CRAN package is required to execute the function. Please install using the following code: ",
                 expression(install.packages("imager")))
        )
      }
      if(Strategy == "Autothreshold"){
        if(!requireNamespace("autothresholdr", quietly = FALSE)) stop(
          paste0("autothresholdr CRAN package is required to execute the function. Please install using the following code: ",
                 expression(install.packages("autothresholdr")))
        )
      }
      if(any(Strategy == "TriClass_Otsu", Strategy == "Multi_level")){
        if(!requireNamespace("imagerExtra", quietly = FALSE)) stop(
          paste0("imagerExtra CRAN package is required to execute the function. Please install using the following code: ",
                 expression(install.packages("imagerExtra")))
        )
      }
      if(Strategy == "Arbitrary"){
        if(!requireNamespace("mmand", quietly = FALSE)) stop(
          paste0("Arbitrary CRAN package is required to execute the function. Please install using the following code: ",
                 expression(install.packages("mmand")))
        )
      }
    }

    if(!(Strategy %in% c("EBI_Otsu", "Kmeans", "Kmeans_Otsu", "Autothreshold", "TriClass_Otsu", "Mean", "Quantile", "Arbitrary", "Multi_level"))) {
      stop("Strategy not correctly specified, please choose between EBI_Otsu, Kmeans, Kmeans_Otsu, Autothreshold, TriClass_Otsu, Mean, Quantile, Arbitrary, Multi_level")
    }

    else if(!identical(names(DATA)[1:4], c("Cell_no", "X", "Y", "Subject_Names"))) {
      stop("DATA not correctly specified, please format appropiatetly (see Step 0)")
    }
    if(!is.logical(Local_thresholding)) stop("Local_thresholding must be a logical value")

    else if(Local_thresholding == T) {
      warning("Local thresholding should be avoided. Take into account that global thresholding is the preferred method.")

      if(Strategy == "EBI_Otsu"){
        dplyr::bind_cols(DATA[(1:4)],
                         purrr::map_dfr(purrr::map(unique(DATA$Subject_Names), function(x) {
                           Tibble <- DATA %>% dplyr::filter(Subject_Names == x)
                           purrr::map2_df(.x = Tibble[-c(1:4)],
                                          .y =purrr::map_dbl(Tibble[-c(1:4)], function(z){
                                            EBImage::otsu(array(z, dim = c(1, length(z))), range = c(min(z), max(z)), levels = length(unique(z)))
                                          }),
                                          function(.x, .y) .x >= .y)
                         }),dplyr::bind_rows, .progress = list(clear = F,
                                                               name = "Calculating thresholds",
                                                               show_after = 2,
                                                               type = "iterator"))
        )
      }

      else if(Strategy == "Kmeans"){
        dplyr::bind_cols(DATA[(1:4)],
                         purrr::map_dfr(unique(DATA$Subject_Names), function(x){
                           Tibble <- DATA %>% dplyr::filter(Subject_Names == x)
                           purrr::map_dfc(Tibble[-c(1:4)], function(z){
                             as.logical(mmand::threshold(z, method = "kmeans", binarise = T))
                           })
                         }, .progress = list(clear = F,
                                             name = "Calculating thresholds",
                                             show_after = 2,
                                             type = "iterator"))
        )
      }

      else if(Strategy == "Kmeans_Otsu"){
        Interim_list <-purrr::map(unique(DATA$Subject_Names), function(Individual) DATA %>% dplyr::filter(Subject_Names == Individual) %>% dplyr::select(-c(1:4)))
        dplyr::bind_cols(DATA[(1:4)],
                         purrr::map_dfr(Interim_list,
                                        function(Individual){
                                          purrr::map_dfc(Individual, function(feature){
                                            if(length(unique(feature))>1){ #Check if vector has more than one unique value (if not no threshold can be calculated)
                                              if(!berryFunctions::is.error(as.vector(imager::threshold(imager::as.cimg(feature, dim = c(x = 1, y = length(feature), z = 1, cc = 1)),thr = "auto")))){
                                                #Check if thresholding by otsu K means yields an error anyway (due to histogram shape)
                                                as.vector(imager::threshold(
                                                  imager::as.cimg(feature, dim = c(x = 1, y = length(feature), z = 1, cc = 1)),
                                                  thr = "auto", approx = F))
                                              }
                                              else(NA)
                                            }
                                            else(NA)
                                          })
                                        }, .progress = list(clear = F,
                                                            name = "Calculating thresholds",
                                                            show_after = 2,
                                                            type = "iterator")))

      }

      else if(Strategy == "Autothreshold"){
        #Check autothreshold argument
        if(!Method_autothreshold %in% c("IJDefault", "Huang", "Huang2", "Intermodes", "IsoData", "Li", "MaxEntropy", "Mean",
                                        "MinErrorI", "Minimum", "Moments", "Otsu", "RenyiEntropy", "Shanbhag", "Triangle", "Yen")) {
          stop("Method_autothreshold must be one of the following: IJDefault, Huang, Huang2, Intermodes, IsoData, Li, MaxEntropy, Mean,
                                         MinErrorI, Minimum, Moments, Otsu, RenyiEntropy, Shanbhag, Triangle or Yen")
        }

        dplyr::bind_cols(DATA[(1:4)],

                         purrr::map_dfr(unique(DATA$Subject_Names), function(x) {

                           Tibble <- DATA %>% dplyr::filter(Subject_Names == x)

                           purrr::map_dfc(
                             Tibble[-c(1:4)], function(z) {
                               if(length(unique(z))>1){
                                 as.logical(autothresholdr::auto_thresh_mask(round((z/max(z))*1000, digits = 0),#Requires integer, just in case multiply by large number and round
                                                                             method = Method_autothreshold, #Specify the method
                                                                             ignore_black = F,
                                                                             ignore_na = F))
                               } else(NA)
                             })
                         }, .progress = list(clear = F,
                                             name = "Calculating thresholds",
                                             show_after = 2,
                                             type = "iterator"))
        )

      }

      else if(Strategy == "TriClass_Otsu"){
        #Check iteration argument
        if(number_iterations_TriClass%%1 != 0) stop("number_iterations_TriClass must be an integer value")

        dplyr::bind_cols(DATA[(1:4)],
                         purrr::map_dfr(purrr::map(unique(DATA$Subject_Names), function(x) {
                           Tibble <- DATA %>% dplyr::filter(Subject_Names == x)
                           purrr::map_df(Tibble[-c(1:4)],
                                         function(z){
                                           if(length(unique(z))>1){ #requires at least 2 levels
                                             as.logical(as.double(imagerExtra::ThresholdTriclass(imager::cimg(array(z, dim = c(1, length(z), 1, 1))),
                                                                                                 repeatnum = number_iterations_TriClass)))#How many times Otsu method should be refined
                                           }else(NA)
                                         })
                         }),dplyr::bind_rows, .progress = list(clear = F,
                                                               name = "Calculating thresholds",
                                                               show_after = 2,
                                                               type = "iterator"))
        )

      }

      else if(Strategy == "Mean"){
        Interim_list <-purrr::map(unique(DATA$Subject_Names), function(Individual) DATA %>% dplyr::filter(Subject_Names == Individual) %>% dplyr::select(-c(1:4)))
        dplyr::bind_cols(DATA[(1:4)],
                         purrr::map_dfr(Interim_list,
                                        function(Individual){
                                          purrr::map_dfc(Individual, function(z){
                                            Threshold <- mean(z)
                                            z >= Threshold
                                          })
                                        }, .progress = list(clear = F,
                                                            name = "Calculating thresholds",
                                                            show_after = 2,
                                                            type = "iterator"))
        )
      }

      else if(Strategy == "Quantile"){
        #Check percentile argument
        if(Percentile < 0.01 | Percentile > 0.99) stop("Percentile must be between 0.01 and 0.99")
        Interim_list <-purrr::map(unique(DATA$Subject_Names), function(Individual) DATA %>% dplyr::filter(Subject_Names == Individual) %>% dplyr::select(-c(1:4)))
        dplyr::bind_cols(DATA[(1:4)],
                         purrr::map_dfr(Interim_list,
                                        function(Individual){
                                          purrr::map_dfc(Individual, function(z){
                                            Threshold <- stats::quantile(z, Percentile)
                                            z >= Threshold
                                          })
                                        }, .progress = list(clear = F,
                                                            name = "Calculating thresholds",
                                                            show_after = 2,
                                                            type = "iterator"))
        )
      }

      else if(Strategy == "Arbitrary"){
        #Check argument
        if(!is.numeric(Defined_threshold)) stop("Defined_threshold must be a numeric value")

        dplyr::bind_cols(DATA[(1:4)],
                         purrr::map_dfr(unique(DATA$Subject_Names), function(x){
                           Tibble <- DATA %>% dplyr::filter(Subject_Names == x)
                           purrr::map_dfc(Tibble[-c(1:4)], function(z){
                             as.logical(mmand::threshold(z, level = Defined_threshold, method = "literal", binarise = T))
                           })
                         }, .progress = list(clear = F,
                                             name = "Calculating thresholds",
                                             show_after = 2,
                                             type = "iterator"))
        )
      }

      else if(Strategy == "Multi_level"){
        #Check level argument
        if(Levels%%1 != 0) stop("Levels must be an integer value")

        dplyr::bind_cols(DATA[(1:4)],
                         purrr::map_dfr(purrr::map(unique(DATA$Subject_Names), function(x) {
                           Tibble <- DATA %>% dplyr::filter(Subject_Names == x)
                           purrr::map_df(Tibble[-c(1:4)],
                                         function(z){
                                           if(length(unique(z))>Levels){ #requires at least n levels levels
                                             as.double(imagerExtra::ThresholdML(imager::cimg(array(z, dim = c(1, length(z), 1, 1))), k = (Levels-1))) #K to specify the amount of cut-off points
                                           }else(NA)
                                         })
                         }),dplyr::bind_rows, .progress = list(clear = F,
                                                               name = "Calculating thresholds",
                                                               show_after = 2,
                                                               type = "iterator"))
        )
      }

    }

    else{
      if(Strategy == "EBI_Otsu"){
        dplyr::bind_cols(DATA[(1:4)],
                         purrr::map2_df(.x = DATA[-c(1:4)],
                                        .y =purrr::map_dbl(DATA[-c(1:4)], function(z){
                                          EBImage::otsu(array(z, dim = c(1, length(z))), range = c(min(z), max(z)), levels = length(unique(z)))
                                        }),
                                        function(.x, .y) .x >= .y,
                                        .progress = list(clear = F,
                                                         name = "Calculating thresholds",
                                                         show_after = 2,
                                                         type = "iterator")))
      }

      else if(Strategy == "Kmeans"){
        dplyr::bind_cols(DATA[(1:4)],
                         purrr::map_dfc(DATA[-c(1:4)],
                                        function(z){
                                          as.logical(mmand::threshold(z, method = "kmeans", binarise = T))
                                        },
                                        .progress = list(clear = F,
                                                         name = "Calculating thresholds",
                                                         show_after = 2,
                                                         type = "iterator"))
        )
      }

      else if(Strategy == "Kmeans_Otsu"){
        dplyr::bind_cols(DATA[(1:4)],
                         purrr::map_df(DATA[-c(1:4)],
                                       function(feature){
                                         if(length(unique(feature))>1){ #Check if vector has more than one unique value (if not no threshold can be calculated)
                                           if(!berryFunctions::is.error(as.vector(imager::threshold(imager::as.cimg(feature, dim = c(x = 1, y = length(feature), z = 1, cc = 1)),thr = "auto")))){
                                             #Check if thresholding by otsu K means yields an error anyway (due to histogram shape)
                                             as.vector(imager::threshold(
                                               imager::as.cimg(feature, dim = c(x = 1, y = length(feature), z = 1, cc = 1)),
                                               thr = "auto", approx = F))
                                           }
                                           else(NA)
                                         }
                                         else(NA)
                                       },
                                       .progress = list(clear = F,
                                                        name = "Calculating thresholds",
                                                        show_after = 2,
                                                        type = "iterator")))

      }

      else if(Strategy == "Autothreshold"){
        #Check autothreshold argument
        if(!Method_autothreshold %in% c("IJDefault", "Huang", "Huang2", "Intermodes", "IsoData", "Li", "MaxEntropy", "Mean",
                                        "MinErrorI", "Minimum", "Moments", "Otsu", "RenyiEntropy", "Shanbhag", "Triangle", "Yen")) {
          stop("Method_autothreshold must be one of the following: IJDefault, Huang, Huang2, Intermodes, IsoData, Li, MaxEntropy, Mean,
                                         MinErrorI, Minimum, Moments, Otsu, RenyiEntropy, Shanbhag, Triangle or Yen")
        }

        dplyr::bind_cols(DATA[(1:4)],
                         purrr::map_df(purrr::map_df(DATA[-c(1:4)], function(x) {
                           round((x/max(x))*1000, digits = 0)
                         }),
                         function(z){
                           if(length(unique(z))>1){
                             as.logical(autothresholdr::auto_thresh_mask(z,#Requires integer, just in case multiply by large number and round
                                                                         method = Method_autothreshold, #Specify the method
                                                                         ignore_black = F,
                                                                         ignore_na = F))
                           } else(NA)
                         },
                         .progress = list(clear = F,
                                          name = "Calculating thresholds",
                                          show_after = 2,
                                          type = "iterator"))
        )

      }

      else if(Strategy == "TriClass_Otsu"){
        #Check iteration argument
        if(number_iterations_TriClass%%1 != 0) stop("number_iterations_TriClass must be an integer value")

        dplyr::bind_cols(DATA[(1:4)],
                         purrr::map_df(DATA[-c(1:4)],
                                       function(z){
                                         if(length(unique(z))>1){ #requires at least 3 levels
                                           as.logical(as.double(imagerExtra::ThresholdTriclass(imager::cimg(array(z, dim = c(1, length(z), 1, 1))),
                                                                                               repeatnum = number_iterations_TriClass)))#How many times Otsu method should be refined
                                         }else(NA)
                                       },
                                       .progress = list(clear = F,
                                                        name = "Calculating thresholds",
                                                        show_after = 2,
                                                        type = "iterator"))
        )

      }

      else if(Strategy == "Mean"){
        dplyr::bind_cols(DATA[(1:4)],
                         purrr::map_df(DATA[-c(1:4)],
                                       function(z){
                                         Threshold <- mean(z)
                                         z >= Threshold
                                       },
                                       .progress = list(clear = F,
                                                        name = "Calculating thresholds",
                                                        show_after = 2,
                                                        type = "iterator"))
        )
      }

      else if(Strategy == "Quantile"){
        #Check percentile argument
        if(Percentile < 0.01 | Percentile > 0.99) stop("Percentile must be between 0.01 and 0.99")

        dplyr::bind_cols(DATA[(1:4)],
                         purrr::map_df(DATA[-c(1:4)],
                                       function(z){
                                         Threshold <- quantile(z, Percentile)
                                         z >= Threshold
                                       },
                                       .progress = list(clear = F,
                                                        name = "Calculating thresholds",
                                                        show_after = 2,
                                                        type = "iterator"))
        )
      }

      else if(Strategy == "Arbitrary"){
        #Check argument
        if(!is.numeric(Defined_threshold)) stop("Defined_threshold must be a numeric value")

        dplyr::bind_cols(DATA[(1:4)],
                         purrr::map_dfc(DATA[-c(1:4)],
                                        function(z){
                                          as.logical(mmand::threshold(z, level = Defined_threshold, method = "literal", binarise = T))
                                        },
                                        .progress = list(clear = F,
                                                         name = "Calculating thresholds",
                                                         show_after = 2,
                                                         type = "iterator"))
        )
      }

      else if(Strategy == "Multi_level"){
        #Check level argument
        if(Levels%%1 != 0) stop("Levels must be an integer value")

        dplyr::bind_cols(DATA[(1:4)],
                         purrr::map_df(DATA[-c(1:4)],
                                       function(z){
                                         if(length(unique(z))>Levels){ #requires at least n Levels to be calculated
                                           as.double(imagerExtra::ThresholdML(imager::cimg(array(z, dim = c(1, length(z), 1, 1))), k = (Levels-1))) #K to specify the amount of cut-off points
                                         }else(NA)
                                       },
                                       .progress = list(clear = F,
                                                        name = "Calculating thresholds",
                                                        show_after = 2,
                                                        type = "iterator"))
        )
      }
    }


  }
