#' Normalizes cell feature expression data
#'
#' Normalization of cell feature data following two potential strategies included in mxnorm and simpleSeg packages.
#'
#' @param DATA A dataframe or tibble containing cell feature data.
#' @param Strategy 'mxnorm' or 'simpleSeg'.
#' @param Parameters List of parameters to be used in the normalization process (see details).
#'
#' @returns Returns a tibble containing normalized cell feature data.
#'
#' @details For mxnorm normalization, the list must contain the following elements:
#'  slide_id (The column name of the DATA that contains the slide ID. this can be the Subject_Names or a user supplied vector containing the slide identity of these cells),
#'  image_id (The column name of the DATA that contains the image ID. Can be the same as slide_id),
#'  marker_cols (The name of the columns to be normalized),
#'  transform ("None", "log10", "mean_divide","log10_mean_divide"),
#'  method ("None", "ComBat","Registration"),
#'  method_override (OPTIONAL) ("specify a custom normalization method. See mxnorm package documentation for more details),
#'  method_override_name (OPTIONAL) (specify the name of the user defined method).
#'
#'  For simpleSeg approach, the list must contain the following elements:
#'  cells (The tibble containing the cell feature data),
#'  markers (The name of the columns to be normalized),
#'  imageID (The column name of the DATA that contains the image ID),
#'  transformation (Choose between:  NULL, "asinh", "sqrt"),
#'  method ("trim99", "minMax"),
#'  cores (Irrelevant argument. If parallelization is required use N_cores argument and [Normalization_function_parallel()] function).
#'
#' @seealso [Normalization_function_parallel()], [Normalization_diagnostics()]
#'
#' @examples
#' \dontrun{
#' #Set normalization parameters-----------
#'mxnorm_Parameters <-
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
#'mxnorm_Parameters_slide <-
#' list(
#'   slide_id = c(rep("Slide1", times = 229), rep("Slide2", times = 818)),
#'   image_id = "Subject_Names", #The column name of the DATA that contains the image ID
#'   marker_cols = names(CSM_Arrangedcellfeaturedata_test[-c(1:4)]),
#'   transform = "None",
#'   method = "ComBat",
#'   method_override = NULL,
#'   method_override_name = NULL
#')
#'
#'simpleSeg_Parameters <-
#' list(
#'   cells = CSM_Arrangedcellfeaturedata_test,
#'   markers = names(CSM_Arrangedcellfeaturedata_test[-c(1:4)]),
#'   imageID = "Subject_Names",
#'   transformation = "sqrt",
#'   method = c("trim99", "minMax"),
#'   cores = 1
#')
#'
#'#Run normalization----------------------
#'Normalization_function(
#'   CSM_Arrangedcellfeaturedata_test,
#'   Strategy = "mxnorm",
#'   Parameters = mxnorm_Parameters
#' )
#'
#' Normalization_function(
#'   CSM_Arrangedcellfeaturedata_test,
#'   Strategy = "simpleSeg",
#'   Parameters = simpleSeg_Parameters
#' )
#' }
#'
#' @export

Normalization_function <-
  function(DATA = NULL,
           Strategy = NULL,
           Parameters = NULL) {

    #Check suggested packages
    {
      if(Strategy == "mxnorm"){
        if(!requireNamespace("mxnorm", quietly = FALSE)) stop(
          paste0("mxnorm CRAN package is required to execute the function. Please install using the following code: ",
                 expression(install.packages("mxnorm")))
        )
      }
      if(Strategy == "simpleSeg"){
        if(!requireNamespace("simpleSeg", quietly = TRUE)) stop(
          paste0("simpleSeg Bioconductor package is required to execute the function. Please install using the following code: ",
                 expression({
                   if (!require("BiocManager", quietly = TRUE))
                     install.packages("BiocManager")

                   BiocManager::install("simpleSeg")
                 })
          )
        )
      }
    }

    #Check strategy
    if(!Strategy %in% c("mxnorm", "simpleSeg")){
      stop("Strategy must be one of the following: mxnorm, simpleSeg")
    }

    #MXNORM CODE
    if(Strategy == "mxnorm") {
      #check that argument list is correct and contains adequate data
      if(!all(c("slide_id", "image_id", "marker_cols", "transform", "method", "method_override", "method_override_name") %in% names(Parameters))){
        Missing_arguments <- c("slide_id", "image_id", "marker_cols", "transform", "method", "method_override", "method_override_name")[
          !c("slide_id", "image_id", "marker_cols", "transform", "method", "method_override", "method_override_name") %in% names(Parameters)
        ]
        stop(paste0(stringr::str_c(Missing_arguments, collapse = ", "), " not present in Paramters list"))
      }
      #check if arguments are well specified
      Argument_checker <- c(data_OK = identical(names(DATA)[1:4], c("Cell_no", "X", "Y", "Subject_Names")),
                            image_id_OK = all(is.character(Parameters[["image_id"]]), Parameters[["image_id"]] %in% names(DATA)),
                            slide_id_OK = any(identical(Parameters[["slide_id"]], Parameters[["image_id"]]),
                                              all(is.vector(Parameters[["slide_id"]], mode = "character"),
                                                  length(Parameters[["slide_id"]]) == nrow(DATA))
                            ),
                            marker_cols_OK = all(Parameters[["marker_cols"]] %in% names(DATA)),
                            transform_OK = Parameters[["transform"]] %in% c("None", "log10", "mean_divide","log10_mean_divide"),
                            method_OK = Parameters[["method"]] %in% c("None", "ComBat","Registration")
      )
      Stop_messages <- c(data_OK = "DATA must be adequately arrange (using Data_arrange_function from step 0)",
                         image_id_OK = stringr::str_c(Parameters[["image_id"]], " must be present in DATA"),
                         slide_id_OK = "If slide_id is not identical to image_id, it must be a vector indicating the slide origin of every cell",
                         marker_cols_OK = "All marker_cols must be present in the DATA provided",
                         transform_OK = "Transform method must be one of the following: None, log10, mean_divide, log10_mean_divide",
                         method_OK = "Method must be one of the following: None, ComBat, Registration")
      #Check arguments and stop if necessary
      if(!all(Argument_checker)){
        stop(cat(Stop_messages[!Argument_checker],
                 fill = sum(!Argument_checker)))
      }

      #Now we generate the argument list again with the ones provided
      if(identical(Parameters[["slide_id"]], Parameters[["image_id"]])){
        Parameters <- list(
          mx_data = mxnorm::mx_dataset(data = DATA,
                                       slide_id = Parameters[["slide_id"]],
                                       image_id = Parameters[["image_id"]],
                                       marker_cols = Parameters[["marker_cols"]]
          ),
          transform = Parameters[["transform"]],
          method = Parameters[["method"]],
          method_override = Parameters[["method_override"]],
          method_override_name = Parameters[["method_override_name"]]
        )
      }
      #If slide_id is not the same as image_id generate another parameter list
      if(!identical(Parameters[["slide_id"]], Parameters[["image_id"]])) {
        DATA <-dplyr::bind_cols(DATA,
                                tibble(slide_id = Parameters[["slide_id"]]))
        Parameters <- list(
          mx_data = mxnorm::mx_dataset(data = DATA,
                                       slide_id = "slide_id",
                                       image_id = Parameters[["image_id"]],
                                       marker_cols = Parameters[["marker_cols"]]
          ),
          transform = Parameters[["transform"]],
          method = Parameters[["method"]],
          method_override = Parameters[["method_override"]],
          method_override_name = Parameters[["method_override_name"]]
        )
      }

      #Execute normalization
      DATA_normalized <- base::do.call(mxnorm::mx_normalize, args = Parameters)
      DATA_normalized <- as_tibble(DATA_normalized$norm_data)[-c(1:2)]
      #Append normalization results to DATA
      return(bind_cols(DATA[1:4], DATA_normalized))
    }

    #simpleSeg CODE
    else if(Strategy == "simpleSeg") {
      #check that argument list is correct and contains adequate data
      if(!all(c("cells", "markers", "imageID", "transformation", "method", "cores") %in% names(Parameters))){
        Missing_arguments <- c("cells", "markers", "imageID", "transformation", "method", "cores")[
          !c("cells", "markers", "imageID", "transformation", "method", "cores") %in% names(Parameters)
        ]
        stop(paste0(stringr::str_c(Missing_arguments, collapse = ", "), " not present in Paramters list"))
      }
      #check if arguments are well specified
      Argument_checker <- c(markers_OK = all(Parameters[["markers"]] %in% names(Parameters[["cells"]])),
                            imageID_OK = Parameters[["imageID"]] == "Subject_Names",

                            transform_OK = if(!is.null(Parameters[["transformation"]])){
                              Parameters[["transformation"]] %in% c("asinh", "sqrt")
                            } else(TRUE),
                            method_OK = if(!is.null(Parameters[["method"]])){
                              all(Parameters[["method"]] %in% c('mean', 'minMax', 'trim99', 'PC1'))
                            } else(TRUE),
                            cores_OK = all(Parameters[["cores"]]%%1 == 0, Parameters[["cores"]] > 0)
      )
      Stop_messages <- c(markers_OK = "Marker names not present in the cell Data",
                         imageID_OK = "Image ID must be Subject_Names",
                         transform_OK = "Transformation method must be NULL or one of the following: asinh, sqrt",
                         method_OK = "Method must be NULL or one of the following: NULL, mean, minMax, trim99, PC1",
                         cores_OK = "Core number must be a positive integer value"
      )
      #Check arguments and stop if necessary
      if(!all(Argument_checker)){
        stop(cat(Stop_messages[!Argument_checker],
                 fill = sum(!Argument_checker)))
      }

      #Execute the code
      base::do.call(simpleSeg::normalizeCells, args = Parameters)
    }
  }
