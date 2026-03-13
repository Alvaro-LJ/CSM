#' Computes sample summary metrics of barrier elements
#'
#' Givena dataset obtained using the [Barrier_effect_calculator()], the function computes sample summary metrics for every element analyzed.
#'
#'
#' @param DATA A tibble containing barrier cell data obtained using the [Barrier_effect_calculator()] function.
#' @param Analysis_type A character value indicating the type of metric to be summarized. One of the following: 'Area' (default), 'Distance', 'N_targets'.
#'
#' @returns A list containing tibbles with by sample summary metrics.
#'
#' @seealso [Barrier_effect_calculator()], [Barrier_graph_maker()]
#'
#' @examples
#' \dontrun{
#' #Create temporary input and output directories------------------------------
#' Input_Dir <- tempfile(pattern = "tempdir1_Input")
#' dir.create(Input_Dir, recursive = TRUE)
#'
#' #Save images in Input directory
#' purrr::map(1:2,
#'           function(Image){
#'             EBImage::writeImage(CSM_MiniMultiTiff_test[[Image]],
#'                                file.path(Input_Dir, names(CSM_MiniMultiTiff_test)[Image]))
#'           })
#'
#' #Generate the Image_parameter_list------------------------------
#' Image_Parameters <- list(
#'  PDL1 = list(
#'    Directory = Input_Dir,
#'    Channel = 2,
#'    Image_rotate = NULL,
#'    Image_x_flip = FALSE,
#'    Image_y_flip = TRUE),
#'  GZMB = list(
#'    Directory = Input_Dir,
#'    Channel = 3,
#'    Image_rotate = NULL,
#'    Image_x_flip = FALSE,
#'    Image_y_flip = TRUE))
#'
#' #Compute the cells and pixels in the barrier------------------------------
#' Barrier_DATA <- Barrier_effect_calculator(
#'  N_cores = 1,
#'
#'  DATA = CSM_Phenotypecell_test,
#'
#'
#'  Cell_Of_Origin = "CD8_GZMBneg",
#'  Target_Cell = "TUMOR",
#'  Barrier_Cell = c("OTHER", "CD8_GZMBpos"),
#'  Image_parameter_list = Image_Parameters,
#'
#'  Distance_sampled = 30,
#'
#'  Perform_edge_correction = FALSE,
#'
#'  Polygon_angle = 30,
#'  Target_cell_Angle_tolerance = 30
#' )
#'
#' #Summarize the results by sample------------------------------
#' Barrier_effect_analyzer(
#'  DATA = Barrier_DATA,
#'  Analysis_type = "Area"
#'  )
#'
#' #Plot the results------------------------------
#' Barrier_graph_maker(
#'   DATA_Phenotypes = CSM_Phenotypecell_test,
#'   Phenotypes_included = unique(CSM_Phenotypecell_test$Phenotype),
#'
#'   Barrier_DATA = Barrier_DATA,
#'   Image_name = "ABCM22001_B14_MiniCrop.tif",
#'
#'   Color_by = "GZMB_PixelValueSum_Area_Ratio",
#'   Point_size = 2,
#'
#'   Image_directory = dir(Input_Dir, full.names = TRUE)[[2]],
#'   Pixel_distance_ratio = NULL,
#'   Channel_to_display = 3,
#'   Image_rotate = 0,
#'   Image_x_flip = FALSE,
#'   Image_y_flip = TRUE
#' )
#'}
#'
#' @export

Barrier_effect_analyzer <-
  function(DATA = NULL,
           Analysis_type = "Area"){

    #Check arguments
    if(!identical(names(DATA)[1:5],  c("Cell_no", "COO_X", "COO_Y", "Subject_Names", "Phenotype"))) { #Check if Data is correctly formatted
      stop("DATA provided should have been generated using the Barrier_effect_calculator function")
    }
    if(!all(c("Original_N_Targets", "Analyzed_N_Targets", "Cumulative_distance", "Analyzed_Area", "Polygon_objects") %in% names(DATA))){
      stop("DATA provided should have been generated using the Barrier_effect_calculator function")
    }
    if(!any(Analysis_type %in% c("Area", "Distance", "N_targets"))) stop("Analysis_type must be one of the following: 'Area', 'Distance', 'N_targets")

    #Obtain the requested data
    if(Analysis_type == "N_targets") DATA_analysis <- DATA %>% dplyr::select(1:5, "Analyzed_N_Targets", dplyr::ends_with("_TargetCell_Ratio"))
    if(Analysis_type == "Distance") DATA_analysis <- DATA %>% dplyr::select(1:5, "Cumulative_distance", dplyr::ends_with("_Distance_Ratio"))
    if(Analysis_type == "Area") DATA_analysis <- DATA %>% dplyr::select(1:5, "Analyzed_Area", dplyr::ends_with("_Area_Ratio"))

    #Generate the list with results

    Result_list <-
      purrr::map(DATA_analysis[-c(1:6)],
          function(Feature){
            Interim <- tibble(Subject_Names = DATA_analysis$Subject_Names,
                              Normalization = DATA_analysis[[6]],
                              Feature = Feature)

            if(Analysis_type == "N_targets"){
              return(
                Interim %>% dplyr::group_by(Subject_Names) %>%
                  dplyr::summarise(N_COO = n(),
                                   Average_targets = mean(Normalization),
                                   Mean = mean(Feature),
                                   SD = sd(Feature),
                                   SE = SD/sqrt(N_COO),
                                   Sample_mean_05CI = Mean - 1.96*SE,
                                   Sample_mean_95CI = Mean + 1.96*SE,
                                   min = min(Feature),
                                   p25 = quantile(Feature, 0.25),
                                   p50 = quantile(Feature, 0.5),
                                   p75 = quantile(Feature, 0.75),
                                   max = max(Feature)) %>%
                  dplyr::ungroup() %>%
                  mutate(Sample_mean_05CI = dplyr::case_when(Sample_mean_05CI < 0 ~ 0,
                                                             TRUE ~ Sample_mean_05CI))
              )
            }

            if(Analysis_type == "Distance"){
              return(
              Interim %>% dplyr::group_by(Subject_Names) %>%
                dplyr::summarise(N_COO = n(),
                                 Average_distance = mean(Normalization),
                                 Mean = mean(Feature),
                                 SD = sd(Feature),
                                 SE = SD/sqrt(N_COO),
                                 Sample_mean_05CI = Mean - 1.96*SE,
                                 Sample_mean_95CI = Mean + 1.96*SE,
                                 min = min(Feature),
                                 p25 = quantile(Feature, 0.25),
                                 p50 = quantile(Feature, 0.5),
                                 p75 = quantile(Feature, 0.75),
                                 max = max(Feature)) %>%
                dplyr::ungroup() %>%
                mutate(Sample_mean_05CI = dplyr::case_when(Sample_mean_05CI < 0 ~ 0,
                                                           TRUE ~ Sample_mean_05CI))
              )
            }

            if(Analysis_type == "Area"){
              return(
                Interim %>% dplyr::group_by(Subject_Names) %>%
                  dplyr::summarise(N_COO = n(),
                                   Average_area = mean(Normalization),
                                   Mean = mean(Feature),
                                   SD = sd(Feature),
                                   SE = SD/sqrt(N_COO),
                                   Sample_mean_05CI = Mean - 1.96*SE,
                                   Sample_mean_95CI = Mean + 1.96*SE,
                                   min = min(Feature),
                                   p25 = quantile(Feature, 0.25),
                                   p50 = quantile(Feature, 0.5),
                                   p75 = quantile(Feature, 0.75),
                                   max = max(Feature)) %>%
                  dplyr::ungroup() %>%
                  mutate(Sample_mean_05CI = dplyr::case_when(Sample_mean_05CI < 0 ~ 0,
                                                             TRUE ~ Sample_mean_05CI))
              )
            }
          })

    #Return the final result
    return(Result_list)
  }

