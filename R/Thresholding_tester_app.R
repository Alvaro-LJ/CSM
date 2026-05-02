#' Launches a shiny APP to explore feature thresholding methods
#'
#' `Thresholding_tester_app()` launches an APP to interactively explore thresholding methods. Parameters can then be used in the [Thresholding_function()] or [Thresholding_function_tailored()] functions.
#'
#' @param DATA A dataframe or tibble containing cell feature data.
#' @param Directory Character specifying the path to the folder where images are present.
#' @param Ordered_Channels Character vector specifying image channels in their exact order.
#' @param Max_Gb_cache A single numeric value indicating the memory size of the image cache (see details). 10 Gb by default.
#'
#' @seealso [Thresholding_function()], [Thresholding_function_tailored()]
#'
#' @details
#' In order to deal with large images and speed up image toggling, the shiny APP works with a memoised version of the [Smart_image_importer()] function. Loaded images will
#' be stored in a temporary cache until the APP is closed or the cache reaches it's limit.
#'
#' Although the APP has been designed to work with images, it can be executed without providing an image directory. If no image directory is provided, only cell information
#' will be loaded into the APP.
#'
#' Image settings in the control panel allow the user to control the image channel display.
#'
#' User can use the 'Color', 'Add channel', 'Remove channel' and 'Reset' buttons to merge different channels into a single image. The add channel button
#' will color the current channel according to the color provided and will be merged to the current image.
#'
#' Thresholding methods and type of thresholding (GLOBAL or LOCAL) can be toggled using the control panel.
#'
#' The reset button removes any currently selected cell.
#'
#' Active results are shown in the lower area of the control panel.
#'
#' 4 results panel
#' \itemize{
#' \item{Upper left: Displays the image (use it to zoom in and out).}
#' \item{Upper right: Feature expression by cell. Cells can be selected to explore the % of positive cells from selection.}
#' \item{Lower left: Cells above threshold for the currently selected image. Cells can be selected to explore the % of positive cells from selection.}
#' \item{Lower right: Sample feature expression histogram with the active threshold depicted.}
#' }
#'
#'
#' @examples
#' \dontrun{
#' #Create temporary directory------------------------------
#' Input_Dir <- tempfile(pattern = "tempdir1_Input")
#' dir.create(Input_Dir, recursive = TRUE)
#'
#' #Save images in Input directory
#' purrr::map(1:2,
#' function(Image){
#'    EBImage::writeImage(CSM_MiniMultiTiff_test[[Image]],
#'    file.path(Input_Dir, names(CSM_MiniMultiTiff_test)[Image]))
#' })
#'
#' #Deploy app----------------------------------------------
#'Thresholding_tester_app(
#'    DATA = CSM_Arrangedcellfeaturedata_test,
#'    Directory = Input_Dir,
#'    Ordered_Channels = c("DAPI", "PDL1", "GZMB", "PD1", "CK-EPCAM", "CD8a", "FOXP3")
#')
#'
#'#Remove directories---------------------------------------------------------
#' unlink(Input_Dir, recursive = TRUE)
#' }
#'
#' @export

Thresholding_tester_app <-
  function(DATA,
           Directory = NULL,
           Ordered_Channels = NULL,
           Max_Gb_cache = 10){

    ####Check required packages####
    {
      if(!requireNamespace("magick", quietly = FALSE)) stop(
        paste0("magick CRAN package is required to execute the function. Please install using the following code: ",
               expression(install.packages("magick"))
        )
      )
      if(!requireNamespace("RBioFormats", quietly = TRUE)) stop(
        paste0("RBioFormats Bioconductor package is required to execute the function. Please install using the following code: ",
               expression({
                 if (!require("BiocManager", quietly = TRUE))
                   install.packages("BiocManager")

                 BiocManager::install("RBioFormats")
               })
        )
      )
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

    ####CHECK ARGUMENTS####
    DATA <- DATA

    #Test that DATA provided is adequate
    if(!identical(names(DATA)[1:4], c("Cell_no", "X", "Y", "Subject_Names"))) {
      stop("Please generate an appropiate data object using the Data_arrange_function")
    }
    #Check that ordered channels and channels to keep are unique
    if(length(Ordered_Channels) != length(unique(Ordered_Channels))) stop("Ordered_Channels must contain non-duplicated character values")

    #Check that the Max_Gb_cache is OK
    if(!all(is.numeric(Max_Gb_cache), Max_Gb_cache > 0, length(Max_Gb_cache) == 1)) stop("Max_Gb_cache must be a numeric value > 0")

    #Generate a memoised version of the the smart image importer
    Memoised_importer <- memoise::memoise(
      #Function to be memorized
      Smart_image_importer,

      #Argument controlling memory use
      cache = cachem::cache_mem(
        max_size = Max_Gb_cache * 1024 * 1024 * 1024, #10 Gb is the max amount of bytes
        max_age = Inf,
        max_n = Inf,
        evict = c("lru", "fifo"),
        missing = cachem::key_missing(),
        logfile = NULL
      )
    )


    ####BUILD Look-up table####
    #Try to build it, if no images provided it will return an error and then the app will be executed without images
    try({
      Image_names <- dir(Directory, full.names = FALSE)
      Directory_path <- dir(Directory, full.names = TRUE)
      Image_directory_tibble <- tibble::tibble(Image_name = Image_names,
                                               Image_path = Directory_path)
      print(paste0(length(Image_names), " files found in Directory. Finding matches between file names and Subject_Names in DATA..."))

      #Generate the look-up table
      Look_up_table <-
        purrr::map_dfr(unique(DATA$Subject_Names), function(Name){
          #Compute the distance from each Subject_Name to its closest image
          Distance_vector <- as.double(adist(Name, Image_names, fixed = TRUE, ignore.case = TRUE))
          names(Distance_vector) <- Image_names
          Distance_vector <- sort(Distance_vector)

          #Generate a tibble
          tibble::tibble(Subject_Names = Name,
                         Image_name = names(Distance_vector)[1],
                         Match = Distance_vector[1] == 0,
                         Distance = Distance_vector[1])
        })
      #Bind the paths
      Look_up_table <- dplyr::left_join(Look_up_table, Image_directory_tibble, by = "Image_name")


      #Find duplicates in Image_names and solve conflicts
      Duplicated_images <- duplicated(Look_up_table$Image_name)
      if(any(Duplicated_images)){
        print("Multiple matches found for certain images. Solving conflicts by closest match...")

        Conflictive_images <- unique(Look_up_table$Image_name[duplicated(Look_up_table$Image_name)])

        Losser_Subject_Names <-
          unique(
            unlist(
              purrr::map(Conflictive_images, function(Duplicated_image){
                Subject_Names_in_conflict <- Look_up_table %>% dplyr::filter(Image_name == Duplicated_image) %>% dplyr::arrange(Distance)
                Subject_Names_in_conflict$Subject_Names[-1] #Get all the subject names except the winner
              })
            )
          )

        Look_up_table$Image_name[Look_up_table$Subject_Names %in% Losser_Subject_Names] <- NA
        Look_up_table$Image_path[Look_up_table$Subject_Names %in% Losser_Subject_Names] <- NA

        N_conflictive_subject_names <- sum(Look_up_table$Subject_Names %in% Losser_Subject_Names)

        print(paste0(N_conflictive_subject_names, " out of ", nrow(Look_up_table), " samples without an adequate image."))
      }

      Images_in_Data <- unique(Look_up_table$Subject_Names)
    })



    Channels_in_Data <- names(DATA)[-c(1:4)]
    Channels_in_images <- Ordered_Channels
    Thresholding_methods <- c("EBI_Otsu", "Kmeans", "Kmeans_Otsu", "Autothreshold", "TriClass_Otsu", "Mean", "Quantile", "Arbitrary", "Multi_level")
    Autothreshold_methods <- c("IJDefault", "Huang", "Huang2", "Intermodes", "IsoData", "Li", "MaxEntropy", "Mean",
                               "MinErrorI", "Minimum", "Moments", "Otsu", "RenyiEntropy", "Shanbhag", "Triangle", "Yen")

    #####BUILD THE USER INTERFACE####
    {
      user_interface <- shiny::fluidPage(

        #To make the sidebar collapsable
        shinyjs::useShinyjs(),

        #Set the title and action button
        shiny::fluidRow(
          shiny::column(width = 3, shiny::h3("Thresholding exploration APP")),
          shiny::column(width = 2, shiny::actionButton("toggleSidebar", "Toggle", icon = shiny::icon(name = "square-caret-up", lib = "font-awesome")))
        ),

        #To make the tags work
        shiny::tags$hr(),

        #We want a two panel layout, one in the left containing the input parameters and the output in the right
        shiny::sidebarLayout(
          #Set the first column (which contains the user defined parameters)
          shiny::sidebarPanel(
            #ID and width
            id="sidebar",

            #Button to toggle different rows
            shiny::fluidRow(
              shiny::actionButton("toggle_Image", "Image settings", icon = shiny::icon(name = "square-caret-up", lib = "font-awesome")),
              shiny::actionButton("toggle_Thresh", "Thresholding settings", icon = shiny::icon(name = "square-caret-up", lib = "font-awesome"))
            ),

            #IMAGE PARAMETERS
            shiny::h4("Image settings"),
            #Collapsable rows
            shiny::tags$div(
              class = "Image_settings",

              #Image to be analyzed and basic controls
              shiny::fluidRow(
                shiny::column(6, shiny::selectInput("Image_name", "Image", sort(unique(DATA$Subject_Names)), multiple = FALSE)),
                #Select the image resolution
                shiny::column(2, shiny::selectInput("Res", "Res", c(Low = 5, Mid = 5.7, High = 6, Ultra = 6.3), selected = 5.7, multiple = FALSE)),
                #Pixel adjustment
                shiny::column(2, shinyWidgets::materialSwitch("Change_coords", "Pixel/dist", value = FALSE)),
                shiny::column(2, shiny::conditionalPanel(condition = "input.Change_coords == '1'",
                                                         shiny::textInput("Ratio", "pixel size", value = "1")
                ))
              ),

              #Image channel, cell info and flipping
              shiny::fluidRow(
                shiny::column(4, shiny::selectInput("Channel", "Channel", Channels_in_images, multiple = FALSE)),
                shiny::column(4, shiny::selectInput("Data_Marker_name", "Marker from data", sort(Channels_in_Data), multiple = FALSE)),
                #Flip image
                shiny::column(2, shiny::selectInput("X_flip", "Flip X", c(YES = TRUE, NO = FALSE), selected = FALSE, multiple = FALSE)),
                shiny::column(2, shiny::selectInput("Y_flip", "Flip Y", c(YES = TRUE, NO = FALSE), selected = FALSE, multiple = FALSE))

              ),

              #image display controls
              shiny::fluidRow(
                #Select the equalization
                shiny::column(2, shiny::selectInput("Equalize", "Equalize", c(YES = TRUE, NO = FALSE), selected = FALSE, multiple = FALSE)),

                shiny::column(2, shiny::sliderInput("Degrees", "Rotate", min = -90, max = 90, value = 0, step = 90, ticks = FALSE)),
                #Select the min point of the image
                shiny::column(3, shiny::sliderInput("Min_Image", "Min", value = 0, min = 0, max = 100, step = 1, ticks = FALSE)),
                #Select the max point of the image
                shiny::column(3, shiny::sliderInput("Max_Image", "Max", value = 100, min = 0, max = 100, step = 1, ticks = FALSE)),
                #Select the Gamma
                shiny::column(2, shiny::sliderInput("Gamma", "Gamma", value = 0, min = -3, max = +3, step = 0.01, ticks = FALSE)),
              ),

              #Color merging parameters
              shiny::fluidRow(
                shiny::column(3, shiny::textInput("Color_channel", "Color", value = "blue")),
                shiny::column(3, shiny::actionButton("Add_color", "Add channel", shiny::icon("plus", library = "font-awesome"))),
                shiny::column(3, shiny::actionButton("Remove_color", "Remove channel", shiny::icon("minus", library = "font-awesome"))),
                shiny::column(3, shiny::actionButton("Reset_image", "Reset", icon = shiny::icon("redo")))
              )

            ),

            #IMAGE PARAMETERS
            shiny::h4("Thresholding settings"),
            #Collapsable rows
            shiny::tags$div(
              class = "Thresholding_settings",

              shiny::fluidRow(
                #Select the thresholding method
                shiny::column(8, shiny::selectInput("Threshold_method", "Thresholding method", Thresholding_methods, multiple = FALSE)),
                #Select the thresholding method
                shiny::column(4, shiny::selectInput("Local", "Type of threshold", c(LOCAL = TRUE, GLOBAL = FALSE), selected = FALSE, multiple = FALSE))
              ),
              #Add thresholding parameters according to the selected threshold method
              shiny::conditionalPanel(condition = "input.Threshold_method == 'Autothreshold'",
                                      shiny::selectInput("Autothreshold_method", "Autothreshold method", Autothreshold_methods, selected = "Otsu", multiple = FALSE)
              ),
              shiny::conditionalPanel(condition = "input.Threshold_method == 'TriClass_Otsu'",
                                      shiny::numericInput("TriClass_Iters", "TriClass iterations", value = 10, min = 2, max = 30)
              ),
              shiny::conditionalPanel(condition = "input.Threshold_method == 'Quantile'",
                                      shiny::sliderInput("Percentile", "Quantile selected", value = 0.5, min = 0.01, max = 0.99, step = 0.01)
              ),
              shiny::conditionalPanel(condition = "input.Threshold_method == 'Arbitrary'",
                                      shiny::textInput("Arbitrary", "User selected threshold", placeholder = "0.9")
              ),
              shiny::conditionalPanel(condition = "input.Threshold_method == 'Multi_level'",
                                      shiny::numericInput("Levels", "Number of Levels", value = 3, min = 2, max = 10)
              )
            ),

            #The UI will be completed with summary tables of the sample
            shiny::fluidRow(
              shiny::column(6, htmltools::p("Final threshold/s is: ", shiny::tableOutput("Final_threshold"))),
              shiny::column(6, htmltools::p("Sample Summary: ", shiny::tableOutput("Summary")))
            ),
            #Also it will include a summary of selected cells
            shiny::fluidRow(
              shiny::column(6, htmltools::p("Selected cells: ", shiny::tableOutput("Cell_selection"))),
              shiny::column(3, shiny::actionButton("reset", shiny::icon("redo"), label = "Reset selection"))
            )
          ),

          #Set the outcome columns
          shiny::mainPanel(

            #Give the main panel an ID
            id = "mainPanel",

            #First row will have the Photo and the overview of marker intensity by cell
            shiny::fluidRow(
              shiny::column(5, shiny::plotOutput("Photo",
                                                 #Controls for zoom in
                                                 dblclick = "Photo_dblclick",
                                                 brush = shiny::brushOpts(id = "Photo_brush",
                                                                          resetOnNew = TRUE)
              )
              ),
              shiny::column(5, ggiraph::girafeOutput("Cell_by_intensity"))
            ),
            #Second row will contain the positive cells and the histogram
            shiny::fluidRow(
              shiny::column(5, ggiraph::girafeOutput("Positive_cells")),
              shiny::column(5, shiny::plotOutput("Histogram"))
            )
          )
        ),
        shiny::tags$head(shiny::tags$style(
          htmltools::HTML('
         #sidebar {
            background-color: #d2dbfa;
        }

        body, label, input, button, select {
          font-family: "Arial";
        }'))),

        shiny::tags$script(htmltools::HTML("
    $(document).on('click', '#toggle_Image', function() {
      $('.Image_settings').toggle();
    });
    $(document).on('click', '#toggle_Thresh', function() {
      $('.Thresholding_settings').toggle();
    });
  "))
      )
    }
    ####BUILD THE SERVER####
    server <- function(input, output, session){

      # JS function to switch classes when togglin sidebar
      toggle_script <- "
    if ($('#sidebar').is(':visible')) {
      // hide sidebar + expand main panel
      $('#sidebar').hide();
      $('#mainPanel').removeClass('col-sm-8').addClass('col-sm-12');
    } else {
      // show sidebar + restore width
      $('#sidebar').show();
      $('#mainPanel').removeClass('col-sm-12').addClass('col-sm-8');
    }
  "
      #What to do if the user hits the toggle button
      shiny::observeEvent(input$toggleSidebar, {
        shinyjs::runjs(toggle_script)
      })

      #All the reactives to be used
      #Generate a reactive with the real Image name and the channel number
      Photo_name <- shiny::reactive(try(as.character((Look_up_table %>% dplyr::filter(Subject_Names == input$Image_name))[1,5])))
      Channel_index <- shiny::reactive(which(input$Channel == Ordered_Channels))
      Photo_min <- shiny::reactive(input$Min_Image)
      Photo_max <- shiny::reactive(input$Max_Image)
      Photo_gamma <- shiny::reactive(10^input$Gamma)
      Equalize <- shiny::reactive(input$Equalize)
      #Reactive expression to control the graphs
      Case_id <- shiny::reactive(input$Image_name)
      Variable <- shiny::reactive(input$Data_Marker_name)
      Local_threshold <- shiny::reactive(input$Local)
      Threshold_method <- shiny::reactive(input$Threshold_method)
      Autothreshold_method <- shiny::reactive(input$Autothreshold_method)
      TriClass_Iters <- shiny::reactive(input$TriClass_Iters)
      Percentile <- shiny::reactive(input$Percentile)
      Arbitrary <- shiny::reactive(input$Arbitrary)
      Levels <- shiny::reactive(input$Levels)
      X_flip <- shiny::reactive(input$X_flip)
      Y_flip <- shiny::reactive(input$Y_flip)
      Degrees_rotate <- shiny::reactive(input$Degrees)
      Pixel_dist_conversion <- shiny::reactive(input$Change_coords)
      Pixel_dist_ratio <- shiny::reactive(input$Ratio)
      ranges <- shiny::reactiveValues(x = NULL, y = NULL)#Control the user selected ranges
      #Generates a reactive that stores the original image coordinates and cropping parameters
      Original_ranges <-
        shiny::reactiveValues(
          #Image dimension Original
          Original_Dims_Width = NULL,
          Original_Dims_Height = NULL
        )
      Resolution <- shiny::reactive(input$Res)
      #Control the data source
      Source_DATA <- shiny::reactive({
        Final_DATA <- DATA %>% dplyr::select(1:4, all_of(Variable()))
        names(Final_DATA)[5] <- "Marker"
        #Modify pixel values if required
        if(as.logical(Pixel_dist_conversion())){
          Final_DATA$X <- Final_DATA$X * as.numeric(Pixel_dist_ratio())
          Final_DATA$Y <- Final_DATA$Y * as.numeric(Pixel_dist_ratio())
        }
        #Return the final data
        return(Final_DATA)
      })

      #Generate the thresholded DATA
      Thresholded_DATA <- shiny::reactive({
        DATA_thresholded <- Thresholding_function(
          DATA = Source_DATA(),
          Strategy = as.character(Threshold_method()),
          Local_thresholding = as.logical(Local_threshold()),
          Method_autothreshold = as.character(Autothreshold_method()),
          number_iterations_TriClass = as.numeric(TriClass_Iters()),
          Percentile = as.numeric(Percentile()),
          Defined_threshold = as.numeric(Arbitrary()),
          Levels = as.numeric(Levels())
        )
        return(DATA_thresholded)
      })
      #Obtain a list with the histogram and the final thresholds using the dedicated function
      Histo_list <- shiny::reactive(SHINY_Thresholding_summary_function(DATA = Source_DATA(),
                                                                        DATA_Thresholded = Thresholded_DATA(),
                                                                        LOCAL = Local_threshold(),
                                                                        CASE = Case_id())
      )

      #The reactive that controls the parameter list for image coloring (Null by default at start)
      Image_coloring_list <- shiny::reactiveValues(Parameter_list = NULL)

      #Control the buttons for the color list
      #The ADD CHANNEL BUTTON
      shiny::observeEvent(
        input$Add_color,
        {

          #First generate the element_list
          Current_Channel_list <-
            list(
              channel_index = which(input$Channel == Ordered_Channels),
              color = input$Color_channel,
              gamma = as.numeric(10^input$Gamma),
              min = as.numeric(input$Min_Image),
              max = as.numeric(input$Max_Image),
              equalize = as.logical(input$Equalize)
            )

          if(berryFunctions::is.error(grDevices::col2rgb(Current_Channel_list$color))){
            shiny::showModal(modalDialog(
              paste0(Current_Channel_list$color, " is not an adequate color name or HEX code, please review"),
              easyClose = TRUE,
              footer = NULL)
            )
          }

          #If the Image_coloring_list is NULL then generate this single output
          else if(is.null(Image_coloring_list$Parameter_list)){
            #Generate the list
            Image_coloring_list$Parameter_list <- list(Current_Channel_list)
            #Add the names
            names(Image_coloring_list$Parameter_list) <- input$Channel
          }

          #If the Image_coloring_list is not NULL, check if the color or the input names are already present
          else{

            #Generate logical values if any of both are present
            Channel_present_in_list <- input$Channel %in% names(Image_coloring_list$Parameter_list)
            Color_present_in_list <- input$Color_channel %in% purrr::map_chr(Image_coloring_list$Parameter_list, ~.[["color"]])

            #If not present, we need to create a new element in the list
            if(!any(Channel_present_in_list, Color_present_in_list)){
              Image_coloring_list$Parameter_list <- append(Image_coloring_list$Parameter_list, list(Current_Channel_list))
              names(Image_coloring_list$Parameter_list)[length(Image_coloring_list$Parameter_list)] <- input$Channel
            }

            #If channel name is already present and color is not present
            if(Channel_present_in_list & !Color_present_in_list){
              #print a dialogue box
              shiny::showModal(shiny::modalDialog(
                paste0(input$Channel, " already present, it will be overwritten"),
                easyClose = TRUE,
                footer = NULL)
              )
              #Proceed to replace the channel
              Image_coloring_list$Parameter_list[[input$Channel]] <- Current_Channel_list
            }

            #If channel name is not present and color is repeated
            if(!Channel_present_in_list & Color_present_in_list){
              #print a dialogue box
              shiny::showModal(modalDialog(
                paste0(input$Color_channel, " already in use, it will be overwritten"),
                easyClose = TRUE,
                footer = NULL)
              )
              #Proceed to replace the channel
              Conflictive_element <- which(purrr::map_chr(Image_coloring_list$Parameter_list, ~.[["color"]]) == input$Color_channel)
              Image_coloring_list$Parameter_list[[Conflictive_element]] <- Current_Channel_list
              names(Image_coloring_list$Parameter_list)[[Conflictive_element]] <- input$Channel
            }

            #If channel name is present and the color is in use delete the color and replace the channel
            if(Channel_present_in_list & Color_present_in_list){
              #print a dialogue box
              shiny::showModal(modalDialog(
                paste0(input$Color_channel, " already in use and ", input$Channel, " already present, they will be overwritten"),
                easyClose = TRUE,
                footer = NULL)
              )
              #If they are the same, then only do one action
              Conflictive_channel <- which(names(Image_coloring_list$Parameter_list) == input$Channel)
              Conflictive_color <- which(purrr::map_chr(Image_coloring_list$Parameter_list, ~.[["color"]]) == input$Color_channel)
              if(Conflictive_channel == Conflictive_color) Image_coloring_list$Parameter_list[[input$Channel]] <- Current_Channel_list
              #Else remove the color and add then replace the channel
              else{
                #Remove the conflictive color
                Image_coloring_list$Parameter_list <- Image_coloring_list$Parameter_list[-Conflictive_color]
                #Replace the conflictive channel
                Image_coloring_list$Parameter_list[[input$Channel]] <- Current_Channel_list
              }
            }

          }
        })

      #REMOVE BUTTON
      shiny::observeEvent(
        input$Remove_color,
        {
          #If channel is not present print dialogue and don't do anything
          if(!input$Channel %in% names(Image_coloring_list$Parameter_list)){
            shiny::showModal(shiny::modalDialog(
              paste0(input$Channel, " not being used in image, therefore it cannot be removed"),
              easyClose = TRUE,
              footer = NULL)
            )
          }
          else{
            #print a dialogue box
            shiny::showModal(shiny::modalDialog(
              paste0(input$Channel, " will be removed"),
              easyClose = TRUE,
              footer = NULL)
            )
            #Find conflictive channel
            Conflictive_channel <- which(names(Image_coloring_list$Parameter_list) == input$Channel)
            Image_coloring_list$Parameter_list <- Image_coloring_list$Parameter_list[-Conflictive_channel]

            #If the length of Image_coloring_list$Parameter_list is 0 then turn it to null again
            if(length(Image_coloring_list$Parameter_list) == 0) Image_coloring_list$Parameter_list <- NULL
          }
        })

      #RESET BUTTON
      shiny::observeEvent(
        input$Reset_image,
        {
          #Show an alert
          shinyalert::shinyalert(title = "WARNING!",
                                 text = "Merged image will be removed",
                                 type = "warning",
                                 closeOnEsc = TRUE,
                                 closeOnClickOutside = TRUE,
                                 showCancelButton = TRUE,
                                 showConfirmButton = TRUE,
                                 confirmButtonText = "Proceed",
                                 cancelButtonText = "Cancel",
                                 callbackR = function() Image_coloring_list$Parameter_list <- NULL
          )

        })

      #Reactive that imports the photograph and does the pre-processing
      Photo_reactive <- shiny::reactive({
        ####MOD CROPPING COORDINATES####
        Final_X_crop_coordinates <- NULL
        Final_Y_crop_coordinates <- NULL

        #If cropping is required and image is rotated or flipped modify accordingly the cropping coordinates
        if(any(!is.null(ranges$x), !is.null(ranges$y))){
          #Obtain the preliminary image cropping coordinates
          Final_X_crop_coordinates <- sort(ceiling(ranges$x))
          Final_Y_crop_coordinates <- sort(ceiling(ranges$y))

          #Account for X_flip
          if(as.logical(X_flip())){
            Final_X_crop_coordinates <- sort(ceiling(Original_ranges$Original_Dims_Width - Final_X_crop_coordinates))
          }
          #Account for Y flip
          if(as.logical(Y_flip())){
            Final_Y_crop_coordinates <- sort(ceiling(Original_ranges$Original_Dims_Height - Final_Y_crop_coordinates))
          }

          #Account for rotation
          if(as.numeric(Degrees_rotate()) == 90){
            Old_X_coordinates <- Final_X_crop_coordinates
            Old_Y_coordinates <- Final_Y_crop_coordinates
            #Y is now X
            Final_Y_crop_coordinates <- Old_X_coordinates
            #X is now Y
            Final_X_crop_coordinates <- Old_Y_coordinates
          }


          if(as.numeric(Degrees_rotate()) == -90){
            Old_X_coordinates <- Final_X_crop_coordinates
            Old_Y_coordinates <- Final_Y_crop_coordinates

            #X is now the width minus the active Y coordinates
            Final_X_crop_coordinates <- sort(ceiling(Original_ranges$Original_Dims_Width - Old_Y_coordinates))

            #Y is active X
            Final_Y_crop_coordinates <- Old_X_coordinates
          }
        }

        ####OBTAIN THE IMAGE####
        Image <- Memoised_importer(
          N_cores = 2,
          Image_directory = Photo_name(),
          Log10_pixel_output = as.numeric(Resolution()),
          X_crop_coordinates = Final_X_crop_coordinates,
          Y_crop_coordinates = Final_Y_crop_coordinates
        )

        ####Update Image and current Image dimension paramteres(will be used if crop is required by user)####
        Original_ranges$Original_Dims_Width <- Image$Original_Dims[1]
        Original_ranges$Original_Dims_Height <- Image$Original_Dims[2]

        ####COLOR MERGING NOT REQUIRED####
        if(is.null(Image_coloring_list$Parameter_list)){
          #Get the single channel to be obtained
          Image$Image <- Image$Image[as.numeric(Channel_index())]

          #Modify min max and gamma
          Image$Image <- Image$Image %>%
            magick::image_level(black_point = as.numeric(Photo_min()),
                                white_point = as.numeric(Photo_max()),
                                mid_point = as.numeric(Photo_gamma()))

          #Equalize if necessary
          if(as.logical(Equalize())) Image$Image <- Image$Image %>% magick::image_equalize()
        }

        ####COLOR MERGING REQUIRED####
        if(!is.null(Image_coloring_list$Parameter_list)){

          #Perform RGB color merging
          Image$Image <- Image_channel_color_merger(Image = Image$Image,
                                                    Parameter_list = Image_coloring_list$Parameter_list)

        }

        ####ROTATE (ALSO ROTATE ORIGINAL DIMS, CURRENT DIMS and CROP COORDS if angle is 90 or -90)####
        if(as.numeric(Degrees_rotate()) != 0){
          #Rotate the image
          Image$Image <- Image$Image %>% magick::image_rotate(degrees = as.numeric(Degrees_rotate()))

          #Change the associated dimensions (for adequate plotting)
          Image$Original_Dims <- rev(Image$Original_Dims)
          Image$Current_Dims <- rev(Image$Current_Dims)

          #If the image has been cropped the switch X for Y
          if(any(!is.null(ranges$x), !is.null(ranges$y))){
            #Get the old crops in a separate object
            Old_X_crop_coords <- Image$X_crop_coords
            Old_Y_crop_coords <- Image$Y_crop_coords

            #Make the switch
            Image$X_crop_coords <- Old_Y_crop_coords
            Image$Y_crop_coords <- Old_X_crop_coords
          }
        }

        ####FLIP AND FLOP####
        if(as.logical(X_flip())) Image$Image <- Image$Image %>% magick::image_flop()
        if(!as.logical(Y_flip())) Image$Image <- Image$Image %>% magick::image_flip()#Opposite due to ggplot2 graph plotting for images


        ####RETURN FINAL PRODUCT####
        return(Image)

      })

      #PHOTO AS GGPLOT OBJECT (with adequate axis)
      Photo_plot_reactive <- shiny::reactive({

        #Obtain the Image
        Photo <- Photo_reactive()$Image

        #Obtain the final axis limits
        #If no cropping has been performed plot according to image dimension
        if(all(is.null(ranges$x), is.null(ranges$y))){
          Axis_Width_Min <- 1
          Axis_Width_Max <- Photo_reactive()$Original_Dims[1]
          Axis_Height_Min <- 1
          Axis_Height_Max <- Photo_reactive()$Original_Dims[2]
        }
        #If cropping has been performed then get the dimensions from the Cropping dimensions
        if(any(!is.null(ranges$x), !is.null(ranges$y))){
          Axis_Width_Min <- min(Photo_reactive()$X_crop_coords)
          Axis_Width_Max <- max(Photo_reactive()$X_crop_coords)
          Axis_Height_Min <- min(Photo_reactive()$Y_crop_coords)
          Axis_Height_Max <- max(Photo_reactive()$Y_crop_coords)
        }

        #Obtain the size in pixels to compute the aspect ratio
        Aspect_ratio_photo <- abs(Axis_Height_Max - Axis_Height_Min) / abs(Axis_Width_Max - Axis_Width_Min)

        #Return the result as a scaffold ggplot_object
        Photo_plot <- ggplot() +
          annotation_raster(Photo, xmin = Axis_Width_Min, xmax = Axis_Width_Max, ymin = Axis_Height_Min, ymax = Axis_Height_Max, interpolate = TRUE)+
          scale_x_continuous(limits = c(Axis_Width_Min, Axis_Width_Max)) +
          scale_y_continuous(limits = c(Axis_Height_Min, Axis_Height_Max)) +
          coord_cartesian(xlim = c(Axis_Width_Min, Axis_Width_Max), ylim = c(Axis_Height_Min, Axis_Height_Max), expand = FALSE)+
          theme(axis.title = element_blank(),
                axis.text = element_blank(),
                axis.ticks = element_blank(),
                axis.line = element_blank(),
                panel.background = element_rect(fill = "black"),
                panel.grid = element_blank(),
                aspect.ratio = Aspect_ratio_photo)
        return(Photo_plot)
      })

      #Print the photo
      output$Photo <- shiny::renderPlot({
        #Try to Plot the result
        try(Photo <- Photo_plot_reactive())

        #Obtain the data to plot
        DATA_plot <- Source_DATA() %>% dplyr::filter(Subject_Names == Case_id())
        #Modify pixel values if required
        if(as.logical(Pixel_dist_conversion())){
          DATA_plot$X <- DATA_plot$X * as.numeric(Pixel_dist_ratio())
          DATA_plot$Y <- DATA_plot$Y * as.numeric(Pixel_dist_ratio())
        }

        #Remove any cell not being plotted
        if(any(!is.null(ranges$x), !is.null(ranges$y))){
          DATA_plot <- DATA_plot %>% dplyr::filter(X >= min(ranges$x), X <= max(ranges$x),
                                                   Y >= min(ranges$y), Y <= max(ranges$y))
        }

        #If the photo returns an error return the point image
        if(berryFunctions::is.error(Photo)){
          return(
            DATA_plot %>% ggplot() +
              geom_point(aes(x = X, y = Y),
                         color = "white",
                         size = 2.5) +
              scale_x_continuous(limits = c(min(DATA_plot$X), max(DATA_plot$X))) +
              scale_y_continuous(limits = c(min(DATA_plot$Y), max(DATA_plot$Y))) +
              coord_cartesian(xlim = ranges$x, ylim = ranges$y, expand = FALSE)+
              theme(axis.title = element_blank(),
                    axis.text = element_blank(),
                    axis.ticks = element_blank(),
                    axis.line = element_blank(),
                    panel.background = element_rect(fill = "black"),
                    panel.grid = element_blank()) +
              annotate("text", x = quantile(DATA_plot$X, 0.5), y = quantile(DATA_plot$Y, 0.5),
                       color = "red", size = 10, hjust = 0.5,
                       label = "UNABLE TO RENDER PHOTO\nUSE ME TO ZOOM IN")
          )

        }
        #else Generate the plot with the photo
        else{
          #If no legend required
          if(is.null(Image_coloring_list$Parameter_list)) return(Photo)

          #If legend is required due to channel merging
          if(!is.null(Image_coloring_list$Parameter_list)){
            #Generate the color tibble to generate the legend
            color_tibble <-
              tibble(Channel_index = 1:length(Image_coloring_list$Parameter_list),
                     Channel_name = names(Image_coloring_list$Parameter_list),
                     Color_names = purrr::map_chr(Image_coloring_list$Parameter_list, ~.[["color"]])
              ) %>% arrange(Channel_index)
            color_tibble$Channel_name = factor(color_tibble$Channel_name, levels = color_tibble$Channel_name)

            return(
              Photo +
                geom_point(aes(x = 1, y = 1, color = Channel_name), size = 0,
                           data = color_tibble) +
                scale_color_manual("", values = color_tibble$Color_names) +
                guides(color = guide_legend(ncol = 4,
                                            byrow = TRUE,
                                            override.aes = list(size = 10, shape = 15),
                                            position = "bottom",
                                            theme = theme(legend.text = element_text(color = "white", size = 10, face = "bold"),
                                                          legend.key = element_rect(fill = "black"),
                                                          legend.background = element_rect(fill = "black"),
                                                          legend.spacing = unit(-10, "pt"),
                                                          legend.margin = margin(-0, -0, -0, -0),
                                                          legend.box.margin = margin(-10, -10, -10, -10))))
            )
          }
        }
      })

      #Control the zoom in of the Photo
      shiny::observeEvent(input$Photo_dblclick, {
        brush <- input$Photo_brush
        if (!is.null(brush)) {
          ranges$x <- c(brush$xmin, brush$xmax)
          ranges$y <- c(brush$ymin, brush$ymax)

        } else {
          ranges$x <- NULL
          ranges$y <- NULL
        }
      })

      #All cells marker expression
      #Create a reactive that will generate the very basic PLOT
      Cell_intensity_plot <-
        shiny::reactive({
          #Get data and change the axis to coordinate all three plots displayed
          DATA_plot <- Source_DATA() %>% dplyr::filter(Subject_Names == Case_id())

          #Try to import the photo plot
          try(Photo <- Photo_plot_reactive())

          #Generate the color
          color_fun <-
            circlize::colorRamp2(breaks = c(min(DATA_plot$Marker), quantile(DATA_plot$Marker, 0.99)),
                                 colors = c(alpha("black", 0),
                                            alpha("red", 1))
            )
          DATA_plot$Color <- color_fun(DATA_plot$Marker)

          #Remove any cell not being plotted
          if(any(!is.null(ranges$x), !is.null(ranges$y))){
            DATA_plot <- DATA_plot %>% dplyr::filter(X >= min(ranges$x), X <= max(ranges$x),
                                                     Y >= min(ranges$y), Y <= max(ranges$y))
          }

          #Generate the final plot
          #If the Photo returns an error do not plot it
          if(berryFunctions::is.error(Photo)){
            return(
              DATA_plot %>% ggplot() +
                ggiraph::geom_point_interactive(aes(x = X, y = Y, color = Color,
                                                    data_id = Cell_no,
                                                    tooltip = stringr::str_c(as.character(Cell_no), " ", "Value = ", as.character(round(Marker, 4)))),
                                                hover_nearest = FALSE,
                                                size = 2.5) +
                scale_color_identity() +
                scale_x_continuous(limits = c(min(DATA_plot$X), max(DATA_plot$X))) +
                scale_y_continuous(limits = c(min(DATA_plot$Y), max(DATA_plot$Y))) +
                coord_cartesian(xlim = ranges$x, ylim = ranges$y, expand = FALSE)+
                cowplot::theme_cowplot()+
                guides(colour = "none")+
                ggtitle("Feature expression by cell") +
                theme(axis.line = element_blank(),
                      axis.ticks = element_blank(),
                      axis.text = element_blank(),
                      axis.title = element_blank(),
                      panel.background = element_rect(fill = "black"),
                      legend.text = element_text(size = 10),
                      plot.title = element_text(size = 12, hjust = 0.5))
            )
          }
          #Else plot as usual
          else{
            Final_plot <-
              Photo +
              ggiraph::geom_point_interactive(aes(x = X, y = Y, color = Color,
                                                  data_id = Cell_no,
                                                  tooltip = stringr::str_c(as.character(Cell_no), " ", "Value = ", as.character(round(Marker, 4)))),
                                              hover_nearest = FALSE,
                                              size = 2.5,
                                              data = DATA_plot) +
              scale_color_identity() +
              ggtitle("Feature expression by cell") +
              guides(colour = "none") +
              theme(plot.title = element_text(size = 12, hjust = 0.5))
            return(Final_plot)
          }

        })
      #Send the plot to the UI
      output$Cell_by_intensity <- ggiraph::renderGirafe({
        plot <- ggiraph::girafe(code = print(Cell_intensity_plot()),
                                options = list(
                                  ggiraph::opts_hover(css = "stroke:black;cursor:pointer;", reactive = TRUE),
                                  ggiraph::opts_selection(type = "multiple", css = "fill:#FF3333;stroke:black;")
                                )
        )
        return(plot)
      })

      #Positive cells
      #Create a reactive that will generate the very basic PLOT
      Positive_cells_plot <- shiny::reactive({
        DATA_threshold <- Thresholded_DATA() %>% dplyr::filter(Subject_Names == Case_id())
        return(DATA_threshold)
      })
      #Send the plot to the UI according to the thresholding method
      output$Positive_cells <- ggiraph::renderGirafe({

        Positive_cells <- Positive_cells_plot()
        #Remove any cell not being plotted
        if(any(!is.null(ranges$x), !is.null(ranges$y))){
          Positive_cells <- Positive_cells %>% dplyr::filter(X >= min(ranges$x), X <= max(ranges$x),
                                                             Y >= min(ranges$y), Y <= max(ranges$y))
        }

        #Try to Import the photo
        try(Photo <- Photo_plot_reactive())

        #If there is an error with the Photo still execute the graph
        if(berryFunctions::is.error(Photo)){
          #Define behavior for 2 levels
          if(Threshold_method() != "Multi_level"){
            #Add specific code for non Multi_level
            Plot_code <- Positive_cells %>%
              ggplot() +
              cowplot::theme_cowplot()+
              theme(axis.line = element_blank(),
                    axis.ticks = element_blank(),
                    axis.text = element_blank(),
                    axis.title = element_blank(),
                    panel.background = element_rect(fill = "black"),
                    legend.text = element_text(size = 10),
                    plot.title = element_text(size = 12, hjust = 0.5)) +
              scale_colour_manual("", values = c(alpha("black", 0), "red")) +
              guides(color = "none") +
              scale_x_continuous(limits = c(min(Positive_cells$X), max(Positive_cells$X))) +
              scale_y_continuous(limits = c(min(Positive_cells$Y), max(Positive_cells$Y))) +
              coord_cartesian(xlim = ranges$x, ylim = ranges$y, expand = FALSE) +
              scale_fill_identity() +
              ggtitle("Cells above threshold") +
              ggiraph::geom_point_interactive(aes(x = X, y = Y, color = Marker,
                                                  data_id = Cell_no,
                                                  tooltip = as.character(Cell_no)),
                                              size = 2.5)
            #Send it to plot
            plot <- ggiraph::girafe(code = print(Plot_code),
                                    options = list(
                                      ggiraph::opts_hover(css = "stroke:black;cursor:pointer;", reactive = TRUE),
                                      ggiraph::opts_selection(type = "multiple", css = "fill:#FF3333;stroke:black;")
                                    ))
            return(plot)
          }
          #Define behavior for Multi-level
          else{
            #Add specific code for multi_level
            Plot_code <- Positive_cells %>%
              ggplot() +
              cowplot::theme_cowplot()+
              theme(axis.line = element_blank(),
                    axis.ticks = element_blank(),
                    axis.text = element_blank(),
                    axis.title = element_blank(),
                    panel.background = element_rect(fill = "black"),
                    legend.text = element_text(size = 10),
                    plot.title = element_text(size = 12, hjust = 0.5)) +
              scale_colour_manual("", values = c(alpha("black", 0), RColorBrewer::brewer.pal(n = Levels()-1, "Reds")))+
              scale_x_continuous(limits = c(min(Positive_cells$X), max(Positive_cells$X))) +
              scale_y_continuous(limits = c(min(Positive_cells$Y), max(Positive_cells$Y))) +
              coord_cartesian(xlim = ranges$x, ylim = ranges$y, expand = FALSE) +
              scale_fill_identity() +
              ggtitle("Cells above threshold") +
              ggiraph::geom_point_interactive(aes(x = X, y = Y, color = as.factor(Marker),
                                                  data_id = Cell_no,
                                                  tooltip = stringr::str_c("Value = ", as.character(round(Marker, 0)))),
                                              size = 2.5)
            #Send it to plot
            plot <- ggiraph::girafe(code = print(Plot_code),
                                    options = list(
                                      ggiraph::opts_hover(css = "stroke:black;cursor:pointer;", reactive = TRUE),
                                      ggiraph::opts_selection(type = "multiple", css = "fill:#FF3333;stroke:black;")
                                    ))
            return(plot)
          }
        }
        #If no error then continue with normal execution
        else{
          #Define behavior for 2 levels
          if(Threshold_method() != "Multi_level"){
            #Add specific code for non Multi_level
            Plot_code <- Photo +
              scale_colour_manual("", values = c(alpha("black", 0), "red")) +
              guides(color = "none") +
              ggiraph::geom_point_interactive(aes(x = X, y = Y, color = Marker,
                                                  data_id = Cell_no,
                                                  tooltip = as.character(Cell_no)),
                                              size = 2.5,
                                              data = Positive_cells) +
              ggtitle("Cells above threshold") +
              theme(plot.title = element_text(size = 12, hjust = 0.5))

            #Send it to plot
            plot <- ggiraph::girafe(code = print(Plot_code),
                                    options = list(
                                      ggiraph::opts_hover(css = "stroke:black;cursor:pointer;", reactive = TRUE),
                                      ggiraph::opts_selection(type = "multiple", css = "fill:#FF3333;stroke:black;")
                                    ))
            return(plot)
          }
          #Define behavior for Multi-level
          else{
            #Add specific code for multi_level
            Plot_code <- Photo +
              theme(legend.text = element_text(size = 10)) +
              scale_colour_manual("", values = c(alpha("black", 0), RColorBrewer::brewer.pal(n = Levels()-1, "Reds")))+
              ggiraph::geom_point_interactive(aes(x = X, y = Y, color = as.factor(Marker),
                                                  data_id = Cell_no,
                                                  tooltip =stringr::str_c("Value = ", as.character(round(Marker, 0)))),
                                              size = 2.5,
                                              data = Positive_cells) +
              ggtitle("Cells above threshold") +
              theme(plot.title = element_text(size = 12, hjust = 0.5))
            #Send it to plot
            plot <- ggiraph::girafe(code = print(Plot_code),
                                    options = list(
                                      ggiraph::opts_hover(css = "stroke:black;cursor:pointer;", reactive = TRUE),
                                      ggiraph::opts_selection(type = "multiple", css = "fill:#FF3333;stroke:black;")
                                    ))
            return(plot)
          }
        }
      })

      #Histogram
      output$Histogram <- shiny::renderPlot(
        plot(Histo_list()[[1]])
      )
      #Threshold text
      output$Final_threshold <- shiny::renderTable(Histo_list()[[2]])

      #Threshold sample summary
      output$Summary <- shiny::renderTable(Histo_list()[[3]])

      #Selected cells and reset button
      selected_cells_intensity <- shiny::reactive(input$Cell_by_intensity_selected)
      selected_cells_positive <- shiny::reactive(input$Positive_cells_selected)
      #What to do in case the user hits the reset button
      shiny::observeEvent(input$reset, {
        session$sendCustomMessage(type = 'Cell_by_intensity_set', message = character(0))
        session$sendCustomMessage(type = 'Positive_cells_set', message = character(0))
      })
      #Generate the output tibble
      output$Cell_selection <- shiny::renderTable({
        #Get cell_no in both plots and remove duplicates
        Cells <- unique(c(selected_cells_intensity(), selected_cells_positive()))
        #Get Marker data and threshold data an start processing
        Marker_data <- Source_DATA() %>% dplyr::filter(Subject_Names == Case_id(), Cell_no %in% Cells) %>% dplyr::select(-Subject_Names)
        Threshold_data <- Thresholded_DATA() %>% dplyr::filter(Subject_Names == Case_id(), Cell_no %in% Cells) %>% dplyr::select(-Subject_Names, -X, -Y)
        Final <- Threshold_data %>% dplyr::count(Marker)
        names(Final) <- c("Cell", "number")
        Final
      })

      #If browser is closed end the app
      session$onSessionEnded(function() {
        future::plan("future::sequential")
        memoise::forget(Memoised_importer)
        gc()
        shiny::stopApp()
      })
    }

    #Run the server
    message("Always stop current R execution if you want to continue with your R session")
    shiny::shinyApp(user_interface, server)
  }

