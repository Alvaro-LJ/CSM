#' Launches a shiny APP to explore feature thresholding methods
#'
#' `Thresholding_tester_app()` launches an APP to interactively explore thresholding methods. Parameters can then be used in the [Thresholding_function()] or [Thresholding_function_tailored()] functions.
#'
#' @param DATA A dataframe or tibble containing cell feature data.
#' @param Directory Character specifying the path to the folder where images are present.
#' @param Ordered_Channels Character vector specifying image channels in their exact order.
#'
#' @seealso [Thresholding_function()], [Thresholding_function_tailored()]
#'
#' @details
#' Image settings in the control panel allow the user to control the image channel display.
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
           Directory,
           Ordered_Channels){

    DATA <- DATA

    #Test that DATA provided is adequate
    if(!identical(names(DATA)[1:4], c("Cell_no", "X", "Y", "Subject_Names"))) {
      stop("Please generate an appropiate data object using the Data_arrange_function")
    }

    Images_in_Data <- unique(DATA$Subject_Names)
    Channels_in_Data <- names(DATA)[-c(1:4)]
    Real_Images <- dir(Directory, full.names = FALSE)
    Channels_in_images <- Ordered_Channels
    Thresholding_methods <- c("EBI_Otsu", "Kmeans", "Kmeans_Otsu", "Autothreshold", "TriClass_Otsu", "Mean", "Quantile", "Arbitrary", "Multi_level")
    Autothreshold_methods <- c("IJDefault", "Huang", "Huang2", "Intermodes", "IsoData", "Li", "MaxEntropy", "Mean",
                               "MinErrorI", "Minimum", "Moments", "Otsu", "RenyiEntropy", "Shanbhag", "Triangle", "Yen")

    #BUILD THE USER INTERFACE
    {
      user_interface <- shiny::fluidPage(

        #Set the title
        shiny::titlePanel("Thresholding exploration APP"),

        #We want a two panel layout, one in the left containing the input parameters and the output in the right
        shiny::sidebarLayout(
          #Set the first column (which contains the user defined parameters)
          shiny::sidebarPanel(
            #ID and width
            id="sidebar",

            shiny::fluidRow(
              #Select Image to be analyzed
              shiny::column(7, shiny::selectInput("Data_Image_name", "Image from data", sort(Images_in_Data), multiple = FALSE)),
              #Select the marker to be analyzed
              shiny::column(5, shiny::selectInput("Data_Marker_name", "Marker from data", sort(Channels_in_Data), multiple = FALSE))
            ),
            shiny::fluidRow(
              #Select the real image to be displayed
              shiny::column(7, shiny::selectInput("Real_Image_name", "Image to display", sort(Real_Images), multiple = FALSE)),
              #Select the channel to be displayed
              shiny::column(5, shiny::selectInput("Channel", "Channel to display", Channels_in_images, multiple = FALSE))
            ),
            shiny::fluidRow(
              #Select the min point of the image
              shiny::column(4, shiny::sliderInput("Min_Image", "Absolute Black", value = 0, min = 0, max = 100, step = 1)),
              #Select the max point of the image
              shiny::column(4, shiny::sliderInput("Max_Image", "Absolute White", value = 100, min = 0, max = 100, step = 1)),
              shiny::column(2, shinyWidgets::materialSwitch("Change_coords", "Pixel/dist", value = FALSE)),
              shiny::column(2, shiny::conditionalPanel(condition = "input.Change_coords == '1'",
                                                       shiny::textInput("Ratio", "pixel size", value = "1")
              ))
            ),

            shiny::fluidRow(
              #Select the Gamma
              shiny::column(6, shiny::sliderInput("Gamma", "Gamma", value = 0, min = -3, max = +3, step = 0.01)),
              #Select the image resolution
              shiny::column(3, shiny::selectInput("Res", "Image Res", c('Very Low' = 300, Low = 500, Mid = 750, High = 1000, Original = 1400), selected = 500, multiple = FALSE)),
              #Select the equalization
              shiny::column(3, shiny::selectInput("Equalize", "Equalize", c(YES = TRUE, NO = FALSE), selected = FALSE, multiple = FALSE))
            ),

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
            ),
            #Finally add a couple of rows more with extra options and the final result
            shiny::fluidRow(
              shiny::column(3, shiny::selectInput("X_flip", "Flip X image", c(YES = TRUE, NO = FALSE), selected = FALSE, multiple = FALSE)),
              shiny::column(3, shiny::selectInput("Y_flip", "Flip Y image", c(YES = TRUE, NO = FALSE), selected = FALSE, multiple = FALSE)),
              shiny::column(3, shiny::sliderInput("Degrees", "Rotate", min = -180, max = 180, value = 0, step = 90)),
              shiny::column(3, shiny::actionButton("reset", shiny::icon("redo"), label = "Reset selection"))
            ),
            #The UI will be completed with summary tables of the sample
            shiny::fluidRow(
              shiny::column(6, htmltools::p("Final threshold/s is: ", shiny::tableOutput("Final_threshold"))),
              shiny::column(6, htmltools::p("Sample Summary: ", shiny::tableOutput("Summary")))
            ),
            #Also it will include a summary of selected cells
            shiny::fluidRow(
              shiny::column(12, htmltools::p("Selected cells: ", shiny::tableOutput("Cell_selection")))
            )
          ),

          #Set the outcome columns
          shiny::mainPanel(
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
        }')))
      )
    }
    #BUILD THE SERVER
    server <- function(input, output, session){
      #All the reactives to be used
      #Generate a reactive with the real Image name and the channel number
      Photo_name <- shiny::reactive(stringr::str_c(Directory, "/", input$Real_Image_name))
      Channel_index <- shiny::reactive(which(input$Channel == Ordered_Channels))
      Photo_min <- shiny::reactive(input$Min_Image)
      Photo_max <- shiny::reactive(input$Max_Image)
      Photo_gamma <- shiny::reactive(10^input$Gamma)
      Equalize <- shiny::reactive(input$Equalize)
      #Reactive expression to control the graphs
      Case_id <- shiny::reactive(input$Data_Image_name)
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
      ranges <- shiny::reactiveValues(x = NULL, y = NULL)
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
        DATA_thresholded
      })
      #Obtain a list with the histogram and the final thresholds using the dedicated function
      Histo_list <- shiny::reactive(SHINY_Thresholding_summary_function(DATA = Source_DATA(),
                                                                        DATA_Thresholded = Thresholded_DATA(),
                                                                        LOCAL = Local_threshold(),
                                                                        CASE = Case_id())
      )

      #Reactive that imports the photograph (will only be updated when data or photo parameters are updated)
      Photo_reactive <- shiny::reactive({
        #Import the Photo
        Photo <- magick::image_read(Photo_name())[Channel_index()]
        #Perform flip and flop if required
        if(as.logical(X_flip())) Photo <- Photo %>% magick::image_flop()
        if(as.logical(Y_flip())) Photo <- Photo %>% magick::image_flip()
        if(as.numeric(Degrees_rotate() != 0)) Photo <- Photo %>% magick::image_rotate(degrees = as.numeric(Degrees_rotate()))
        #Perform image equalization as requested by user
        if(as.logical(Equalize())) Photo <- Photo %>% magick::image_equalize()
        #Perform image white adjustment
        Photo <- Photo %>%
          magick::image_level(black_point = Photo_min(),
                              white_point = Photo_max(),
                              mid_point = Photo_gamma())


        #Obtain dimensions (these will set the axis limits) BEFORE the resolution is modified
        Photo_Dim <- magick::image_info(Photo)
        #Change resolution if appropiate
        if(as.numeric(Resolution() != 1400)){
          Image_Resolution <-stringr::str_c("X", Resolution())
          Photo <- magick::image_scale(Photo, Image_Resolution)
        }
        return(list(Photo = Photo,
                    Dims = Photo_Dim))
      })

      #Generate the photo plot reactive
      Photo_plot_reactive <- shiny::reactive({
        Photo <- Photo_reactive()[["Photo"]]
        Photo_Dim <- Photo_reactive()[["Dims"]]
        #Return the result as a scaffold ggplot_object
        Photo_plot <- ggplot() +
          annotation_raster(Photo, xmin = 0, xmax = Photo_Dim$width, ymin = 0, ymax = Photo_Dim$height, interpolate = TRUE)+
          scale_x_continuous(limits = c(0, Photo_Dim$width)) +
          scale_y_continuous(limits = c(0, Photo_Dim$height)) +
          coord_cartesian(xlim = ranges$x, ylim = ranges$y, expand = FALSE)+
          theme(axis.title = element_blank(),
                axis.text = element_blank(),
                axis.ticks = element_blank(),
                axis.line = element_blank(),
                panel.background = element_rect(fill = "black"),
                panel.grid = element_blank())
        return(Photo_plot)
      })

      #Print the photo
      output$Photo <- shiny::renderPlot({
        #Try to Plot the result
        try(Photo <- Photo_plot_reactive())

        #If the photo returns an error return the point image
        if(berryFunctions::is.error(Photo)){
          DATA_plot <- Source_DATA() %>% dplyr::filter(Subject_Names == Case_id())
          #Modify pixel values if required
          if(as.logical(Pixel_dist_conversion())){
            DATA_plot$X <- DATA_plot$X * as.numeric(Pixel_dist_ratio())
            DATA_plot$Y <- DATA_plot$Y * as.numeric(Pixel_dist_ratio())
          }

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
                       color = "red", size = 2, hjust = 0.5,
                       label = "UNABLE TO RENDER PHOTO\nUSE ME TO ZOOM IN")
          )

        }
        #else Generate the plot with the photo
        else{
          return(Photo)
        }
      }, res = 300)

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
                theme(axis.line = element_blank(),
                      axis.ticks = element_blank(),
                      axis.text = element_blank(),
                      axis.title = element_blank(),
                      panel.background = element_rect(fill = "black"),
                      legend.text = element_text(size = 10))
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
              guides(colour = "none")
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
        #Try to Import the photo
        try(Photo <- Photo_plot_reactive())

        #If there is an error with the Photo still execute the graph
        if(berryFunctions::is.error(Photo)){
          #Define behavior for 2 levels
          if(Threshold_method() != "Multi_level"){
            #Add specific code for non Multi_level
            Plot_code <- Positive_cells_plot() %>%
              ggplot() +
              cowplot::theme_cowplot()+
              theme(axis.line = element_blank(),
                    axis.ticks = element_blank(),
                    axis.text = element_blank(),
                    axis.title = element_blank(),
                    panel.background = element_rect(fill = "black"),
                    legend.text = element_text(size = 10)) +
              scale_colour_manual("", values = c(alpha("black", 0), "red")) +
              guides(color = "none") +
              scale_x_continuous(limits = c(min(Positive_cells_plot()$X), max(Positive_cells_plot()$X))) +
              scale_y_continuous(limits = c(min(Positive_cells_plot()$Y), max(Positive_cells_plot()$Y))) +
              coord_cartesian(xlim = ranges$x, ylim = ranges$y, expand = FALSE) +
              scale_fill_identity() +
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
            Plot_code <- Positive_cells_plot() %>%
              ggplot() +
              cowplot::theme_cowplot()+
              theme(axis.line = element_blank(),
                    axis.ticks = element_blank(),
                    axis.text = element_blank(),
                    axis.title = element_blank(),
                    panel.background = element_rect(fill = "black"),
                    legend.text = element_text(size = 10)) +
              scale_colour_manual("", values = c(alpha("black", 0), RColorBrewer::brewer.pal(n = Levels()-1, "Reds")))+
              scale_x_continuous(limits = c(min(Positive_cells_plot()$X), max(Positive_cells_plot()$X))) +
              scale_y_continuous(limits = c(min(Positive_cells_plot()$Y), max(Positive_cells_plot()$Y))) +
              coord_cartesian(xlim = ranges$x, ylim = ranges$y, expand = FALSE) +
              scale_fill_identity() +
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
                                              data = Positive_cells_plot())
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
                                              data = Positive_cells_plot())
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
      session$onSessionEnded(function() { shiny::stopApp() })
    }

    #Run the server
    message("Always stop current R execution if you want to continue with your R session")
    shiny::shinyApp(user_interface, server)
  }
