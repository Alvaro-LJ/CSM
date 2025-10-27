#' Launches a shiny App to explore cell phenotyping results
#'
#' After phenotyping, cell phenotypes can be explored interactively proyected on original images.
#'
#'
#' @param DATA A dataframe or tibble containing cell feature data and a column named 'Phenotype' containing cell phenotype labels.
#' @param Directory Character specifying the path to the folder where images are present.
#' @param Ordered_Channels Character vector specifying image channels in their exact order.
#'
#' @details
#' Image parameters control the image and channel to be displayed.
#' Sample summary table is depicted at the bottom of the control panel.
#'
#'Image Panels:
#'\itemize{
#' \item{Upper left: Image being explored. Use it to zoom in.}
#' \item{Upper right: Cell phenotype labels overlayed. Cells of interest can be selected by clicking or by the area selection tools.}
#' \item{Lower left: Cell selection summary. If any cells are selected, the resulting cell phenotype counts be displayed here.}
#' \item{Lower right: Heatmap. If any cells are selected, a heatmap of the relative feature expression will be generated.}
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
#'    EBImage::writeImage(CSM_MiniMultiTiff_test[[Image]], file.path(Input_Dir, names(CSM_MiniMultiTiff_test)[Image]))
#' })
#'
#' #Deploy app----------------------------------------------
#' Phenotyping_evaluator_shiny_app_launcher(
#'     DATA = CSM_Phenotypecell_test,
#'     Directory = Input_Dir,
#'     Ordered_Channels = c("DAPI", "PDL1", "GZMB", "PD1", "CK-EPCAM", "CD8a", "FOXP3")
#')
#'
#' #Remove directories---------------------------------------------------------
#' unlink(Input_Dir, recursive = TRUE)
#' }
#'
#'
#' @export

Phenotyping_evaluator_shiny_app_launcher <-
  function(DATA = NULL,
           Directory = NULL,
           Ordered_Channels = NULL){

    #Check required packages
    {
      if(!requireNamespace("ComplexHeatmap", quietly = TRUE)) stop(
        paste0("ComplexHeatmap Bioconductor package is required to execute the function. Please install using the following code: ",
               expression({
                 if (!require("BiocManager", quietly = TRUE))
                   install.packages("BiocManager")

                 BiocManager::install("ComplexHeatmap")
               })
        )
      )
    }

    DATA <- DATA

    #Test that DATA provided is adequate
    if(!identical(names(DATA)[1:4], c("Cell_no", "X", "Y", "Subject_Names"))) {
      stop("Please generate an appropiate data object using the Data_arrange_function")
    }
    if(!"Phenotype" %in% names(DATA)) stop("DATA must contain a column specifying cell phenotypes")

    #Generate vectors that will guide the sliders and selectors of UI
    Images_in_Data <- unique(DATA$Subject_Names)
    Phenotypes <- unique(DATA$Phenotype)
    Real_Images <- dir(Directory, full.names = FALSE)
    Channels_in_images <- Ordered_Channels
    N_colors <- length(unique(DATA$Phenotype))

    #Set the colors for the plots
    if(N_colors > 22) warning("Currently the color only supports displaying 22 phenotype simulataneously. The same color will be assigned to several phenotypes")
    Color_tibble <- tibble(Phenotype = unique(DATA$Phenotype), Color_code = unname(pals::trubetskoy(n = length(unique(DATA$Phenotype))))) %>%
      dplyr::mutate(Color_code = case_when(is.na(Color_code) ~ "white",
                                           TRUE ~ Color_code))
    DATA <-dplyr::left_join(DATA, Color_tibble, by = "Phenotype")

    #Generate a scaled version of the data if it comes from numeric values
    Variables <- names(DATA)[which(!names(DATA) %in% c("Cell_no", "X", "Y", "Subject_Names", "Phenotype", "Color_code"))]
    if(is.logical(unlist(DATA[Variables]))){
      print("DATA has been thresholded. Heatmap may yield misleading results")
      DATA[Variables] <- DATA[Variables] %>% scale()
    }
    else{DATA[Variables] <- DATA[Variables] %>% scale()}

    #Define quantiles 95 and 5 values of all the dataset for the heatmap
    Min_HEATMAP <- quantile(unlist(DATA[Variables]), 0.025)
    Max_HEATMAP <- quantile(unlist(DATA[Variables]), 0.975)

    #BUILD THE USER INTERFACE
    {
      user_interface <- shiny::fluidPage(

        #Set the title
        shiny::titlePanel("Phenotyping exploration APP"),

        #We want a two panel layout, one in the left containing the input parameters and the output in the right
        shiny::sidebarLayout(
          #Set the first column (which contains the user defined parameters)
          shiny::sidebarPanel(
            #ID and width
            id="sidebar",

            shiny::fluidRow(
              #Select Image to be analyzed
              shiny::column(7, shiny::selectInput("Data_Image_name", "Image from data", sort(Images_in_Data), multiple = FALSE)),
              shiny::column(2, shinyWidgets::materialSwitch("Change_coords", "Pixel/dist", value = FALSE)),
              shiny::column(2, shiny::conditionalPanel(condition = "input.Change_coords == '1'",
                                                       shiny::textInput("Ratio", "pixel size", value = "1")
              ))),
            #Real image to be displayed
            shiny::fluidRow(
              #Select the real image to be displayed
              shiny::column(7, shiny::selectInput("Real_Image_name", "Image to display", sort(Real_Images), multiple = FALSE)),
              #Select the channel to be displayed
              shiny::column(5, shiny::selectInput("Channel", "Channel to display", Channels_in_images, multiple = FALSE))
            ),

            #Add the check box
            shiny::fluidRow(
              shinyWidgets::virtualSelectInput("Checkbox", label = "Phenotypes to display",
                                               choices = sort(Phenotypes),
                                               selected = sort(Phenotypes),
                                               search = TRUE,
                                               multiple = TRUE)
            ),

            #Image parameters
            shiny::fluidRow(
              #Select the min point of the image
              shiny::column(4, shiny::sliderInput("Min_Image", "Absolute Black", value = 0, min = 0, max = 100, step = 1)),
              #Select the max point of the image
              shiny::column(4, shiny::sliderInput("Max_Image", "Absolute White", value = 100, min = 0, max = 100, step = 1)),
              #Select the rotation
              shiny::column(4, shiny::sliderInput("Degrees", "Rotate", min = -180, max = 180, value = 0, step = 90)),

            ),
            shiny::fluidRow(
              #Select the Gamma
              shiny::column(6, shiny::sliderInput("Gamma", "Gamma", value = 0, min = -3, max = +3, step = 0.01)),
              #Select the image resolution
              shiny::column(3, shiny::selectInput("Res", "Image Res", c('Very Low' = 300, Low = 500, Mid = 750, High = 1000, Original = 1400), selected = 500, multiple = FALSE)),
              #Select the equalization
              shiny::column(3, shiny::selectInput("Equalize", "Equalize", c(YES = TRUE, NO = FALSE), selected = FALSE, multiple = FALSE))
            ),

            #Finally add a couple of rows more with extra options and the final result
            shiny::fluidRow(
              shiny::column(4, shiny::selectInput("X_flip", "Flip X image", c(YES = TRUE, NO = FALSE), selected = FALSE, multiple = FALSE)),
              shiny::column(4, shiny::selectInput("Y_flip", "Flip Y image", c(YES = TRUE, NO = FALSE), selected = FALSE, multiple = FALSE)),
              shiny::column(4, shiny::actionButton("reset", shiny::icon("redo"), label = "Reset selection"))
            ),
            #The UI will be completed with summary tables of the sample
            shiny::fluidRow(
              shiny::column(12, htmltools::p("Sample Summary: ", shiny::tableOutput("Summary")))
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
              shiny::column(5, ggiraph::girafeOutput("All_phenotypes"))
            ),
            #Second row will contain the positive cells table and the heatmap
            shiny::fluidRow(
              shiny::column(4, htmltools::p("Selected cells: ", shiny::tableOutput("Cell_selection"))),
              shiny::column(5, shiny::plotOutput("Heatmap"))
            )
          )
        ),
        shiny::tags$head(shiny::tags$style(
          htmltools::HTML('
         #sidebar {
            background-color: #fad2d2;
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

      X_flip <- shiny::reactive(input$X_flip)
      Y_flip <- shiny::reactive(input$Y_flip)
      Degrees_rotate <- shiny::reactive(input$Degrees)
      ranges <- shiny::reactiveValues(x = NULL, y = NULL) #Controls the zoom in
      Pixel_dist_conversion <- shiny::reactive(input$Change_coords)
      Pixel_dist_ratio <- shiny::reactive(input$Ratio)
      Resolution <- shiny::reactive(input$Res)

      #Generate the subject specific DATA that will be used in all plots
      #Control the data source
      Source_DATA <- shiny::reactive({
        Final_DATA <- DATA %>% dplyr::filter(Subject_Names == Case_id())
        #Modify pixel values if required
        if(as.logical(Pixel_dist_conversion())){
          Final_DATA$X <- Final_DATA$X * as.numeric(Pixel_dist_ratio())
          Final_DATA$Y <- Final_DATA$Y * as.numeric(Pixel_dist_ratio())
        }
        #Return the final data
        return(Final_DATA)
      })

      #Control the checkbox output
      Checkbox_output <- shiny::reactive(input$Checkbox)

      #Reactive that imports the photograph
      Photo_reactive <- shiny::reactive({
        #Import the Photo
        Photo <- magick::image_read(Photo_name())[Channel_index()]
        #Perform flip and flop and rotate if required
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


        ##Obtain dimensions (these will set the axis limits) BEFORE the resolution is modified
        Photo_Dim <- magick::image_info(Photo)
        #Change resolution if appropiate
        if(as.numeric(Resolution() != 1400)){
          Image_Resolution <-stringr::str_c("X", Resolution())
          Photo <- magick::image_scale(Photo, Image_Resolution)
        }
        return(list(Photo = Photo,
                    Dims = Photo_Dim))
      })

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
        #Plot the result
        try(Photo <- Photo_plot_reactive())

        #If the photo returns an error return the point image
        if(berryFunctions::is.error(Photo)){
          Final_DATA <- Source_DATA()

          Photo_plot <-
            Final_DATA %>%
            ggplot() +
            geom_point(aes(x = X, y = Y),
                       color = "white",
                       size = 2.5) +
            scale_x_continuous(limits = c(min(Final_DATA$X), max(Final_DATA$X))) +
            scale_y_continuous(limits = c(min(Final_DATA$Y), max(Final_DATA$Y))) +
            coord_cartesian(xlim = ranges$x, ylim = ranges$y, expand = FALSE)+
            theme(axis.title = element_blank(),
                  axis.text = element_blank(),
                  axis.ticks = element_blank(),
                  axis.line = element_blank(),
                  panel.background = element_rect(fill = "black"),
                  panel.grid = element_blank()) +
            annotate("text", x = quantile(Final_DATA$X, 0.5), y = quantile(Final_DATA$Y, 0.5),
                     color = "red", size = 2, hjust = 0.5,
                     label = "UNABLE TO RENDER PHOTO\nUSE ME TO ZOOM IN")
          return(Photo_plot)
        }

        else{
          #Generate the plot
          return(Photo)
        }

      }, res = 300)

      #Control the zoom of the Photo and the rest of the graphs
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

      #All cells phenotype
      Cell_Phenotype_plot <-
        shiny::reactive({
          #Get sample data and change the axis to coordinate all three plots displayed
          Final_DATA <- Source_DATA()
          #Get the min and max of the photo before selecting the required phenotypes
          plot_x_min <- min(Final_DATA$X)
          plot_x_max <- max(Final_DATA$X)
          plot_y_min <- min(Final_DATA$Y)
          plot_y_max <- max(Final_DATA$Y)

          #Filter phenotypes in checkbox
          Final_DATA <- Final_DATA %>% dplyr::filter(Phenotype %in% Checkbox_output())

          #Import the photo
          try(Photo <- Photo_plot_reactive())

          #If image not available plot without image
          if(berryFunctions::is.error(Photo)){
            return(
              #Generate the final plot
              Final_DATA %>%
                ggplot() +
                scale_color_identity()+
                ggiraph::geom_point_interactive(aes(x = X, y = Y, group = Phenotype, color = Color_code,
                                                    data_id = Cell_no,
                                                    tooltip = stringr::str_c(as.character(Cell_no)," Type = ", as.character(Phenotype)),
                                                    hover_nearest = FALSE),
                                                size = 2) +
                cowplot::theme_cowplot()+
                guides(color = "none") +
                theme(axis.line = element_blank(),
                      axis.ticks = element_blank(),
                      axis.text = element_blank(),
                      axis.title = element_blank(),
                      panel.background = element_rect(fill = "black"),
                      legend.position = "bottom",
                      legend.text = element_text(size = 10)) +
                scale_x_continuous(limits = c(plot_x_min, plot_x_max)) +
                scale_y_continuous(limits = c(plot_y_min,  plot_y_max)) +
                coord_cartesian(xlim = ranges$x, ylim = ranges$y, expand = FALSE)
            )
          }
          #If not proceed as usual
          else{
            return(
              #Generate the final plot
              Photo +
                scale_color_identity()+
                ggiraph::geom_point_interactive(aes(x = X, y = Y, group = Phenotype, color = Color_code,
                                                    data_id = Cell_no,
                                                    tooltip = stringr::str_c(as.character(Cell_no)," Type = ", as.character(Phenotype)),
                                                    hover_nearest = FALSE),
                                                size = 2,
                                                data = Final_DATA) +
                cowplot::theme_cowplot()+
                guides(color = "none") +
                theme(axis.title = element_blank(),
                      axis.text = element_blank(),
                      axis.ticks = element_blank(),
                      axis.line = element_blank(),
                      panel.background = element_rect(fill = "black"),
                      panel.grid = element_blank())
            )
          }
        })
      #Send plot to UI
      output$All_phenotypes <- ggiraph::renderGirafe({
        plot <- ggiraph::girafe(code = print(Cell_Phenotype_plot()),
                                options = list(
                                  ggiraph::opts_hover(css = "stroke:black;cursor:pointer;", reactive = TRUE),
                                  ggiraph::opts_selection(type = "multiple", css = "fill:#FF3333;stroke:black;")
                                )
        )
        plot
      })

      #The sample summary (phenotype count table)
      output$Summary <- gt::render_gt({
        Sample_DATA <- Source_DATA()
        Sample_table <- Sample_DATA %>% dplyr::count(Phenotype)

        Sample_table <- Sample_table %>%dplyr::mutate('%' = round(n/nrow(Source_DATA())*100, 2))
        Sample_table <-dplyr::left_join(Sample_table, Color_tibble, by = "Phenotype")
        Sample_table <- Sample_table %>% dplyr::arrange(desc(n))

        Sample_table %>% gt::gt() %>%
          gt::opt_table_outline(style = "solid", color = "black") %>%
          gt::tab_style(
            style = list(
              gt::cell_text(color = "black", size = gt::px(20), font = "Calibri", align = "center"),
              gt::cell_fill(color = gt::from_column(column = "Color_code"))
            ),
            locations = gt::cells_body(columns = c(1:4))
          ) %>%
          gt::tab_style(
            style = list(
              gt::cell_text(color = "black", size = gt::px(22), weight = "bold", font = "Calibri", align = "center")
            ),
            locations = gt::cells_column_labels()
          ) %>%
          gt::cols_width(
            1 ~ gt::pct(50),
            2	~ gt::pct(30),
            3	~ gt::pct(20),
            4 ~ gt::pct(0)
          )
      })

      #Selected cells and reset button
      selected_cells_phenotype <- shiny::reactive(input$All_phenotypes_selected)
      #What to do in case the user hits the reset button
      shiny::observeEvent(input$reset, {
        session$sendCustomMessage(type = 'All_phenotypes_set', message = character(0))
      })
      #Generate the output tibble for the selected cells
      output$Cell_selection <- shiny::renderTable({
        #Get cell_no in both plots and remove duplicates
        Cells <- unique(selected_cells_phenotype())
        #Get Marker data and threshold data an start processing
        Selected_Cells <- Source_DATA() %>% dplyr::filter(Cell_no %in% Cells) %>% dplyr::select(Cell_no, X, Y, Phenotype)
        Final_cells <- Selected_Cells %>% dplyr::count(Phenotype)
        dplyr::bind_rows(Final_cells, tibble(Phenotype = "TOTAL", n = sum(Final_cells$n)))
      })

      #Generate the heatmap
      output$Heatmap <- shiny::renderPlot({
        #Get cell_no in both plots and remove duplicates
        Cells <- unique(selected_cells_phenotype())

        #Get Marker data and threshold data an start processing
        Selected_Cells <- Source_DATA() %>% dplyr::filter(Cell_no %in% Cells)

        #If user has not selected cells then print a fixed value matrix
        if(nrow(Selected_Cells) == 0){
          Matrix <- matrix(0, 1, 1)
          ComplexHeatmap::Heatmap(Matrix)
        }
        #If cells have been selected print the actual matrix
        else{
          #Prepare the matrix
          HEATMAP_MATRIX <- Selected_Cells %>% dplyr::select(-c(1:4), -Phenotype, -Color_code)
          HEATMAP_MATRIX <- as.matrix(HEATMAP_MATRIX)
          row.names(HEATMAP_MATRIX) <- Selected_Cells[["Cell_no"]]

          #Prepare the color function for the matrix itself and the rowside annotation
          col_fun <- circlize::colorRamp2(c(Min_HEATMAP, 0, Max_HEATMAP), c("#0000ff", "white", "#ff0000"))
          row_color_code <- Color_tibble[[2]]
          names(row_color_code) <- Color_tibble[[1]]


          #Prepare the side annotation
          Side_annotation <- ComplexHeatmap::rowAnnotation(
            Phenotype = Selected_Cells[["Phenotype"]],
            col = list(Phenotype = row_color_code),
            show_annotation_name = F,
            show_legend = F)

          ComplexHeatmap::Heatmap(HEATMAP_MATRIX,
                                  col = col_fun,
                                  cluster_rows = T,
                                  show_row_names = FALSE,
                                  show_heatmap_legend = F,
                                  left_annotation = Side_annotation,
                                  show_row_dend = F,
                                  show_column_dend = F,
                                  column_names_side = "top",
                                  column_names_rot = 0,
                                  column_names_centered = T,
                                  column_names_gp = grid::gpar(fontsize = 10),
                                  border = T
          )
        }
      })

      #If browser is closed end the app
      session$onSessionEnded(function() { shiny::stopApp() })
    }

    #Run the server
    message("Always stop current R execution if you want to continue with your R session")
    shiny::shinyApp(user_interface, server)
  }
