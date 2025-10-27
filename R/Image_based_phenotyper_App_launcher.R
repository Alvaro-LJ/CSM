#' Launches a shiny APP that can generate cell phenotyping models based on user expertise
#'
#' Launches an APP to interactively explore cell segmented images. The user can manually assign cell phenotype labels to cells and train a model.
#' The model can then be used in [Model_cell_phenotyper()] function to generalize it to all the cells in a dataset.
#'
#' @param DATA A dataframe or tibble containing cell feature data.
#' @param Directory Character specifying the path to the folder where images are present.
#' @param Ordered_Channels Character vector specifying image channels in their exact order.
#' @param Channels_to_keep Character vector indicating the channels to be kept in the analysis.
#'
#' @seealso [Model_cell_phenotyper()]
#'
#' @details
#' Image setting allow the user to control the image channel display
#'
#' Cell label assigner has several buttons:
#' \itemize{
#' \item{Number of selected cells: Number of currently selected cells is shown. Cells can be selected or deselected using the upper or bottom right panels.}
#' \item{-Remove selected cells: Eliminates currently selected cells from the training dataset.}
#' \item{-Remove all cells: Eliminates all cells from the dataset.}
#' \item{Cell label: Entry box that allows the user to type a cell label.}
#' \item{+Assign labels: Assigns the label in the "Cell label" box to the currently selected cells and sends the cells to the training dataset.}
#' \item{Remove labels: Removes any cells from the training dataset matching the label in the "Cell label" box.}
#' }
#'
#' Model settings allows the user to control how the model will be calculated. If spatial interaction features are required the user must
#' specify how neighbors are defined.
#'
#'Relevan buttons in the lower area of the control panel
#'\itemize{
#'\item{Fit button computes the model with the current model settings and test dataset}
#'\item{Test button tests the model in the current displayed image}
#'\item{Download button saves the current active model in the Global environment}
#'}
#'
#' Upper left panel: Displays the image (use it to zoom in and out)
#'
#' Upper right panel: Allows assigning cells for the training dataset. Select cells by clicking or use the area selection tools provided. Cell labels will be updated as they are assigned.
#'
#' Lower left: Shows features of the training dataset and the results of the test dataset if tests have been performed.
#'
#' Lower right: If a test has been performed, it shows the model results. Missclassified cells can be selected and labels can be re-assigned to refine the model.
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
#'Image_based_phenotyper_App_launcher(
#'    DATA = CSM_Arrangedcellfeaturedata_test,
#'    Directory = Input_Dir,
#'    Ordered_Channels = c("DAPI", "PDL1", "GZMB", "PD1", "CK-EPCAM", "CD8a", "FOXP3"),
#'    Channels_to_keep = c("DAPI", "PDL1", "GZMB", "PD1", "CK-EPCAM", "CD8a", "FOXP3")
#')
#'
#'#Remove directories---------------------------------------------------------
#' unlink(Input_Dir, recursive = TRUE)
#' }
#'
#'
#' @export

Image_based_phenotyper_App_launcher <-
  function(DATA = NULL,
           Directory = NULL,
           Ordered_Channels = NULL,
           Channels_to_keep = NULL
  ){

    #Check suggested packages
    {
      if(!requireNamespace("Matrix", quietly = FALSE)) stop(
        paste0("Matrix CRAN package is required to execute the function. Please install using the following code: ",
               expression(install.packages("Matrix")))
      )
      if(!requireNamespace("rtree", quietly = FALSE)) stop(
        paste0("rtree GitHub package is required to execute the function. Please install using the following code: ",
               expression(remotes::install_github("akoyabio/rtree")))
      )
      if(!requireNamespace("magick", quietly = FALSE)) stop(
        paste0("magick CRAN package is required to execute the function. Please install using the following code: ",
               expression(install.packages("magick")))
      )
      if(!requireNamespace("tidymodels", quietly = FALSE)) stop(
        paste0("tidymodels CRAN package is required to execute the function. Please install using the following code: ",
               expression(install.packages("tidymodels")))
      )
      if(!requireNamespace("randomForest", quietly = FALSE)) stop(
        paste0("randomForest CRAN package is required to execute the function. Please install using the following code: ",
               expression(install.packages("randomForest")))
      )
      if(!requireNamespace("xgboost", quietly = FALSE)) stop(
        paste0("xgboost CRAN package is required to execute the function. Please install using the following code: ",
               expression(install.packages("xgboost")))
      )
      if(!requireNamespace("brulee", quietly = FALSE)) stop(
        paste0("brulee CRAN package is required to execute the function. Please install using the following code: ",
               expression(install.packages("brulee")))
      )
    }

    on.exit({
      future::plan("future::sequential")
      gc()
    })

    #Basic argument checker
    Argument_checker <- c(Empty_directory = length(dir(Directory)) >= 1,
                          Channels_OK = all(Channels_to_keep %in% Ordered_Channels)
    )

    Stop_messages <- c(Empty_directory = "No files found at the directory provided. Please check out the path.",
                       Channels_OK =stringr::str_c(
                         "The following channels are not present the channel names provided: ",
                         stringr::str_c(Channels_to_keep[!(Channels_to_keep %in% Ordered_Channels)], collapse = ", "),
                         sep = "")
    )
    if(!all(Argument_checker)){
      stop(cat(Stop_messages[!Argument_checker],
               fill = sum(!Argument_checker)))
    }
    #Check data is adequately formatted
    if(!identical(names(DATA)[1:4],  c("Cell_no", "X", "Y", "Subject_Names"))) stop("DATA provided should have an adecuate format")

    #Generate look-up table
    Image_names <- dir(Directory, full.names = FALSE)
    Directory_path <- dir(Directory, full.names = TRUE)
    DATA <- DATA

    #Remove any features in the data that are not numeric
    Numeric_features <-purrr::map_lgl(DATA[-c(1:4)], ~is.numeric(.))
    if(any(!Numeric_features)){
      print(paste0("The following non-numeric features will be removed from DATA: ", stringr::str_c(names(DATA)[-c(1:4)][!Numeric_features], collapse = ", ")))
      DATA <- DATA %>% dplyr::select(-all_of(stringr::str_c(names(DATA)[-c(1:4)][!Numeric_features])))
    }

    print(paste0(length(Image_names), " files found in Directory. Finding matches between file names and Subject_Names in DATA..."))

    Look_up_table <-
      purrr::map_dfr(Image_names, function(Name){
        Distance_vector <- as.double(adist(Name, unique(DATA$Subject_Names), fixed = TRUE, ignore.case = TRUE))
        names(Distance_vector) <-  unique(DATA$Subject_Names)
        Distance_vector <- sort(Distance_vector)

        tibble(Image_name = Name,
               Subject_Names = names(Distance_vector)[1],
               Match = Distance_vector[1] == 0,
               Distance = Distance_vector[1])

      })
    Look_up_table$Path <- Directory_path

    if(any(!Look_up_table$Match)){
      print("The following file names have not been matched with Subject_Names in DATA. Please select a distance to perform partial matching")
      View(Look_up_table %>% dplyr::filter(!Match))
      Answer <- menu(1:max(Look_up_table$Distance))


      Look_up_table <- Look_up_table %>% dplyr::filter(Distance <= Answer)
      print(paste0(length(Image_names) - nrow(Look_up_table), " out of ", length(Image_names), " images have been removed from the phenotyping process"))
    }

    ##APP UI
    #BUILD THE USER INTERFACE
    {
      user_interface <- shiny::fluidPage(

        #Set the title
        shiny::titlePanel("Supervised Phenotyper"),

        #We want a two panel layout, one in the left containing the input parameters and the output in the right
        shiny::sidebarLayout(
          #Set the first column (which contains the user defined parameters)
          shiny::sidebarPanel(
            #ID and width
            id="sidebar",

            #IMAGE PARAMETERS
            shiny::h4("Image settings"),
            shiny::fluidRow(
              shiny::column(6, shiny::selectInput("Image_name", "Image", sort(Look_up_table$Subject_Names), multiple = FALSE)),
              shiny::column(2, shiny::selectInput("Res", "Image Res", c('Very Low' = 300, Low = 500, Mid = 750, High = 1000, Original = 1400), selected = 500, multiple = FALSE)),
              shiny::column(2, shinyWidgets::materialSwitch("Change_coords", "Pixel/dist", value = FALSE)),
              shiny::column(2, shiny::conditionalPanel(condition = "input.Change_coords == '1'",
                                                       shiny::textInput("Ratio", "pixel size", value = "1")
              ))),
            shiny::fluidRow(
              shiny::column(4, shiny::selectInput("Channel", "Channel to display", Channels_to_keep, multiple = FALSE)),
              shiny::column(3, shiny::selectInput("X_flip", "Flip X image", c(YES = TRUE, NO = FALSE), selected = FALSE, multiple = FALSE)),
              shiny::column(3, shiny::selectInput("Y_flip", "Flip Y image", c(YES = TRUE, NO = FALSE), selected = FALSE, multiple = FALSE)),
              shiny::column(2, shiny::selectInput("Equalize", "Equalize", c(YES = TRUE, NO = FALSE), selected = FALSE, multiple = FALSE))

            ),

            shiny::fluidRow(
              #Rotate if necessary
              shiny::column(3, shiny::sliderInput("Degrees", "Rotate", min = -180, max = 180, value = 0, step = 90)),
              #Select the min point of the image
              shiny::column(3, shiny::sliderInput("Min_Image", "Absolute Black", value = 0, min = 0, max = 100, step = 1)),
              shiny::column(3, shiny::sliderInput("Max_Image", "Absolute White", value = 100, min = 0, max = 100, step = 1)),
              shiny::column(3, shiny::sliderInput("Gamma", "Gamma", value = 0, min = -3, max = +3, step = 0.01))

            ),

            #IMAGE PARAMETERS
            shiny::h4("Cell label assigner"),
            shiny::fluidRow(
              shiny::column(1, shiny::textOutput("N_cells_selected")),
              shiny::column(3, shiny::actionButton("reset", shiny::icon("redo"), label = "Reset selection")),
              shiny::column(4, shiny::actionButton("Remove_Cells", shiny::icon("minus", library = "font-awesome"), label = "Remove selected cells")),
              shiny::column(4, shiny::actionButton("Remove_All", shiny::icon("minus", library = "font-awesome"), label = "Remove all cells"))
            ),

            shiny::fluidRow(
              shiny::column(6, shiny::textInput("Label_name", label = "Cell label", value = "Cell_type_1")),
              shiny::column(3, shiny::actionButton("Add_label", shiny::icon("plus", library = "font-awesome"), label = "Assign labels")),
              shiny::column(3, shiny::actionButton("Remove_label", shiny::icon("minus", library = "font-awesome"), label = "Remove labels"))
            ),

            #Model PARAMETERS
            shiny::h4("Model settings"),
            shiny::fluidRow(
              shiny::column(3, shiny::selectInput("Method", "Method", c("Random forest", "XG boost", "NNET"), multiple = FALSE)),
              shiny::column(5, shinyWidgets::virtualSelectInput("Model_vars", label = "Variables included",
                                                                choices = sort(names(DATA)[-c(1:4)]),
                                                                selected = sort(names(DATA)[-c(1:4)]),
                                                                search = TRUE,
                                                                multiple = TRUE)),
              shiny::column(2, shiny::numericInput("Threshold", "Thresh", value = 0.5, min = 0.01, max = 1, step = 0.01)),
              shiny::column(2, shinyWidgets::materialSwitch("GO_spatial", "Spatial", value = FALSE))
            ),

            #Random Forest
            shiny::conditionalPanel(
              condition = "input.Method == 'Random forest'",
              shiny::fluidRow(
                shiny::column(3, shiny::numericInput("RF_mtry", "% Features", value = 1, min = 0.01, max = 1, step = 0.01)),
                shiny::column(3, shiny::numericInput("RF_trees", "N Trees", value = 100, min = 1, max = NA, step = 1))
              )
            ),

            #XGB
            shiny::conditionalPanel(
              condition = "input.Method == 'XG boost'",
              shiny::fluidRow(
                shiny::column(3, shiny::numericInput("XG_mtry", "% Features", value = 1, min = 0.01, max = 1, step = 0.01)),
                shiny::column(3, shiny::numericInput("XG_sample_size", "% Cells", value = 1, min = 0.01, max = 1, step = 0.01)),
                shiny::column(3, shiny::numericInput("XG_trees", "N Trees", value = 100, min = 1, max = NA, step = 1)),
                shiny::column(3, shiny::numericInput("XG_tree_depth", "Depth", value = 10, min = 1, max = NA, step = 1))
              )
            ),

            #NNET
            shiny::conditionalPanel(
              condition = "input.Method == 'NNET'",
              shiny::fluidRow(
                shiny::column(3, shiny::numericInput("Hidden", "Hidden units", value = 1, min = 1, max = NA, step = 1)),
                shiny::column(3, shiny::numericInput("Layers", "Hidden layers", value = 1, min = 1, max = NA, step = 1)),
                shiny::column(3, shiny::numericInput("Epochs", "Epochs", value = 1, min = 1, max = NA, step = 1)),
                shiny::column(3, shiny::numericInput("Penalty", "Penalty", value = 0.001, min = 0.001, max = NA, step = 0.001))
              )
            ),

            #SPATIAL PARAMETERS
            shiny::conditionalPanel(
              condition = "input.GO_spatial == '1'",
              shiny::h4("Spatial settings"),
              shiny::fluidRow(
                shiny::column(4, shiny::selectInput("Neighbor_strategy", "Strategy", c("Number", "Distance", "Both"), multiple = FALSE)),
                shiny::column(4, shiny::selectInput("Message_strategy", "Summary Strategy", c("Averaging", "Sum"), multiple = FALSE)),
                shiny::column(4, shiny::selectInput("Weighting_Strategy", "Weighting", c("None", "Proximity"), multiple = FALSE))
              ),
              shiny::fluidRow(
                shiny::column(4, shiny::numericInput("N_neighbors", "N Neighbors", value = 1, min = 1, max = NA, step = 1)),
                shiny::column(4, shiny::numericInput("Max_dist_allowed", "Max Dist", value = 100, min = 0.001, max = NA, step = 1)),
                shiny::column(4, shiny::numericInput("N_cores", "Cores", value = 1, min = 1, max = NA, step = 1))
              )
            ),

            shiny::fluidRow(
              shiny::column(2, shiny::actionButton("Fit_model", shiny::icon("bolt-lightning"), label = "Fit")),
              shiny::column(2, shiny::actionButton("Test_model", shiny::icon("square-check", library = "font-awesome"), label = "Test")),
              shiny::column(2, shiny::actionButton("Download_model", shiny::icon("download"), label = "Download"))
            ),

            #Active model
            shiny::fluidRow(
              shiny::column(12, shiny::verbatimTextOutput("Active_Model"))
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
              shiny::column(5, ggiraph::girafeOutput("Cells_unassigned"))
            ),
            #Second row will contain the positive cells table and the heatmap
            shiny::fluidRow(
              shiny::column(5, shiny::verbatimTextOutput("Training_selection_summary")),
              shiny::column(5, ggiraph::girafeOutput("Cells_assigned"))
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

    server <- function(input, output, session){
      #All the reactives to be used
      #Generate a reactive with the real Image name and the channel number
      Photo_name <- shiny::reactive(as.character((Look_up_table %>% dplyr::filter(Subject_Names == input$Image_name))[1,5]))
      Channel_index <- shiny::reactive(which(input$Channel == Ordered_Channels))
      Photo_min <- shiny::reactive(input$Min_Image)
      Photo_max <- shiny::reactive(input$Max_Image)
      Photo_gamma <- shiny::reactive(10^input$Gamma)
      Equalize <- shiny::reactive(input$Equalize)
      #Reactive expression to control the graphs
      Case_id <- shiny::reactive(input$Image_name)
      X_flip <- shiny::reactive(input$X_flip)
      Y_flip <- shiny::reactive(input$Y_flip)
      Degrees_rotate <- shiny::reactive(input$Degrees)
      ranges <- shiny::reactiveValues(x = NULL, y = NULL) #Controls the zoom in
      Pixel_dist_conversion <- shiny::reactive(input$Change_coords)
      Pixel_dist_ratio <- shiny::reactive(input$Ratio)
      Resolution <- shiny::reactive(input$Res)

      #The reactive that controls the final output
      Result_list <- shiny::reactiveValues(Training_Dataset = NULL,
                                           Model_Param = NULL,
                                           Model = NULL,
                                           Test_Dataset = NULL)

      #Generate a dynamic look up table to match labels and colors
      Color_lookuptable <- shiny::reactive({
        tibble(Label = c("Unassigned", unique(Result_list$Training_Dataset$Label)),
               Color = c("black", unname(pals::polychrome(n = length(unique(Result_list$Training_Dataset$Label))))))
      })

      #Selected cells using the unassigned plot
      selected_cells_unassigned <- shiny::reactive(input$Cells_unassigned_selected)
      selected_cells_assigned <- shiny::reactive(input$Cells_assigned_selected)


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

      #Reactive that imports the photograph
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
        return(Photo)
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

      #Generate the cell unassigned plot
      Cell_Unassigned_plot <-
        shiny::reactive({
          #Get the data
          Final_DATA <- Source_DATA()
          Final_DATA <- Final_DATA %>% dplyr::select(1:4)

          #Add the Unassigned label
          Final_DATA$Label <- "Unassigned"

          #Modify the label according to the training dataset
          if(!is.null(Result_list$Training_Dataset)){
            Final_DATA <- Final_DATA %>% dplyr::filter(!Cell_no %in% Result_list$Training_Dataset$Cell_no)
            Cells_to_bind <- Result_list$Training_Dataset %>% dplyr::filter(Subject_Names == Case_id())
            Final_DATA <- dplyr::bind_rows(Final_DATA, Cells_to_bind)
          }

          #Add the color code
          Final_DATA <-dplyr::left_join(Final_DATA, Color_lookuptable(), by = "Label")

          #Import the photo
          try(Photo <- Photo_plot_reactive())
          return(
            #Generate the final plot
            Photo +
              ggiraph::geom_point_interactive(aes(x = X, y = Y, color = Color,
                                                  data_id = Cell_no,
                                                  tooltip = stringr::str_c(Cell_no, " ", Label, sep = ""),
                                                  hover_nearest = FALSE),
                                              size = 3,
                                              data = Final_DATA) +
              scale_color_identity() +
              cowplot::theme_cowplot()+
              guides(color = "none") +
              theme(axis.title = element_blank(),
                    axis.text = element_blank(),
                    axis.ticks = element_blank(),
                    axis.line = element_blank(),
                    panel.background = element_rect(fill = "black"),
                    panel.grid = element_blank())
          )

        })
      #Send plot to UI
      output$Cells_unassigned <- ggiraph::renderGirafe({
        plot <- ggiraph::girafe(code = print(Cell_Unassigned_plot()),
                                options = list(
                                  ggiraph::opts_hover(css = "stroke:black;cursor:pointer;", reactive = TRUE),
                                  ggiraph::opts_selection(type = "multiple", css = "fill:#FF3333;stroke:black;")
                                )
        )
        return(plot)
      })

      #Generate the test dataset plot
      Test_Dataset_plot <-
        shiny::reactive({
          #If no test dataset or the test dataset does not match the case id, then abort
          if(any(is.null(Result_list$Test_Dataset), unique(Result_list$Test_Dataset$Subject_Names) != Case_id())){
            return(ggplot())
          }

          #If test dataset is equal to the case id then proceed
          Final_DATA <- Result_list$Test_Dataset

          #Add the Unassigned label
          Final_DATA$Label[Final_DATA$Probability < Result_list$Model_Param$Model_threshold] <- "Unassigned"

          #Add the color code
          Final_DATA <-dplyr::left_join(Final_DATA, Color_lookuptable(), by = "Label")

          #Import the photo
          try(Photo <- Photo_plot_reactive())
          return(
            #Generate the final plot
            Photo +
              ggiraph::geom_point_interactive(aes(x = X, y = Y, color = Color,
                                                  data_id = Cell_no,
                                                  tooltip = stringr::str_c(Cell_no, " ", Label, ". Prob = ", round(Probability, 2), sep = ""),
                                                  hover_nearest = FALSE),
                                              size = 3,
                                              data = Final_DATA) +
              scale_color_identity() +
              cowplot::theme_cowplot()+
              guides(color = "none") +
              theme(axis.title = element_blank(),
                    axis.text = element_blank(),
                    axis.ticks = element_blank(),
                    axis.line = element_blank(),
                    panel.background = element_rect(fill = "black"),
                    panel.grid = element_blank())
          )

        })
      #Send plot to UI
      output$Cells_assigned <- ggiraph::renderGirafe({
        plot <- ggiraph::girafe(code = print(Test_Dataset_plot()),
                                options = list(
                                  ggiraph::opts_hover(css = "stroke:black;cursor:pointer;", reactive = TRUE),
                                  ggiraph::opts_selection(type = "multiple", css = "fill:#FF3333;stroke:black;")
                                )
        )
        return(plot)
      })

      #The total number of cells
      output$N_cells_selected <- shiny::renderText({
        length(unique(c(selected_cells_unassigned(), selected_cells_assigned())))
      })

      #The summary of the Training dataset and the test dataset
      output$Training_selection_summary <- shiny::renderPrint({
        if(is.null(Result_list$Training_Dataset)) print("No cells selected for Training Dataset")
        else{
          #If is null the test dataset print only training datatset
          if(is.null(Result_list$Test_Dataset)){
            print(list('Training Cells by label' = Result_list$Training_Dataset %>% dplyr::count(Label),
                       'Training Cells by image' = Result_list$Training_Dataset %>% dplyr::count(Subject_Names))
            )
          }
          else{
            print(list('Training Cells by label' = Result_list$Training_Dataset %>% dplyr::count(Label),
                       'Training Cells by image' = Result_list$Training_Dataset %>% dplyr::count(Subject_Names),
                       'Test dataset cells by label' = Result_list$Test_Dataset %>% dplyr::count(Label)
            )
            )
          }

        }
      })

      #Print the active model if available
      output$Active_Model <- shiny::renderPrint({
        #If no model is present print a message
        if(is.null(Result_list$Model)) print("No model has been created")

        #Else print the model parameters of the active model
        else{
          print(Result_list$Model_Param)
        }
      })

      #What to do when the user hits the following buttons:
      #What to do in case the user hits the reset button
      shiny::observeEvent(input$reset, {
        session$sendCustomMessage(type = 'Cells_unassigned_set', message = character(0))
        session$sendCustomMessage(type = 'Cells_assigned_set', message = character(0))
      })

      #If the user hits assign cells button
      shiny::observeEvent(input$Add_label,
                          {
                            #IF label is not provided
                            if(input$Label_name == "") shiny::showModal(modalDialog(
                              paste0("Cell label must be provided"),
                              easyClose = TRUE,
                              footer = NULL)
                            )

                            else if(input$Label_name == "Unassigned") shiny::showModal(modalDialog(
                              paste0("Unassigned is not a valid cell label"),
                              easyClose = TRUE,
                              footer = NULL)
                            )

                            #IF no cells have been selected
                            else if(length(unique(c(selected_cells_unassigned(), selected_cells_assigned()))) == 0) shiny::showModal(modalDialog(
                              paste0("Cells must have been selected"),
                              easyClose = TRUE,
                              footer = NULL)
                            )

                            #If everything OK proceed with label assignment
                            else{
                              if(all(!is.null(Result_list$Training_Dataset), any(unique(c(selected_cells_unassigned(), selected_cells_assigned())) %in% Result_list$Training_Dataset$Cell_no))){
                                shiny::showModal(modalDialog(
                                  paste0("Some cells are already present in the Training dataset. Their labels will be overwritten"),
                                  easyClose = TRUE,
                                  footer = NULL))
                              }

                              #If training dataset is NULL then create a new one
                              if(is.null(Result_list$Training_Dataset)){
                                Interim <- DATA %>% dplyr::select(1:4) %>% dplyr::filter(Cell_no %in% unique(c(selected_cells_unassigned(), selected_cells_assigned()))) %>% dplyr::mutate(Label = input$Label_name)
                                Result_list$Training_Dataset <- Interim
                              }

                              #If training dataset has already been created bind the new rows
                              if(!is.null(Result_list$Training_Dataset)){
                                #First we will remove cells that are already present in the Training dataset and need to be overwritten
                                Interim <- Result_list$Training_Dataset %>% dplyr::filter(!Cell_no %in% unique(c(selected_cells_unassigned(), selected_cells_assigned())))
                                To_bind <- DATA %>% dplyr::select(1:4) %>% dplyr::filter(Cell_no %in% unique(c(selected_cells_unassigned(), selected_cells_assigned()))) %>% dplyr::mutate(Label = input$Label_name)
                                Result_list$Training_Dataset <- dplyr::bind_rows(Interim, To_bind)
                              }
                            }
                          },
                          ignoreInit = TRUE)

      #If the user hits the remove selected cells button
      shiny::observeEvent(input$Remove_Cells,
                          {
                            #IF no cells have been selected
                            if(length(unique(c(selected_cells_unassigned(), selected_cells_assigned()))) == 0) shiny::showModal(modalDialog(
                              paste0("Cells must have been selected"),
                              easyClose = TRUE,
                              footer = NULL)
                            )

                            #IF Training dataset has not been created
                            else if(is.null(Result_list$Training_Dataset)) shiny::showModal(modalDialog(
                              paste0("Training dataset has not been created yet"),
                              easyClose = TRUE,
                              footer = NULL)
                            )

                            #If everything OK proceed with removal of cells from the dataset
                            else{
                              Result_list$Training_Dataset <- Result_list$Training_Dataset %>% dplyr::filter(!Cell_no %in% unique(c(selected_cells_unassigned(), selected_cells_assigned())))
                            }
                          },
                          ignoreInit = TRUE)

      #If the user hits the remove all cells button
      shiny::observeEvent(input$Remove_All,
                          {
                            #IF Training dataset has not been created
                            if(is.null(Result_list$Training_Dataset)) shiny::showModal(modalDialog(
                              paste0("Training dataset has not been created yet"),
                              easyClose = TRUE,
                              footer = NULL)
                            )

                            else{
                              #Show an alert
                              shinyalert::shinyalert(title = "WARNING!",
                                                     text = "All cells from the training dataset will be removed. Proceed?",
                                                     type = "warning",
                                                     closeOnEsc = TRUE,
                                                     closeOnClickOutside = TRUE,
                                                     showCancelButton = TRUE,
                                                     showConfirmButton = TRUE,
                                                     confirmButtonText = "Proceed",
                                                     cancelButtonText = "Cancel",
                                                     callbackR = function() Result_list$Training_Dataset <- NULL
                              )
                            }
                          },
                          ignoreInit = TRUE)

      #If the user hits the remove all labels button
      shiny::observeEvent(input$Remove_label,
                          {
                            #IF Training dataset has not been created
                            if(is.null(Result_list$Training_Dataset)) shiny::showModal(modalDialog(
                              paste0("Training dataset has not been created yet"),
                              easyClose = TRUE,
                              footer = NULL)
                            )

                            #IF label is not provided
                            else if(input$Label_name == "") shiny::showModal(modalDialog(
                              paste0("Cell label must be provided"),
                              easyClose = TRUE,
                              footer = NULL)
                            )

                            #If label is not present in training dataset
                            else if(!input$Label_name %in% unique(Result_list$Training_Dataset$Label)) shiny::showModal(modalDialog(
                              paste0(input$Label_name, " is not present in the training dataset"),
                              easyClose = TRUE,
                              footer = NULL)
                            )

                            else{
                              #Show an alert
                              shinyalert::shinyalert(title = "WARNING!",
                                                     text = paste0(input$Label_name, " cells will be removed from the training datasets"),
                                                     type = "warning",
                                                     closeOnEsc = TRUE,
                                                     closeOnClickOutside = TRUE,
                                                     showCancelButton = TRUE,
                                                     showConfirmButton = TRUE,
                                                     confirmButtonText = "Proceed",
                                                     cancelButtonText = "Cancel",
                                                     callbackR = function() Result_list$Training_Dataset <- Result_list$Training_Dataset %>% dplyr::filter(Label != input$Label_name)
                              )
                            }
                          },
                          ignoreInit = TRUE)

      #If the user hits fit model button
      shiny::observeEvent(input$Fit_model,
                          {
                            #IF Training dataset has not been created
                            if(is.null(Result_list$Training_Dataset)) shiny::showModal(modalDialog(
                              paste0("A model cannot be created without an active Training dataset"),
                              easyClose = TRUE,
                              footer = NULL)
                            )

                            else{
                              #IF NO SPATIAL INFORMATION IS REQUIRED
                              if(!as.logical(input$GO_spatial)){
                                #Show an alert
                                shinyalert::shinyalert(title = "WARNING!",
                                                       text = paste0("The model training process will be executed using ", nrow(Result_list$Training_Dataset), " cells"),
                                                       type = "warning",
                                                       closeOnEsc = TRUE,
                                                       closeOnClickOutside = TRUE,
                                                       showCancelButton = TRUE,
                                                       showConfirmButton = TRUE,
                                                       confirmButtonText = "Proceed",
                                                       cancelButtonText = "Cancel",
                                                       callbackR = function(){

                                                         #GENERATE THE FEATURE DATA
                                                         shiny::showModal(modalDialog("Obtaining the features from the training dataset", footer=NULL))
                                                         Feature_DATA <- DATA %>% dplyr::filter(Cell_no %in% Result_list$Training_Dataset$Cell_no) %>% dplyr::select(Cell_no, dplyr::all_of(as.character(input$Model_vars)))
                                                         Final_Training_Data <- dplyr::left_join(Result_list$Training_Dataset, Feature_DATA, by = "Cell_no")
                                                         shiny::removeModal()

                                                         #FIT THE MODEL
                                                         shiny::showModal(modalDialog("Fitting the model", footer=NULL))


                                                         #Random Forest pathway
                                                         if(input$Method == "Random forest"){
                                                           Model_recipe <-
                                                             Final_Training_Data %>% dplyr::select(-c(1:4)) %>%
                                                             recipes::recipe(Label ~ .) %>% recipes::step_normalize(recipes::all_predictors()) %>%
                                                             workflows::workflow(
                                                               parsnip::rand_forest(mode = "classification",
                                                                                    engine = "randomForest",
                                                                                    mtry = ceiling(as.numeric(input$RF_mtry)* ncol(Final_Training_Data)-5),
                                                                                    trees = input$RF_trees,
                                                                                    min_n = 1)

                                                             )
                                                           Model <- parsnip::fit(Model_recipe , Final_Training_Data)

                                                           #Define the model and the parameters in the result list
                                                           Result_list$Model <- Model
                                                           Result_list$Model_Param <- list(Model_type = "Random forest",
                                                                                           Model_features = as.character(input$Model_vars),
                                                                                           Model_threshold = as.numeric(input$Threshold),
                                                                                           Spatial_context = as.logical(input$GO_spatial),
                                                                                           Per_features = as.numeric(input$RF_mtry),
                                                                                           Trees = as.numeric(input$RF_trees)
                                                           )
                                                         }
                                                         #XG BOOST pathway
                                                         if(input$Method == "XG boost"){
                                                           Model_recipe <-
                                                             Final_Training_Data %>% dplyr::select(-c(1:4)) %>%
                                                             recipes::recipe(Label ~ .) %>% recipes::step_normalize(recipes::all_predictors()) %>%
                                                             workflows::workflow(
                                                               parsnip::boost_tree(mode = "classification",
                                                                                   engine = "xgboost",
                                                                                   mtry = ceiling(as.numeric(input$XG_mtry) * ncol(Final_Training_Data)-5),
                                                                                   trees = as.numeric(input$XG_trees),
                                                                                   min_n = 1,
                                                                                   tree_depth = as.numeric(input$XG_tree_depth),
                                                                                   learn_rate = 0.1,
                                                                                   loss_reduction = 0,
                                                                                   sample_size = as.numeric(input$XG_sample_size),
                                                                                   stop_iter = Inf)
                                                             )
                                                           Model <- parsnip::fit(Model_recipe , Final_Training_Data)

                                                           #Define the model and the parameters in the result list
                                                           Result_list$Model <- Model
                                                           Result_list$Model_Param <- list(Model_type = "XG boost",
                                                                                           Model_features = as.character(input$Model_vars),
                                                                                           Model_threshold = as.numeric(input$Threshold),
                                                                                           Spatial_context = as.logical(input$GO_spatial),
                                                                                           Per_features = as.numeric(input$XG_mtry),
                                                                                           Per_cells = as.numeric(input$XG_sample_size),
                                                                                           Trees = as.numeric(input$XG_trees),
                                                                                           Tree_depth = as.numeric(input$XG_tree_depth)
                                                           )
                                                         }
                                                         #NNET pathway
                                                         if(input$Method == "NNET"){
                                                           Model_recipe <-
                                                             Final_Training_Data %>% dplyr::select(-c(1:4)) %>%
                                                             recipes::recipe(Label ~ .) %>% recipes::step_normalize(recipes::all_predictors()) %>%
                                                             workflows::workflow(
                                                               parsnip::mlp(mode = "classification",
                                                                            hidden_units = rep(as.numeric(input$Hidden), as.numeric(input$Layers)),
                                                                            epochs = as.numeric(input$Epochs),
                                                                            engine = "brulee",
                                                                            penalty = as.numeric(input$Penalty))
                                                             )
                                                           Model <- parsnip::fit(Model_recipe , Final_Training_Data)

                                                           #Define the model and the parameters in the result list
                                                           Result_list$Model <- Model
                                                           Result_list$Model_Param <- list(Model_type = "NNET",
                                                                                           Model_features = as.character(input$Model_vars),
                                                                                           Model_threshold = as.numeric(input$Threshold),
                                                                                           Spatial_context = as.logical(input$GO_spatial),
                                                                                           Hidden_units = as.numeric(input$Hidden),
                                                                                           Hidden_layers = as.numeric(input$Layers),
                                                                                           Epochs = as.numeric(input$Epochs),
                                                                                           Penalty = as.numeric(input$Penalty)
                                                           )
                                                         }
                                                         shiny::removeModal()

                                                       }
                                )
                              }

                              #IF SPATIAL
                              if(as.logical(input$GO_spatial)){
                                #Show an alert that spatial data will be used
                                shinyalert::shinyalert(title = "WARNING!",
                                                       text = paste0("The model training process will be executed using ", nrow(Result_list$Training_Dataset), " cells. ",
                                                                     "Spatial information will be used. This can be computationally intensive"),
                                                       type = "warning",
                                                       closeOnEsc = TRUE,
                                                       closeOnClickOutside = TRUE,
                                                       showCancelButton = TRUE,
                                                       showConfirmButton = TRUE,
                                                       confirmButtonText = "Proceed",
                                                       cancelButtonText = "Cancel",
                                                       callbackR = function(){

                                                         #GENERATE THE FEATURE DATA
                                                         shiny::showModal(modalDialog("Obtaining the features from the training dataset and performing spatial information sharing...", footer=NULL))

                                                         #Obtain all the data and select only the features of interest
                                                         Neighbor_Feature_DATA <- DATA %>% dplyr::select(1:4, dplyr::all_of(as.character(input$Model_vars)))
                                                         #Perform a sort of message passing based on the customized version of UTAG message passing function
                                                         Neighbor_Feature_DATA <- UTAG_message_passing_Image_based_phenotyper(
                                                           DATA = Neighbor_Feature_DATA,
                                                           COO_to_visit = Neighbor_Feature_DATA$Cell_no %in% Result_list$Training_Dataset$Cell_no,
                                                           Neighbor_strategy = input$Neighbor_strategy,
                                                           Message_strategy = input$Message_strategy,
                                                           N_neighbors = input$N_neighbors,
                                                           Max_dist_allowed = input$Max_dist_allowed,
                                                           Weighting_Strategy = input$Weighting_Strategy,
                                                           N_cores = input$N_cores
                                                         )
                                                         #Retain only the cell_no and the features of the neighbors
                                                         Neighbor_Feature_DATA <- Neighbor_Feature_DATA %>% dplyr::select(-X, -Y, -Subject_Names, -mean_DIST, -max_DIST, -N_neighbors)
                                                         #Change the names to reflect that this info comes from the neighbors
                                                         names(Neighbor_Feature_DATA)[-1] <-stringr::str_c("Neighbor_", names(Neighbor_Feature_DATA)[-1], sep = "")
                                                         shiny::removeModal()

                                                         #Obtaining the traditional features
                                                         shiny::showModal(modalDialog("Generating the final training datasets", footer=NULL))
                                                         Feature_DATA <- DATA %>% dplyr::filter(Cell_no %in% Result_list$Training_Dataset$Cell_no) %>% dplyr::select(Cell_no, dplyr::all_of(as.character(input$Model_vars)))
                                                         Feature_DATA <-dplyr::left_join(Feature_DATA, Neighbor_Feature_DATA, by = "Cell_no")
                                                         Final_Training_Data <- dplyr::left_join(Result_list$Training_Dataset, Feature_DATA, by = "Cell_no")
                                                         shiny::removeModal()

                                                         #FIT THE MODEL
                                                         shiny::showModal(modalDialog("Fitting the model", footer=NULL))

                                                         #Random Forest pathway
                                                         if(input$Method == "Random forest"){
                                                           Model_recipe <-
                                                             Final_Training_Data %>% na.omit() %>% dplyr::select(-c(1:4)) %>%
                                                             recipes::recipe(Label ~ .) %>% recipes::step_normalize(recipes::all_predictors()) %>%
                                                             workflows::workflow(
                                                               parsnip::rand_forest(mode = "classification",
                                                                                    engine = "randomForest",
                                                                                    mtry = ceiling(as.numeric(input$RF_mtry)* ncol(Final_Training_Data)-5),
                                                                                    trees = input$RF_trees,
                                                                                    min_n = 1)
                                                             )
                                                           Model <- parsnip::fit(Model_recipe , Final_Training_Data)

                                                           #Define the model and the parameters in the result list
                                                           Result_list$Model <- Model
                                                           Result_list$Model_Param <- list(Model_type = "Random forest",
                                                                                           Model_features = as.character(input$Model_vars),
                                                                                           Model_threshold = as.numeric(input$Threshold),
                                                                                           Spatial_context = as.logical(input$GO_spatial),
                                                                                           Neighbor_strategy = input$Neighbor_strategy,
                                                                                           Message_strategy = input$Message_strategy,
                                                                                           N_neighbors = input$N_neighbors,
                                                                                           Max_dist_allowed = input$Max_dist_allowed,
                                                                                           Weighting_Strategy = input$Weighting_Strategy,
                                                                                           N_cores = input$N_cores,
                                                                                           Per_features = as.numeric(input$RF_mtry),
                                                                                           Trees = as.numeric(input$RF_trees)
                                                           )
                                                         }
                                                         #XG BOOST pathway
                                                         if(input$Method == "XG boost"){
                                                           Model_recipe <-
                                                             Final_Training_Data %>% na.omit() %>% dplyr::select(-c(1:4)) %>%
                                                             recipes::recipe(Label ~ .) %>% recipes::step_normalize(recipes::all_predictors()) %>%
                                                             workflows::workflow(
                                                               parsnip::boost_tree(mode = "classification",
                                                                                   engine = "xgboost",
                                                                                   mtry = ceiling(as.numeric(input$XG_mtry) * ncol(Final_Training_Data)-5),
                                                                                   trees = as.numeric(input$XG_trees),
                                                                                   min_n = 1,
                                                                                   tree_depth = as.numeric(input$XG_tree_depth),
                                                                                   learn_rate = 0.1,
                                                                                   loss_reduction = 0,
                                                                                   sample_size = as.numeric(input$XG_sample_size),
                                                                                   stop_iter = Inf)
                                                             )
                                                           Model <- parsnip::fit(Model_recipe , Final_Training_Data)

                                                           #Define the model and the parameters in the result list
                                                           Result_list$Model <- Model
                                                           Result_list$Model_Param <- list(Model_type = "XG boost",
                                                                                           Model_features = as.character(input$Model_vars),
                                                                                           Model_threshold = as.numeric(input$Threshold),
                                                                                           Spatial_context = as.logical(input$GO_spatial),
                                                                                           Neighbor_strategy = input$Neighbor_strategy,
                                                                                           Message_strategy = input$Message_strategy,
                                                                                           N_neighbors = input$N_neighbors,
                                                                                           Max_dist_allowed = input$Max_dist_allowed,
                                                                                           Weighting_Strategy = input$Weighting_Strategy,
                                                                                           N_cores = input$N_cores,
                                                                                           Per_features = as.numeric(input$XG_mtry),
                                                                                           Per_cells = as.numeric(input$XG_sample_size),
                                                                                           Trees = as.numeric(input$XG_trees),
                                                                                           Tree_depth = as.numeric(input$XG_tree_depth)
                                                           )
                                                         }
                                                         #NNET pathway
                                                         if(input$Method == "NNET"){
                                                           Model_recipe <-
                                                             Final_Training_Data %>% na.omit() %>% dplyr::select(-c(1:4)) %>%
                                                             recipes::recipe(Label ~ .) %>% recipes::step_normalize(recipes::all_predictors()) %>%
                                                             workflows::workflow(
                                                               parsnip::mlp(mode = "classification",
                                                                            hidden_units = rep(as.numeric(input$Hidden), as.numeric(input$Layers)),
                                                                            epochs = as.numeric(input$Epochs),
                                                                            engine = "brulee",
                                                                            penalty = as.numeric(input$Penalty))
                                                             )
                                                           Model <- parsnip::fit(Model_recipe , Final_Training_Data)

                                                           #Define the model and the parameters in the result list
                                                           Result_list$Model <- Model
                                                           Result_list$Model_Param <- list(Model_type = "NNET",
                                                                                           Model_features = as.character(input$Model_vars),
                                                                                           Model_threshold = as.numeric(input$Threshold),
                                                                                           Spatial_context = as.logical(input$GO_spatial),
                                                                                           Neighbor_strategy = input$Neighbor_strategy,
                                                                                           Message_strategy = input$Message_strategy,
                                                                                           N_neighbors = input$N_neighbors,
                                                                                           Max_dist_allowed = input$Max_dist_allowed,
                                                                                           Weighting_Strategy = input$Weighting_Strategy,
                                                                                           N_cores = input$N_cores,
                                                                                           Hidden_units = as.numeric(input$Hidden),
                                                                                           Hidden_layers = as.numeric(input$Layers),
                                                                                           Epochs = as.numeric(input$Epochs),
                                                                                           Penalty = as.numeric(input$Penalty)
                                                           )
                                                         }
                                                         shiny::removeModal()
                                                       }
                                )
                              }
                            }
                          },
                          ignoreInit = TRUE)

      #If the user hits the test button
      shiny::observeEvent(input$Test_model,
                          {
                            #IF no model STOP
                            if(is.null(Result_list$Model)) shiny::showModal(modalDialog(
                              paste0("There is no active model"),
                              easyClose = TRUE,
                              footer = NULL)
                            )

                            #If there is an active model proceed
                            else{
                              #If no spatial context is needed
                              if(!as.logical(Result_list$Model_Param$Spatial_context)){
                                Test_Data <- DATA %>% dplyr::filter(Subject_Names == input$Image_name)
                                Test_DATA_features <- Test_Data %>% dplyr::select(dplyr::all_of(Result_list$Model_Param$Model_features))

                                Predictions <- predict(Result_list$Model, new_data = Test_DATA_features, type = "prob")
                                Col_index <- max.col(Predictions, ties.method = "random")
                                Predictions <- tibble(Label = colnames(Predictions)[Col_index],
                                                      Probability =purrr::map2_dbl(.x = 1:nrow(Predictions), .y = Col_index, function(.x, .y) Predictions[[.x, .y]]))
                                Predictions$Label <- substr(Predictions$Label, start = 7, stop = nchar(Predictions$Label))

                                Result_list$Test_Dataset <-dplyr::bind_cols(Test_Data %>% dplyr::select(1:4),
                                                                            Predictions)
                              }

                              #If Spatial context is needed
                              if(as.logical(Result_list$Model_Param$Spatial_context)){

                                shinyalert::shinyalert(title = "WARNING!",
                                                       text = paste0("Spatial information will be used to test the model. This can be computationally intensive"),
                                                       type = "warning",
                                                       closeOnEsc = TRUE,
                                                       closeOnClickOutside = TRUE,
                                                       showCancelButton = TRUE,
                                                       showConfirmButton = TRUE,
                                                       confirmButtonText = "Proceed",
                                                       cancelButtonText = "Cancel",
                                                       callbackR = function(){
                                                         #Obtain the image data
                                                         Test_Data <- DATA %>% dplyr::filter(Subject_Names == input$Image_name)
                                                         #Generate the neighbor features
                                                         Neighbors_DATA_features <- Test_Data %>% dplyr::select(1:4, dplyr::all_of(Result_list$Model_Param$Model_features))

                                                         #Perform a sort of message passing based on the customized version of UTAG message passing function
                                                         Neighbors_DATA_features <- UTAG_message_passing_Image_based_phenotyper(
                                                           DATA = Neighbors_DATA_features,
                                                           COO_to_visit = NULL,
                                                           Neighbor_strategy = Result_list$Model_Param$Neighbor_strategy,
                                                           Message_strategy = Result_list$Model_Param$Message_strategy,
                                                           N_neighbors = Result_list$Model_Param$N_neighbors,
                                                           Max_dist_allowed = Result_list$Model_Param$Max_dist_allowed,
                                                           Weighting_Strategy = Result_list$Model_Param$Weighting_Strategy,
                                                           N_cores = Result_list$Model_Param$N_cores
                                                         )
                                                         #Retain only the cell_no and the features of the neighbors
                                                         Neighbors_DATA_features <- Neighbors_DATA_features %>% dplyr::select(-X, -Y, -Subject_Names, -mean_DIST, -max_DIST, -N_neighbors)
                                                         #Change the names to reflect that this info comes from the neighbors
                                                         names(Neighbors_DATA_features)[-1] <-stringr::str_c("Neighbor_", names(Neighbors_DATA_features)[-1], sep = "")

                                                         #Get the actual data
                                                         Test_DATA_features <- Test_Data %>% dplyr::select(Cell_no, dplyr::all_of(Result_list$Model_Param$Model_features))

                                                         #Bind both datasets
                                                         Test_DATA_features <-dplyr::left_join(Test_DATA_features, Neighbors_DATA_features, by = "Cell_no") %>% dplyr::select(-Cell_no)

                                                         #Proceed with the predictions
                                                         Predictions <- predict(Result_list$Model, new_data = Test_DATA_features, type = "prob")
                                                         Col_index <- max.col(Predictions, ties.method = "random")
                                                         Predictions <- tibble(Label = colnames(Predictions)[Col_index],
                                                                               Probability =purrr::map2_dbl(.x = 1:nrow(Predictions), .y = Col_index, function(.x, .y) Predictions[[.x, .y]]))
                                                         Predictions$Label <- substr(Predictions$Label, start = 7, stop = nchar(Predictions$Label))

                                                         Result_list$Test_Dataset <-dplyr::bind_cols(Test_Data %>% dplyr::select(1:4),
                                                                                                     Predictions)
                                                       }
                                )

                              }
                            }
                          },
                          ignoreInit = TRUE)

      #If the user hits the Download button
      shiny::observeEvent(input$Download_model,
                          {
                            if(is.null(Result_list$Model)){
                              shiny::showModal(modalDialog(
                                paste0("There is no active model"),
                                easyClose = TRUE,
                                footer = NULL)
                              )
                            }
                            else{
                              Cell_phenotyping_model <<- list(Training_dataset = Result_list$Training_Dataset,
                                                              Model_parameters = Result_list$Model_Param ,
                                                              Model = Result_list$Model)

                              shiny::showModal(modalDialog(
                                paste0("An object called 'Cell_phenotyping_model' has been created in the GlobalEnvironment. Use it in the Model_cell_phenotyper function to predict cell phenotypes"),
                                easyClose = TRUE,
                                footer = NULL)
                              )
                            }
                          },
                          ignoreInit = TRUE)

      #If browser is closed end the app
      session$onSessionEnded(function() { shiny::stopApp() })
    }

    #Run the server
    message("Always stop current R execution if you want to continue with your R session")
    shiny::shinyApp(user_interface, server)
  }
