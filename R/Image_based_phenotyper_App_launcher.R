#' Launches a shiny APP that can generate cell phenotyping models based on user expertise
#'
#' Launches an APP to interactively explore cell segmented images. The user can manually assign cell phenotype labels to cells and train a model.
#' The model can then be used in [Model_cell_phenotyper()] function to generalize it to all the cells in a dataset.
#'
#' @param DATA A dataframe or tibble containing cell feature data.
#' @param Directory Character specifying the path to the folder where images are present.
#' @param Ordered_Channels Character vector specifying image channels in their exact order.
#' @param Channels_to_keep Character vector indicating the channels to be kept in the analysis.
#' @param Max_Gb_cache A single numeric value indicating the memory size of the image cache (see details). 10 Gb by default.
#' 
#' @seealso [Model_cell_phenotyper()]
#'
#' @details
#' In order to deal with large images and speed up image toggling, the shiny APP works with a memoised version of the [Smart_image_importer()] function. Loaded images will
#' be stored in a temporary cache until the APP is closed or the cache reaches it's limit. 
#' 
#' Image setting allow the user to control the image channel display.
#' 
#' User can use the 'Color', 'Add channel', 'Remove channel' and 'Reset' buttons to merge different channels into a single image. The add channel button
#' will color the current channel according to the color provided and will be merged to the current image.
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
#' Model settings allows the user to control how the model will be calculated. If spatial interaction features are required, the user must
#' specify how neighbors are defined. The 'thresh' (threshold) parameter defines the minimum confidence of the prediction required to
#' assign a cell a label. If threshold is not reached, the cell will remain unnasigned.
#'
#'Relevan buttons in the lower area of the control panel
#'\itemize{
#'\item{Fit button computes the model with the current model settings and test dataset}
#'\item{Test button tests the model in the current displayed image}
#'\item{Download button saves the current active model in the Global environment}
#'}
#'
#' Image Panels:
#'\itemize{
#' \item{Upper left panel: Displays the image (use it to zoom in and out).}
#' \item{Upper right panel: Allows assigning cells for the training dataset. Select cells by clicking or use the area selection tools provided. Cell labels will be updated as they are assigned.}
#' \item{Lower left: Shows features of the training dataset and the results of the test dataset if tests have been performed.}
#' \item{Lower right: If a test has been performed, it shows the model results. Misclassified cells can be selected and labels can be re-assigned to refine the model.}
#' }
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

Image_based_phenotyper_App_launcher2 <-
  function(DATA,
           Directory,
           Ordered_Channels,
           Channels_to_keep,
           Max_Gb_cache = 10
  ){
    
    ####Check suggested packages####
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
    
    ####What to do on exit####
    on.exit({
      future::plan("future::sequential")
      gc()
    })
    
    ####Check arguments and generate memoised version####
    #Check that ordered channels and channels to keep are unique
    if(length(Ordered_Channels) != length(unique(Ordered_Channels))) stop("Ordered_Channels must contain non-duplicated character values")
    if(length(Channels_to_keep) != length(unique(Channels_to_keep))) stop("Channels_to_keep must contain non-duplicated character values")
    
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
    
    ####Look-up table####
    #Generate look-up table
    Image_names <- dir(Directory, full.names = FALSE)
    Directory_path <- dir(Directory, full.names = TRUE)
    DATA <- DATA
    
    #Remove any features in the data that are not numeric
    Numeric_features <- purrr::map_lgl(DATA[-c(1:4)], ~is.numeric(.))
    if(any(!Numeric_features)){
      print(paste0("The following non-numeric features will be removed from DATA: ", stringr::str_c(names(DATA)[-c(1:4)][!Numeric_features], collapse = ", ")))
      DATA <- DATA %>% dplyr::select(-all_of(stringr::str_c(names(DATA)[-c(1:4)][!Numeric_features])))
    }
    
    print(paste0(length(Image_names), " files found in Directory. Finding matches between file names and Subject_Names in DATA..."))
    
    Look_up_table <-
      purrr::map_dfr(Image_names, function(Name){
        Distance_vector <- as.double(adist(Name, unique(DATA$Subject_Names), fixed = FALSE, ignore.case = TRUE))
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
    
    ####APP UI####
    #BUILD THE USER INTERFACE
    {
      user_interface <- shiny::fluidPage(
        
        #To make the sidebar collapsable
        shinyjs::useShinyjs(),
        
        #Set the title and action button
        shiny::fluidRow(
          shiny::column(width = 3, shiny::h3("Supervised Phenotyper")),
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
              shiny::actionButton("toggle_cells", "Cell labels", icon = shiny::icon(name = "square-caret-up", lib = "font-awesome")),
              shiny::actionButton("toggle_params", "Model settings", icon = shiny::icon(name = "square-caret-up", lib = "font-awesome"))
            ),
            
            
            #IMAGE PARAMETERS
            shiny::h4("Image settings"),
            #Collapsable rows
            shiny::tags$div(
              class = "Image_settings",
              
              shiny::fluidRow(
                shiny::column(6, shiny::selectInput("Image_name", "Image", sort(Look_up_table$Subject_Names), multiple = FALSE)),
                shiny::column(2, shiny::selectInput("Res", "Res", c(Low = 5, Mid = 5.7, High = 6, Ultra = 6.3), selected = 5.7, multiple = FALSE)),
                shiny::column(2, shinyWidgets::materialSwitch("Change_coords", "Pixel/dist", value = FALSE)),
                shiny::column(2, shiny::conditionalPanel(condition = "input.Change_coords == '1'",
                                                         shiny::textInput("Ratio", "pixel size", value = "1")
                ))),
              shiny::fluidRow(
                shiny::column(3, shiny::sliderInput("Degrees", "Rotate", min = -90, max = 90, value = 0, step = 90, ticks = FALSE)),
                shiny::column(3, shiny::selectInput("X_flip", "Flip X image", c(YES = TRUE, NO = FALSE), selected = FALSE, multiple = FALSE)),
                shiny::column(3, shiny::selectInput("Y_flip", "Flip Y image", c(YES = TRUE, NO = FALSE), selected = FALSE, multiple = FALSE)),
                shiny::column(3, shiny::selectInput("Equalize", "Equalize", c(YES = TRUE, NO = FALSE), selected = FALSE, multiple = FALSE))
              ),
              
              shiny::fluidRow(
                shiny::column(3, shiny::selectInput("Channel", "Channel", Channels_to_keep, multiple = FALSE)),
                shiny::column(3, shiny::sliderInput("Min_Image", "Min", value = 0, min = 0, max = 100, step = 1, ticks = FALSE)),
                shiny::column(3, shiny::sliderInput("Max_Image", "Max", value = 100, min = 0, max = 100, step = 1, ticks = FALSE)),
                shiny::column(3, shiny::sliderInput("Gamma", "Gamma", value = 0, min = -3, max = +3, step = 0.01, ticks = FALSE))
              ),
              
              shiny::fluidRow(
                shiny::column(3, shiny::textInput("Color_channel", "Color", value = "blue")),
                shiny::column(3, shiny::actionButton("Add_color", "Add channel", shiny::icon("plus", library = "font-awesome"))),
                shiny::column(3, shiny::actionButton("Remove_color", "Remove channel", shiny::icon("minus", library = "font-awesome"))),
                shiny::column(3, shiny::actionButton("Reset_image", "Reset", icon = shiny::icon("redo")))
              )
              
            ),
            
            #Cell label parameters
            shiny::h4("Cell label assigner"),
            
            #Collapsable rows
            shiny::tags$div(
              class = "Cell_label",
              
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
              )
            ),
            
            #Model PARAMETERS
            shiny::h4("Model settings"),
            shiny::tags$div(
              class = "Model_settings",
              
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
              )
            ),
            
            #BUTTONS
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
          }'))),
        
        
        shiny::tags$script(htmltools::HTML("
    $(document).on('click', '#toggle_Image', function() {
      $('.Image_settings').toggle();
    });
    $(document).on('click', '#toggle_cells', function() {
      $('.Cell_label').toggle();
    });
    $(document).on('click', '#toggle_params', function() {
      $('.Model_settings').toggle();
    });
  "))
      )
    }
    
    
    
    ####SERVER####
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
      #Generates a reactive that stores the original image coordinates and cropping parameters
      Original_ranges <- 
        shiny::reactiveValues(
          #Image dimension Original
          Original_Dims_Width = NULL, 
          Original_Dims_Height = NULL
        )
      Pixel_dist_conversion <- shiny::reactive(input$Change_coords)
      Pixel_dist_ratio <- shiny::reactive(input$Ratio)
      Resolution <- shiny::reactive(input$Res)
      
      #The reactive that controls the final output
      Result_list <- shiny::reactiveValues(Training_Dataset = NULL,
                                           Model_Param = NULL,
                                           Model = NULL,
                                           Test_Dataset = NULL)
      
      #The reactive that controls the parameter list for image coloring (Null by default at start)
      Image_coloring_list <- shiny::reactiveValues(Parameter_list = NULL)
      
      #Control the buttons for the color list
      #The ADD CHANNEL BUTTON
      shiny::observeEvent(
        input$Add_color, {
          
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
        
        #Remove any cell not being plotted
        if(any(!is.null(ranges$x), !is.null(ranges$y))){
          Final_DATA <- Final_DATA %>% dplyr::filter(X >= min(ranges$x), X <= max(ranges$x),
                                                     Y >= min(ranges$y), Y <= max(ranges$y))
        }
        
        #Return the final data
        return(Final_DATA)
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
        if(is.null(Image_coloring_list$Parameter_list)){
          #Plot the result
          try(Photo <- 
                Photo_plot_reactive() +
                ggtitle("Original") +
                theme(plot.title = element_text(size = 12, hjust = 0.5)))
        }
        if(!is.null(Image_coloring_list$Parameter_list)){
          #Generate the color tibble to generate the legend
          color_tibble <- 
            tibble(Channel_index = 1:length(Image_coloring_list$Parameter_list),
                   Channel_name = names(Image_coloring_list$Parameter_list),
                   Color_names = purrr::map_chr(Image_coloring_list$Parameter_list, ~.[["color"]])
            ) %>% arrange(Channel_index)
          color_tibble$Channel_name = factor(color_tibble$Channel_name, levels = color_tibble$Channel_name)
          
          try(Photo <- 
                Photo_plot_reactive() +
                ggtitle("Original") +
                theme(plot.title = element_text(size = 12, hjust = 0.5)) +
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
                                                          legend.box.margin = margin(-10, -10, -10, -10)))))
          
        }
        
        #Return the final result
        return(Photo)
      })
      
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
          Final_DATA <- dplyr::left_join(Final_DATA, Color_lookuptable(), by = "Label")
          
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
              guides(color = "none") +
              ggtitle("Assign cells") +
              theme(plot.title = element_text(size = 12, hjust = 0.5))
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
          Final_DATA <- dplyr::left_join(Final_DATA, Color_lookuptable(), by = "Label")
          
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
              guides(color = "none") +
              ggtitle("Model performance") +
              theme(plot.title = element_text(size = 12, hjust = 0.5))
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






  