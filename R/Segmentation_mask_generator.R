#' Launches a shiny APP to generate or refine cell segmentation masks
#'
#' The function launches an app to generate from scratch cell segmentation masks. If masks have already generated, the app allows the user to refine them by
#' removing foulty masks and re-drawing them
#'
#' @param Directory Character specifying the path to the folder where images to be segmented are stored.
#' @param Output_Directory Character specifying the path to the folder where cell mask images need to be written or where the already existing masks are present.
#' @param Ordered_Channels Character vector specifying image channels in their exact order.
#' @param Max_Gb_cache A single numeric value indicating the memory size of the image cache (see details). 10 Gb by default.
#'
#' @seealso [Watershed_cell_mask_generator()]
#'
#' @details
#' 
#' The function will first scope the image names in the Directory and in the Output_Directory. Cell mask image names in the Output_Directory are expected to
#' end with the appendix "_CellMask.tif". After matching images in the Directory with their corresponding cell masks in the Output_Directory the function
#' will organize the cell segmentation work starting with the first image without a corresponding cell mask. 
#' 
#' In order to deal with large images and speed up image toggling, the shiny APP works with a memoised version of the `magick::image_read()` function. Loaded images will
#' be stored in a temporary cache until the APP is closed or the cache reaches it's limit.
#'
#' User can use the 'Color', 'Add channel', 'Remove channel' and 'Reset' buttons to merge different channels into a single image. The add channel button
#' will color the current channel according to the color provided and will be merged to the current image.
#' 
#' WARNING!!! The app has veen developed to work with small tiled images. It is not a matter of the total pixel size but the number of cell masks in the image can
#' easily saturate the app. Every cell mask polygon needs to be computed every time a single change in the mask is made. This is not troublesome for low small images with
#' few cells but in can slow things down for large images with many cells.
#' 
#' 2nd WARNING!!!! The app works by generating sf objects containing polygons that are then turned into an EBImage object containing the image mask. When a cell mask
#' is already available it will load it and turn it into an sf object containing the polygons. Multi-polygon cell masks are not compatible with the app (cell masks
#' that consist of non-spatially-conected tissue areas). In addition, polygon cell masks fully containing other polygon cell masks are not allowed. The user needs to be
#' aware that their original cell masks generated with algorithms that allow these specific cell masks may be inadvertently severely modified when using the app.
#' 
#'
#' Relevant buttons:
#' \itemize{
#' \item{Previous and next button can be used to navigate through the images.}
#' \item{The Done! button will save the current changes by writing a cell mask image in the Output_Directory.}
#' \item{The Draw and Remove buttons will activate the drawing and erasing modes, respectively. While in drawing mode individual cells can be
#'outlined with the mouse cursor. While in erasing mode the user can click on the center of the desired cell to remove the corresponding mask. The erase all button
#'can be used to remove all currently drawn cell outlines.}
#' \item{The navigation panel in the bottom provides the user with information regarding the number of images in the Directory that have a cell mask in the Output_Directoyr.
#' The user can click on any of the images to navigate to the selected position.}
#' 
#' }
#'
#'
#' @details
#' For more information please see the vignette 'Performing cell segmentation'.
#'
#'
#' @examples
#' \dontrun{
#'#Generate our images
#'Input_Dir <- tempfile(pattern = "tempdir1_Input")
#'dir.create(Input_Dir, recursive = TRUE)
#'
#'Output_Dir <- tempfile(pattern = "tempdir2_Output")
#'dir.create(Output_Dir, recursive = TRUE)
#'
#'#Save images in Input directory
#'purrr::map(1:2,
#'          function(Image){
#'             EBImage::writeImage(CSM_MiniMultiTiff_test[[Image]],
#'                                 file.path(Input_Dir, names(CSM_MiniMultiTiff_test)[Image]))
#'           })
#'
#'#Use the app to segment the cells
#'Segmentation_mask_generator(Directory = Input_Dir, 
#'                            Output_Directory = Output_Dir, 
#'                           Ordered_Channels = CSM_SegmentParams_test$Ordered_Channels
#')
#'
#'#Regenerate the output dir
#'Output_Dir <- tempfile(pattern = "tempdir2_Output")
#'dir.create(Output_Dir, recursive = TRUE)
#'
#'#Check a segmentation parameters list obtained using the dedicated function----------
#'print(CSM_SegmentParams_test)
#'
#'
#'#Generate the whole cell masks using watershed---------------------------------------
#'Watershed_cell_mask_generator(
#'  Directory = Input_Dir,
#'  Output_directory = Output_Dir,
#'  Parameter_list = CSM_SegmentParams_test
#')
#'
#'#Use the app to refine cell masks
#'Segmentation_mask_generator(Directory = Input_Dir, 
#'                            Output_Directory = Output_Dir, 
#'                            Ordered_Channels = CSM_SegmentParams_test$Ordered_Channels
#')
#'}
#'
#' @export

Segmentation_mask_generator <-
  function(Directory,
           Output_Directory,
           Ordered_Channels,
           Max_Gb_cache = 10){
    
    ####Check suggested packages####
    {
      if(!requireNamespace("EBImage", quietly = TRUE)) stop(
        paste0("EBImage Bioconductor package is required to execute the function. Please install using the following code: ",
               expression({
                 if (!require("BiocManager", quietly = TRUE))
                   install.packages("BiocManager")
                 
                 BiocManager::install("EBImage")
               })
        )
      )
      if(!requireNamespace("magick", quietly = FALSE)) stop(
        paste0("magick CRAN package is required to execute the function. Please install using the following code: ",
               expression(install.packages("magick")))
      )
      if(!requireNamespace("later", quietly = FALSE)) stop(
        paste0("later CRAN package is required to execute the function. Please install using the following code: ",
               expression(install.packages("later")))
      )
      if(!requireNamespace("plotly", quietly = FALSE)) stop(
        paste0("plotly CRAN package is required to execute the function. Please install using the following code: ",
               expression(install.packages("plotly")))
      )
      if(!requireNamespace("terra", quietly = FALSE)) stop(
        paste0("terra CRAN package is required to execute the function. Please install using the following code: ",
               expression(install.packages("terra")))
      )
    }
    
    #Check arguments provided
    if(length(dir(Directory)) <1 ) stop("No files found in the Directory provided")
    if(is.null(Ordered_Channels)) stop("Ordered_Channels must not be NULL")
    #Check that ordered channels and channels to keep are unique
    if(length(Ordered_Channels) != length(unique(Ordered_Channels))) stop("Ordered_Channels must contain non-duplicated character values")
    #check that Output_Directory exists
    if(!dir.exists(Output_Directory)) stop("Issue with Ouput_Directory")
    
    #Check that the Max_Gb_cache is OK
    if(!all(is.numeric(Max_Gb_cache), Max_Gb_cache > 0, length(Max_Gb_cache) == 1)) stop("Max_Gb_cache must be a numeric value > 0")
    
    #Generate a memoised version of the the magick image reader
    Memoised_magick_reader <- memoise::memoise(
      #Function to be memorized
      magick::image_read,
      
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
    
    #Build the look up table to match input and output images
    Image_names <- dir(Directory, full.names = FALSE)
    Image_paths <- dir(Directory, full.names = TRUE)
    Look_up_tibble <- tibble::tibble(Image_name = Image_names,
                                     Image_path = Image_paths,
                                     Mask_name = stringr::str_c(Image_names, "_CellMask.tif") #Predict the cell mask name
    )
    Look_up_tibble$Mask_paths <- paste0(Output_Directory, "/", Look_up_tibble$Mask_name)
    Look_up_tibble$Mask_DONE <- Look_up_tibble$Mask_name %in% list.files(Output_Directory)
    
    #Arrange the Look_up_tibble to contain first the ones completed and then the other ones
    Look_up_tibble <- Look_up_tibble %>% dplyr::arrange(dplyr::desc(Mask_DONE), Image_name)
    #Set the first image to be observed
    Primitive_image_index <- min(which(Look_up_tibble$Mask_DONE == FALSE))
    #If Primitive index is Inf (no min found) then set to 1
    if(is.infinite(Primitive_image_index)) Primitive_image_index <- 1
    

    #Check if primitive index has an associated image to start the app
    if(file.exists(Look_up_tibble$Mask_paths[Primitive_image_index])){
      Existing_mask <- EBImage::readImage(Look_up_tibble$Mask_paths[Primitive_image_index])
      
      #If the mask has no cells turn to NULL
      if(max(Existing_mask) == 0){
        Primitive_polygons <- NULL
        Primitive_centroids <- NULL
      }
      
      #If cells present compute the cell polygons
      else{
        #If the dims are not equal then pad
        if(dim(Existing_mask)[1] != dim(Existing_mask)[2]){
          #If width is smaller than height
          if(dim(Existing_mask)[1] < dim(Existing_mask)[2]){
            #Compute the difference
            Difference_width <- dim(Existing_mask)[2] - dim(Existing_mask)[1]
            #Pad
            Existing_mask <- Pad_image(img = Existing_mask, Right_padding = Difference_width, Bottom_padding = 0)
          }
          #If width is larger than height
          else{
            #Compute difference
            Difference_height <- dim(Existing_mask)[1] - dim(Existing_mask)[2]
            #Pad
            Existing_mask <- Pad_image(img = Existing_mask, Right_padding = 0, Bottom_padding = Difference_height)
          }
        }
        
        #Now extract the polygons from the image with the dedicated function
        Polygons_list <- Mask_to_sf(Existing_mask)
        
        #Add the polygons to the reactive
        Primitive_polygons <- Polygons_list$sf_object
        Primitive_centroids <- Polygons_list$sf_centroids
      }
      
      
    } else{
      Primitive_polygons <- NULL
      Primitive_centroids <- NULL
    }
    
    #Add batch index and rows and columns for the navigation plot
    #We will have batches of 100 images
    N_Batches <- ceiling(nrow(Look_up_tibble)/100)
    Look_up_tibble$Batch_index <- rep(1:N_Batches, each = 100)[1:nrow(Look_up_tibble)]
    
    #25 columns and 5 rows per batch
    Look_up_tibble$Col_ID <- rep(1:25, times = N_Batches*4)[1:nrow(Look_up_tibble)]
    Look_up_tibble$Row_ID <- rep(rep(4:1, each = 25), times = N_Batches)[1:nrow(Look_up_tibble)]
    Look_up_tibble$Border_color <- "#c2c0c0" #Start without any selected
    Look_up_tibble$Fill_color <- "#f07587" #Start with background according to status
    Look_up_tibble$Fill_color[Look_up_tibble$Mask_DONE] <- "#9df5ad" #Continue with background according to status
    
    #Get the primitive batch
    Primitive_batch <- Look_up_tibble$Batch_index[Primitive_image_index]
    
    
    ####BUILD THE USER INTERFACE####
    {
      user_interface <- shiny::fluidPage(
        
        shinyjs::useShinyjs(),
        
        shiny::tags$head(
          #Style of control panel
          shiny::tags$style(
            htmltools::HTML("
        
        html, body {
          height: 100vh;
          margin: 0;
          overflow: hidden;
          font-family: Arial;
        }
        
        .container-fluid {
          height: 100vh;
        }
        
        #titlePanel {
          height: 70px;
          padding-top: 10px;
          padding-left: 15px;
        }
        
        #mainLayout {
          height: calc(100vh - 80px);
        }
        
        #sidebar {
          background-color: #c2c0c0;
          height: calc(100vh - 100px);
          overflow-y: auto;
        }
        
        #mainPanel {
          height: calc(100vh - 100px);
          padding: 0px;
          overflow: hidden;
        }
        
        #Img_plot {
          width: 100% !important;
          height: 100% !important;
        }
        
        .selected-mode {
        background-color: #00cc44 !important;
        color: white !important;
        font-weight: bold !important;
      }

      .unselected-mode {
        background-color: white !important;
        color: black !important;
        font-weight: normal !important;
      }
        
      ")
          )
        ),
        
        #--------------------------------------------------
        # HEADER
        #--------------------------------------------------
        
        shiny::tags$div(
          id = "titlePanel",
          
          shiny::h3("Cell mask generation APP"),
          
          shiny::tags$hr(style = "
        margin-top:5px;
        margin-bottom:5px;
      ")
        ), 
        
        #--------------------------------------------------
        # MAIN LAYOUT
        #--------------------------------------------------
        
        shiny::tags$div(
          id = "mainLayout",
          
          shiny::sidebarLayout(
            
            #------------------------------------------------
            # SIDEBAR
            #------------------------------------------------
            
            sidebarPanel = shiny::sidebarPanel(
              
              id = "sidebar",
              
              # IMAGE PARAMETERS
              shiny::h4("Image settings"),
              
              shiny::fluidRow(
                shiny::column(5,  shiny::selectInput("Channel", "Channel to display", Ordered_Channels, multiple = FALSE)),
                shiny::column(3, shiny::selectInput("Equalize", "Equalize", c(YES = TRUE, NO = FALSE), selected = FALSE, multiple = FALSE)),
                shiny::column(3, shiny::textInput("Color_channel", "Color", value = "blue")),
              ),
              
              shiny::fluidRow(
                shiny::column(4, shiny::sliderInput("Min_Image", "Black", value = 0, min = 0, max = 100, step = 1, ticks = FALSE)),
                
                shiny::column(4,shiny::sliderInput("Max_Image", "White", value = 100, min = 0, max = 100, step = 1, ticks = FALSE)),
                
                shiny::column(4, shiny::sliderInput("Gamma", "Gamma", value = 0, min = -3, max = 3, step = 0.01, ticks = FALSE))
              ),
              
              shiny::fluidRow(
                
                shiny::column(4,shiny::actionButton("Add_color", "Add channel", icon = shiny::icon("plus"))),
                
                shiny::column(4, shiny::actionButton("Remove_color", "Remove channel", icon = shiny::icon("minus"))),
                
                shiny::column(3, shiny::actionButton("Reset_image", "Reset", icon = shiny::icon("redo")))
              ),
              
              #Add line
              shiny::hr(),
              
              #Image Status
              shiny::fluidRow(
                shiny::column(9, shiny::uiOutput("Current_image_name_text")),
                shiny::column(3, shiny::uiOutput("Current_image_status_text"))
              ),
              
              
              #Add line
              shiny::hr(),
              
              
              # NAVIGATION
              shiny::fluidRow(
                shiny::column(4,shiny::actionButton("Previous_button", label = "Previous",icon = shiny::icon("left-long"))),
                
                shiny::column(4, shiny::actionButton("Next_button", label = "Next", icon = shiny::icon("right-long"))),
                
                shiny::column(4, shiny::actionButton("Done_button", label = "Done!", icon = shiny::icon("square-check")))
              ),
              
              shiny::hr(),
              
              shiny::fluidRow(
                shiny::column(4, shiny::actionButton("Draw_button", label = "Draw", icon = shiny::icon("pencil"), class = "selected-mode")),
                
                shiny::column(4, shiny::actionButton("Erase_button", label = "Remove",icon = shiny::icon("eraser"), class = "unselected-mode")),
                
                shiny::column(4, shiny::actionButton("Erase_all_button",label = "Erase all", icon = shiny::icon("trash")))
              ),
              
              shiny::hr(),
              
              #The polygon status
              shiny::fluidRow(
                shiny::column(12, shiny::uiOutput("status"))
              ),
              
              shiny::hr(),
              
              #The navigation plot
              shiny::fluidRow(
                style = "margin-left:0px;margin-right:0px;",
                shiny::column(12, style = "padding-left:0px;padding-right:0px;",
                              shiny::plotOutput("Navigation_plot", 
                                                width = "100%",
                                                height = "100px",
                                                hover = shiny::hoverOpts("nav_hover",
                                                                         delay = 100,
                                                                         delayType = "debounce",
                                                                         clip = TRUE,
                                                                         nullOutside = TRUE),
                                                click = shiny::clickOpts("nav_click", clip = TRUE)))
              ),
              
              #Buttons to up and down page and global progress
              shiny::fluidRow(
                style = "margin-left:0px;margin-right:0px;",
                shiny::column(1, style = "padding-left:0px;padding-right:0px;",
                              shiny::actionButton("Previous_Batch", label = "", icon = shiny::icon("chevron-left"))
                              ),
                shiny::column(2, style = "padding-left:0px;padding-right:0px;",
                              shiny::uiOutput("Batch_status")
                              ),
                shiny::column(1, style = "padding-left:0px;padding-right:0px;",
                              shiny::actionButton("Next_Batch", label = "", icon = shiny::icon("chevron-right"))
                              ),
                shiny::column(7, shiny::uiOutput("Progress_status")
                ),
              ),
              #The hovered tile and images display
              shiny::fluidRow(
                shiny::column(12, shiny::uiOutput("hovered_tile"))
                )
             
            ),
            
            
            
            #------------------------------------------------
            # MAIN PANEL
            #------------------------------------------------
            
            mainPanel = shiny::mainPanel(
              
              id = "mainPanel",
              
              plotly::plotlyOutput(
                "img",
                width = "100%",
                height = "100%"
              )
            )
          )
        )
      )
      
    }
    
    
    ####BUILD THE SERVER#####
    server <- function(input, output, session){
      
      # -------------------------------------------------------
      # BASIC REACTIVES
      # -------------------------------------------------------
      
      #All the basic reactives to be used
      Photo_name <- shiny::reactive(Look_up_tibble$Image_path[Reactive_value$Current_Image_index])
      Channel_index <- shiny::reactive(which(input$Channel == Ordered_Channels))
      Photo_min <- shiny::reactive(input$Min_Image)
      Photo_max <- shiny::reactive(input$Max_Image)
      Photo_gamma <- shiny::reactive(10^input$Gamma)
      Equalize <- shiny::reactive(input$Equalize)
      
      #Reactive value list containing polygons and drawmode
      Reactive_value <- shiny::reactiveValues(
        
        sf_object = Primitive_polygons, #No sf objects initially
        sf_centroids = Primitive_centroids, #No sf centroids initially
        
        mode = "draw", #Starting mode is always "draw"
        
        Current_Image_index = Primitive_image_index, #The initializing value is the Primitive_image_index
        Look_up_tibble = Look_up_tibble, #The look up tibble will also be used
        Batch_index = Primitive_batch, #The primitive batch will be the initiator
        
        Current_image_width = NULL, #The width and the height
        Current_image_height = NULL, #The width and the height
        
        Padded_image = FALSE, #Add padding if required
        Original_image_width = NULL,
        Original_image_height = NULL
      )
      
      
      # -------------------------------------------------------
      # PHOTO REACTIVE
      # -------------------------------------------------------
      
      #Photo displayed
      Photo_reactive <- shiny::reactive({
        
        #Import the Photo
        Photo <- Memoised_magick_reader(Photo_name())
        
        #If the Image_coloring_list$Parameter_list is null then proceed as usual
        if(is.null(Image_coloring_list$Parameter_list)){
          Photo <- Photo[Channel_index()]
          #Perform image equalization as requested by user
          if(as.logical(Equalize())) Photo <- Photo %>% magick::image_equalize()
          #Perform image white adjustment
          Photo <- Photo %>%
            magick::image_level(black_point = Photo_min(),
                                white_point = Photo_max(),
                                mid_point = Photo_gamma())
        }
        
        #If the Image_coloring_list$Parameter_list is present then use the Channel merger
        if(!is.null(Image_coloring_list$Parameter_list)){
          Photo <- Image_channel_color_merger(Image = Photo,
                                              Parameter_list = Image_coloring_list$Parameter_list)
        }
        
        #Change to an EBImage object and return
        Photo <- Photo %>% magick::as_EBImage()
        
        #Set the Padded_image to FALSE
        Reactive_value$Padded_image <- FALSE
        
        #If the image is not square pad the image on the right or on the bottom
        if(dim(Photo)[1] != dim(Photo)[2]){
          
          #Keep the original dims and set padding to true
          Reactive_value$Padded_image <- TRUE
          Reactive_value$Original_image_width <- dim(Photo)[1]
          Reactive_value$Original_image_height <- dim(Photo)[2]
          
          #If width is smaller than height
          if(dim(Photo)[1] < dim(Photo)[2]){
            #Compute the difference
            Difference_width <- dim(Photo)[2] - dim(Photo)[1]
            #Pad
            Photo <- Pad_image(img = Photo, Right_padding = Difference_width, Bottom_padding = 0)
          }
          #If width is larger than height
          else{
            #Compute difference
            Difference_height <- dim(Photo)[1] - dim(Photo)[2]
            #Pad
            Photo <- Pad_image(img = Photo, Right_padding = 0, Bottom_padding = Difference_height)
          }
        }
        
        #Add the width and the height to the reactive values
        Reactive_value$Current_image_width <- dim(Photo)[1]
        Reactive_value$Current_image_height <- dim(Photo)[2]
        
        #Return the photo
        return(Photo)
      })
      
      
      # -------------------------------------------------------
      # NEXT, PREVIOUS AND CLICK on navigation
      # -------------------------------------------------------
      
      #NEXT BUTTON BEHAVIOR
      shiny::observeEvent(
        input$Next_button, {
          
          #Place an alert
          shinyalert::shinyalert(title = "Proceed to next image",
                                 text = paste0("Any unsaved cell mask data will be lost"),
                                 type = "warning",
                                 closeOnEsc = TRUE,
                                 closeOnClickOutside = TRUE,
                                 showCancelButton = TRUE,
                                 showConfirmButton = TRUE,
                                 confirmButtonText = "Proceed",
                                 cancelButtonText = "Cancel",
                                 callbackR = function(){
                                   
                                   #Set the sf objects and centroids back to NULL
                                   Reactive_value$sf_object <- NULL
                                   Reactive_value$sf_centroids <- NULL
                                   
                                   #Get the current Image_index
                                   Current_index <- Reactive_value$Current_Image_index
                                   #If no more image available print a message and return to Primitive_image_index
                                   if(Current_index >= nrow(Reactive_value$Look_up_tibble)){
                                     shiny::showModal(modalDialog(
                                       paste0("Returning to initial image in the session"),
                                       easyClose = TRUE,
                                       footer = NULL)
                                     )
                                     Reactive_value$Current_Image_index <- Primitive_image_index
                                   }
                                   #Else add +1 to the current index
                                   else{
                                     Current_index <- Current_index + 1
                                     Reactive_value$Current_Image_index <- Current_index
                                   }
                                   
                                   #Update the current batch index to show in the navigation plot
                                   Reactive_value$Batch_index <- Reactive_value$Look_up_tibble$Batch_index[Reactive_value$Current_Image_index]
                                   
                                   #Check if the next current image has available cell mask
                                   New_index <- Reactive_value$Current_Image_index
                                   New_mask_path <- Reactive_value$Look_up_tibble$Mask_paths[New_index]
                                   
                                   #If the mask exists then proceed to import it
                                   if(file.exists(New_mask_path)){
                                     Existing_mask <- EBImage::readImage(New_mask_path)
                                     
                                     #If the mask has no cells turn to NULL
                                     if(max(Existing_mask) == 0){
                                       Reactive_value$sf_object <- NULL
                                       Reactive_value$sf_centroids <- NULL
                                     }
                                     
                                     else{
                                       #If the dims are not equal then pad
                                       if(dim(Existing_mask)[1] != dim(Existing_mask)[2]){
                                         #If width is smaller than height
                                         if(dim(Existing_mask)[1] < dim(Existing_mask)[2]){
                                           #Compute the difference
                                           Difference_width <- dim(Existing_mask)[2] - dim(Existing_mask)[1]
                                           #Pad
                                           Existing_mask <- Pad_image(img = Existing_mask, Right_padding = Difference_width, Bottom_padding = 0)
                                         }
                                         #If width is larger than height
                                         else{
                                           #Compute difference
                                           Difference_height <- dim(Existing_mask)[1] - dim(Existing_mask)[2]
                                           #Pad
                                           Existing_mask <- Pad_image(img = Existing_mask, Right_padding = 0, Bottom_padding = Difference_height)
                                         }
                                       }
                                       
                                       #Now extract the polygons from the image with the dedicated function
                                       Polygons_list <- Mask_to_sf(Existing_mask)
                                       
                                       #Add the polygons to the reactive
                                       Reactive_value$sf_object <- Polygons_list$sf_object
                                       Reactive_value$sf_centroids <- Polygons_list$sf_centroids
                                     }
                                     
                                     
                                   }
                                   
                                 }
          )
        }
      )
      
      #PREVIOUS BUTTON BEHAVIOR
      shiny::observeEvent(
        input$Previous_button, {
          
          #Place an alert
          shinyalert::shinyalert(title = "Proceed to previous image",
                                 text = paste0("Any unsaved cell mask data will be lost"),
                                 type = "warning",
                                 closeOnEsc = TRUE,
                                 closeOnClickOutside = TRUE,
                                 showCancelButton = TRUE,
                                 showConfirmButton = TRUE,
                                 confirmButtonText = "Proceed",
                                 cancelButtonText = "Cancel",
                                 callbackR = function(){
                                   
                                   #Set the sf objects and centroids back to NULL
                                   Reactive_value$sf_object <- NULL
                                   Reactive_value$sf_centroids <- NULL
                                   
                                   #Get the current Image_index
                                   Current_index <- Reactive_value$Current_Image_index
                                   #If no more image available print a message and return to Primitive_image_index
                                   if(Current_index == 1){
                                     shiny::showModal(modalDialog(
                                       paste0("Returning to initial image in the session"),
                                       easyClose = TRUE,
                                       footer = NULL)
                                     )
                                     Reactive_value$Current_Image_index <- Primitive_image_index
                                   }
                                   #Else substract -1 to the current index
                                   else{
                                     Current_index <- Current_index - 1
                                     Reactive_value$Current_Image_index <- Current_index
                                   }
                                   
                                   #Update the current batch index to show in the navigation plot
                                   Reactive_value$Batch_index <- Reactive_value$Look_up_tibble$Batch_index[Reactive_value$Current_Image_index]
                                   
                                   #Check if the previous current image has available cell mask
                                   New_index <- Reactive_value$Current_Image_index
                                   New_mask_path <- Reactive_value$Look_up_tibble$Mask_paths[New_index]
                                   
                                   #If the mask exists then proceed to import it
                                   if(file.exists(New_mask_path)){
                                     Existing_mask <- EBImage::readImage(New_mask_path)
                                     
                                     #If the mask has no cells turn to NULL
                                     if(max(Existing_mask) == 0){
                                       Reactive_value$sf_object <- NULL
                                       Reactive_value$sf_centroids <- NULL
                                     }
                                     else{
                                       #If the dims are not equal then pad
                                       if(dim(Existing_mask)[1] != dim(Existing_mask)[2]){
                                         #If width is smaller than height
                                         if(dim(Existing_mask)[1] < dim(Existing_mask)[2]){
                                           #Compute the difference
                                           Difference_width <- dim(Existing_mask)[2] - dim(Existing_mask)[1]
                                           #Pad
                                           Existing_mask <- Pad_image(img = Existing_mask, Right_padding = Difference_width, Bottom_padding = 0)
                                         }
                                         #If width is larger than height
                                         else{
                                           #Compute difference
                                           Difference_height <- dim(Existing_mask)[1] - dim(Existing_mask)[2]
                                           #Pad
                                           Existing_mask <- Pad_image(img = Existing_mask, Right_padding = 0, Bottom_padding = Difference_height)
                                         }
                                       }
                                       
                                       #Now extract the polygons from the image with the dedicated function
                                       Polygons_list <- Mask_to_sf(Existing_mask)
                                       
                                       #Add the polygons to the reactive
                                       Reactive_value$sf_object <- Polygons_list$sf_object
                                       Reactive_value$sf_centroids <- Polygons_list$sf_centroids
                                     }
                                     }
                                 }
          )
        }
      )
      
      #Observe event if click
      shiny::observeEvent(input$nav_click,
                          {
                            #Get the current look up tibble for navigation and get the current batch
                            Current_look_up_tibble <- Reactive_value$Look_up_tibble %>% dplyr::filter(Batch_index == Reactive_value$Batch_index)
                            
                            #Get the navigation data
                            Navigation_point <- 
                              shiny::nearPoints(
                                Current_look_up_tibble,
                                coordinfo = input$nav_click,
                                xvar = "Col_ID",
                                yvar = "Row_ID",
                                maxpoints = 1
                              )
                            #If no navigation data available then return empty
                            if(any(is.null(Navigation_point), nrow(Navigation_point) == 0)) return()
                            
                            #Else get the name of the image
                            New_image_name <- Navigation_point$Image_name[1]
                            New_index <- which(Reactive_value$Look_up_tibble$Image_name == New_image_name)
                            
                            #If the new index is identical to the current index return empty
                            if(New_index == Reactive_value$Current_Image_index) return()
                            
                            #Place an alert
                            shinyalert::shinyalert(title = paste0("Proceed to ", New_image_name),
                                                   text = paste0("Any unsaved cell mask data will be lost"),
                                                   type = "warning",
                                                   closeOnEsc = TRUE,
                                                   closeOnClickOutside = TRUE,
                                                   showCancelButton = TRUE,
                                                   showConfirmButton = TRUE,
                                                   confirmButtonText = "Proceed",
                                                   cancelButtonText = "Cancel",
                                                   callbackR = function(){
                                                     
                                                     #Set the sf objects and centroids back to NULL
                                                     Reactive_value$sf_object <- NULL
                                                     Reactive_value$sf_centroids <- NULL
                                                     
                                                     #Update the reactive value of the current index
                                                     Reactive_value$Current_Image_index <- New_index
                                                     
                                                     #Update the current batch index to show in the navigation plot
                                                     Reactive_value$Batch_index <- Reactive_value$Look_up_tibble$Batch_index[Reactive_value$Current_Image_index]
                                                     
                                                     #Check if the next current image has available cell mask
                                                     New_index <- Reactive_value$Current_Image_index
                                                     New_mask_path <- Reactive_value$Look_up_tibble$Mask_paths[New_index]
                                                     
                                                     #If the mask exists then proceed to import it
                                                     if(file.exists(New_mask_path)){
                                                       Existing_mask <- EBImage::readImage(New_mask_path)
                                                       
                                                       #If the mask has no cells turn to NULL
                                                       if(max(Existing_mask) == 0){
                                                         Reactive_value$sf_object <- NULL
                                                         Reactive_value$sf_centroids <- NULL
                                                       }
                                                       
                                                       else{
                                                         #If the dims are not equal then pad
                                                         if(dim(Existing_mask)[1] != dim(Existing_mask)[2]){
                                                           #If width is smaller than height
                                                           if(dim(Existing_mask)[1] < dim(Existing_mask)[2]){
                                                             #Compute the difference
                                                             Difference_width <- dim(Existing_mask)[2] - dim(Existing_mask)[1]
                                                             #Pad
                                                             Existing_mask <- Pad_image(img = Existing_mask, Right_padding = Difference_width, Bottom_padding = 0)
                                                           }
                                                           #If width is larger than height
                                                           else{
                                                             #Compute difference
                                                             Difference_height <- dim(Existing_mask)[1] - dim(Existing_mask)[2]
                                                             #Pad
                                                             Existing_mask <- Pad_image(img = Existing_mask, Right_padding = 0, Bottom_padding = Difference_height)
                                                           }
                                                         }
                                                         
                                                         #Now extract the polygons from the image with the dedicated function
                                                         Polygons_list <- Mask_to_sf(Existing_mask)
                                                         
                                                         #Add the polygons to the reactive
                                                         Reactive_value$sf_object <- Polygons_list$sf_object
                                                         Reactive_value$sf_centroids <- Polygons_list$sf_centroids
                                                       }
                                                     }
                                                     
                                                   }
                            )
                          })
      
      # -------------------------------------------------------
      # DONE button 
      # -------------------------------------------------------
      
      shiny::observeEvent(
        input$Done_button, {
          
          #Get the current look_up_tibble
          Current_look_up_tibble <- Reactive_value$Look_up_tibble
          
          #If no polygons then return a warning
          if(any(is.null(Reactive_value$sf_object), nrow(Reactive_value$sf_object) == 0)){
            
            #Write a blank cell mask image
            Cell_mask <- EBImage::Image(array(0, dim = c(Reactive_value$Current_image_width, Reactive_value$Current_image_height)))
            
            #If the image has padding, remove it before saving it
            if(Reactive_value$Padded_image){
              Cell_mask <- Cell_mask[1:Reactive_value$Original_image_width, 1:Reactive_value$Original_image_height]
            }
            
            #Get the mask name
            Mask_path <- Current_look_up_tibble$Mask_paths[Reactive_value$Current_Image_index]
            
            #Write the mask
            EBImage::writeImage(Cell_mask, Mask_path, bits.per.sample = 16, compression = "LZW")
            
            shiny::showModal(modalDialog(
              paste0("No cell masks have been identified. Writing empty mask"),
              easyClose = TRUE,
              footer = NULL)
            )
          }
          
          else{
            
            #Turn the SF object to an image
            Cell_mask <- Polygon_to_EBImage_mask(sf_polygons = Reactive_value$sf_object, 
                                                 Image_width = Reactive_value$Current_image_width,
                                                 Image_height = Reactive_value$Current_image_height)
            
            #If the image has padding, remove it before saving it
            if(Reactive_value$Padded_image){
              Cell_mask <- Cell_mask[1:Reactive_value$Original_image_width, 1:Reactive_value$Original_image_height]
            }
            
            #Get the mask name
            Mask_path <- Current_look_up_tibble$Mask_paths[Reactive_value$Current_Image_index]
            
            #Write the mask
            EBImage::writeImage(Cell_mask, Mask_path, bits.per.sample = 16, compression = "LZW")
            
            #Show a message
            Mask_image_name <- Current_look_up_tibble$Mask_name[Reactive_value$Current_Image_index]
            
            shiny::showModal(modalDialog(
              paste0(Mask_image_name, " saved in output directory"),
              easyClose = TRUE,
              footer = NULL)
            )
          }
          
          #Turn mask status to TRUE
          #Turn the current Image status to TRUE
          Current_look_up_tibble$Mask_DONE[Reactive_value$Current_Image_index] <- TRUE
          
          #Replace the old Look_up_tibble with the new one
          Reactive_value$Look_up_tibble <- Current_look_up_tibble
          
          
        }
      )
      
      # -------------------------------------------------------
      # CHANNEL COLOR BEHAVIOR
      # -------------------------------------------------------
      
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
      
      #REMOVE CHANNEL BUTTON
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
      
      #RESET CHANNELS BUTTON
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
      
      
      # -------------------------------------------------------
      # DRAW, ERASE and ERASE ALL BEHAVIOR
      # -------------------------------------------------------
      
      #MODIFY BUTTONS FOR DRAWING AND ERASING IF SELECTED
      observeEvent(input$Draw_button, {
        
        shinyjs::removeClass("Draw_button", "unselected-mode")
        shinyjs::addClass("Draw_button", "selected-mode")
        
        shinyjs::removeClass("Erase_button", "selected-mode")
        shinyjs::addClass("Erase_button", "unselected-mode")
      })
      
      observeEvent(input$Erase_button, {
        
        shinyjs::removeClass("Erase_button", "unselected-mode")
        shinyjs::addClass("Erase_button", "selected-mode")
        
        shinyjs::removeClass("Draw_button", "selected-mode")
        shinyjs::addClass("Draw_button", "unselected-mode")
      })
      
      #Change drawing behavior if draw button is selected
      observeEvent(
        input$Draw_button,
        {
          
          Reactive_value$mode <- "draw"
          
          session$sendCustomMessage(
            "setDrawMode",
            TRUE
          )
          
        }
      )
      
      #Change drawing behavior if erase button is selected button is selected
      observeEvent(
        input$Erase_button,
        {
          
          Reactive_value$mode <- "delete"
          
          session$sendCustomMessage(
            "setDrawMode",
            FALSE
          )
          
        }
      )
      
      #If the Erase all button is selected
      observeEvent(
        input$Erase_all_button,
        {
          
          n_pol <- if (
            is.null(Reactive_value$sf_object)
          ) {
            0
          } else {
            nrow(Reactive_value$sf_object)
          }
          
          shiny::showModal(
            
            shiny::modalDialog(
              
              title = "Delete all cells",
              
              paste(
                n_pol,
                " cells will be deleted. Proceed?"
              ),
              
              easyClose = TRUE,
              
              footer = shiny::tagList(
                
                shiny::modalButton(
                  "Cancel"
                ),
                
                shiny::actionButton(
                  "confirm_clear",
                  "Delete all",
                  class = "btn-danger"
                )
                
              )
              
            )
            
          )
          
        }
      )
      
      #If erase all button is confirmed
      observeEvent(
        input$confirm_clear,
        {
          
          Reactive_value$sf_object <- NULL
          
          shiny::removeModal()
          
        }
      )
      
      
      # -------------------------------------------------------
      # IMAGE VISUALIZATION
      # -------------------------------------------------------
      
      output$img <- plotly::renderPlotly({
        
        #Get the image and dimensions
        img <- Photo_reactive()
        h <- dim(img)[1]
        w <- dim(img)[2]
        
        #Write temporary file
        tmp <- tempfile(fileext = ".png")
        EBImage::writeImage(img, tmp)
        
        
        #Write the basic plot
        p <- plotly::plot_ly(
          source = "img_source"
        ) |>
          
          plotly::layout(
            
            dragmode = FALSE,
            
            images = list(
              list(
                
                source = base64enc::dataURI(
                  file = tmp
                ),
                
                xref = "x",
                yref = "y",
                
                x = 0,
                y = h,
                
                sizex = w,
                sizey = h,
                
                sizing = "stretch",
                
                layer = "below"
                
              )
            ),
            
            xaxis = list(
              range = c(0, w),
              visible = FALSE
            ),
            
            yaxis = list(
              range = c(0, h),
              visible = FALSE,
              scaleanchor = "x"
            )
            
          ) |>
          
          plotly::config(
            
            editable = FALSE,
            
            modeBarButtonsToAdd = list(
              "drawclosedpath"
            )
            
          )
        
        #If polygons are present then add the polygons
        if (!is.null(Reactive_value$sf_object)) {
          
          for (i in seq_len(nrow(Reactive_value$sf_object))) {
            
            coords <- sf::st_coordinates(
              Reactive_value$sf_object[i, ]
            )
            
            p <- p |>
              
              plotly::add_trace(
                
                x = coords[, 1],
                y = coords[, 2],
                
                type = "scatter",
                
                mode = "lines",
                
                fill = "toself",
                
                fillcolor = "rgba(0,0,0,0)",  # fully transparent
                
                hoveron = "fills",
                
                line = list(
                  color = "red",
                  width = 2
                ),
                
                hoverinfo = "skip",
                
                showlegend = FALSE,
                
                inherit = FALSE
                
              )
            
          }
          
          #If erase mode is on
          if (Reactive_value$mode == "delete") {
            
            p <- p |>
              plotly::add_markers(
                data = Reactive_value$sf_centroids,
                x = ~X,
                y = ~Y,
                key = ~id,
                text = "Click to delete",
                hoverinfo = "text",
                marker = list(
                  color = "red",
                  size = 14,
                  symbol = "x"
                )
              )
            
            #Register the clicks
            p <- plotly::event_register(p, event = 'plotly_click')
            
            
            #Return it without the rendering options
            return(p)
            
          }
        }
        
        
        #Modify the rendering behaviour
        p <- htmlwidgets::onRender(
          p,
          "
  function(el,x){

    var gd = document.getElementById(el.id);

    window.plotlyDrawMode = true;

    //-------------------------------------------------------
    // Only once handler installation
    //-------------------------------------------------------

    if(!window.plotlyHandlersInstalled){

      window.plotlyHandlersInstalled = true;

      Shiny.addCustomMessageHandler(
        'setDrawMode',
        function(enabled){

          window.plotlyDrawMode = enabled;

          if(enabled){

            Plotly.relayout(
              gd,
              {
                dragmode:'drawclosedpath'
              }
            );

          }else{

            Plotly.relayout(
              gd,
              {
                dragmode:false
              }
            );

          }

        }
      );

    }

    //-------------------------------------------------------
    // Initial state
    //-------------------------------------------------------

    if(window.plotlyDrawMode){

      Plotly.relayout(
        gd,
        {
          dragmode:'drawclosedpath'
        }
      );

    }

  }
  "
        )
        
        p
        
      })
      
      # -------------------------------------------------------
      # POLYGON DRAWING
      # -------------------------------------------------------
      
      observeEvent(
        #Observe an event of drawing
        plotly::event_data(
          "plotly_relayout",
          source = "img_source"
        ),
        
        #If detected
        {
          #Try it safely
          tryCatch({
            #Get the polygon
            ev <- plotly::event_data(
              "plotly_relayout",
              source = "img_source"
            )
            
            if (is.null(ev$shapes))
              return()
            #Generate the polygon outline coords from the shape
            paths <- extract_plotly_paths(
              ev$shapes
            )
            
            if (length(paths) == 0)
              return()
            
            #Generate the sf object
            polygons <- lapply(
              paths,
              function(path) {
                
                coords <- svg_path_to_matrix(
                  path
                )
                
                sf::st_polygon(
                  list(coords)
                )
                
              }
            )
            
            #generate an ID 
            if(any(is.null(Reactive_value$sf_object), nrow(Reactive_value$sf_object) == 0)) ID <- 1 else ID <- nrow(Reactive_value$sf_object) + 1 
            
            #Transform to a sf polygon object
            new_sf <- sf::st_sf(
              
              id = ID,
              
              geometry = sf::st_sfc(polygons)
              
            )
            
            #Check that the new Sf is between the image margins
            Photo_square <- 
              sf::st_as_sfc(sf::st_bbox(c(xmin = 0, ymin = 0, xmax = Reactive_value$Current_image_width, ymax = Reactive_value$Current_image_height))) %>%
              sf::st_sf()
            new_sf <- sf::st_intersection(new_sf, Photo_square)
            
            
            
            #Check that is a polygon, if not don´t do anything
            if(!all(sf::st_geometry_type(new_sf) == "POLYGON")) stop()
            
            #Check that the polygon is valid
            if(any(!sf::st_is_valid(new_sf), is.na(sf::st_is_valid(new_sf)))) stop()
            
            
            #If no previous polygons in list generate a new one
            if(any(is.null(Reactive_value$sf_object), nrow(Reactive_value$sf_object) == 0)) {
              
              Reactive_value$sf_object <- new_sf
              
              #Generate a tibble with the polygon centroids
              Centroid_tibble <- dplyr::bind_cols(id = ID,
                                                  tibble::as_tibble(sf::st_coordinates(sf::st_centroid(new_sf)))
              )
              Reactive_value$sf_centroids <- Centroid_tibble
              
              
            } else {
              #If already present add it safely (removing intersecting polygon area)
              Reactive_value$sf_object <- Add_polygon_no_overlap(sf_polygons = Reactive_value$sf_object, 
                                                                 new_polygon = new_sf)
              
              #If the new polygon generates a conflict remove it and return
              if(any(!sf::st_is_valid(Reactive_value$sf_object))){
                
                Reactive_value$sf_object <- Reactive_value$sf_object[-ID,]
                
                stop()
              }
              
              #Get the final centroid coordinates of the sf (taking into account potential erosion)
              To_bind_tibble <- dplyr::bind_cols(id = ID,
                                                 as_tibble(sf::st_coordinates(sf::st_centroid(Reactive_value$sf_object[Reactive_value$sf_object$id == ID, ]$geometry)))
              )
              
              #bind it to the already present tibble
              Reactive_value$sf_centroids <- dplyr::bind_rows(Reactive_value$sf_centroids,
                                                              To_bind_tibble)
              
              
              
            }
            
            #Return shapes to an empty list before proceeding
            plotly::plotlyProxy(
              "img",
              session
            ) |>
              plotly::plotlyProxyInvoke(
                "relayout",
                list(
                  shapes = list()
                )
              )
            
            #Continue with drawmode 
            if (Reactive_value$mode == "draw") {
              
              later::later(
                
                function() {
                  
                  session$sendCustomMessage(
                    "setDrawMode",
                    TRUE
                  )
                  
                },
                
                delay = 0.15
                
              )
              
            }
          
        },
        error = function(e) {
          #Return shapes to an empty list before proceeding
          plotly::plotlyProxy(
            "img",
            session
          ) |>
            plotly::plotlyProxyInvoke(
              "relayout",
              list(
                shapes = list()
              )
            )
          
          #Continue with drawmode 
          if (Reactive_value$mode == "draw") {
            
            later::later(
              
              function() {
                
                session$sendCustomMessage(
                  "setDrawMode",
                  TRUE
                )
                
              },
              
              delay = 0.15
              
            )
            
          }
          
          shiny::showModal(shiny::modalDialog(title = "Invalid cell mask",  easyClose = TRUE))
        })
          
          
        }
      )
      
      # -------------------------------------------------------
      # POLYGON ERASING
      # -------------------------------------------------------
      
      observeEvent(
        plotly::event_data("plotly_click", 
                           source = "img_source"),
        {
          
          #Require that the mode is delete and that the sf_object exists
          req(Reactive_value$mode == "delete")
          req(Reactive_value$sf_object)
          
          #If nothing in the sf_object list
          if(any(is.null(Reactive_value$sf_object), nrow(Reactive_value$sf_object) == 0)){
            shiny::showModal(shiny::modalDialog(title = "No cells detected",  easyClose = TRUE))
          }
          
          #else proceed
          else{
            #Get the click site
            click <- plotly::event_data(
              "plotly_click",
              source = "img_source"
            )
            #Get the point ID
            clicked_id <- click$key
            
            #Remove sf_object
            Reactive_value$sf_object <- Reactive_value$sf_object[Reactive_value$sf_object$id != clicked_id, ]
            #Remove the centroid object
            Reactive_value$sf_centroids <- Reactive_value$sf_centroids[Reactive_value$sf_centroids$id != clicked_id, ]
            
            if (nrow(Reactive_value$sf_object) > 0) {
              #Change ID for sf objects
              Reactive_value$sf_object$id <- seq_len(nrow(Reactive_value$sf_object))
              
              #Change ID for centroid objects
              Reactive_value$sf_centroids$id <- seq_len(nrow(Reactive_value$sf_object))
              
            }
          }
          
        }
      )
      
      
      # -------------------------------------------------------
      # CURRENT IMAGE INFO and STATUS
      # -------------------------------------------------------
      #Output current image
      output$Current_image_name_text <- shiny::renderUI({
        #Get the current image name
        Current_name <- Look_up_tibble$Image_name[Reactive_value$Current_Image_index]
        
        #Return the current name
        shiny::div(
          style = paste0(
            "background-color:", "#fcffff", ";",
            "color:black;",
            "font-weight:bold;",
            "padding:10px;",
            "text-align:center;",
            "border-radius:5px;"
          ),
          Current_name
        )
      })
      
      #Ouput current image status
      output$Current_image_status_text <- shiny::renderUI({
        #Get the current look_up_tibble
        Current_look_up_tibble <- Reactive_value$Look_up_tibble
        
        #Get the current image status
        Current_status <- Current_look_up_tibble$Mask_DONE[Reactive_value$Current_Image_index]
        
        bg_color <- if(Current_status) "#9df5ad" else "#f07587"
        printed_text <- if(Current_status) "DONE!" else "PENDING"
        
        
        shiny::div(
          style = paste0(
            "background-color:", bg_color, ";",
            "color:black;",
            "font-weight:bold;",
            "padding:2px;",
            "text-align:center;",
            "border-radius:2px;"
          ),
          printed_text
        )
      })
      
      #Output polygon status
      output$status <- shiny::renderUI({
        if(is.null(nrow(Reactive_value$sf_object))) N_cells <- 0 else N_cells <- nrow(Reactive_value$sf_object)
        #Return the current number of cells
        shiny::div(
          style = paste0(
            "background-color:", "#fcffff", ";",
            "color:black;",
            "font-weight:bold;",
            "padding:10px;",
            "text-align:center;",
            "border-radius:5px;"
          ),
          paste0(N_cells, " cells present in image")
        )
      })
      
      # -------------------------------------------------------
      # NAVIGATION AND WHOLE PROCESS STATUS
      # -------------------------------------------------------

      #Batch buttons
      shiny::observeEvent(input$Previous_Batch, 
                          {
                            #Get current batch
                            Current_batch <- Reactive_value$Batch_index
                            
                            #If Current_bath is 1 don't do anything
                            if(Current_batch == 1) return()
                            
                            #Else update batch
                            Reactive_value$Batch_index <- Current_batch - 1
                          })
      
      shiny::observeEvent(input$Next_Batch, 
                          {
                            #Get current batch
                            Current_batch <- Reactive_value$Batch_index
                            
                            #If Current_bath is equal to the max batches don't do anything
                            if(Current_batch == N_Batches) return()
                            
                            #Else update batch
                            Reactive_value$Batch_index <- Current_batch + 1
                          })
      
      #Batch page status
      output$Batch_status <- shiny::renderUI({
        Current_batch_index <- Reactive_value$Batch_index
        Total_batches <- N_Batches
        
        #Print message
        shiny::div(
          style = paste0(
            "background-color:", "#ffffff", ";",
            "color:black;",
            "font-weight:plain;",
            "padding:1px;",
            "text-align:center;",
            "border-radius:1px;"
          ),
          paste0("Page ", Current_batch_index, " / ", Total_batches)
        )
      })
      
      #Total progress status
      output$Progress_status <- shiny::renderUI({
        Current_done_images <- sum(Reactive_value$Look_up_tibble$Mask_DONE)
        Total_images <- nrow(Reactive_value$Look_up_tibble)
        
        Progress <- round(Current_done_images/Total_images, digits = 2) * 100
        
        #Print message
        shiny::div(
          style = paste0(
            "
      padding: 12px;
      border: 1px solid #ccc;
      border-radius: 5px;
      text-align:center;
      font-weight:bold;
      background: linear-gradient(
        to right,
        #9df5ad 0%,
        #9df5ad ", Progress, "%,
        #f07587 ", Progress, "%,
        #f07587 100%
      );
      "
          ),
          paste0(Current_done_images, " / ", Total_images, " images processed")
        )
        
      })
      
      #Output for the nagivation plot
      output$Navigation_plot <- shiny::renderPlot({
        #Get the current look_up_tibble
        Current_look_up_tibble <- Reactive_value$Look_up_tibble
        
        #Turn the current image border to black
        Current_look_up_tibble$Border_color[Reactive_value$Current_Image_index] <- "black"
        
        #Select the tiles present in the current batch
        Current_look_up_tibble <- Current_look_up_tibble %>% dplyr::filter(Batch_index == Reactive_value$Batch_index)
        #Turn the ones that are done to green
        Current_look_up_tibble$Fill_color[Current_look_up_tibble$Mask_DONE] <- "#9df5ad" 
        
        #Plot the result
        return(
          Current_look_up_tibble %>%
            ggplot(aes(x = Col_ID, y = Row_ID, fill = Fill_color, color = Border_color)) +
            geom_tile(width = 0.9, height = 0.9, linewidth = 2) +
            ggplot2::geom_point(
              alpha = 1,
              size = 1) +
            ggplot2::coord_fixed() +
            scale_x_continuous(limits = c(0, 26)) +
            scale_y_continuous(-1, max(Current_look_up_tibble$Row_ID)+1) +
            scale_fill_identity() +
            scale_color_identity() +
            ggplot2::theme(
              # Remove axis text and titles
              axis.text.x  = ggplot2::element_blank(),
              axis.text.y  = ggplot2::element_blank(),
              axis.title.x = ggplot2::element_blank(),
              axis.title.y = ggplot2::element_blank(),
              
              # Remove ticks
              axis.ticks = ggplot2::element_blank(),
              
              # Remove legend
              legend.position = "none",
              
              # White background
              panel.background = ggplot2::element_rect(
                fill = "#c2c0c0",
                colour = "#c2c0c0"
              ),
              plot.background = ggplot2::element_rect(
                fill = "#c2c0c0",
                colour = "#c2c0c0"
              ),
              
              # Remove borders and grid
              panel.border = ggplot2::element_blank(),
              panel.grid.major = ggplot2::element_blank(),
              panel.grid.minor = ggplot2::element_blank(),
              
              # Minimal margins
              plot.margin = ggplot2::margin(
                t = 0, r = 0, b = 0, l = 0,
                unit = "pt"
              )
            )
        )
      },
      bg = "#c2c0c0")
      
      #Reactive for the hovering of the navigation plot 
      hovered <- shiny::reactive({
        #Get the current look up tibble and filter in the batch
        Current_look_up_tibble <- Reactive_value$Look_up_tibble %>% dplyr::filter(Batch_index == Reactive_value$Batch_index)
        
        
        shiny::nearPoints(
          Current_look_up_tibble,
          coordinfo = input$nav_hover,
          xvar = "Col_ID",
          yvar = "Row_ID",
          maxpoints = 1
        )
      })
      
      #Output of the hovering
      output$hovered_tile <- shiny::renderUI({
        
        #If no hovering
        if(any(is.null(hovered()), nrow(hovered()) == 0)){
          shiny::div(
            style = paste0(
              "background-color:", "#c2c0c0", ";",
              "color:black;",
              "font-weight:plain;",
              "padding:5px;",
              "text-align:left;",
              "border-radius:1px;"
            ),
            paste0("")
          )
        }
        else{
          #Get image name and current status
          Image_name <- hovered()$Mask_name[1]
          Current_status <- hovered()$Mask_DONE[1]
          
          #Obtain backgroun color and message
          bg_color <- if(Current_status) "#9df5ad" else "#f07587"
          if(Current_status) Current_status <- "DONE" else Current_status <- "PENDING"
          
          #Print message
          shiny::div(
            style = paste0(
              "background-color:", bg_color, ";",
              "color:black;",
              "font-weight:plain;",
              "padding:5px;",
              "text-align:center;",
              "border-radius:1px;"
            ),
            paste0(Image_name, " - ", Current_status)
          )
        }

      })
      
      # -------------------------------------------------------
      # END OF SESSION
      # -------------------------------------------------------
      
      #If browser is closed end the app
      session$onSessionEnded(function() {
        memoise::forget(Memoised_magick_reader)
        gc()
        shiny::stopApp()
      })
      
      
    }
    
    #Run the server
    message("Always stop current R execution if you want to continue with your R session")
    shiny::shinyApp(user_interface, server)
  }
