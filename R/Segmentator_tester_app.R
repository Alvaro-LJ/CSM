#' Launches a shiny APP to explore cell segmentation parameters
#'
#' `Segmentator_tester_app()` launches an APP to interactively explore cell segmentation parameters. Parameters can then be used to feed the [Cell_segmentator_quantificator()] function.
#'
#' @param Directory Character specifying the path to the folder where images to be segmented are stored.
#' @param Ordered_Channels Character vector specifying image channels in their exact order.
#' @param Max_Gb_cache A single numeric value indicating the memory size of the image cache (see details). 10 Gb by default.
#'
#' @seealso [Cell_segmentator_quantificator()]
#'
#' @details
#' In order to deal with large images and speed up image toggling, the shiny APP works with a memoised version of the [magick::image_read()] function. Loaded images will
#' be stored in a temporary cache until the APP is closed or the cache reaches it's limit.
#'
#' Control panel controls the image display settings. If Pre process is active (default) image pre-processing will be applied to the nuclear channels before running segmentation (this can enhance performance in some scenarios).
#'
#' User can use the 'Color', 'Add channel', 'Remove channel' and 'Reset' buttons to merge different channels into a single image. The add channel button
#' will color the current channel according to the color provided and will be merged to the current image.
#'
#' Relevant buttons:
#' \itemize{
#' \item{GO!: Runs the cell segmentation and feature extraction algorithm with the current setting using the current image.}
#' \item{Download Parameters: Saves parameters in the Global environment.}
#' }
#'
#' The lower panel of the control panel shows the current active parameters.
#'
#' 4 main panels:
#' \itemize{
#' \item{Upper left: Displays the image (use it to zoom in and out).}
#' \item{Upper right: Cell boundaries are displayed.}
#' \item{Lower left: Nucleus are depicted}
#' \item{Lower right: Heatmap of marker expression by cell. If image is zoomed in, only cells in the field of view are shown.}
#' }
#'
#'
#' @examples
#' \dontrun{
#' #Create temporary input directory----------------------------------------
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
#' #Launch the app------------------------------------------------------------
#' Segmentator_tester_app(
#'     Directory = Input_Dir,
#'     Ordered_Channels = c("DAPI", "PDL1", "GZMB", "PD1", "CK-EPCAM", "CD8a", "FOXP3")
#'     )
#'
#'#Remove directories---------------------------------------------------------
#'unlink(Input_Dir, recursive = TRUE)
#'}
#' @export

Segmentator_tester_app <-
  function(Directory,
           Ordered_Channels,
           Max_Gb_cache = 10){

    ####Check suggested packages####
    {
      if(!requireNamespace("simpleSeg", quietly = TRUE)) stop(
        paste0("simpleSeg Bioconductor package is required to execute the function. Please install using the following code: ",
               expression({
                 if (!require("BiocManager", quietly = TRUE))
                   install.packages("BiocManager")

                 BiocManager::install("simpleSeg")
               })
        )
      )
      if(!requireNamespace("S4Vectors", quietly = TRUE)) stop(
        paste0("S4Vectors Bioconductor package is required to execute the function. Please install using the following code: ",
               expression({
                 if (!require("BiocManager", quietly = TRUE))
                   install.packages("BiocManager")

                 BiocManager::install("S4Vectors")
               })
        )
      )
      if(!requireNamespace("cytomapper", quietly = TRUE)) stop(
        paste0("cytomapper Bioconductor package is required to execute the function. Please install using the following code: ",
               expression({
                 if (!require("BiocManager", quietly = TRUE))
                   install.packages("BiocManager")

                 BiocManager::install("cytomapper")
               })
        )
      )
      if(!requireNamespace("ComplexHeatmap", quietly = TRUE)) stop(
        paste0("ComplexHeatmap Bioconductor package is required to execute the function. Please install using the following code: ",
               expression({
                 if (!require("BiocManager", quietly = TRUE))
                   install.packages("BiocManager")

                 BiocManager::install("ComplexHeatmap")
               })
        )
      )
      if(!requireNamespace("magick", quietly = FALSE)) stop(
        paste0("magick CRAN package is required to execute the function. Please install using the following code: ",
               expression(install.packages("benchmarkme")))
      )
    }

    ####CHECK ARGUMENTS and generate memoise magick::read_image####
    #check that the directory provided contains at least one file
    if(length(dir(Directory)) <1 ) stop("No files found in the Directory provided")
    if(is.null(Ordered_Channels)) stop("Ordered_Channels must not be NULL")
    #Check that ordered channels and channels to keep are unique
    if(length(Ordered_Channels) != length(unique(Ordered_Channels))) stop("Ordered_Channels must contain non-duplicated character values")

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


    #Generate a memoised version of the cell segmentation function
    Shiny_segmentator_function <-
      function(Image,
               Nuclear_pre_processing,

               Equalize = FALSE,
               Photo_min = 0,
               Photo_max = 100,
               Photo_gamma = 0,
               Opening = 1,
               Closing = 1,

               Ordered_Channels,
               Channels_to_keep,
               Nuclear_marker,
               Cell_body_method,
               Min_pixel,
               Smooth_amount,
               Normalization,
               Watershed_type,
               Tolerance_value = NULL,
               Neighborhood_distance,
               Disc_size,
               Tissue_mask_markers,
               Perform_PCA = FALSE){
        #First we need to import the photograph as EBI, select the required channels and transform it to a cytomapper object
        shiny::showModal(modalDialog("Importing image", footer=NULL))
        Image <- Image
        shiny::removeModal()

        #Mopdify nuclear channels if required by user
        if(as.logical(Nuclear_pre_processing)){
          shiny::showModal(modalDialog("Performing Nuclear channels image Pre-processing", footer=NULL))
          Nuclear_channels_number <- match(Nuclear_marker, Ordered_Channels)

          #Apply changes to every nuclear channel
          for(index in Nuclear_channels_number){
            Image_Modified <- Image[index]

            if(Equalize) Image_Modified <- Image_Modified %>% magick::image_equalize() #Equalize if necessary
            Image_Modified <- Image_Modified %>% magick::image_level(black_point = Photo_min,
                                                                     white_point = Photo_max,
                                                                     mid_point = 10^Photo_gamma) #Change white, black and gamma
            Image_Modified <- Image_Modified %>% magick::as_EBImage() #turn to EBImage object
            Image_Modified  <-
              Image_Modified %>% EBImage::opening(EBImage::makeBrush(size = Opening, shape = "disc")) %>%
              EBImage::closing(EBImage::makeBrush(size = Closing, shape = "disc")) #opening and closing

            Image_Modified <- magick::image_read(Image_Modified) #again as Magick image

            Image[index] <- Image_Modified
          }
          shiny::removeModal()
        }

        #Produce an EBImage object
        Image <- Image %>% magick::as_EBImage()


        Image <- cytomapper::CytoImageList(Image)#Transform it to cytoImage object
        cytomapper::channelNames(Image) <- Ordered_Channels #define channel names
        S4Vectors::mcols(Image)$img_id <- as.character(" ")#Modify name
        Image <- cytomapper::getChannels(Image, Channels_to_keep) #Keep only user defined channels

        #Perform cell segmentation
        shiny::showModal(modalDialog("Generating segmentation mask. This can take some time. Please wait", footer=NULL))
        Tolerance <- if(is.na(Tolerance_value)) NULL else(Tolerance_value)
        Seg_results <- simpleSeg::simpleSeg(Image,
                                            nucleus = Nuclear_marker,
                                            cellBody = Cell_body_method,
                                            sizeSelection = Min_pixel,
                                            smooth = Smooth_amount,
                                            transform = Normalization,
                                            watershed = Watershed_type,
                                            tolerance = Tolerance_value,
                                            ext = Neighborhood_distance,
                                            discSize = Disc_size,
                                            tissue = Tissue_mask_markers,
                                            pca = Perform_PCA,
                                            cores = 1
        )
        S4Vectors::mcols(Seg_results)$img_id <- as.character(" ")
        shiny::removeModal()

        #Generate a list with the final parameters
        Parameter_list <- list(
          Ordered_Channels = Ordered_Channels,
          Channels_to_keep = Channels_to_keep,
          Nuclear_marker = Nuclear_marker,
          Cell_body_method = Cell_body_method,
          Watershed_type = Watershed_type,
          Normalization = Normalization,
          Tissue_mask_markers = Tissue_mask_markers,
          Smooth_amount = Smooth_amount,
          Min_pixel = Min_pixel,
          Neighborhood_distance = Neighborhood_distance,
          Disc_size = Disc_size,
          Tolerance_value = Tolerance_value,
          Perform_PCA = Perform_PCA,
          Perform_nuclear_channel_processing = Nuclear_pre_processing,
          Black_level = Photo_min,
          White_level = Photo_max,
          Gamma_level = Photo_gamma,
          Equalize = Equalize,
          Opening_kernel_size = Opening,
          Closing_kernel_size = Closing
        )

        #Obtain the first nuclear channel of the image (for the cell limit plot)
        Nuclear_channel_image <- cytomapper::getChannels(Image, Nuclear_marker[1])
        S4Vectors::mcols(Nuclear_channel_image)$img_id <- as.character(" ")
        cytomapper::channelNames(Nuclear_channel_image) <- Nuclear_marker[1]

        #Return all the results in a list
        return(list(Nuclear_image = Nuclear_channel_image,
                    Mask = Seg_results,
                    Segmentation_Parameters = Parameter_list))
      }



    #Obtain image names and channel names
    Real_Images <- dir(Directory, full.names = FALSE)
    Channels_in_images <- Ordered_Channels


    ####BUILD THE USER INTERFACE####
    {
      user_interface <- shiny::fluidPage(

        #To make the sidebar collapsable
        shinyjs::useShinyjs(),

        #Set the title and action button
        shiny::fluidRow(
          shiny::column(width = 3, shiny::h3("Segmentation exploration APP")),
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
              shiny::actionButton("toggle_preprocess", "Nuclear Preprocessing", icon = shiny::icon(name = "square-caret-up", lib = "font-awesome")),
              shiny::actionButton("toggle_segmentation", "Segmentation parameters", icon = shiny::icon(name = "square-caret-up", lib = "font-awesome"))
            ),


            #IMAGE PARAMETERS
            shiny::h4("Image settings"),
            #Collapsable rows
            shiny::tags$div(
              class = "Image_settings",

              #Image to display
              shiny::fluidRow(
                #Select the real image to be displayed
                shiny::column(8, shiny::selectInput("Real_Image_name", "Image to display", sort(Real_Images), multiple = FALSE)),
                #Select the channel to be displayed
                shiny::column(4, shiny::selectInput("Channel", "Channel to display", Channels_in_images, multiple = FALSE))

              ),

              #Color merging
              shiny::fluidRow(
                shiny::column(3, shiny::textInput("Color_channel", "Color", value = "blue")),
                shiny::column(3, shiny::actionButton("Add_color", "Add channel", shiny::icon("plus", library = "font-awesome"))),
                shiny::column(3, shiny::actionButton("Remove_color", "Remove channel", shiny::icon("minus", library = "font-awesome"))),
                shiny::column(3, shiny::actionButton("Reset_image", "Reset", icon = shiny::icon("redo")))
              )
            ),

            #Nuclear pre processing
            #Optional nuclear-prepocessing parameters
            shiny::h4("Nuclear pre-processing"),
            #Collapsable rows
            shiny::tags$div(
              class = "Nuclear_settings",
              shiny::fluidRow(
                #Should this be performed or not
                shiny::column(3, shinyWidgets::materialSwitch("Pre_processing", "Pre process", value = TRUE)),
                #Select the kernels
                shiny::column(3, shiny::conditionalPanel(condition = "input.Pre_processing == '1'",
                                                         shiny::numericInput("Opening", "Opening", value = 1, min = 1, max = NA, step = 1)
                )),
                shiny::column(3, shiny::conditionalPanel(condition = "input.Pre_processing == '1'",
                                                         shiny::numericInput("Closing", "Closing", value = 1, min = 1, max = NA, step = 1)
                )),
                #Select the equalization
                shiny::column(3, shiny::selectInput("Equalize", "Equalize", c(YES = TRUE, NO = FALSE), selected = FALSE, multiple = FALSE))
              ),

              #Image control
              shiny::fluidRow(
                #Select the min point of the image
                shiny::column(4, shiny::sliderInput("Min_Image", "Black", value = 0, min = 0, max = 100, step = 1, ticks = FALSE)),
                #Select the max point of the image
                shiny::column(4, shiny::sliderInput("Max_Image", "White", value = 100, min = 0, max = 100, step = 1, ticks = FALSE)),
                #Select the Gamma
                shiny::column(4, shiny::sliderInput("Gamma", "Gamma", value = 0, min = -3, max = +3, step = 0.01, ticks = FALSE))
              )
            ),


            #Segmentation parameters
            shiny::h4("Segmentation parameters"),
            #Basic segmentation parameters
            #Collapsable rows
            shiny::tags$div(
              class = "Segmentation_settings",
              shiny::fluidRow(
                #Select the channels to keep
                shiny::column(4, shinyWidgets::virtualSelectInput("Keep_channels", label = "Channels to keep",
                                                                  choices = Channels_in_images,
                                                                  search = TRUE,
                                                                  multiple = TRUE
                )),
                #Select the nuclear marker
                shiny::column(4, shinyWidgets::virtualSelectInput("Nuclear", "Nuclear Marker", choices = Channels_in_images, multiple = TRUE, search = TRUE)),
                #Select the cell body method
                shiny::column(4, shiny::selectInput("Cell_body", "Cell-body method", c("none", "dilate", "discModel"), selected = "discModel", multiple = FALSE))
              ),
              #Image pre-processing parameters
              shiny::fluidRow(
                #Select the watershed type
                shiny::column(4, shiny::selectInput("Watershed", "Watershed type", c("intensity", "distance", "combine"), selected = "combine", multiple = FALSE)),
                #Select the image normalization steps
                shiny::column(4, shinyWidgets::pickerInput("Normalization", label = "Normalization",
                                                           choices = c("sqrt", "asinh", "norm99", "maxThresh", "tissueMask"),
                                                           multiple = TRUE, width = "fit",
                                                           options = shinyWidgets::pickerOptions(
                                                             selectedTextFormat = 'count',
                                                             showSubtext = FALSE)
                )),
                #Select the markers for tissueMask
                shiny::column(4, shinyWidgets::pickerInput("Tissue_mask", label = "Tissue-mask markers",
                                                           choices = Channels_in_images,
                                                           multiple = TRUE, width = "fit",
                                                           options = shinyWidgets::pickerOptions(
                                                             selectedTextFormat = 'count',
                                                             showSubtext = FALSE)
                ))

              ),
              #Cell identification parameters
              shiny::fluidRow(
                #Min pixel
                shiny::column(4, shiny::numericInput("Min_pix", "Min pixel value", value = 10, min = 1, max = NA, step = 1)),
                #Neighbor distance
                shiny::column(4, shiny::numericInput("Neigh_dist", "Neighbor distance", value = 1, min = 1, max = NA, step = 1)),
                #Select the disc size
                shiny::column(4, shiny::numericInput("Disc_size", "Disc size", value = 1, min = 1, max = NA, step = 1))
              ),
              #Minor parameters
              shiny::fluidRow(
                #Tolerance
                shiny::column(4, shiny::numericInput("Tolerance", "Tolerance", value = NULL, min = 1, max = NA, step = 1)),
                #Perform PCA
                shiny::column(4, shiny::selectInput("PCA", "Nuclear PCA", c(FALSE, TRUE), selected = FALSE, multiple = FALSE)),
                #Select the smooth amount
                shiny::column(4, shiny::numericInput("Gaussian", "Smoothening", value = 1, min = 0, max = NA, step = 0.1))
              )
            ),

            #BUTTONS
            shiny::fluidRow(
              #Action buttons
              shiny::column(2, shiny::actionButton("GO_button", shiny::icon("bolt-lightning"), label = "GO!")),
              shiny::column(4, shiny::actionButton("Download_Param", shiny::icon("download"), label = "Download Parameters"))
            )
          ),

          #Set the outcome columns
          shiny::mainPanel(
            #Give the main panel an ID
            id = "mainPanel",

            #First row will have the Photo and the overview of marker intensity by cell
            shiny::fluidRow(
              shiny::column(6, shiny::plotOutput("Photo",
                                                 width = "auto",
                                                 #Controls for zoom in
                                                 dblclick = "Photo_dblclick",
                                                 brush = shiny::brushOpts(id = "Photo_brush",
                                                                          resetOnNew = TRUE)
              )
              ),
              shiny::column(6, shiny::plotOutput("Cell_borders",
                                                 width = "auto"))
            ),
            #Second row will contain the positive cells and the histogram
            shiny::fluidRow(
              shiny::column(6, shiny::plotOutput("Cell_surface",
                                                 width = "auto")),

              #Parameters
              shiny::column(6, shiny::verbatimTextOutput("Parameters"))

            )
          )
        ),
        shiny::tags$head(shiny::tags$style(
          htmltools::HTML('
         #sidebar {
            background-color: #acf2ae;
        }

        body, label, input, button, select {
          font-family: "Arial";
        }'))),

        shiny::tags$script(htmltools::HTML("
    $(document).on('click', '#toggle_Image', function() {
      $('.Image_settings').toggle();
    });
    $(document).on('click', '#toggle_preprocess', function() {
      $('.Nuclear_settings').toggle();
    });
    $(document).on('click', '#toggle_segmentation', function() {
      $('.Segmentation_settings').toggle();
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
      #Generate a reactive with the real Image name and the channel number and pre-processing steps
      Photo_name <- shiny::reactive(stringr::str_c(Directory, "/", input$Real_Image_name))
      Channel_index <- shiny::reactive(which(input$Channel == Ordered_Channels))
      Photo_min <- shiny::reactive(input$Min_Image)
      Photo_max <- shiny::reactive(input$Max_Image)
      Photo_gamma <- shiny::reactive(10^input$Gamma)
      Equalize <- shiny::reactive(input$Equalize)
      Overlay <- shiny::reactive(input$Overlay)

      Pre_process <- shiny::reactive(input$Pre_processing)
      Opening_kernel <- shiny::reactive(input$Opening)
      Closing_kernel <- shiny::reactive(input$Closing)


      #Generate a reactivevalue that controls the image
      Segmentation_results <- shiny::reactiveValues(Nuclear_image = NULL,
                                                    Mask = NULL,
                                                    Segmentation_Parameters = NULL)

      #Generate a reactive that controls the zoom in
      ranges <- shiny::reactiveValues(x = NULL, y = NULL)

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

      #Reactive that imports the photograph and will serve for all the images
      Photo_reactive <- shiny::reactive({

        #Import the Photo
        Photo <- Memoised_magick_reader(Photo_name())
        Photo_list <- list(Image = Photo,
                           Original_Dims = c(magick::image_info(Photo)$width, magick::image_info(Photo)$height))

        #Crop the photo if required
        if(any(!is.null(ranges$x), !is.null(ranges$y))){


          Geometry_crop <- stringr::str_c(max(ranges$x) - min(ranges$x),
                                          "x",
                                          max(ranges$y) - min(ranges$y),
                                          "+",
                                          min(ranges$x),
                                          "+",
                                          min(ranges$y))

          Photo_list$Image <- magick::image_crop(Photo_list$Image,
                                                 geometry = Geometry_crop,
                                                 gravity = "SouthWest")
          Photo_list$X_crop_coords <- ranges$x
          Photo_list$Y_crop_coords <- ranges$y
        }


        #Return the photo result in a list

        return(Photo_list)

      })

      #PHOTO AS GGPLOT OBJECT (with adequate axis) that will serve as basis for all the plots with a photo on the background
      Photo_plot_reactive <- shiny::reactive({

        #Obtain the Image
        Photo <- Photo_reactive()$Image

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
          if(as.logical(Pre_process())){
            Photo <- Photo %>% magick::as_EBImage() %>%
              EBImage::opening(EBImage::makeBrush(size = Opening_kernel(), shape = "disc")) %>%
              EBImage::closing(EBImage::makeBrush(size = Closing_kernel(), shape = "disc"))
            Photo <- magick::image_read(Photo)
          }
        }

        #If the Image_coloring_list$Parameter_list is present then use the Channel merger
        if(!is.null(Image_coloring_list$Parameter_list)){
          if(as.logical(Pre_process())){
            Photo <- Photo %>% magick::as_EBImage() %>%
              EBImage::opening(EBImage::makeBrush(size = Opening_kernel(), shape = "disc")) %>%
              EBImage::closing(EBImage::makeBrush(size = Closing_kernel(), shape = "disc"))
            Photo <- magick::image_read(Photo)
          }
          #Perform RGB color merging
          Photo <- Image_channel_color_merger(Image = Photo,
                                              Parameter_list = Image_coloring_list$Parameter_list)
        }

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
        #Obtain the photo
        Photo_plot <- Photo_plot_reactive()

        #If no channel merging required
        if(is.null(Image_coloring_list$Parameter_list)) return(Photo_plot)

        #If channel merging required print the legend
        if(!is.null(Image_coloring_list$Parameter_list)){
          #Generate the color tibble to generate the legend
          color_tibble <-
            tibble(Channel_index = 1:length(Image_coloring_list$Parameter_list),
                   Channel_name = names(Image_coloring_list$Parameter_list),
                   Color_names = purrr::map_chr(Image_coloring_list$Parameter_list, ~.[["color"]])
            ) %>% arrange(Channel_index)
          color_tibble$Channel_name = factor(color_tibble$Channel_name, levels = color_tibble$Channel_name)

          return(
            Photo_plot +
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
      })

      #Control the zoom in of the Photo and the other plots
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

      #What to do when user hits the go button
      shiny::observeEvent(input$GO_button, {
        Photo <- Photo_reactive()$Image

        Segmentation_result_list <- Shiny_segmentator_function(Image = Photo,
                                                               Nuclear_pre_processing = as.logical(input$Pre_processing),

                                                               Equalize = as.logical(input$Equalize),
                                                               Photo_min = input$Min_Image,
                                                               Photo_max = input$Max_Image,
                                                               Photo_gamma = input$Gamma,
                                                               Opening = input$Opening,
                                                               Closing = input$Closing,

                                                               Ordered_Channels = Channels_in_images,
                                                               Channels_to_keep = input$Keep_channels,
                                                               Nuclear_marker = input$Nuclear,
                                                               Cell_body_method = input$Cell_body,
                                                               Min_pixel = input$Min_pix,
                                                               Smooth_amount = input$Gaussian,
                                                               Normalization = input$Normalization,
                                                               Watershed_type = input$Watershed,
                                                               Tolerance_value = as.numeric(input$Tolerance),
                                                               Neighborhood_distance = input$Neigh_dist,
                                                               Disc_size = input$Disc_size,
                                                               Tissue_mask_markers = input$Tissue_mask,
                                                               Perform_PCA = as.logical(input$PCA)
        )

        #Export the results to the reactive values object
        Segmentation_results$Nuclear_image <- Segmentation_result_list$Nuclear_image
        Segmentation_results$Mask <- Segmentation_result_list$Mask
        Segmentation_results$Segmentation_Parameters <- Segmentation_result_list$Segmentation_Parameters
      })


      #First output cell borders
      output$Cell_borders <- shiny::renderPlot({
        if(is.null(Segmentation_results[["Mask"]])) ggplot()
        else{
          #Load the nuclear image
          Image <- Segmentation_results$Nuclear_image

          #Load the mask
          Seg_results <- Segmentation_results$Mask[[1]]
          Seg_results <- cytomapper::CytoImageList(Seg_results)
          S4Vectors::mcols(Seg_results)$img_id <- as.character(" ")

          #Plot the actual plot
          color_list <- list(A = c(scales::alpha("black", 0), "#0549fc"))
          names(color_list) <- "DAPI"
          cytomapper::plotPixels(image = Image,
                                 mask = Seg_results,
                                 img_id = "img_id",
                                 colour_by = input$Nuclear[1],
                                 colour = color_list,
                                 display = "single",
                                 legend = NULL)
        }
      })
      #Second the cell surface
      output$Cell_surface <- shiny::renderPlot({
        if(is.null(Segmentation_results[["Mask"]])) ggplot()
        else{
          EBImage::image(EBImage::flip(EBImage::colorLabels(Segmentation_results$Mask[[1]])))
        }
      })
      #Finally plot the parameters
      output$Parameters <- shiny::renderPrint({
        if(is.null(Segmentation_results$Segmentation_Parameters)) "Awaiting initial test"
        else{Segmentation_results$Segmentation_Parameters}
      })

      #Download segmentation parameters to R session if required
      observeEvent(input$Download_Param, {
        if(is.null(Segmentation_results$Segmentation_Parameters)){
          showModal(modalDialog(
            "Segmentation Parameters not found. Please run a test before re-trying",
            easyClose = TRUE,
            footer = NULL
          )
          )
        }
        else{
          Segmentation_Parameters <<- Segmentation_results$Segmentation_Parameters
          showModal(modalDialog(
            "An object called 'Segmentation_Parameters' has been created in the Global environment.",
            easyClose = TRUE,
            footer = NULL
          )
          )
        }

      })


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




