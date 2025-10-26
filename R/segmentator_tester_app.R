#' Launches a shiny APP to explore cell segmentation parameters
#'
#' `Segmentator_tester_app()` launches an APP to interactively explore cell segmentation parameters interactively
#' Parameters can then be used to feed the [Cell_segmentator_quantificator()].
#' @param Directory Character specifying the path to the folder where images to be segmented are present.
#' @param Ordered_Channels Character vector specifying image channels in their exact order.
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
#'    EBImage::writeImage(CSM_MiniMultiTiff_test[[Image]], file.path(Input_Dir, names(CSM_MiniMultiTiff_test)[Image]))
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
  function(Directory = NULL,
           Ordered_Channels = NULL){

    #Check suggested packages
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

    #check that the directory provided contains at least one file
    if(length(dir(Directory)) <1 ) stop("No files found in the Directory provided")
    if(is.null(Ordered_Channels)) stop("Ordered_Channels must not be NULL")

    #Obtain image names and channel names
    Real_Images <- dir(Directory, full.names = FALSE)
    Channels_in_images <- Ordered_Channels


    #BUILD THE USER INTERFACE
    user_interface <- shiny::fluidPage(

      #Set the title
      shiny::titlePanel("Segmentation exploration APP"),

      #We want a two panel layout, one in the left containing the input parameters and the output in the right
      shiny::sidebarLayout(
        #Set the first column (which contains the user defined parameters)
        shiny::sidebarPanel(
          #ID and width
          id="sidebar",

          #Image to display
          shiny::fluidRow(
            #Select the real image to be displayed
            shiny::column(5, shiny::selectInput("Real_Image_name", "Image to display", sort(Real_Images), multiple = FALSE)),
            #Select the channel to be displayed
            shiny::column(4, shiny::selectInput("Channel", "Channel to display", Channels_in_images, multiple = FALSE)),
            #Select the channels to keep
            shiny::column(3, shinyWidgets::virtualSelectInput("Keep_channels", label = "Channels to keep",
                                                              choices = Channels_in_images,
                                                              search = TRUE,
                                                              multiple = TRUE
            ))
          ),
          #Image control
          shiny::fluidRow(
            #Select the min point of the image
            shiny::column(3, shiny::sliderInput("Min_Image", "Black", value = 0, min = 0, max = 100, step = 1)),
            #Select the max point of the image
            shiny::column(3, shiny::sliderInput("Max_Image", "White", value = 100, min = 0, max = 100, step = 1)),
            #Select the Gamma
            shiny::column(3, shiny::sliderInput("Gamma", "Gamma", value = 0, min = -3, max = +3, step = 0.01)),
            #Select the equalization
            shiny::column(3, shiny::selectInput("Equalize", "Equalize", c(YES = TRUE, NO = FALSE), selected = FALSE, multiple = FALSE))
          ),

          #Optional nuclear-prepocessing parameters
          shiny::fluidRow(
            #Should this be performed or not
            shiny::column(3, shinyWidgets::materialSwitch("Pre_processing", "Pre process", value = TRUE)),
            #Select the kernels
            shiny::column(3, shiny::conditionalPanel(condition = "input.Pre_processing == '1'",
                                                     shiny::numericInput("Opening", "Opening", value = 1, min = 1, max = NA, step = 1)
            )),
            shiny::column(3, shiny::conditionalPanel(condition = "input.Pre_processing == '1'",
                                                     shiny::numericInput("Closing", "Closing", value = 1, min = 1, max = NA, step = 1)
            ))
          ),

          #Basic segmentation parameters
          shiny::fluidRow(
            #Select the nuclear marker
            shiny::column(4, shinyWidgets::virtualSelectInput("Nuclear", "Nuclear Marker", choices = Channels_in_images, multiple = TRUE, search = TRUE)),
            #Select the cell body method
            shiny::column(4, shiny::selectInput("Cell_body", "Cell-body method", c("none", "dilate", "discModel"), selected = "discModel", multiple = FALSE)),
            #Select the watershed type
            shiny::column(4, shiny::selectInput("Watershed", "Watershed type", c("intensity", "distance", "combine"), selected = "combine", multiple = FALSE))
          ),
          #Image pre-processing parameters
          shiny::fluidRow(
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
            )),
            #Select the smooth amount
            shiny::column(4, shiny::numericInput("Gaussian", "Smoothening", value = 1, min = 0, max = NA, step = 0.1))
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
            shiny::column(4, shiny::numericInput("Tolerance", "Tolerance ", value = NULL, min = 1, max = NA, step = 1)),
            #Perform PCA
            shiny::column(4, shiny::selectInput("PCA", "Nuclear PCA", c(FALSE, TRUE), selected = FALSE, multiple = FALSE)),
            #Switch for cell overlay
            shiny::column(4, shinyWidgets::materialSwitch("Overlay", "Cell overlay", value = TRUE))
          ),

          shiny::fluidRow(
            #Action buttons
            shiny::column(2, shiny::actionButton("GO_button", shiny::icon("bolt-lightning"), label = "GO!")),
            shiny::column(4, shiny::actionButton("Download_Param", shiny::icon("download"), label = "Download Parameters"))
          ),
          shiny::fluidRow(
            #Parameters
            shiny::column(12, shiny::verbatimTextOutput("Parameters"))
          )
        ),

        #Set the outcome columns
        shiny::mainPanel(
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
            shiny::column(6, shiny::plotOutput("Heatmap"))
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
        }')))
    )

    #BUILD THE SERVER
    server <- function(input, output, session){
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
      Segmentation_results <- shiny::reactiveValues(Image = NULL,
                                                    Mask = NULL,
                                                    Cell_data = NULL,
                                                    Segmentation_Parameters = NULL)

      #Generate a reactive that controls the zoom in
      ranges <- shiny::reactiveValues(x = NULL, y = NULL)

      #Reactive that imports the photograph
      Photo_reactive <- shiny::reactive({
        #Import the Photo
        Photo <- magick::image_read(Photo_name())[Channel_index()]
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

        #Transform to raster and generate plot
        Photo <- Photo %>% magick::image_raster()
        #Return the result
        Photo
      })

      #Print the photo
      output$Photo <- shiny::renderPlot({
        #Obtain the photo
        Photo <- Photo_reactive()

        #If no cell data are available just plot the photo
        if(is.null(Segmentation_results$Cell_data) | !as.logical(Overlay())){
          #Generate the plot if no data are available
          Photo_plot <- ggplot() + geom_raster(aes(x = x, y = y, fill = col), data = Photo) +
            scale_fill_identity() +
            coord_cartesian(xlim = ranges$x, ylim = ranges$y, expand = FALSE)+
            theme(axis.title = element_blank(),
                  axis.text = element_blank(),
                  axis.ticks = element_blank(),
                  axis.line = element_blank(),
                  panel.background = element_rect(fill = "black"),
                  panel.grid = element_blank())
        }

        else{
          #Generate the plot if data is available
          Photo_plot <- ggplot() + geom_raster(aes(x = x, y = y, fill = col), data = Photo) +
            scale_fill_identity() +
            coord_cartesian(xlim = ranges$x, ylim = ranges$y, expand = FALSE)+
            geom_point(aes(x = m.cx, y = m.cy), size = 0.1, color = "red", data = Segmentation_results$Cell_data) +
            theme(axis.title = element_blank(),
                  axis.text = element_blank(),
                  axis.ticks = element_blank(),
                  axis.line = element_blank(),
                  panel.background = element_rect(fill = "black"),
                  panel.grid = element_blank())
        }
        Photo_plot
      }, res = 300)

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
        #First we need to import the photograph as EBI, select the required channels and transform it to a cytomapper object
        shiny::showModal(modalDialog("Importing image", footer=NULL))
        Image <- EBImage::readImage(stringr::str_c(Directory, "/", input$Real_Image_name))
        shiny::removeModal()

        #Mopdify nuclear channels if required by user
        if(as.logical(input$Pre_processing)){
          shiny::showModal(modalDialog("Performing Nuclear channels image Pre-processing", footer=NULL))
          Image <- magick::image_read(Image)
          Nuclear_channels_number <- match(input$Nuclear, Ordered_Channels)

          #Apply changes to every nuclear channel
          for(index in Nuclear_channels_number){
            Image_Modified <- Image[index]

            if(as.logical(input$Equalize)) Image_Modified <- Image_Modified %>% magick::image_equalize() #Equalize if necesary
            Image_Modified <- Image_Modified %>% magick::image_level(black_point = input$Min_Image,
                                                                     white_point = input$Max_Image,
                                                                     mid_point = 10^input$Gamma) #Chane withe, black and gamma
            Image_Modified <- Image_Modified %>% magick::as_EBImage() #turn to EBImage object
            Image_Modified  <-
              Image_Modified %>% EBImage::opening(EBImage::makeBrush(size = input$Opening, shape = "disc")) %>%
              EBImage::closing(EBImage::makeBrush(size = input$Closing, shape = "disc")) #opening and closing

            Image_Modified <- magick::image_read(Image_Modified) #again as Magick image

            Image[index] <- Image_Modified
          }
          #Returna a EBImage object
          Image <- Image %>% magick::as_EBImage()
          shiny::removeModal()
        }


        Image <- cytomapper::CytoImageList(Image)#Transform it to cytoImage object
        cytomapper::channelNames(Image) <- Ordered_Channels #define channel names
        S4Vectors::mcols(Image)$img_id <- as.character(" ")#Modify name
        Image <- cytomapper::getChannels(Image, input$Keep_channels) #Keep only user defined channels
        Segmentation_results$Image <- Image

        #Perform cell segmentation
        shiny::showModal(modalDialog("Generating segmentation mask. This can take some time. Please wait", footer=NULL))
        Tolerance <- if(is.na(input$Tolerance)) NULL else(input$Tolerance)
        Seg_results <- simpleSeg::simpleSeg(Image,
                                            nucleus = input$Nuclear,
                                            cellBody = input$Cell_body,
                                            sizeSelection = input$Min_pix,
                                            smooth = input$Gaussian,
                                            transform = input$Normalization,
                                            watershed = input$Watershed,
                                            tolerance = Tolerance,
                                            ext = input$Neigh_dist,
                                            discSize = input$Disc_size,
                                            tissue = input$Tissue_mask,
                                            pca = as.logical(input$PCA),
                                            cores = 1
        )
        S4Vectors::mcols(Seg_results)$img_id <- as.character(" ")
        Segmentation_results$Mask <- Seg_results
        shiny::removeModal()

        #Obtain cell data
        shiny::showModal(modalDialog("Retrieving cell data", footer=NULL))
        Cells <- cytomapper::measureObjects(mask = Seg_results,
                                            image = Image,
                                            img_id = "img_id",
                                            feature_types = c("basic", "moment"),
                                            moment_feature = c("cx", "cy"),
                                            basic_feature = "mean")

        Cells_tibble <- as_tibble(cbind(as_tibble(SummarizedExperiment::colData(Cells)),
                                        as_tibble(t(SummarizedExperiment::assays(Cells)[[1]])))) %>% dplyr::select(-objectNum)
        Segmentation_results$Cell_data <- Cells_tibble
        shiny::removeModal()

        #Generate a list with the final parameters
        Parameter_list <- list(
          Ordered_Channels = Channels_in_images,
          Channels_to_keep = input$Keep_channels,
          Nuclear_marker = input$Nuclear,
          Cell_body_method = input$Cell_body,
          Watershed_type = input$Watershed,
          Normalization = input$Normalization,
          Tissue_mask_markers = input$Tissue_mask,
          Smooth_amount = input$Gaussian,
          Min_pixel = input$Min_pix,
          Neighborhood_distance = input$Neigh_dist,
          Disc_size = input$Disc_size,
          Tolerance_value = Tolerance,
          Perform_PCA = as.logical(input$PCA),

          Perform_nuclear_channel_processing = as.logical(input$Pre_processing),
          Black_level = input$Min_Image,
          White_level = input$Max_Image,
          Gamma_level = input$Gamma,
          Equalize = as.logical(input$Equalize),
          Opening_kernel_size = input$Opening,
          Closing_kernel_size = input$Closing
        )
        Segmentation_results$Segmentation_Parameters <- Parameter_list
      })

      #First output cell borders
      output$Cell_borders <- shiny::renderPlot({
        if(is.null(Segmentation_results[["Image"]])) ggplot()
        else{
          #If no zoom in required
          if(is.null(ranges$x) | is.null(ranges$y)){
            Image <- EBImage::flip(Segmentation_results$Image[[1]])
            Image <- cytomapper::CytoImageList(Image)
            Image <- cytomapper::getChannels(Image, input$Nuclear[1])
            S4Vectors::mcols(Image)$img_id <- as.character(" ")
            cytomapper::channelNames(Image) <- input$Nuclear[1]

            Seg_results <- EBImage::flip(Segmentation_results$Mask[[1]])
            Seg_results <- cytomapper::CytoImageList(Seg_results)
            S4Vectors::mcols(Seg_results)$img_id <- as.character(" ")
          }
          #If user zooms in
          else{
            Image <- Segmentation_results$Image[[1]][ranges$x[1]:ranges$x[2], ranges$y[1]:ranges$y[2], 1]
            Image <- EBImage::flip(Image)
            Image <- cytomapper::CytoImageList(Image)
            S4Vectors::mcols(Image)$img_id <- as.character(" ")
            cytomapper::channelNames(Image) <- input$Nuclear[1]

            Seg_results <- Segmentation_results$Mask[[1]][ranges$x[1]:ranges$x[2], ranges$y[1]:ranges$y[2]]
            Seg_results <- EBImage::flip(Seg_results)
            Seg_results <- cytomapper::CytoImageList(Seg_results)
            S4Vectors::mcols(Seg_results)$img_id <- as.character(" ")
          }



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
        if(is.null(Segmentation_results[["Image"]])) ggplot()
        else{
          #If no zoom in required
          if(is.null(ranges$x) | is.null(ranges$y)){
            Seg_results <- Segmentation_results$Mask[[1]]
          }
          #If user zooms in
          else{
            Seg_results <- Segmentation_results$Mask[[1]][ranges$x[1]:ranges$x[2], ranges$y[1]:ranges$y[2]]
          }
          EBImage::image(EBImage::colorLabels(Seg_results))
        }
      })
      #Finally plot the parameters
      output$Parameters <- shiny::renderPrint({
        if(is.null(Segmentation_results$Image)) "Awaiting initial test"
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

      #Lets plot the Heatmap of selected cells
      output$Heatmap <- shiny::renderPlot({
        if(is.null(Segmentation_results[["Cell_data"]])) tibble(a = NA, b = NA)
        else{
          Data <- Segmentation_results$Cell_data[1:ncol(Segmentation_results$Cell_data)]

          #Scale our data and bind it to the original data
          Data_scaled <- Data %>% dplyr::select(-c(1:5)) %>% scale()
          Data <- dplyr::bind_cols(Data["m.cx"], Data["m.cy"], Data_scaled)

          #Prepare the min and max for our heatmap
          HEATMAP_MATRIX <- Data %>%dplyr::select(-c(1:2)) %>% unlist()
          Min_HEATMAP <- quantile(HEATMAP_MATRIX, 0.05)
          Max_HEATMAP <- quantile(HEATMAP_MATRIX, 0.95)

          #Now we will select cells contained within area
          Min_x <- if(is.null(ranges$x)) min(Data$m.cx) else(ranges$x[1])
          Max_x <- if(is.null(ranges$x)) max(Data$m.cx) else(ranges$x[2])
          Min_y <- if(is.null(ranges$y)) min(Data$m.cy) else(ranges$y[1])
          Max_y <- if(is.null(ranges$y)) max(Data$m.cy) else(ranges$y[2])

          HEATMAP_MATRIX <-
            Data %>%
            dplyr::filter(m.cx >= Min_x,
                          m.cx <= Max_x,
                          m.cy >= Min_y,
                          m.cy <= Max_y) %>%
            dplyr::select(-c(1:2)) %>% as.matrix()

          #Prepare color function
          col_fun <- circlize::colorRamp2(c(Min_HEATMAP, 0, Max_HEATMAP), c("#0000ff", "white", "#ff0000"))

          #Plot our heatmap
          ComplexHeatmap::Heatmap(HEATMAP_MATRIX,
                                  col = col_fun,
                                  cluster_rows = TRUE,
                                  cluster_columns = FALSE,
                                  show_row_names = FALSE,
                                  show_heatmap_legend = FALSE,
                                  show_row_dend = FALSE,
                                  show_column_dend = FALSE,
                                  column_names_side = "top",
                                  column_names_rot = 0,
                                  column_names_centered = TRUE,
                                  column_names_gp = grid::gpar(fontsize = 10),
                                  border = TRUE
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
