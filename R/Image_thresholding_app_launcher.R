#' Launches a shiny APP to explore image pixel thresholding parameters
#'
#' `Image_thresholding_app_launcher()` launches an APP to interactively explore image pixel segmentation parameters interactively.
#' @param Directory Character specifying the path to the folder where images to be thresholded are present.
#' @param Ordered_Channels Character vector specifying image channels in their exact order.
#'
#' @seealso [Pixel_Threshold_calculator()], [Binary_threshold_image_combinator()], [MFI_Experimet_Calculator()]
#'
#' @details
#' Image parameters control the image and channel to be displayed.
#'
#'Tissue mask parameters allow customizing how the tissue mask will be computed (decide foreground / background pixels).
#'
#'Tissue thresholding parameters control the tresholding process. If arbitrary, the "value" box is used to write the desired threshold value.
#'If Multilevel, the "value" box can be used to write a vector with the thresholds (for example "c(0.1, 0.5, 0.9)") or the number of levels
#'can be specified with the Levels box. If "value" is not provided, the app will calculate threshold values using the imagerExtra::ThresholdML function.
#'
#'The GO! button executes the thresholding according to active parameters.
#'
#'Image Panels:
#'\itemize{
#' \item{Upper left: Image being thresholded. Use it to zoom in.}
#' \item{Upper right: Tissue mask. Shows the tissue mask generated.}
#' \item{Lower left: Pixels above the threshold included in the tissue mask.}
#' \item{Lower right: Final result after applying blurring to the target (if required).}
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
#'    EBImage::writeImage(CSM_MiniMultiTiff_test[[Image]], file.path(Input_Dir, names(CSM_MiniMultiTiff_test)[Image]))
#' })
#'
#' #Deploy app----------------------------------------------
#' Image_thresholding_app_launcher(
#'     Directory = Input_Dir,
#'     Ordered_Channels = c("DAPI", "PDL1", "GZMB", "PD1", "CK-EPCAM", "CD8a", "FOXP3")
#' )
#'
#' #Remove directories---------------------------------------------------------
#' unlink(Input_Dir, recursive = TRUE)
#' }
#'
#' @export

Image_thresholding_app_launcher <-
  function(Directory = NULL,
           Ordered_Channels = NULL){

    #Check suggested packages
    {
      if(!requireNamespace("magick", quietly = FALSE)) stop(
        paste0("magick CRAN package is required to execute the function. Please install using the following code: ",
               expression(install.packages("magick")))
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


    #check that the directory provided contains at least one file
    if(length(dir(Directory)) <1 ) stop("No files found in the Directory provided")

    #Obtain image names and channel names
    Real_Images <- dir(Directory, full.names = FALSE)
    Complete_names <- dir(Directory, full.names = TRUE)
    Channels_in_images <- Ordered_Channels

    #BUILD THE USER INTERFACE
    {
      user_interface <- shiny::fluidPage(

        #Set the title
        shiny::titlePanel("Pixel thresholding exploration APP"),

        #We want a two panel layout, one in the left containing the input parameters and the output in the right
        shiny::sidebarLayout(
          #Set the first column (which contains the user defined parameters)
          shiny::sidebarPanel(
            #ID and width
            id="sidebar",

            #Image to display and resolution
            shiny::h4("Image parameters"),
            shiny::fluidRow(
              #Select the real image to be displayed
              shiny::column(6, shiny::selectInput("Real_Image_name", "Image to display", sort(Real_Images), multiple = FALSE)),
              #Select the channel to be displayed
              shiny::column(4, shiny::selectInput("Channel", "Channel to display", Channels_in_images, multiple = FALSE)),
              #Select the resolution
              shiny::column(2, shiny::selectInput("Res", "Image Res", c(Low = 150, Mid = 300, high = 600, Original = 1000), selected = 150, multiple = FALSE))
            ),

            #Channels included in tissue threshold
            shiny::h4("Tissue mask parameters"),
            shiny::fluidRow(
              #Select the channels to include in the tissue mask generation
              shiny::column(3, shinyWidgets::virtualSelectInput("Tissue_mask_channels", label = "Channels used",
                                                                choices = Channels_in_images,
                                                                search = TRUE,
                                                                multiple = TRUE
              )),
              #Select the type of threshold being used for tissue mask
              shiny::column(3, shiny::selectInput("Tissue_threshold_type", "Threshold type", c("Otsu", "Arbitrary", "Absolute"), multiple = FALSE)),
              shiny::column(2, shiny::conditionalPanel(condition = "input.Tissue_threshold_type == 'Arbitrary'",shiny::textInput("Tissue_value", "Value", value = "0.01"))),
              #Select tissue mask blurr
              shiny::column(2, shiny::selectInput("Tissue_blurr", "Blurr", c(TRUE,FALSE), selected = FALSE, multiple = FALSE)),
              shiny::column(2, shiny::conditionalPanel(condition = "input.Tissue_blurr == 'TRUE'",shiny::textInput("Tissue_sigma", "Sigma", value = "0.5")))
            ),

            shiny::h4("Tissue thresholding parameters"),
            shiny::fluidRow(
              #Local or global thresholding
              shiny::column(3, shiny::selectInput("Target_threshold_local", "Type", c("Global", "Local"), selected = "Local", multiple = FALSE)),

              #Threshold type
              shiny::column(3, shiny::selectInput("Target_threshold_type", "Type", c("Otsu", "Arbitrary", "Multilevel"), selected = "Otsu", multiple = FALSE)),

              #Select target mask blurr
              shiny::column(3, shiny::selectInput("Target_blurr", "Blurr", c(TRUE,FALSE), selected = FALSE, multiple = FALSE)),
              shiny::column(3, shiny::conditionalPanel(condition = "input.Target_blurr == 'TRUE'",shiny::textInput("Target_sigma", "Sigma", value = "0.5")))
            ),
            #Conditional panels for values
            shiny::fluidRow(
              shiny::column(6, shiny::conditionalPanel(condition = "input.Target_threshold_type == 'Arbitrary' || input.Target_threshold_type == 'Multilevel'",
                                                       shiny::textInput("Target_Value", "Value", value = NULL))),
              shiny::column(6, shiny::conditionalPanel(condition = "input.Target_threshold_type == 'Multilevel'",
                                                       shiny::numericInput("Target_Levels", "Levels", value = 3, min = 3, max = NA, step = 1)))
            ),

            shiny::fluidRow(shiny::column(2, shiny::actionButton("GO_button", shiny::icon("bolt-lightning"), label = "GO!")))

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
              shiny::column(6, shiny::plotOutput("Tissue_mask",
                                                 width = "auto"))
            ),
            #Second row will contain the positive cells and the histogram
            shiny::fluidRow(
              shiny::column(6, shiny::plotOutput("Target_mask",
                                                 width = "auto")),
              shiny::column(6, shiny::verbatimTextOutput("Pixel_summary"))
            )
          )
        ),
        shiny::tags$head(shiny::tags$style(
          htmltools::HTML('
         #sidebar {
            background-color: #fc92fc;
        }

        body, label, input, button, select {
          font-family: "Arial";
        }')))
      )
    }


    #BUILD THE SERVER
    server <- function(input, output, session){
      #First generate the reactives (those that need to be continuously tested for reactivity)
      Photo_name <- shiny::reactive(stringr::str_c(Directory, "/", input$Real_Image_name))
      Channel_index <- shiny::reactive(which(input$Channel == Ordered_Channels))
      Image_resolution <- shiny::reactive(input$Res)
      #Generate a reactive that controls the zoom in
      ranges <- shiny::reactiveValues(x = NULL, y = NULL)
      #Generate the reactive that controls the results in the panel
      Image_results <- shiny::reactiveValues(Tissue_mask = NULL,
                                             Target_mask = NULL,
                                             Pixel_data = NULL)


      #Import the actual photo
      Photo_reactive <- shiny::reactive({
        #Import the Photo
        Photo <- magick::image_read(Photo_name())[Channel_index()]
        return(Photo)
      })

      #Print the photo
      output$Photo <- shiny::renderPlot({
        #Obtain the photo
        Photo <- Photo_reactive()
        #Change resolution if required
        if(as.numeric(Image_resolution()) < 1000){
          Image_Resolution <- stringr::str_c("X", Image_resolution())
          Photo <- magick::image_scale(Photo, geometry = Image_Resolution)
        }
        #Get image info
        Image_information <- magick::image_info(Photo)
        Width <- Image_information$width
        Height <- Image_information$height

        #Generate the plot as annotation_raster
        return(
          ggplot() +
            annotation_raster(Photo, xmin = 0, xmax = Width, ymin = 0, ymax = Height, interpolate = TRUE)+
            scale_x_continuous(limits = c(0, Width)) +
            scale_y_continuous(limits = c(0, Height)) +
            coord_cartesian(xlim = ranges$x, ylim = ranges$y, expand = FALSE)+
            theme(axis.title = element_blank(),
                  axis.text = element_blank(),
                  axis.ticks = element_blank(),
                  axis.line = element_blank(),
                  panel.background = element_rect(fill = "black"),
                  panel.grid = element_blank())
        )
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

      shiny::observeEvent(input$GO_button, {
        #First we need to import the photograph as EBI, select the required channels and transform it to a cytomapper object
        shiny::showModal(modalDialog("Importing image", footer=NULL))
        Image <- magick::image_read(Photo_name())
        Image <- Image[which(Ordered_Channels %in% input$Tissue_mask_channels)]
        Image <- magick::as_EBImage(Image)
        shiny::removeModal()

        #Generate tissue mask
        shiny::showModal(modalDialog("Generating tissue mask", footer=NULL))
        Tissue_mask <- Tissue_mask_generator(Image = Image,
                                             Threshold_type = input$Tissue_threshold_type,
                                             Threshold_value = as.numeric(input$Tissue_value),
                                             Blurr = as.logical(input$Tissue_blurr),
                                             Sigma = as.numeric(input$Tissue_sigma))
        Image_results$Tissue_mask <- Tissue_mask
        shiny::removeModal()

        #Perform tissue thresholding according to user preferences
        shiny::showModal(modalDialog("Thresholding pixels", footer=NULL))
        Image_index <- match(input$Channel, input$Tissue_mask_channels)
        Target_image <- EBImage::getFrame(y = Image, i = Image_index)


        #If Arbitrary (Global or local)
        if(input$Target_threshold_type == "Arbitrary"){
          Target_result <- Pixel_thresholder(Target = Target_image ,
                                             Tissue_mask = Tissue_mask,
                                             Threshold_type = "Arbitrary",
                                             Threshold_value = as.numeric(input$Target_Value),
                                             Blurr = as.logical(input$Target_blurr),
                                             Sigma = as.numeric(input$Target_sigma)
          )
          Image_results$Target_mask <- Target_result
        }

        #If local and Otsu
        if(all(input$Target_threshold_type == "Otsu", input$Target_threshold_local == "Local")){
          Target_result <- Pixel_thresholder(Target = Target_image ,
                                             Tissue_mask = Tissue_mask,
                                             Threshold_type = "Otsu",
                                             Threshold_value = NULL,
                                             Blurr = as.logical(input$Target_blurr),
                                             Sigma = as.numeric(input$Target_sigma)
          )
          Image_results$Target_mask <- Target_result
        }

        #If global and Otsu
        if(all(input$Target_threshold_type == "Otsu", input$Target_threshold_local == "Global")){
          shiny::showModal(modalDialog("Calculating global Otsu images taking into account all images in Directory. This can take a while", footer=NULL))
          #Calculate a Tissue mask list
          Tissue_mask_list <-purrr::map(seq_along(1:length(Complete_names)),
                                        function(Image_index){
                                          Image <- magick::image_read(Complete_names[Image_index])[which(Ordered_Channels %in% input$Tissue_mask_channels)]
                                          Image <- magick::as_EBImage(Image)

                                          #First generate the tissue mask and then remove the original image (no longer required)
                                          Tissue_mask <- Tissue_mask_generator(Image = Image,
                                                                               Threshold_type = input$Tissue_threshold_type,
                                                                               Threshold_value = as.numeric(input$Tissue_value),
                                                                               Blurr = as.logical(input$Tissue_blurr),
                                                                               Sigma = as.numeric(input$Tissue_sigma))
                                          rm(Image)
                                          gc()
                                          return(Tissue_mask)
                                        })
          #Generate a full vector containing all image values with an applied tissue mask
          Composite_Image_list <-purrr::map2(.x = 1:length(Complete_names), .y = Tissue_mask_list, function(.x, .y){
            #Import target image
            Image <- magick::image_read(Complete_names[.x])[which(Ordered_Channels %in% input$Tissue_mask_channels)]
            Target_Image <- Image[which(input$Channel == input$Tissue_mask_channels)]
            Target_Image <- magick::as_EBImage(Target_Image)

            #Turn target image values outside tissue mask to 0
            Target_Image[!.y] <- 0

            #return as vector
            return(as.vector(Target_Image))
          })

          #Generare a unified vector using all images
          Common_vector <- unlist(Composite_Image_list)
          #Apply otsu algorithm
          Threshold_global_otsu <- EBImage::otsu(array(Common_vector, dim = c(1, length(Common_vector))), range = c(min(Common_vector), max(Common_vector)), levels = length(unique(Common_vector)))
          Target_result <- Pixel_thresholder(Target = Target_image ,
                                             Tissue_mask = Tissue_mask,
                                             Threshold_type = "Arbitrary",
                                             Threshold_value = Threshold_global_otsu,
                                             Blurr = as.logical(input$Target_blurr),
                                             Sigma = as.numeric(input$Target_sigma)
          )
          Image_results$Target_mask <- Target_result
        }

        #If Multilevel input$Target_Value != "" (Global or local)
        if(all(input$Target_threshold_type == "Multilevel", input$Target_Value != "")){
          #Generate the adequate parsing of expression contained in value
          values_for_function <- eval(parse(text = as.character(input$Target_Value)))

          Target_result <- Pixel_Multilevel_thresholder(Target = Target_image,
                                                        Tissue_mask = Tissue_mask,
                                                        Threshold_values = values_for_function,
                                                        Blurr = as.logical(input$Target_blurr),
                                                        Sigma = as.numeric(input$Target_sigma)
          )
          Target_result$Image <- Target_result$Image/length(values_for_function)
          Image_results$Target_mask <- Target_result
        }

        #If Multilevel Local and input$Target_Value == ""
        if(all(input$Target_threshold_type == "Multilevel", input$Target_Value == "", input$Target_threshold_local == "Local")){
          Threshold_levels <- imagerExtra::ThresholdML(imager::cimg(array(as.vector(Target_image), dim = c(1, length(as.vector(Target_image)), 1, 1))),
                                                       k = (as.numeric(input$Target_Levels)-1),
                                                       returnvalue = TRUE)

          Target_result <- Pixel_Multilevel_thresholder(Target = Target_image,
                                                        Tissue_mask = Tissue_mask,
                                                        Threshold_values = Threshold_levels,
                                                        Blurr = as.logical(input$Target_blurr),
                                                        Sigma = as.numeric(input$Target_sigma)
          )
          Target_result$Image <- Target_result$Image/length(Threshold_levels)
          Image_results$Target_mask <- Target_result
        }

        #If Multilevel Global and input$Target_Value == ""
        if(all(input$Target_threshold_type == "Multilevel", input$Target_Value == "", input$Target_threshold_local == "Global")){
          shiny::showModal(modalDialog("Calculating global Multilevels images taking into account all images in Directory. This can take a while", footer=NULL))
          #Calculate a Tissue mask list
          Tissue_mask_list <-purrr::map(seq_along(1:length(Complete_names)),
                                        function(Image_index){
                                          Image <- magick::image_read(Complete_names[Image_index])[which(Ordered_Channels %in% input$Tissue_mask_channels)]
                                          Image <- magick::as_EBImage(Image)

                                          #First generate the tissue mask and then remove the original image (no longer required)
                                          Tissue_mask <- Tissue_mask_generator(Image = Image,
                                                                               Threshold_type = input$Tissue_threshold_type,
                                                                               Threshold_value = as.numeric(input$Tissue_value),
                                                                               Blurr = as.logical(input$Tissue_blurr),
                                                                               Sigma = as.numeric(input$Tissue_sigma))
                                          rm(Image)
                                          gc()
                                          return(Tissue_mask)
                                        })
          #Generate a full vector containing all image values with an applied tissue mask
          Composite_Image_list <-purrr::map2(.x = 1:length(Complete_names), .y = Tissue_mask_list, function(.x, .y){
            #Import target image
            Image <- magick::image_read(Complete_names[.x])[which(Ordered_Channels %in% input$Tissue_mask_channels)]
            Target_Image <- Image[which(input$Channel == input$Tissue_mask_channels)]
            Target_Image <- magick::as_EBImage(Target_Image)

            #Turn target image values outside tissue mask to 0
            Target_Image[!.y] <- 0

            #return as vector
            return(as.vector(Target_Image))
          })

          #Generare a unified vector using all images
          Common_vector <- unlist(Composite_Image_list)

          Threshold_levels <- imagerExtra::ThresholdML(imager::cimg(array(as.vector(Common_vector), dim = c(1, length(as.vector(Common_vector)), 1, 1))),
                                                       k = (as.numeric(input$Target_Levels)-1),
                                                       returnvalue = TRUE)
          Target_result <- Pixel_Multilevel_thresholder(Target = Target_image,
                                                        Tissue_mask = Tissue_mask,
                                                        Threshold_values = Threshold_levels,
                                                        Blurr = as.logical(input$Target_blurr),
                                                        Sigma = as.numeric(input$Target_sigma)
          )
          Target_result$Image <- Target_result$Image/length(Threshold_levels)
          Image_results$Target_mask <- Target_result
        }

        #Add final parameters
        Image_results$Target_mask$mask_channels <- stringr::str_c(input$Tissue_mask_channels, collapse = ", ")
        Image_results$Target_mask$Tissue_mask <- c("Threshold_type" = input$Tissue_threshold_type,
                                                   "Value" = as.character(input$Tissue_value),
                                                   "Blurr" = as.character(input$Tissue_blurr),
                                                   "Sigma" = as.character(input$Tissue_sigma))

        Image_results$Target_mask$Target <- c("Threshold_type" = input$Target_threshold_type,
                                              "Value" = as.character(input$Target_Value),
                                              "Blurr" = as.character(input$Target_blurr),
                                              "Sigma" = as.character(input$Target_sigma),
                                              "Type" = as.character(input$Target_threshold_local))
        Image_results$Target_mask$Image_analized <- input$Real_Image_name

        shiny::removeModal()
      })

      #Tissue mask output
      output$Tissue_mask <- shiny::renderPlot({
        #If no plot has been generated the print a plot
        if(is.null(Image_results[["Tissue_mask"]])) ggplot()
        else{
          Photo <- Image_results[["Tissue_mask"]]
          Photo[Photo] <- 1
          Photo <- magick::image_read(Photo)
          #Change resolution if required
          if(as.numeric(Image_resolution()) < 1000){
            Image_Resolution <- stringr::str_c("X", Image_resolution())
            Photo <- magick::image_scale(Photo, geometry = Image_Resolution)
          }
          #Get image info
          Image_information <- magick::image_info(Photo)
          Width <- Image_information$width
          Height <- Image_information$height

          #Generate the plot as annotation_raster
          return(
            ggplot() +
              annotation_raster(Photo, xmin = 0, xmax = Width, ymin = 0, ymax = Height, interpolate = TRUE)+
              scale_x_continuous(limits = c(0, Width)) +
              scale_y_continuous(limits = c(0, Height)) +
              coord_cartesian(xlim = ranges$x, ylim = ranges$y, expand = FALSE)+
              theme(axis.title = element_blank(),
                    axis.text = element_blank(),
                    axis.ticks = element_blank(),
                    axis.line = element_blank(),
                    panel.background = element_rect(fill = "black"),
                    panel.grid = element_blank())
          )

        }
      })

      #Target output
      output$Target_mask <- shiny::renderPlot({
        #If no plot has been generated the print a plot
        if(is.null(Image_results[["Target_mask"]])) ggplot()
        else{
          Photo <- Image_results[["Target_mask"]][["Image"]]
          Photo <- Photo*1
          Photo <- magick::image_read(Photo)
          #Change resolution if required
          if(as.numeric(Image_resolution()) < 1000){
            Image_Resolution <- stringr::str_c("X", Image_resolution())
            Photo <- magick::image_scale(Photo, geometry = Image_Resolution)
          }
          #Get image info
          Image_information <- magick::image_info(Photo)
          Width <- Image_information$width
          Height <- Image_information$height

          #Generate the plot as annotation_raster
          return(
            ggplot() +
              annotation_raster(Photo, xmin = 0, xmax = Width, ymin = 0, ymax = Height, interpolate = TRUE)+
              scale_x_continuous(limits = c(0, Width)) +
              scale_y_continuous(limits = c(0, Height)) +
              coord_cartesian(xlim = ranges$x, ylim = ranges$y, expand = FALSE)+
              theme(axis.title = element_blank(),
                    axis.text = element_blank(),
                    axis.ticks = element_blank(),
                    axis.line = element_blank(),
                    panel.background = element_rect(fill = "black"),
                    panel.grid = element_blank())
          )

        }
      })

      #Summary of results
      output$Pixel_summary <- shiny::renderPrint({
        if(is.null(Image_results[["Target_mask"]])) cat("Awaiting initial test")
        else{
          Result_list <- Image_results[["Target_mask"]]
          Result_list <- Result_list[-1]


          Result_list
        }

      }
      )

      #If browser is closed end the app
      session$onSessionEnded(function() { shiny::stopApp() })
    }

    #Run the server
    message("Always stop current R execution if you want to continue with your R session")
    shiny::shinyApp(user_interface, server)
  }
