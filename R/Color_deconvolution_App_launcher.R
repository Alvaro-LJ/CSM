#' Launches a shiny APP to explore color deconvolution parameters
#'
#' `Color_deconvolution_App_launcher()` launches an APP to interactively explore color deconvolution parameters
#' Parameters can then be used to feed the [Image_deconvolution_function()].
#'
#' @param Directory Character string specifying the path to the folder where RGB color images are present.
#'
#' @details
#' Control panel controls the image display settings. Channel name box allows the user to specify the name of the channel to be extracted.
#' The "+Add" button adds the channel extraction parameters to the parameters.
#' The "-Remove" button removes any channel parameters that match the "Channel name" box.
#'
#' Initial Pre-processing allows the initial modification of the image to improve results. The colors option allows the user to limit the maximum amount of different colors used in the image.
#'
#' Color targeting allows different options:
#'\itemize{
#' \item{Color: The user can select the target color by it's HEX code.}
#' \item{RED GREEN BLUE: The user can select the target color by providing RGB values.}
#' \item{Upper right panel: The user can select the target color by clicking on the desired target pixel color value.}
#' }
#'
#' Tolerance values are key to extract the color adequately:
#'\itemize{
#' \item{Individual tolerance: Allows the user to specify the flexibility of individual color channel extraction.}
#' \item{End tolerance: After the individual RGB channels have been extracted and merged, individual tolerance regulates flexibility of pixel values with respect to the color target.}
#' }
#'
#' Final processing stablishes how the image is going to be modified after color channel extraction.
#'
#' "Download Parameters" button saves parameters in the Global Environment.
#'
#' Image Panels:
#'\itemize{
#' \item{Upper left: Image being color processed. Initial Pre-processing steps will be applied. Use it to zoom in.}
#' \item{Upper right: Pixel selector. Initial Pre-processing steps will be applied.}
#' \item{Lower left: Initial color extracted according to active parameters.}
#' \item{Lower right: Final channel result after Final processing steps.}
#' }
#'
#' @examples
#' \dontrun{
#' #Create temporary input and output directories------------------------------
#' Input_Dir <- tempfile(pattern = "tempdir1_Input")
#' dir.create(Input_Dir, recursive = TRUE)
#'
#' #Save images in Input directory
#' purrr::map(1:2,
#' function(Image){
#'    EBImage::writeImage(CSM_MiniHE_test[[Image]], file.path(Input_Dir, names(CSM_MiniHE_test)[Image]))
#' })
#'
#'#Launch the app------------------------------------------------------------
#'Color_deconvolution_App_launcher(Directory = Input_Dir)
#'
#'#Remove directories---------------------------------------------------------
#'unlink(Input_Dir, recursive = TRUE)
#'
#' }
#'
#'
#' @export

Color_deconvolution_App_launcher <-
  function(Directory = NULL){
    #Check suggested packages
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
               expression(install.packages("benchmarkme")))
      )
    }

    #Run a gc after exiting because it consumes a load of RAM
    on.exit(gc())

    #Obtain image names and channel names
    Real_Images <- dir(Directory, full.names = FALSE)

    #BUILD THE USER INTERFACE
    {
      user_interface <- shiny::fluidPage(

        #Set the title
        shiny::titlePanel("Image color deconvolution APP"),

        #We want a two panel layout, one in the left containing the input parameters and the output in the right
        shiny::sidebarLayout(
          #Set the first column (which contains the user defined parameters)
          shiny::sidebarPanel(
            #ID and width
            id="sidebar",

            #Select the image to be tested
            shiny::fluidRow(
              shiny::column(3, shiny::selectInput("Real_Image_name", "Image to display", sort(Real_Images), multiple = FALSE)),
              shiny::column(3, shiny::textInput("Channel_name", "Channel name", value = "Channel_1")),
              shiny::column(2, shiny::selectInput("Res", "Image Res", c(Low = 150, Mid = 300, high = 600, Original = 1000), selected = 150, multiple = FALSE)),
              shiny::column(2, shiny::actionButton("ADD_button", shiny::icon("plus", library = "font-awesome"), label = "Add")),
              shiny::column(2, shiny::actionButton("Remove", shiny::icon("minus", library = "font-awesome"), label = "Remove"))
            ),

            shiny::h4("Initial Pre-processing"),
            #Select basic image features
            shiny::fluidRow(
              shiny::column(4, shiny::sliderInput("Brightness", "Brightness", value = 100, min = 0, max = 100, step = 1)),
              shiny::column(4, shiny::sliderInput("Saturation", "Saturation", value = 100, min = 0, max = 100, step = 1)),
              shiny::column(4, shiny::sliderInput("Hue", "Hue", value = 100, min = 0, max = 100, step = 1))
            ),

            #Normalize, Equalize, Contrast and color reduction
            shiny::fluidRow(
              shiny::column(2, shinyWidgets::materialSwitch("Normalize", "Normalize", value = FALSE)),
              shiny::column(2, shinyWidgets::materialSwitch("Equalize", "Equalize", value = FALSE)),
              shiny::column(2, shinyWidgets::materialSwitch("Contrast", "Contrast", value = FALSE)),
              shiny::column(2, shiny::conditionalPanel(condition = "input.Contrast == '1'",
                                                       shiny::numericInput("Sharpen", "Sharpen", value = 1, min = 1, max = NA, step = 1)
              )),
              shiny::column(2, shinyWidgets::materialSwitch("Quantize", "Colors", value = FALSE)),
              shiny::column(2, shiny::conditionalPanel(condition = "input.Quantize == '1'",
                                                       shiny::numericInput("Max_colors", "Max", value = 100, min = 1, max = NA, step = 1)
              ))

            ),

            #Target color
            shiny::h4("Color targeting"),
            shiny::fluidRow(
              shiny::column(3, shiny::textInput("Color_text", "Color", value = "#ffffff")),
              shiny::column(2, shiny::numericInput("RED", "RED", value = 255, min = 0, max = 255, step = 1)),
              shiny::column(2, shiny::numericInput("GREEN", "GREEN", value = 255, min = 0, max = 255, step = 1)),
              shiny::column(2, shiny::numericInput("BLUE", "BLUE", value = 255, min = 0, max = 255, step = 1)),
              shiny::column(1, shiny::plotOutput("Target", width = 45, height = 45, inline = FALSE))
            ),
            shiny::fluidRow(
              shiny::column(3, shiny::numericInput("Tolerance_RED", "Red Tol", value = 0.1, min = 0, max = 1, step = 0.01)),
              shiny::column(3, shiny::numericInput("Tolerance_GREEN", "Green Tol", value = 0.1, min = 0, max = 1, step = 0.01)),
              shiny::column(3, shiny::numericInput("Tolerance_BLUE", "Blue Tol", value = 0.1, min = 0, max = 1, step = 0.01)),
              shiny::column(3, shiny::numericInput("Tolerance_Final", "End Tol", value = 0.1, min = 0, max = 1, step = 0.01))
            ),

            #Opening and closing
            shiny::h4("Final Processing"),
            #Normalize and equalize
            shiny::fluidRow(
              shiny::column(2, shinyWidgets::materialSwitch("Post_Normalize", "Normalize", value = FALSE)),
              shiny::column(2, shinyWidgets::materialSwitch("Post_Equalize", "Equalize", value = FALSE)),
            ),

            #Erode dilate
            shiny::fluidRow(
              shiny::column(2, shinyWidgets::materialSwitch("Erode", "Erode", value = FALSE)),
              shiny::column(2, shiny::conditionalPanel(condition = "input.Erode == '1'",
                                                       shiny::numericInput("Erode_kern", "Kern", value = 1, min = 1, max = NA, step = 1)
              )),
              shiny::column(2, shiny::conditionalPanel(condition = "input.Erode == '1'",
                                                       shiny::numericInput("Erode_round", "Rounds", value = 1, min = 1, max = NA, step = 1)
              )),
              shiny::column(2, shinyWidgets::materialSwitch("Dilate", "Dilate", value = FALSE)),
              shiny::column(2, shiny::conditionalPanel(condition = "input.Dilate == '1'",
                                                       shiny::numericInput("Dilate_kern", "Kern", value = 1, min = 1, max = NA, step = 1)
              )),
              shiny::column(2, shiny::conditionalPanel(condition = "input.Dilate == '1'",
                                                       shiny::numericInput("Dilate_round", "Rounds", value = 1, min = 1, max = NA, step = 1)
              ))
            ),

            shiny::h4("Parameter results"),
            #Buttons to add or remove
            shiny::fluidRow(
              shiny::column(3, shiny::actionButton("Download_Param", shiny::icon("download"), label = "Download Parameters"))
            ),

            #The current results
            shiny::fluidRow(
              shiny::column(12, shiny::verbatimTextOutput("Result_list"))
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
              shiny::column(5, ggiraph::girafeOutput("Pre_Processed"))
            ),
            #Second row will contain the positive cells and the histogram
            shiny::fluidRow(
              shiny::column(5, shiny::plotOutput("Color")),
              shiny::column(5, shiny::plotOutput("Final_channel"))
            )
          )
        ),
        shiny::tags$head(shiny::tags$style(
          htmltools::HTML('
         #sidebar {
            background-color: #cee02f;
        }

        body, label, input, button, select {
          font-family: "Arial";
        }')))
      )
    }


    #BUILD THE SERVER
    server <- function(input, output, session){
      #BASIC reactives
      Image_dir <- shiny::reactive(stringr::str_c(Directory, "/", input$Real_Image_name))
      Channel_name <- shiny::reactive(input$Channel_name)
      Resolution <- shiny::reactive(input$Res)

      Brightness <- shiny::reactive(input$Brightness)
      Saturation <- shiny::reactive(input$Saturation)
      Hue <- shiny::reactive(input$Hue)
      Normalize <- shiny::reactive(input$Normalize)
      Equalize <- shiny::reactive(input$Equalize)
      Contrast <- shiny::reactive(input$Contrast)
      Sharpen <- shiny::reactive(input$Sharpen)
      Reduce_color <- shiny::reactive(input$Quantize)
      Max_color <- shiny::reactive(input$Max_colors)

      Color_name <- shiny::reactive(input$Color_text)
      RED_value <- shiny::reactive(input$RED)
      GREEN_value <- shiny::reactive(input$GREEN)
      BLUE_value <- shiny::reactive(input$BLUE)
      Tolerance_value_RED <- shiny::reactive(input$Tolerance_RED)
      Tolerance_value_GREEN <- shiny::reactive(input$Tolerance_GREEN)
      Tolerance_value_BLUE <- shiny::reactive(input$Tolerance_BLUE)
      Tolerance_value_FINAL <- shiny::reactive(input$Tolerance_Final)

      Normalize_post <- shiny::reactive(input$Post_Normalize)
      Equalize_post <- shiny::reactive(input$Post_Equalize)
      Erode <- shiny::reactive(input$Erode)
      Erode_kern <- shiny::reactive(input$Erode_kern)
      Erode_rounds <- shiny::reactive(input$Erode_round)
      Dilate <- shiny::reactive(input$Dilate)
      Dilate_kern <- shiny::reactive(input$Dilate_kern)
      Dilate_rounds <- shiny::reactive(input$Dilate_round)

      #SPECIAL reactives
      #4 chained reactives to process the photo
      #Original photo
      Photo_original <- shiny::reactive({
        Photo <- magick::image_read(Image_dir())
        return(Photo)
      })
      Pre_processed_photo <- shiny::reactive({
        Photo <- Photo_original()

        #Apply required processing
        Photo <- Photo %>% magick::image_modulate(brightness = Brightness(), saturation = Saturation(), hue = Hue())
        if(as.logical(Equalize())) Photo <- Photo %>% magick::image_equalize()
        if(as.logical(Normalize())) Photo <- Photo %>% magick::image_normalize()
        if(as.logical(Contrast())) Photo <- Photo %>% magick::image_contrast(sharpen = Sharpen())
        if(as.logical(Reduce_color())) Photo <- Photo %>% magick::image_quantize(max = Max_color(), dither = TRUE)

        return(Photo)

      })
      Color_extracted_photo <- shiny::reactive({
        Photo <- Pre_processed_photo()

        #Transform to EBI object
        Photo <- Photo %>% magick::image_convert(type = NULL,
                                                 colorspace = "RGB",
                                                 depth = 8,
                                                 antialias = NULL,
                                                 matte = FALSE) %>%
          magick::as_EBImage() %>%
          EBImage::toRGB()
        #Extract Image
        Photo <- EBImage_color_thresholder(Image = Photo,
                                           Target = c(RED_value(), GREEN_value(), BLUE_value()),
                                           Color_Tolerance = c(Tolerance_value_RED(), Tolerance_value_GREEN(), Tolerance_value_BLUE()),
                                           Final_Tolerance = Tolerance_value_FINAL())
        #Return to magick object
        Photo <- magick::image_read(Photo)
        return(Photo)
      })
      Color_processed_photo <- shiny::reactive({
        Photo <- Color_extracted_photo()
        #Apply transformations before erosion and dilation
        if(as.logical(Normalize_post())) Photo <- Photo %>% magick::image_normalize()
        if(as.logical(Equalize_post())) Photo <- Photo %>% magick::image_equalize()

        #Transform to EBI object
        Photo <- Photo %>% magick::as_EBImage()

        #if erosion is required then execute the rounds of erosion the user requires
        if(as.logical(Erode())){
          Photo <- purrr::reduce(
            .x = seq_along(1:Erode_rounds()),
            .f = function(img, ...) EBImage::erode(img, kern = EBImage::makeBrush(Erode_kern(), "disc")),
            .init = Photo
          )
        }

        #If dilation is required then execute the dilation rounds
        if(as.logical(Dilate())){
          Photo <- purrr::reduce(
            .x = seq_along(1:Dilate_rounds()),
            .f = function(img, ...) EBImage::dilate(img, kern = EBImage::makeBrush(Dilate_kern(), "disc")),
            .init = Photo
          )
        }

        #Return to magick object
        Photo <- magick::image_read(Photo)
        return(Photo)
      })

      #Zoom reactive
      ranges <- shiny::reactiveValues(x = NULL, y = NULL)


      #Bi-directional reactives (Color name and Channel values)
      #First add RGB values if color is specified
      shiny::observeEvent({
        input$Color_text #Observe if text is modified by writing
      },
      {
        shiny::req(input$Color_text)

        #If color is introduced then get the color and rgb vector
        if(!berryFunctions::is.error(grDevices::col2rgb(as.character(Color_name()), alpha = FALSE)[,1])){
          RGB_vector <- grDevices::col2rgb(as.character(Color_name()), alpha = FALSE)[,1]

          shiny::updateNumericInput(session, "RED", value = RGB_vector[["red"]])
          shiny::updateNumericInput(session, "GREEN", value = RGB_vector[["green"]])
          shiny::updateNumericInput(session, "BLUE", value = RGB_vector[["blue"]])
        }
      },
      ignoreInit = TRUE)
      #in addition specify what to do if the user clicks on the pre-processed image
      shiny::observeEvent({
        input$Pre_Processed_selected #Observe if text is modified by writing
      },
      {
        shiny::req(input$Pre_Processed_selected)
        Color_name <- stringr::str_split_i(input$Pre_Processed_selected, "_", i = 1)

        #If color is introduced then get the color and rgb vector
        if(!berryFunctions::is.error(grDevices::col2rgb(Color_name, alpha = FALSE)[,1])){
          RGB_vector <- grDevices::col2rgb(Color_name, alpha = FALSE)[,1]

          shiny::updateNumericInput(session, "RED", value = RGB_vector[["red"]])
          shiny::updateNumericInput(session, "GREEN", value = RGB_vector[["green"]])
          shiny::updateNumericInput(session, "BLUE", value = RGB_vector[["blue"]])
        }
      },
      ignoreInit = TRUE)
      #Second add color name if RGB colors are modified
      shiny::observeEvent(
        { #Observe any modification in RGB values
          input$RED
          input$GREEN
          input$BLUE
        },
        {

          #If color is introduced then get the color and rgb vector
          RED_value <- if(RED_value()/255 > 1) 1 else RED_value()/255
          GREEN_value <- if(GREEN_value()/255 > 1) 1 else GREEN_value()/255
          BLUE_value <- if(BLUE_value()/255 > 1) 1 else BLUE_value()/255

          Color_name <- grDevices::rgb(RED_value, GREEN_value, BLUE_value)
          session$sendCustomMessage(type = 'Pre_Processed_set', message = as.character(Color_name[1]))
          shiny::updateTextInput(session, "Color_text", value = Color_name[1])

        },
        ignoreInit = TRUE)
      #Color target mini-rectangle
      output$Target <- shiny::renderPlot({
        ggplot() +
          geom_rect(
            aes(xmin = 0, xmax = 1, ymin = 0, ymax = 1),
            fill = Color_name(),
            color = "black"
          ) +
          scale_x_continuous("", labels = NULL) +
          scale_y_continuous("", labels = NULL) +
          theme_minimal() +
          theme(axis.ticks = element_blank(),
                axis.line = element_blank(),
                plot.margin = margin(-100, -100, -100, -100, "pt"))
      })


      #PARAMETER LIST
      #generate an object that will be modified by the user when the Add button is clicked
      Parameters_object <- shiny::reactiveValues(Parameters_list = list()) #Parameter list will be a list of lists
      #If the user hits the add button, place all parameters into a list and add a new list item to the Parameters list ()
      shiny::observeEvent(input$ADD_button,
                          {
                            #If name is present print a warning message
                            if(input$Channel_name %in% names(Parameters_object$Parameters_list)){
                              showModal(modalDialog(
                                paste0(input$Channel_name, " already present in Parameter list. ", input$Channel_name, " will be overwritten" ),
                                easyClose = TRUE,
                                footer = NULL
                              )
                              )
                            }

                            #Generate parameter list
                            Deconv_parameters <-
                              list(Brightness = input$Brightness,
                                   Saturation = input$Saturation,
                                   Hue = input$Hue,
                                   Normalize = as.logical(input$Normalize),
                                   Equalize = as.logical(input$Equalize),
                                   Contrast = as.logical(input$Contrast),
                                   Sharpen = input$Sharpen,
                                   Reduce_color = as.logical(input$Quantize),
                                   Max_color = input$Max_colors,

                                   RED_value = input$RED,
                                   GREEN_value = input$GREEN,
                                   BLUE_value = input$BLUE,
                                   Color_Tolerance = c(input$Tolerance_RED, input$Tolerance_GREEN, input$Tolerance_BLUE),
                                   Final_Tolerance = input$Tolerance_Final,

                                   Post_normalize = as.logical(input$Post_Normalize),
                                   Post_equalize = as.logical(input$Post_Equalize),
                                   Erode = as.logical(input$Erode),
                                   Erode_kern = input$Erode_kern,
                                   Erode_rounds = input$Erode_round,
                                   Dilate = as.logical(input$Dilate),
                                   Dilate_kern = input$Dilate_kern,
                                   Dilate_rounds = input$Dilate_round
                              )

                            #If name is present replace the exact element of the list
                            if(input$Channel_name %in% names(Parameters_object$Parameters_list)){
                              Name_index <- match(input$Channel_name, names(Parameters_object$Parameters_list))
                              Parameters_object$Parameters_list[[Name_index]] <- Deconv_parameters
                            }

                            else{
                              Parameters_object$Parameters_list[[length(Parameters_object$Parameters_list)+1]] <- Deconv_parameters
                              names(Parameters_object$Parameters_list)[length(Parameters_object$Parameters_list)] <- as.character(input$Channel_name)
                            }
                          })

      #If the user hits the remove button then look for that element in the Parameter_list and remove it
      shiny::observeEvent(input$Remove,
                          {
                            Name_index <- match(input$Channel_name, names(Parameters_object$Parameters_list))
                            if(is.na(Name_index)){
                              showModal(modalDialog(
                                paste0(input$Channel_name, " not present in Parameter list. Please check channel name provided"),
                                easyClose = TRUE,
                                footer = NULL
                              )
                              )
                            }
                            else{
                              Parameters_object$Parameters_list <- Parameters_object$Parameters_list[-Name_index]
                            }
                          })
      #Generate the output
      output$Result_list <- shiny::renderPrint({
        if(length(Parameters_object$Parameters_list) == 0) "Awaiting initial parameters"
        else{
          Parameters_object$Parameters_list
        }
      })
      #If download parameter button is clicked, download paramater list
      shiny::observeEvent(input$Download_Param,
                          {
                            if(length(Parameters_object$Parameters_list) == 0){
                              showModal(modalDialog(
                                "Channel deconvolution parameters not specified. Please add parameters to the list and retry",
                                easyClose = TRUE,
                                footer = NULL
                              )
                              )
                            }
                            else{
                              Deconvolution_Parameters <<- Parameters_object$Parameters_list
                              showModal(modalDialog(
                                "An object called 'Deconvolution_Parameters' has been created in the Global environment.",
                                easyClose = TRUE,
                                footer = NULL
                              )
                              )
                            }
                          })



      #Original photo (top left)
      #Print the photo
      output$Photo <- shiny::renderPlot({
        #Obtain the photo
        Photo <- Photo_original()

        #Modify resolution before printing
        if(as.numeric(Resolution()) != 1000){
          Image_Resolution <- stringr::str_c("X", Resolution())
          Photo <- magick::image_scale(Photo, Image_Resolution)
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

      #Pre-processed photo (top right)
      #Generate the plot
      Pre_processed_plot <-
        shiny::reactive({
          Photo <- Pre_processed_photo() #Get processed photo

          #If no zoom in required then low down the resolution of the image(always as geom_tile is very computationally slow)
          if(is.null(ranges$x)){
            Image_Resolution <- stringr::str_c("X", Resolution())
            Photo <- magick::image_scale(Photo, "X100")
            Photo <- Photo %>% magick::image_raster()
            #Reverse pixels to match annotate_raster system
            Photo$y <- rev(Photo$y)
          }

          #If zoom is required then keep original user defined resolution (required to match zoomed in areas)
          if(!is.null(ranges$x)){
            if(as.numeric(Resolution()) == 1000){
              Photo <- Photo %>% magick::image_raster()
              #Reverse pixels to match annotate_raster system
              Photo$y <- rev(Photo$y)
            }
            else{
              Image_Resolution <- stringr::str_c("X", Resolution())
              Photo <- magick::image_scale(Photo, Image_Resolution)
              Photo <- Photo %>% magick::image_raster()
              #Reverse pixels to match annotate_raster system
              Photo$y <- rev(Photo$y)
            }

            #remove pixels outside the scale
            Photo <- Photo %>% dplyr::filter(x >= ranges$x[[1]],
                                             x <= ranges$x[[2]],
                                             y >= ranges$y[[1]],
                                             y <= ranges$y[[2]])

          }

          #Generate an ID that contains the color code (used for tooltip data)
          Photo$Pixel_id <- stringr::str_c(Photo$col, 1:nrow(Photo), sep = "_")


          return(
            Photo %>%
              ggplot() +
              ggiraph::geom_tile_interactive(aes(x = x, y = y, fill = col,
                                                 data_id = Pixel_id,
                                                 tooltip = as.character(col),
                                                 hover_nearest = FALSE
              ),
              linewidth = 0) +
              scale_fill_identity() +
              cowplot::theme_cowplot()+
              guides(color = "none") +
              theme(axis.line = element_blank(),
                    axis.ticks = element_blank(),
                    axis.text = element_blank(),
                    axis.title = element_blank(),
                    panel.background = element_rect(fill = "black")) +
              coord_cartesian(xlim = ranges$x, ylim = ranges$y, expand = FALSE)
          )
        })
      #send the graph to the UI
      output$Pre_Processed <- ggiraph::renderGirafe({
        plot <- ggiraph::girafe(code = print(Pre_processed_plot()),
                                options = list(
                                  ggiraph::opts_hover(css = "stroke:black;cursor:pointer;", reactive = TRUE),
                                  ggiraph::opts_selection(type = "single", css = "fill:#FF3333;stroke:black;")
                                )
        )
        return(plot)
      })


      #Color extraction (bottom left)
      output$Color <- shiny::renderPlot({
        #Obtain the photo
        Photo <- Color_extracted_photo()

        #Modify resolution before printing
        if(as.numeric(Resolution()) != 1000){
          Image_Resolution <- stringr::str_c("X", Resolution())
          Photo <- magick::image_scale(Photo, Image_Resolution)
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

      #Final color channel (bottom right)
      output$Final_channel <- shiny::renderPlot({
        #Obtain the photo
        Photo <- Color_processed_photo()

        #Modify resolution before printing
        if(as.numeric(Resolution()) != 1000){
          Image_Resolution <- stringr::str_c("X", Resolution())
          Photo <- magick::image_scale(Photo, Image_Resolution)
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

      #If browser is closed end the app
      session$onSessionEnded(function() { shiny::stopApp() })
    }

    #Run the server
    message("Always stop current R execution if you want to continue with your R session")
    shiny::shinyApp(user_interface, server)
  }
