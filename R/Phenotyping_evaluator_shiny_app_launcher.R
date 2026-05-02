#' Launches a shiny App to explore cell phenotyping results
#'
#' After phenotyping, cell phenotypes can be explored interactively proyected on original images.
#'
#'
#' @param DATA A dataframe or tibble containing cell feature data and a column named 'Phenotype' containing cell phenotype labels.
#' @param Directory Character specifying the path to the folder where images are present.
#' @param Ordered_Channels Character vector specifying image channels in their exact order.
#' @param Max_Gb_cache A single numeric value indicating the memory size of the image cache (see details). 10 Gb by default.
#'
#' @details
#' In order to deal with large images and speed up image toggling, the shiny APP works with a memoised version of the [Smart_image_importer()] function. Loaded images will
#' be stored in a temporary cache until the APP is closed or the cache reaches it's limit.
#'
#' Image parameters control the image and channel to be displayed.
#'
#' User can use the 'Color', 'Add channel', 'Remove channel' and 'Reset' buttons to merge different channels into a single image. The add channel button
#' will color the current channel according to the color provided and will be merged to the current image.
#'
#' Sample summary table is depicted at the bottom of the control panel.
#'
#'Image Panels:
#'\itemize{
#' \item{Upper left: Image being explored. Use it to zoom in.}
#' \item{Upper right: Cell phenotype labels overlayed. Cells of interest can be selected by clicking or by the area selection tools.}
#' \item{Lower left: Cell selection summary.}
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
#'    EBImage::writeImage(CSM_MiniMultiTiff_test[[Image]],
#'    file.path(Input_Dir, names(CSM_MiniMultiTiff_test)[Image]))
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
  function(DATA,
           Directory,
           Ordered_Channels,
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

    ####CHECK ARGUMENTS AND GENERATE MEMOISE FUNCTION####
    DATA <- DATA

    if(!identical(names(DATA)[1:4], c("Cell_no", "X", "Y", "Subject_Names"))) {
      stop("Please generate an appropiate data object using the Data_arrange_function")
    }
    if(!"Phenotype" %in% names(DATA)) stop("DATA must contain a column specifying cell phenotypes")

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

    #####PRE-PROCESSING OF DATA####

    ####Build the look-up table####
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

    })

    #Generate vectors that will guide the sliders and selectors of UI
    Images_in_Data <- unique(DATA$Subject_Names)
    Phenotypes <- unique(DATA$Phenotype)
    Channels_in_images <- Ordered_Channels
    N_colors <- length(unique(DATA$Phenotype))

    #Set the colors for the plots
    if(N_colors > 22) warning("Currently the color only supports displaying 22 phenotype simulataneously. The same color will be assigned to several phenotypes")
    Color_tibble <- tibble(Phenotype = unique(DATA$Phenotype), Color_code = unname(pals::trubetskoy(n = length(unique(DATA$Phenotype))))) %>%
      dplyr::mutate(Color_code = case_when(is.na(Color_code) ~ "white",
                                           TRUE ~ Color_code))
    DATA <- dplyr::left_join(DATA, Color_tibble, by = "Phenotype")

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

    ####BUILD THE USER INTERFACE####
    {
      user_interface <- shiny::fluidPage(

        #To make the sidebar collapsable
        shinyjs::useShinyjs(),

        #Set the title and action button
        shiny::fluidRow(
          shiny::column(width = 3, shiny::h3("Phenotyping exploration APP")),
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

            #Image basic controls
            shiny::fluidRow(
              shiny::column(6, shiny::selectInput("Image_name", "Image", sort(Images_in_Data), multiple = FALSE)),
              shiny::column(2, shiny::selectInput("Res", "Image Res", c(Low = 5, Mid = 5.7, High = 6, Ultra = 6.3), selected = 5.7, multiple = FALSE)),
              shiny::column(2, shinyWidgets::materialSwitch("Change_coords", "Pixel/dist", value = FALSE)),
              shiny::column(2, shiny::conditionalPanel(condition = "input.Change_coords == '1'",
                                                       shiny::textInput("Ratio", "pixel size", value = "1")
              ))),

            #Image rotation and equalize
            shiny::fluidRow(
              shiny::column(3, shiny::sliderInput("Degrees", "Rotate", min = -90, max = 90, value = 0, step = 90, ticks = FALSE)),
              shiny::column(3, shiny::selectInput("X_flip", "Flip X image", c(YES = TRUE, NO = FALSE), selected = FALSE, multiple = FALSE)),
              shiny::column(3, shiny::selectInput("Y_flip", "Flip Y image", c(YES = TRUE, NO = FALSE), selected = FALSE, multiple = FALSE)),
              shiny::column(3, shiny::selectInput("Equalize", "Equalize", c(YES = TRUE, NO = FALSE), selected = FALSE, multiple = FALSE))
            ),

            #Channel settings
            shiny::fluidRow(
              shiny::column(3, shiny::selectInput("Channel", "Channel", Channels_in_images, multiple = FALSE)),
              shiny::column(3, shiny::sliderInput("Min_Image", "Min", value = 0, min = 0, max = 100, step = 1, ticks = FALSE)),
              shiny::column(3, shiny::sliderInput("Max_Image", "Max", value = 100, min = 0, max = 100, step = 1, ticks = FALSE)),
              shiny::column(3, shiny::sliderInput("Gamma", "Gamma", value = 0, min = -3, max = +3, step = 0.01, ticks = FALSE))
            ),

            #Color merging parameters
            shiny::fluidRow(
              shiny::column(3, shiny::textInput("Color_channel", "Color", value = "blue")),
              shiny::column(3, shiny::actionButton("Add_color", "Add channel", shiny::icon("plus", library = "font-awesome"))),
              shiny::column(3, shiny::actionButton("Remove_color", "Remove channel", shiny::icon("minus", library = "font-awesome"))),
              shiny::column(3, shiny::actionButton("Reset_image", "Reset", icon = shiny::icon("redo")))
            ),

            #Add the check box
            shiny::fluidRow(
              shinyWidgets::virtualSelectInput("Checkbox", label = "Phenotypes to display",
                                               choices = sort(Phenotypes),
                                               selected = sort(Phenotypes),
                                               search = TRUE,
                                               multiple = TRUE)
            ),


            #Finally add a couple of rows more with extra options and the final result
            shiny::fluidRow(
              shiny::column(4, shiny::actionButton("reset", shiny::icon("redo"), label = "Reset selection"))
            ),
            #The UI will be completed with summary tables of the sample
            shiny::fluidRow(
              shiny::column(12, htmltools::p("Sample Summary: ", shiny::tableOutput("Summary")))
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
    #####BUILD THE SERVER####
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

      #Control the checkbox output
      Checkbox_output <- shiny::reactive(input$Checkbox)

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

      #PHOTO AS GGPLOT OBJECT (with adequate axis) that will serve as basis for all the plots with a photo on the background
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
                     color = "red", size = 10, hjust = 0.5,
                     label = "UNABLE TO RENDER PHOTO\nUSE ME TO ZOOM IN")
          return(Photo_plot)
        }

        else{
          if(is.null(Image_coloring_list$Parameter_list)){
            #Generate the plot
            return(Photo)
          }
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

          #Get the coords before filtering for phenotypes (will only be used in case the image returns an error)
          Image_X_min <- min(Final_DATA$X)
          Image_X_max <- max(Final_DATA$X)
          Image_Y_min <- min(Final_DATA$Y)
          Image_Y_max <- max(Final_DATA$Y)


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
                scale_x_continuous(limits = c(Image_X_min, Image_X_max)) +
                scale_y_continuous(limits = c(Image_Y_min, Image_Y_max)) +
                coord_cartesian(xlim = ranges$x, ylim = ranges$y, expand = FALSE)
            )
          }
          #If not proceed as usual
          else{
            return(
              #Generate the final plot
              Photo +
                ggiraph::geom_point_interactive(aes(x = X, y = Y, group = Phenotype, color = Color_code,
                                                    data_id = Cell_no,
                                                    tooltip = stringr::str_c(as.character(Cell_no)," Type = ", as.character(Phenotype)),
                                                    hover_nearest = FALSE),
                                                size = 2,
                                                data = Final_DATA) +
                scale_color_identity()+
                guides(color = "none")
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
