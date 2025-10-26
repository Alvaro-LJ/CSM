#' Separates RGB color images into color channels according to user preferences
#'
#' `Image_deconvolution_function()` generates multi-channel tiff images according to color deconvolution parameters.
#' @param Directory Character string specifying the path to the folder where RGB images are present.
#' @param Output_directory Character string specifying the path to the folder where  output images are written. It must be an empty folder.
#' @param Deconvolution_parameters A list containing parameters. Created using [Color_deconvolution_App_launcher()].
#' @param N_cores Integer. Number of cores to parallelize your computation.
#' @returns The function writes multi-tiff images in the ouput directory.
#' @seealso [Color_deconvolution_App_launcher()]
#'
#' @examples
#' |dontrun{
#' #Create temporary input and output directories------------------------------
#' Input_Dir <- tempfile(pattern = "tempdir1_Input")
#' Output_Dir <- tempfile(pattern = "tempdir2_Output")
#' dir.create(Input_Dir, recursive = TRUE)
#' dir.create(Output_Dir, recursive = TRUE)
#'
#' #Save images in Input directory
#' purrr::map(1:2,
#' function(Image){
#'    EBImage::writeImage(CSM_MiniHE_test[[Image]], file.path(Input_Dir, names(CSM_MiniHE_test)[Image]))
#' })
#'
#' #Print an example of the deconvolution parameters obtained using the App------------
#' print(CSM_colordeconv_test)
#'
#' #Run color channel separation------------------------------------------
#' Image_deconvolution_function(
#'        Directory = Input_Dir,
#'        Output_directory = Output_Dir,
#'        Deconvolution_parameters = CSM_colordeconv_test,
#'        N_cores = 1
#' )
#'
#'#Check that files have been created-----------------------------------------
#'list.files(Output_Dir)
#'
#'#Remove directories---------------------------------------------------------
#'unlink(c(Input_Dir, Output_Dir), recursive = TRUE)
#' }
#'
#' @export

Image_deconvolution_function <-
  function(Directory = NULL,
           Output_directory = NULL,
           Deconvolution_parameters = NULL,
           N_cores = NULL){
    #Check required packages are installed
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


    print("Checking supplied arguments")
    #Check arguments
    if(!length(dir(Directory)) >= 1) stop("Directory does not contain any files") #Check directory files
    if(berryFunctions::is.error(dir(Output_directory))) stop("Invalid output directory") #Check output directory
    if(length(dir(Output_directory)) != 0) stop("Output directory should be an empty folder") #Check that output directory is empty
    if(!all(N_cores >= 1 & N_cores%%1 == 0)) stop("N_cores must be an integer value > 0")

    #The expected names in each element of the parameter list
    Adequate_Names <- c("Brightness","Saturation","Hue", "Normalize", "Equalize", "Contrast", "Sharpen", "Reduce_color",
                        "Max_color", "RED_value", "GREEN_value", "BLUE_value", "Color_Tolerance", "Final_Tolerance",
                        "Post_normalize", "Post_equalize", "Erode", "Erode_kern", "Erode_rounds", "Dilate", "Dilate_kern", "Dilate_rounds")
    Deconvolution_Parameters <- Deconvolution_parameters
    Names_ok <-purrr::map_lgl(Deconvolution_Parameters, function(Channel) identical(names(Channel), Adequate_Names)) #Check that the names are equal
    if(!all(Names_ok)){
      stop(paste0("The following channels do not contain all the required deconvolution parameters: ",
                  stringr::str_c(names(Deconvolution_Parameters)[!Names_ok], collapse = ", "),
                  ". The required parameters should be the following in the following exact order: ",
                  stringr::str_c(Adequate_Names, collapse = ", ")))
    }

    #Check parameters list
    Error_messages <-purrr::map(Deconvolution_Parameters, function(Channel){
      #Conditions to test
      Testing_OK <- c(Brightness_OK = all(is.numeric(Channel[["Brightness"]]), Channel[["Brightness"]] >= 0, Channel[["Brightness"]] <= 100),
                      Saturation_OK = all(is.numeric(Channel[["Saturation"]]), Channel[["Saturation"]] >= 0, Channel[["Saturation"]] <= 100),
                      Hue_OK = all(is.numeric(Channel[["Hue"]]), Channel[["Hue"]] >= 0, Channel[["Hue"]] <= 100),
                      Normalize_OK = is.logical(Channel[["Normalize"]]),
                      Equalize_OK = is.logical(Channel[["Equalize"]]),
                      Contrast_OK = is.logical(Channel[["Contrast"]]),
                      Sharpen_OK = all(is.numeric(Channel[["Sharpen"]]), Channel[["Sharpen"]] > 0),
                      Reduce_color_OK = is.logical(Channel[["Reduce_color"]]),
                      Max_color_OK = all(is.numeric(Channel[["Max_color"]]), Channel[["Max_color"]] >= 1),
                      RED_value_OK = Channel[["RED_value"]] >= 0 & Channel[["RED_value"]] <= 255,
                      GREEN_value_OK = Channel[["GREEN_value"]] >= 0 & Channel[["GREEN_value"]] <= 255,
                      BLUE_value_OK = Channel[["BLUE_value"]] >= 0 & Channel[["BLUE_value"]] <= 255,
                      Color_Tolerance_value_OK = all(is.numeric(Channel[["Color_Tolerance"]]), length(Channel[["Color_Tolerance"]]) == 3, all(Channel[["Color_Tolerance"]] >= 0), all(Channel[["Color_Tolerance"]] <= 1)),
                      Final_Tolerance_value_OK = all(is.numeric(Channel[["Final_Tolerance"]]), Channel[["Final_Tolerance"]] >= 0, Channel[["Final_Tolerance"]] <= 1),
                      Post_normalize_OK = is.logical(Channel[["Post_normalize"]]),
                      Post_equalize_OK = is.logical(Channel[["Post_equalize"]]),
                      Erode_OK = is.logical(Channel[["Erode"]]),
                      Erode_kern_OK = all(is.numeric(Channel[["Erode_kern"]]), Channel[["Erode_kern"]]>0),
                      Erode_rounds_OK = all(is.numeric(Channel[["Erode_rounds"]]), Channel[["Erode_rounds"]]>0),
                      Dilate_OK = is.logical(Channel[["Dilate"]]),
                      Dilate_kern_OK = all(is.numeric(Channel[["Dilate_kern"]]), Channel[["Dilate_kern"]]>0),
                      Dilate_rounds_OK = all(is.numeric(Channel[["Dilate_rounds"]]), Channel[["Erode_rounds"]]>0)
      )

      #Messages to deploy
      Messages_to_return <- c(Brightness_OK = "Brightness must be a positive numeric value between 0 and 100",
                              Saturation_OK = "Saturation must be a positive numeric value between 0 and 100",
                              Hue_OK = "Hue must be a positive numeric value between 0 and 100",
                              Normalize_OK = "Normalize must be a logical value",
                              Equalize_OK = "Equalize must be a logical value",
                              Contrast_OK = "Contrast must be a logical value",
                              Sharpen_OK = "Sharpen must be a numeric value larger than 0",
                              Reduce_color_OK = "Reduce_color must be a logical value",
                              Max_color_OK = "Max_color must be a numeric value equal or larger than 1",
                              RED_value_OK = "RED_value must be a numeric value between 0 and 255",
                              GREEN_value_OK = "GREEN_value must be a numeric value between 0 and 255",
                              BLUE_value_OK = "BLUE_value must be a numeric value between 0 and 255",
                              Color_Tolerance_value_OK = "Tolerance value must be a numeric vector of length = 3, containing values between 0 and 1",
                              Final_Tolerance_value_OK = "Final_Tolerance must be a single numeric value between 0 and 1",
                              Post_normalize_OK = "Post_normalize must be a logical value",
                              Post_equalize_OK = "Post_equalize must be a logical value",
                              Erode_OK = "Erode must be a logical value",
                              Erode_kern_OK = "Erode_kern must be a numeric value equal to or higher than 1",
                              Erode_rounds_OK = "Erode_rounds must be a numeric value equal to or higher than 1",
                              Dilate_OK = "Dilate must be a logical value",
                              Dilate_kern_OK = "Dilate_kern must be a numeric value equal to or higher than 1",
                              Dilate_rounds_OK = "Dilate_rounds must be a numeric value equal to or higher than 1"
      )
      Messages_to_return[!Testing_OK]
    })
    #If badly specified then return an error with adequate messages
    if(any(purrr::map_lgl(Error_messages, ~length(.)>0))) {
      Error_messages <- Error_messages[map_lgl(Error_messages, ~length(.)>0)]
      Messages <-purrr::map_chr(1:length(Error_messages), function(Index){
        paste0(names(Error_messages)[Index], ": ", stringr::str_c(Error_messages[[Index]], sep = ". "))
      })
      stop(cat(Messages, fill = length(Messages)))
    }

    #Before proceeding print a summary and ask the user if they want to proceed
    print(paste0(length(dir(Directory)), " images will be processed"))
    print(paste0(length(Deconvolution_Parameters), " channels will be obtained: ", stringr::str_c(names(Deconvolution_Parameters), collapse = ", ")))
    print(paste0(N_cores, " cores will be used in the computation"))
    answer <- menu(c("Proceed", "Abort"), title = "Should the analysis proceed")

    #If user decides to stop then abort function and return stop message
    if(answer == 2) stop("The function has been stopped.")

    print("Obtaining individual channels for every image and writing results")
    #save exit function if parallelization fails
    on.exit({
      future::plan("future::sequential")
      gc()
    })

    #Now we calculate our distance matrix
    future::plan("future::multisession", workers = N_cores)
    options(future.globals.maxSize = Inf, future.rng.onMisuse = "ignore")
    furrr::furrr_options(scheduling = Inf)


    #Ready to execute. It will be image based (image directory is the input and writing a file is the output)
    furrr::future_map(seq_along(1:length(dir(Directory))), function(Image_index){

      #Get the original image as a magick object
      Image <- magick::image_read(dir(Directory, full.names = TRUE)[Image_index])

      #Obtain a list with each element corresponding to each of the channels
      List_of_image_channel <-purrr::map(Deconvolution_Parameters,
                                         function(Channel_Parameters){
                                           Channel <- Channel_deconvolution_function(Image = Image,
                                                                                     Parameters = Channel_Parameters)
                                           return(Channel)
                                         })
      FINAL_Image <- EBImage::combine(List_of_image_channel)
      FINAL_Image_Name <- stringr::str_c(Output_directory, "/", dir(Directory, full.names = FALSE)[Image_index], "_Channels.tiff")

      EBImage::writeImage(FINAL_Image, FINAL_Image_Name, type = "tiff")
    },
    .progress = TRUE)

    future::plan("future::sequential")
    gc()

    #Finally also print some messages
    print(paste0("Processed images can be found at: ", Output_directory))
    print(tibble(Channel = stringr::str_c("Channel_", 1:length(Deconvolution_Parameters), sep = ""),
                 Name = names(Deconvolution_Parameters)
    )
    )
  }
