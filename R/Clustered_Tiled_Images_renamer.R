#' Modifies neighborhood labels
#'
#' The function allows the user to modify neighborhood labels. The function works with any object created using the [Tiled_Image_Clustering_function()] function.
#'
#' @param DATA An object containing neighborhood labels generated using [Tiled_Image_Clustering_function()] function.
#' @param New_names A character vector indicating the new names. The length must be equal to the number of unique neighborhood labels.
#' @param Old_names (OPTIONAL) A character vector indicating the names to be overwritten by New_names. It must have the same length as New_names.
#'
#' @returns Returns a list with tiles per image with the modified labels.
#'
#' @seealso [Tiled_Image_Clustering_function()]
#'
#' @examples
#' \dontrun{
#' Clustered_Tiled_Images_renamer(
#'     Tiled_images = Clustered_Tiled_Images,
#'     New_names = c("Neighborhood_name_1",
#'                   "Neighborhood_name_2",
#'                   "Neighborhood_name_3",
#'                   "Neighborhood_name_4",
#'                   "Neighborhood_name_5")
#')
#' }
#'
#'
#' @export

Clustered_Tiled_Images_renamer <-
  function(Tiled_images = NULL,
           New_names = NULL,
           Old_names = NULL) {

    #If New names need to be assigned from scratch
    if(is.null(Old_names)){
      #Check if provided names are equal to number of hoods
      if(length(New_names) != length(unique(unlist(purrr::map(Tiled_images, function(Image) Image$Cluster_assignment))))) {
        stop(paste0("Provided New_names should match the number of Neighborhoods in the analysis. Number of neighborhoods: ",
                    length(unique(unlist(purrr::map(Tiled_images, function(Image) Image$Cluster_assignment)))),
                    ". Names provided: ", length(New_names)))
      }

      else{
        #Create a names tibble
        names_tibble <- tibble(Cluster_assignment = factor(1:length(unique(unlist(purrr::map(Tiled_images, function(Image) Image$Cluster_assignment))))),
                               New_names = New_names)

        return(purrr::map(Tiled_images, function(Image){
          dplyr::left_join(Image, names_tibble, by = "Cluster_assignment") %>% dplyr::mutate(Cluster_assignment = New_names) %>% dplyr::select(-New_names)
        }))
      }
    }

    #else
    else{
      if(length(Old_names) != length(New_names)) stop("New_names length and Old_names length are not equal")

      All_old_names <- unique(unlist(map(Tiled_images, ~.$Cluster_assignment)))
      if(!all(Old_names %in% All_old_names)){
        Problematic_names <- Old_names[!Old_names %in% All_old_names]
        stop(paste0("The following Old_names are not present in data: ",
                    stringr::str_c(Problematic_names, collapse = ", ")))
      }

      #If everything OK proceed
      New_results <- map(Tiled_images, function(Image){
        Image <- Image
        for(i in 1:length(New_names)){
          Image$Cluster_assignment[Image$Cluster_assignment == Old_names[i]] <- New_names[i]
        }
        return(Image)
      })
      names(New_results) <- names(Tiled_images)
      return(New_results)
    }
  }



