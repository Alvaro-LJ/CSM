#' Modifies neighborhood labels
#'
#' The function allows the user to modify neighborhood labels. The function works with any object created using the [Tiled_Image_Clustering_function()] function.
#'
#' @param DATA An object containing neighborhood labels generated using [Tiled_Image_Clustering_function()] function.
#' @param New_names A character vector indicating the new names. The length must be equal to the number of unique neighborhood labels.
#'
#' @returns Returns a list with tiles per image with the modified labels.
#'
#' @seealso [Tiled_Image_Clustering_function()]
#'
#' @export

Clustered_Tiled_Images_renamer <-
  function(Tiled_images = NULL,
           New_names = NULL) {

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
