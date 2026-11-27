#' Functions to support the Segmentation_mask_generator function
#'
#' Intended for internal use only
#'
#'
#' @details
#' Used in [Segmentation_mask_generator()]
#'
#' svg_path_to_matrix turns coordiantes captured by plotly mouse movements into a matrix to define shapes.
#'
#' extract_plotly_paths turns shapes into paths that can be ploted by plotly.
#' 
#' Add_polygon_no_overlap adds an sf polygon to a polygon list removing any overlapping area .
#' 
#' Polygon_to_EBImage_mask turns an sf polygon object into an EBImage image containing cell masks.
#' 
#' Mask_to_sf turns an EBImage object containing masks into sf polygons
#'   
#' @returns Several types of objects
#'
#' @keywords Internal

svg_path_to_matrix <- function(path) {
  
  path <- gsub("^M\\s*", "", path)
  path <- gsub("\\s*Z$", "", path)
  
  vertices <- strsplit(
    path,
    "\\s*L\\s*"
  )[[1]]
  
  coords <- do.call(
    rbind,
    lapply(
      vertices,
      function(v) {
        
        as.numeric(
          strsplit(v, ",")[[1]]
        )
        
      }
    )
  )
  
  coords <- rbind(
    coords,
    coords[1, ]
  )
  
  colnames(coords) <- c(
    "x",
    "y"
  )
  
  coords
  
}

#' @rdname svg_path_to_matrix
extract_plotly_paths <- function(shapes) {
  
  if (is.null(shapes)) {
    return(character())
  }
  
  paths <- character()
  
  recurse <- function(x) {
    
    if (is.null(x)) {
      return()
    }
    
    if (
      is.list(x) &&
      !is.null(names(x)) &&
      "path" %in% names(x)
    ) {
      
      paths <<- c(
        paths,
        x$path
      )
      
      return()
      
    }
    
    if (
      is.character(x) &&
      length(x) == 1 &&
      grepl("^M\\s", x)
    ) {
      
      paths <<- c(
        paths,
        x
      )
      
      return()
      
    }
    
    if (is.list(x)) {
      
      lapply(
        x,
        recurse
      )
      
    }
    
  }
  
  recurse(shapes)
  
  unique(paths)
  
}

#' @rdname svg_path_to_matrix
Add_polygon_no_overlap <- function(sf_polygons, new_polygon) {
  
  #Join all current polygons
  existing <- sf::st_union(sf_polygons)
  
  #Check any overlap with new polygon
  if (length(sf::st_intersects(new_polygon, existing)[[1]]) > 0) {
    
    #Remove overlapping areas
    new_polygon$geometry <- sf::st_difference(
      sf::st_geometry(new_polygon),
      existing
    )
  }
  
  #Do not add anything if the polygon is completely overlapping
  if (all(sf::st_is_empty(new_polygon))) {
    return(sf_polygons)
  }
  
  rbind(sf_polygons, new_polygon)
}

#' @rdname svg_path_to_matrix
Polygon_to_EBImage_mask <- function(sf_polygons, Image_width, Image_height){
  
  #Create the canvas of the image
  Image_canvas <- terra::rast(
    xmin = 0,
    xmax = Image_width,
    ymin = 0,
    ymax = Image_height,
    ncols = Image_width,
    nrows = Image_height,
    resolution = 1
  )
  
  #Keep only the polygons, remove any multipolygon
  sf_polygons <- sf_polygons[sf::st_geometry_type(sf_polygons) == "POLYGON", ]
  
  #Add the polygons to the canvas and name them by the polygon id
  Canvas_id <- terra::rasterize(
    terra::vect(sf_polygons),
    Image_canvas,
    field = "id",
    background = 0
  )
  
  #Turn this to a matrix and turn it to EBImage
  Image <- EBImage::Image(
    terra::as.matrix(Canvas_id, wide = TRUE)
  )
  
  #Perform angle correction and flopping
  Image <- EBImage::rotate(Image, angle = 90)
  Image <- EBImage::flop(Image)
  
  #Normalize it to 0-1
  Image <- Image / max(Image)
  
  #Return the final product
  return(Image)
}

#' @rdname svg_path_to_matrix
Mask_to_sf <- function(mask) {
  
  #Turn mask to integer values
  Object_id_tibble <- tibble::tibble(value = unique(sort(as.vector(mask))))
  Object_id_tibble <- Object_id_tibble[Object_id_tibble$value != 0,]
  Object_id_tibble$Object_id <- 1:nrow(Object_id_tibble)
  Match <- match(mask, Object_id_tibble$value)
  mask[!is.na(Match)] <- Object_id_tibble$Object_id[Match[!is.na(Match)]]
  
  #Add third value
  dim(mask)[3] <- 1
  
  #Flip
  mask <- EBImage::flip(mask)
  
  if (inherits(mask, "Image")) {
    mask <- EBImage::imageData(mask)
  }
  
  #Get the unique ID and discard 0 values
  ids <- sort(unique(as.vector(mask)))
  ids <- ids[ids > 0]
  
  geometries <- vector("list", length(ids))
  
  for (i in seq_along(ids)) {
    
    id <- ids[i]
    
    # binary mask for the current cell
    cell <- mask == id
    
    # outer contour
    contour <- EBImage::ocontour(cell)[[1]]
    
    # close ring if necessary
    if (!all(contour[1, ] == contour[nrow(contour), ])) {
      contour <- rbind(contour, contour[1, ])
    }
    
    geometries[[i]] <- sf::st_polygon(
      list(as.matrix(contour))
    )
  }
  
  #Build the sf_object
  sf_object <- sf::st_sf(
    id = ids,
    geometry = sf::st_sfc(geometries)
  )
  
  
  
  #Check that the sf object is error-proof
  #Remove empty objects
  sf_object <- sf_object[!sf::st_is_empty(sf_object), ]
  
  #Remove non-evaluable objects
  NA_in_validation_test <- which(is.na(sf::st_is_valid(sf_object, reason = FALSE)))
  if(length(NA_in_validation_test) > 0) sf_object <- sf_object[-NA_in_validation_test, ]
  
  #Try to make the objects valid
  sf_object <- sf::st_make_valid(sf_object)
  
  #If still there are invalid objects then remove them
  Valid_polygons <- which(sf::st_is_valid(sf_object))
  sf_object <- sf_object[Valid_polygons, ]
  
  #Update the object ID
  sf_object$id <- 1:length(sf_object$geometry)
  
  #Get the centroids
  Centroid_tibble <- dplyr::bind_cols(id = 1:length(sf_object$geometry),
                                      tibble::as_tibble(sf::st_coordinates(sf::st_centroid(sf_object)))
  )
  
  #Return everything in a list
  return(
    list(sf_object = sf_object,
         sf_centroids = Centroid_tibble)
  )
  
}
