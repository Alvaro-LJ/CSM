#' Generates triangles given two point coordinates (a COO and a target point)
#'
#' Intended for internal use only
#'
#' @param DATA A dataframe or tibble containing information of coordinates of two data points.
#' @param Angle A numeric value indicating the angle of the triangle edges in the COO vertex
#'
#' @details
#' Used in [Barrier_effect_calculator()]
#'
#'
#' @returns A tibble containing the coordinates of triangle vertex points. Also returns distance and angular position of COO and target cells.
#' @keywords Internal

Triangle_generator <- 
  function(DATA, Angle = 10){
    
    #check names of data
    if(!identical(names(DATA), c("COO_X", "COO_Y", "Target_X", "Target_Y"))) stop("DATA names must be COO_x, COO_Y, Target_X, Target_Y")
    
    #Generate the main DF
    DATA <- dplyr::bind_cols(DATA, tibble(Angle = Angle))
    
    #Angle into radians
    DATA$angle_rad <- Angle * pi / 180
    
    #Compute triangle height
    DATA <- DATA %>% dplyr::mutate(height = sqrt((COO_X - Target_X)^2 + (COO_Y - Target_Y)^2))
    
    #Compute triangle base length
    DATA$base_length <- 2 * DATA$height * tan(DATA$angle_rad / 2)
    
    #Compute the triangle height vector and the perpendicular height
    DATA <- DATA %>% dplyr::mutate(hvec_X = COO_X - Target_X,
                                   hvec_Y = COO_Y - Target_Y) %>%
      dplyr::mutate(h_unit_X = hvec_X / sqrt(hvec_X^2 + hvec_Y^2),
                    h_unit_Y = hvec_Y / sqrt(hvec_X^2 + hvec_Y^2)) %>%
      dplyr::mutate(Perpend_vec_X = -h_unit_Y,
                    Perpend_vec_Y = h_unit_X)
    
    #Base coordinates
    DATA <- DATA %>% dplyr::mutate(Vertex_1_X = Target_X + (base_length/2) * Perpend_vec_X,
                                   Vertex_1_Y = Target_Y + (base_length/2) * Perpend_vec_Y,
                                   
                                   Vertex_2_X = Target_X - (base_length/2) * Perpend_vec_X,
                                   Vertex_2_Y = Target_Y - (base_length/2) * Perpend_vec_Y)
    
    #Angle between COO and target points
    DATA <- DATA %>% 
      dplyr::mutate(COO_tar_angl_rad = atan2(Target_Y - COO_Y, Target_X - COO_X)) %>%
      dplyr::mutate(COO_tar_angl_deg = COO_tar_angl_rad * 180 / pi)
    
    
    #Return final triangle info
    DATA <- DATA %>% dplyr::select(COO_X, COO_Y, Target_X, Target_Y, Vertex_1_X, Vertex_1_Y, Vertex_2_X, Vertex_2_Y, height, COO_tar_angl_rad, COO_tar_angl_deg)
    return(DATA)
    
  }


