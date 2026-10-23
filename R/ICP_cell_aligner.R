#' Performs alignment based on the Iterative Closest Point algorithm
#'
#' Intended for internal use only
#'
#'
#' @details
#' Used in [Cell_aligner_function()].
#' 
#' Given to set of matrices containing X and Y coordinates, the function will exhaustively identify the combination of rotation, translation and
#' flipping that reduces the average distance between each point in the moving dataset and it's closest neighbor in the target dataset.
#' 
#' @returns A list containing the modified Moving_coords, the error of the model and the final model parameters.
#'
#' @keywords Internal

ICP_cell_aligner <- 
  function(Moving_coords, 
           Target_coords, 
           Iterations = 5){
    #First align the centroids of both dataset using a fast centroid based algorithm and then a first ICP algorithm
    Moving_coords_centered <- Centroid_translate_points(moving = Moving_coords, target = Target_coords)
    Moving_coords_centered <- Morpho::icpmat(x = Moving_coords_centered, y = Target_coords, iterations = Iterations, type = "rigid")
    
    #Then generate all possible combinations
    Combination_tibble <- 
      tidyr::expand_grid(Flip_action = c("none", "x", "y"),
                         Angle_degrees = 1:360)
    Combination_tibble$Angle_radians <- Combination_tibble$Angle_degrees * pi / 180
    
    #Now compute the model error for every combination
    Error_values <- 
      purrr::map_dbl(1:nrow(Combination_tibble), function(Combination_index){
        #Get parameters
        Flip_action <- Combination_tibble$Flip_action[Combination_index]
        Angle_radians <- Combination_tibble$Angle_radians[Combination_index]
        
        #Perform algorithm
        Moving_coords_centered <- Flip_xy(Moving_coords_centered, mode = Flip_action) #Perform flipping
        Moving_coords_centered <- Rotate_points(Moving_coords_centered, theta = Angle_radians) #Perform angle turn
        Moving_coords_centered <- Morpho::icpmat(x = Moving_coords_centered, y = Target_coords, iterations = Iterations, type = "rigid")
        
        #Compute error for target to moving and from moving to target
        Final_error_sense <- sqrt(sum(Rvcg::vcgKDtree(Target_coords, Moving_coords_centered, k=1)$distance^2))
        Final_error_antisense <- sqrt(sum(Rvcg::vcgKDtree(Moving_coords_centered, Target_coords, k=1)$distance^2))
        
        return(max(c(Final_error_sense, Final_error_antisense))) #return the max error of both
        
      })
    #Add the error values
    Combination_tibble$Error <- Error_values
   
    #Select the combination with the lowest error
    Selected_combination <- Combination_tibble[which.min(Combination_tibble$Error), ]
    
    #Compute the actual final points by replicating the flipping => rotating => ICP with the winning combination
    Final_points <- Flip_xy(Moving_coords_centered, mode = Selected_combination$Flip_action[1])
    Final_points <- Rotate_points(Final_points, theta = Selected_combination$Angle_radians[1]) 
    Final_points <- Morpho::icpmat(x = Final_points, y = Target_coords, iterations = Iterations, type = "rigid")
    colnames(Final_points) <- c("X", "Y")
    
    #Compute the model taking into account that Morpho::icpmat can also rotate the points
    Initial_points <- Flip_xy(Moving_coords, #we select the initial points before centering to compute the actual translation in the final step
                              mode = Selected_combination$Flip_action[1]
                              )
    
    #Compute the final angle rotation
    Angle_rotation <- Recover_rigid_rotation(original = Initial_points, aligned = Final_points) #will recover the initial rotation + the ICP rotation
    Initial_flipped_rotated <- Rotate_points(Initial_points, theta = -Angle_rotation$Rotation_radians)
    
    #Compute the final translation
    X_translation <- min(as.numeric(Final_points[, 1] - Initial_flipped_rotated[, 1]))
    Y_translation <- min(as.numeric(Final_points[, 2] - Initial_flipped_rotated[, 2]))
    
    
    return(list(Coordinates = Final_points,
                Error = Selected_combination$Error[1],
                Model = list(Flip_action = Selected_combination$Flip_action[1],
                             Angle_radians = -Angle_rotation$Rotation_radians,
                             Angle_degrees = -Angle_rotation$Rotation_degrees,
                             Translation = c(X_shift = X_translation, Y_shift = Y_translation))))
    
    
    
    
  }
