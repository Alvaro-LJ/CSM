#' Functions to perform spatial operations over X Y coordinates or images
#'
#' Intended for internal use only
#'
#'
#' @details
#' Used in [ICP_cell_aligner()], [Cell_aligner_function()].
#' 
#' Rotate_points rotates a matrix of X Y spatial coordinates by a provided theta (radians). If non_negative argument is set to TRUE, then negative output coords will not be tolerated.
#' 
#' Flip_xy flips a matrix of X Y spatial coordinates. It can be along the X axis, the Y axis or both.
#' 
#' Centroid_translate_points aligns two matrices of X Y spatial coordinates based on the centroids of the target. The moving matrix is modified.
#'
#' Recover_rigid_rotation provides the approximate angle rotation required to align two matrices of X Y spatial coordinates as long as both have the same length.
#'
#' Translate_image_no_crop translates an image in the X and Y dimension adding extra padding if necessary.
#' 
#' @returns A matrix or an image
#'
#' @keywords Internal

Rotate_points <- function(MATRIX, theta, non_negative = TRUE){
  
  R <- matrix(
    c(
      cos(theta), -sin(theta),
      sin(theta),  cos(theta)
    ),
    2, 2,
    byrow = TRUE
  )
  
  xy_rot <- MATRIX %*% R
  colnames(xy_rot) <- c("X", "Y")
  
  # Enforce X >= 0 and Y >= 0
  if (non_negative) {
    
    min_x <- min(xy_rot[, "X"])
    min_y <- min(xy_rot[, "Y"])
    
    if (min_x < 0)
      xy_rot[, "X"] <- xy_rot[, "X"] - min_x
    
    if (min_y < 0)
      xy_rot[, "Y"] <- xy_rot[, "Y"] - min_y
  }
  
  return(xy_rot)
}

Flip_xy <- function(MATRIX,
                    mode = c("none", "x", "y", "xy")) {
  
  mode <- match.arg(mode)
  
  xy <- as.matrix(MATRIX)
  
  xmin <- min(xy[, 1])
  xmax <- max(xy[, 1])
  
  ymin <- min(xy[, 2])
  ymax <- max(xy[, 2])
  
  if (mode %in% c("x", "xy")) {
    xy[, 1] <- xmax + xmin - xy[, 1]
  }
  
  if (mode %in% c("y", "xy")) {
    xy[, 2] <- ymax + ymin - xy[, 2]
  }
  
  colnames(xy) <- c("X", "Y")
  
  return(xy)
}

Centroid_translate_points <- function(moving, target) {
  
  trans <- colMeans(target) - colMeans(moving)
  
  moving_translated <- sweep(
    moving,
    2,
    trans,
    "+"
  )
  return(moving_translated)
}

Recover_rigid_rotation <- function(original, aligned) {
  
  stopifnot(
    nrow(original) == nrow(aligned),
    ncol(original) == 2,
    ncol(aligned) == 2
  )
  
  #centroids
  centroid_original <- colMeans(original)
  centroid_aligned <- colMeans(aligned)
  
  # center clouds
  X <- sweep(
    original,
    2,
    centroid_original,
    "-"
  )
  
  Y <- sweep(
    aligned,
    2,
    centroid_aligned,
    "-"
  )
  
  # cross-covariance
  H <- t(X) %*% Y
  
  sv <- base::svd(H)
  
  # rotation matrix
  R <- sv$v %*% t(sv$u)
  
  # force proper rotation
  if (det(R) < 0) {
    
    sv$v[,2] <- -sv$v[,2]
    
    R <- sv$v %*% t(sv$u)
  }
  
  # rotation angle
  theta_rad <- atan2(
    R[2,1],
    R[1,1]
  )
  
  theta_deg <- theta_rad * 180 / pi
  
  #If negative turn to positive
  if(theta_deg < 0){
    theta_deg <- 360 + theta_deg
    theta_rad <- theta_deg * pi / 180
  }
  
  #Return the final list
  return(
    list(
      Rotation_radians = theta_rad,
      Rotation_degrees = theta_deg
    )
  )
  
}

Translate_image_no_crop <- function(img, dx, dy, bg = 0) {
  
  d <- dim(img)
  
  h <- d[1]
  w <- d[2]
  
  newdim <- d
  newdim[1] <- h + abs(dy)
  newdim[2] <- w + abs(dx)
  
  out <- EBImage::Image(
    array(bg, dim = newdim),
    colormode = EBImage::colorMode(img)
  )
  
  r0 <- max(1, dy + 1)
  c0 <- max(1, dx + 1)
  
  rows <- r0:(r0 + h - 1)
  cols <- c0:(c0 + w - 1)
  
  lhs <- vector("list", length(newdim))
  lhs[[1]] <- rows
  lhs[[2]] <- cols
  
  if(length(newdim) > 2){
    for(i in 3:length(newdim))
      lhs[[i]] <- seq_len(newdim[i])
  }
  
  do.call(
    "[<-",
    c(
      list(out),
      lhs,
      list(value = img)
    )
  )
}

Pad_image <- function(img, Right_padding, Bottom_padding){
  #Check dimensions of the tile 
  Width <- dim(img)[1]
  Height <- dim(img)[2]
  Channels <- dim(img)[3]
  if(any(is.null(Channels), is.na(Channels))) Channels <- 1
  
  #Compute the white padding for RGB color images
  if(EBImage::colorMode(img) == 2){
    Padding_image <- EBImage::Image(
      array(
        1, #White color
        dim = c(
          Width + Right_padding ,
          Height + Bottom_padding,
          Channels
        )
      ),
      colormode = "Color"
    )
  }
  
  #Compute the black padding for non-RGB color images
  if(EBImage::colorMode(img) == 0){
    Padding_image <- EBImage::Image(
      array(
        0, #Black color
        dim = c(
          Width + Right_padding ,
          Height + Bottom_padding,
          Channels
        )
      ),
      colormode = "Grayscale"
    )
  }
  
  #Place the Tile image in the top-left center of the padded images
  Padding_image[
    seq_len(Width),
    seq_len(Height),
    ] <- img
  
  #Return the padded image
  return(Padding_image)
}