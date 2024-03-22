
#' @export get.rotation.angle
get.rotation.angle <- function(coord_df, group_var, col_names = c("dim1", "dim2"), degree = TRUE){
	if(length(unique(coord_df[[group_var]])) > 2) {
		print("This function only works for two groups")
	} else {
		centr_df <- data.frame(coord_df %>% group_by(get(group_var)) %>% summarise(across(col_names, mean)))
		angle <- atan((centr_df[1,3] - centr_df[2,3]) / (centr_df[1,2] - centr_df[2,2]))
		if(degree) angle <- angle * 180 / pi
	}
	return(angle)
}

#' @export get.rotation.matrix
get.rotation.matrix <- function(theta, degree = TRUE){
  if (degree){
    theta <- theta * pi / 180
  }
  rmat <- matrix(c(cos(theta), -sin(theta),
                   sin(theta), cos(theta)), nrow = 2, byrow = T)
  return(rmat)
}

#' @export get.rotated.coords
get.rotated.coords <- function(coord_df, rotation_matrix, col_names = c("dim1", "dim2")){
	return(t(rotation_matrix %*% t(coord_df[,col_names])))
}

#' @export rotate.coords
rotate.coords <- function(coord_df, group_var, col_names = c("dim1", "dim2"), degree = TRUE){
	angle <- get.rotation.angle(coord_df, group_var, col_names, degree = degree)
	# why negative?!
	rot_mat <- get.rotation.matrix(-angle, degree = degree)
	rot_coords <- get.rotated.coords(coord_df, rot_mat, col_names = col_names)
	out_df <- data.frame(coord_df[[group_var]],rot_coords)
	colnames(out_df) <- c(group_var, col_names)
	return(out_df)
}
