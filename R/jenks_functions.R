#' @export calc.ss
calc.ss <- function(values, group){
	if(is.numeric(group)){
		group <- factor(group)
	}
	fit <- aov(values ~ group)
	between_ss <- anova(fit)["group", "Sum Sq"]
	within_ss <- anova(fit)["Residuals", "Sum Sq"]
	total_ss <- sum((values - mean(values))^2)
	return(c(between = between_ss, within = within_ss))
}

#' @export get.jenks.labels
get.jenks.labels <- function(values, breaks){
	z <- sapply(values, function(x){
		y <- sum(x <= breaks)
		y[y == length(breaks)] <- length(breaks) - 1
		factor(y)
	})
	return(z)
}

#' @export calc.gvf
calc.gvf <- function(values, group){
	ss <- calc.ss(values, group)
	gvf <- ss[1] / sum(ss)
	return(gvf)
}

#' @export run.jenks
run.jenks <- function(values, n = NULL){
	jenks<-classInt::classIntervals(values, style = "jenks", unique=T, n = n, samp_prop = 1)
	labels <- get.jenks.labels(values, jenks$brks)
	ss <- unname(calc.ss(values, labels))
	return(list(values = values, labels = labels, breaks = jenks$brks, 
		between.ss = ss[1], within.ss = ss[2], 
		GVF = unname(calc.gvf(values, labels))))
}

#' @export optimize.jenks
optimize.jenks <- function(values, min_n = 2, max_n = 6, plot = F){
	df <- data.frame(matrix(nrow = (max_n - min_n), ncol = 2))
	colnames(df) <- c("N", "GVF")
	for(i in 1:(max_n - min_n+1)){
		n = seq(min_n, max_n)[i]
		results <- run.jenks(values, n)
		df[i,] <- c(n, results$GVF)
	}
	if(plot){
		gg <- ggplot2::ggplot(df, aes(x = N, y = 1 - GVF)) + ggplot2::geom_point() + ggplot2::theme_classic()
 		print(gg + ggplot2::labs(y = ""))
	}
	return(df)
}
