#' @export ggplot.colors
ggplot.colors <- function(n=6, h=c(0, 360) +15){
  if ((diff(h)%%360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}

#' @export get.matched.lines
get.matched.lines <- function(df, id_var, group_var, groups){
	temp <- data.frame(
		df %>% 
		dplyr::group_by(get(id_var)) %>% 
		dplyr::summarise(try = all(groups %in% get(group_var)))
	)
	ids <- temp[temp[,2],1]
	return(ids)
}

#' @importFrom ggplot2 theme element_text margin
#' @export
theme_nb_discrete <- function(){
	theme(axis.text.x = element_text(size = 9, color = "black"), axis.text.y = element_text(size = 8, color = "black"),
		axis.title = element_text(size = 10, color = "black"), legend.position = "top", plot.margin = margin(2.5,2.5,-10.5,2.5, "pt"))
}


#' @importFrom ggplot2 ggplot aes geom_segment geom_violin geom_point scale_fill_manual theme_classic labs
#' @export plot.matched.sample.divergence
plot.matched.sample.divergence <- function(df, palette = NULL, gene, group_var = "Sample.Type", id_var = "Line",  matched = F, dotted = 0, pt_size = 1, lwd = 1, levels = NULL, ylab = NULL){
	df <- df[order(df[[id_var]]),]
	if(!is.factor(df[[group_var]])){
		df[[group_var]] <- factor(df[[group_var]])
	}
	if(!is.null(levels)){
		df[[group_var]] <- factor(df[[group_var]], levels = levels)
	}
	full_matched_lines <- get.matched.lines(df, id_var, group_var, groups = levels(df[[group_var]]))
	if(matched){
		df <- df[df[[id_var]] %in% full_matched_lines,]		
	}
	levs <- levels(df[[group_var]])
	num_pairs <- length(full_matched_lines)
	if(is.null(palette)){
		palette <- ggplot.colors(length(levs))
	}
	string <- ""
	for(i in 1:(length(levs)-1)){
		if(i %in% dotted){
			string <- paste0(string, "geom_segment(aes(x = rep(levs[", i,"], num_pairs), xend = rep(levs[",i + 1,"], num_pairs), y = df[df[[group_var]] == levs[",i,"] & df[[id_var]] %in% full_matched_lines,][[gene]], yend = df[df[[group_var]] == levs[", i+1,"] & df[[id_var]] %in% full_matched_lines,][[gene]]), lwd = ", lwd,", lty = 2) + ")
		} else {
			string <- paste0(string, "geom_segment(aes(x = rep(levs[", i,"], num_pairs), xend = rep(levs[",i + 1,"], num_pairs), y = df[df[[group_var]] == levs[",i,"] & df[[id_var]] %in% full_matched_lines,][[gene]], yend = df[df[[group_var]] == levs[", i+1,"] & df[[id_var]] %in% full_matched_lines,][[gene]]), lwd = ", lwd,") + ")
		}
	}
	if(!is.null(ylab)){
		eval(parse(text = paste0("gg <- ggplot() + geom_violin(aes(x = df[[group_var]], y = df[[gene]], fill = df[[group_var]])) + geom_point(aes(x = df[[group_var]], y = df[[gene]]), size = ", pt_size,") + ",string, "scale_fill_manual(values = palette) + theme_classic() + theme_nb_discrete() + labs(x = '', y = bquote(.(ylab)), fill = '')")))
	} else {
		eval(parse(text = paste0("gg <- ggplot() + geom_violin(aes(x = df[[group_var]], y = df[[gene]], fill = df[[group_var]])) + geom_point(aes(x = df[[group_var]], y = df[[gene]]), size = ", pt_size,") + ",string, "scale_fill_manual(values = palette) + theme_classic() + theme_nb_discrete() + labs(x = '', y = bquote(.(gene)*' log'[2]*' CPM'), fill = '')")))
	}
	
	return(gg)
}
