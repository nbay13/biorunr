#' @export plot.pca
plot.pca <- function(pca, meta, color_var = NA, fill_var = NA, shape_var = NA, plot_pcs = c("PC1", "PC2"), color_pal = NA, ellipse = T, size = 1){
	if(length(intersect(rownames(pca$x), rownames(meta))) != length(union(rownames(pca$x), rownames(meta)))){
		stop("Your PCA data and metadata rownames do not match :(\n  Please check your inputs")
	}
	if(is.na(color_pal)){
		color_pal <- biorunr::ggplot.colors(length(unique(meta[,color_var])))
	}
	vars <- pca$sdev^2
	prop_vars <- round(vars / sum(vars) * 100,2)
	names(prop_vars) <- paste0("PC", seq(1:length(prop_vars)))
	pca_df <- data.frame(meta,pca$x[,plot_pcs])
	if(ellipse){
		gg <- ggplot(pca_df, aes_string(x = plot_pcs[1], y = plot_pcs[2], color = color_var)) + geom_point(size = size) + theme_classic() + stat_ellipse() +
			labs(x = paste0(plot_pcs[1], " (", prop_vars[plot_pcs[1]], "%)"), y = paste0(plot_pcs[2], " (", prop_vars[plot_pcs[2]], "%)")) +
			scale_color_manual(values = color_pal) + guides(fill = guide_legend(override.aes = list(color = "transparent"))) +
			theme(axis.text = element_text(size = 10, color = "black"))
	} else {
		gg <- ggplot(pca_df, aes_string(x = plot_pcs[1], y = plot_pcs[2], color = color_var)) + geom_point(size = size) + theme_classic() + 
			labs(x = paste0(plot_pcs[1], " (", prop_vars[plot_pcs[1]], "%)"), y = paste0(plot_pcs[2], " (", prop_vars[plot_pcs[2]], "%)")) +
			scale_color_manual(values = color_pal) + guides(fill = guide_legend(override.aes = list(color = "transparent"))) +
			theme(axis.text = element_text(size = 10, color = "black"))
	}

	print(gg)
	return(gg)
}
