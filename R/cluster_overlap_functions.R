#' @export prep.cluster.DE
prep.cluster.DE <- function(de_results, p_val_thresh = 0.01, logfc_thresh = 0.5, p_val_col_name = "p_val_adj", fc_col_name = "avg_log2FC"){
	filt_de_results <- de_results[de_results[[p_val_col_name]] < p_val_thresh & abs(de_results[[fc_col_name]]) > logfc_thresh,]
	filt_de_results$pct.1[filt_de_results$pct.1 == 0] <- 0.001
	filt_de_results$pct.2[filt_de_results$pct.2 == 0] <- 0.001
	filt_de_results$score <- 2^filt_de_results[[fc_col_name]] * (filt_de_results$pct.1 / filt_de_results$pct.2)
	return(filt_de_results)
}

#' @export get.cluster.correlations
get.cluster.correlations <- function(de_results, template_results, cor.method = "spearman", score_col_name = "score"){
	results_list <- list()
	de_clusters <- unique(de_results$cluster)
	template_clusters <- unique(template_results$cluster)
	for(i in 1:length(de_clusters)){
		temp_de_cluster <- de_clusters[i]
		temp_de_results <- de_results[de_results$cluster == temp_de_cluster,]
		# make.unique should not be necessary, but was to read in Matthew's results...could be Excel renaming
		rownames(temp_de_results) <- make.unique(temp_de_results$gene)
		temp_df <- data.frame(matrix(ncol = 7, nrow = length(template_clusters)))
		rownames(temp_df) <- template_clusters
		colnames(temp_df) <- c("nTarget Markers", "nTemplate Markers", "nOverlap", "Correlation", "Correlation pvalue","Hypergeometric -log10 pvalue", "Overlap avgL2FC")
		for(j in 1:length(template_clusters)){
			temp_template_cluster <- template_clusters[j]
			temp_template_de_results <- template_results[template_results$cluster == temp_template_cluster,]
			rownames(temp_template_de_results) <- temp_template_de_results$gene
			gene_list <- intersect(temp_template_de_results$gene, temp_de_results$gene)
			if(length(gene_list) > 2){
				cor <- cor.test(
					temp_template_de_results[gene_list, score_col_name], 
					temp_de_results[gene_list, score_col_name], 
					method = cor.method)
				hyper_p <- phyper(
					length(gene_list)-1, 
					length(unique(temp_de_results$gene)), 
					length(unique(union(de_results$gene, template_results$gene))) - length(unique(temp_de_results$gene)), 
					length(unique(temp_template_de_results$gene)), 
					lower.tail = F)
				temp_df[j,] <- round(c(length(unique(temp_de_results$gene)), length(unique(temp_template_de_results$gene)), length(gene_list), cor$estimate, 
					cor$p.value, -log10(hyper_p), mean(temp_de_results[gene_list,grep("FC$", colnames(temp_de_results))])),2)
			} else {
				temp_df[j,] <- round(c(nrow(temp_de_results), nrow(temp_template_de_results), length(gene_list), NA, NA, NA, mean(temp_de_results[gene_list,grep("FC$", colnames(temp_de_results))])),2)
			}

		}
		results_list[[i]] <- temp_df[order(temp_df[,6], decreasing = T),]
	}
	names(results_list) <- de_clusters
	return(results_list)
}

#' @export combine.cluster.results
combine.cluster.results <- function(results_list){
	all_combined <- lapply(results_list, function(x){
		y <- lapply(1:length(x), function(i){
			x[[i]]$cluster <- rownames(x[[i]])
			x[[i]]$source <- names(x)[i]
			return(x[[i]])
		})
		return(do.call(rbind, y))
	})
}

#' @export write.cluster.results
write.cluster.results <- function(combined_results_list, prefix, suffix = "cluster correlation.tsv"){
	for(i in 1:length(combined_results_list)){
		temp_name <- names(combined_results_list)[i]
		write.table(combined_results_list[[i]], paste(prefix, temp_name, suffix), sep = "\t", quote = F, row.names = F, col.names = T)
	}
}

# # read in analysis-ready template data (pre-filtered, pre-scored)
# test_template <- read.table("Allen_humanareas_clustermarkers_v3_for_celltype_annotation.txt", sep = "\t", header = T, stringsAsFactors = F)
# test_template2 <- read.table("Nowakowski_2017_Science_cell_type_DE_table_for_celltype_annotation.txt", sep = "\t", header = T, stringsAsFactors = F)

# # read in un-annotated cluster DE results
# test_results <- read.table("C:/Users/Nick/Downloads/Find_all_markers_Neftel.txt", sep = "\t",
# 	header = T, stringsAsFactors = F)

# # filter  DE results, adjust pct.1 pct.2
# filt_test_results <- prep.cluster.DE(test_results)

# # get annotation results for each cluster with the given template
# cor_res <- get.cluster.correlations(de_results = filt_test_results, template_results = test_template)
# cor_res2 <- get.cluster.correlations(de_results = filt_test_results, template_results = test_template2)

# # combine results from one template across clusters
# final_res <- combine.cluster.results(list(allen = cor_res))
# # output is a list containing all cluster annotation results for each template

# # combine results from two templates
# final_res <- combine.cluster.results(list(allen = cor_res, nowakowski = cor_res2))

# prefix <- "MGreun Neftel 2019 unsupervised clustering"

# for(i in 1:length(final_res)){
# 	temp_name <- names(final_res)[i]
# 	write.table(final_res[[i]], paste(prefix, temp_name, "cluster correlations.tsv"), sep = "\t", quote = F, row.names = F, col.names = T)
# }

