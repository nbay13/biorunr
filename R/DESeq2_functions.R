#' @export filter.genes.for.DE
filter.genes.for.DE <- function(exp_tbl, by = c("mean", "q", "name"), thresh = 50, q = 0.33, gene_names = NULL){
  if(by[1] == "mean") tbl <- exp_tbl[rowMeans(exp_tbl) > thresh,]
  if(by[1] == "q") tbl <- exp_tbl[rowMeans(exp_tbl) > quantile(rowMeans(exp_tbl), q),]
  if(by[1] == "name") tbl <- exp_tbl[rownames(exp_tbl) %in% gene_names,]
  return(tbl)
}

#' @export get.contrast.results
get.contrast.results <- function(deseq_dataset, contrast, 
  lfc_shrink = T, shrink_type = "ashr", filename = F, log.p = T, ...){
	results <- DESeq2::results(deseq_dataset, contrast = contrast, ...)
	if(lfc_shrink){
		shrink_results <- DESeq2::lfcShrink(deseq_dataset, res = results, contrast = contrast, type = shrink_type)
    	results$log2FoldChange <- shrink_results$log2FoldChange
    	results$lfcSE <- shrink_results$lfcSE
  	}
  	final_results <- data.frame(results) %>%
  	dplyr::filter(!is.na(padj)) %>%
  	dplyr::arrange(dplyr::desc(stat))
    if(log.p) final_results <- add.signed.log.p(final_results)
  	if(is.character(filename)){
    	cat(paste0("\nWriting results to file: ", filename, "\nat: ", getwd(), "\n"))
    	write.table(final_results, filename, sep = "\t", quote = F, row.names = T, col.names = NA)
  	}
  	return(final_results)
}

#' @export add.signed.log.p
add.signed.log.p <- function(de_table, col_name = "padj", sign_name = "stat"){
	de_table$signed.log.p <- -log10(de_table[,col_name]) * sign(de_table[,sign_name])
	return(de_table)
}

#' @export write.rnk.file
write.rnk.file <- function(df, filename, metric = "signed.log.p"){
  out_df <- data.frame(gene = rownames(df), value = df[[metric]])
  ord_df <- out_df[order(out_df$value, decreasing = T),]
  write.table(ord_df, filename, row.names = F, col.names = F, quote = F, sep = "\t")
}

#' @export prep.volcano.plot
prep.volcano.plot <- function(de_table, stat_name = "padj", stat_code = "FDR", effect_name = "log2FoldChange", 
  effect_code = "L2FC", sig_name = "sig", stat_thresh = 0.05, effect_thresh = 1){
  de_table$sig <- ifelse(de_table[[stat_name]] < stat_thresh & abs(de_table[[effect_name]]) > effect_thresh, paste(stat_code, "&", effect_code, sep = " "), "NS")
  de_table$sig[de_table[[stat_name]] >= stat_thresh & abs(de_table[[effect_name]]) > effect_thresh] <- effect_code
  de_table$sig[de_table[[stat_name]] < stat_thresh & abs(de_table[[effect_name]]) <= effect_thresh] <- stat_code
  de_table$sig <- factor(de_table$sig, levels = c(paste(stat_code, "&", effect_code, sep = " "), stat_code, effect_code, "NS"))
  colnames(de_table)[colnames(de_table) = "sig"] <- sig_name
  return(de_table)
}