#' @export get.contrast.results
get.contrast.results <- function(deseq_dataset, contrast, 
  lfc_shrink = T, shrink_type = "ashr", filename = F){
	results <- DESeq2::results(deseq_dataset, contrast = contrast)
	if(lfc_shrink){
		shrink_results <- DESeq2::lfcShrink(deseq_dataset, res = results, contrast = contrast, type = shrink_type)
    	results$log2FoldChange <- shrink_results$log2FoldChange
    	results$lfcSE <- shrink_results$lfcSE
  	}
  	final_results <- data.frame(results) %>%
  	dplyr::filter(!is.na(padj)) %>%
  	dplyr::arrange(dplyr::desc(stat))
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
