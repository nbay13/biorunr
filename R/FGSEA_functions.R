#' @export get.rnk.vector
get.rnk.vector <- function(DE_results, column_name = "signed.log.p"){
	rnk <- DE_results[[column_name]]
	names(rnk) <- rownames(DE_results)
	return(rnk)
}

#' @export load.MSigDB
load.MSigDB <- function(species){
	## R package for MSigDB
	msig_df <- msigdbr::msigdbr(species = species)
	cat(paste0("Printing available geneset categories for: ", species, "\n"))
	print(data.frame(msig_df %>% dplyr::distinct(gs_cat, gs_subcat) %>% dplyr::arrange(gs_cat, gs_subcat)))
	return(msig_df)
}

#' @export get.MSigDB.genesets
get.MSigDB.genesets <- function(msig_df, genesets = c("KEGG", "REACTOME", "H$")){
	# Default gets KEGG, Reactome, and Hallmarks genesets
	reg_exp <- paste(genesets, collapse = "|")
	gsea_df <- msig_df %>% dplyr::filter(grepl(reg_exp, gs_subcat) | grepl(reg_exp, gs_cat))
	## format as list
	gsea_list <- gsea_df %>% split(x = .$gene_symbol, f = .$gs_name)
	return(gsea_list)
}

#' @export run.FGSEA
run.FGSEA <- function(rnk, genesets, nproc = 2, minGenes = 3, maxGenes = 5000, reformat = T, filename = F, minP = 1e-20){
	# make sure package is loaded
	gsea_res <- fgsea::fgsea(stats = rnk, pathways = genesets, minSize = minGenes, maxSize = maxGenes, eps = minP, nproc = nproc)
	## order results based on NES
	ord_gsea_res <- gsea_res %>% dplyr::arrange(desc(NES))
	## reformat leadingEdge column
	new_edge <- apply(ord_gsea_res, 1, function(x){
		paste(unlist(x$leadingEdge), collapse = ",")
	})
	if(reformat){
		ord_gsea_res$leadingEdge <- new_edge
		if(is.character(filename)){
			cat(paste0("\nWriting results to file: ", filename, "\nat: ", getwd(), "\n"))
			write.table(ord_gsea_res, filename, sep = "\t", 
				quote = F, row.names = F, col.names = T)
		}
	} else {
		if(is.character(filename)){
			cat("\nError: If writing to file reformat must be set to TRUE")
		}
	}
	return(ord_gsea_res)
}