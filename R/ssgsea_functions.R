#' @export load.genesets
load.genesets <- function(wd, filename = "GBM Master ssGSEA genesets.txt"){
	genesets <- as.list(read.table(paste(wd, filename, sep = "/"), sep = "\t", stringsAsFactors = F, header = T))
	genesets <- lapply(genesets, function(x){
		x[!is.na(x)]
		})
	return(genesets)
}

#' @export run.ssGSEA
run.ssGSEA <- function(exp_mat, gene_list, method = "ssgsea", norm = F){
	return(GSVA::gsva(exp_mat, gene_list, method = method, ssgsea.norm = norm))
}

#' @export run.ssGSEA2
run.ssGSEA2 <- function(exp_mat, gene_list, method = "ssgsea", norm = F){
	gsvaPar <- GSVA::ssgseaParam(exp_mat, gene_list, normalize = norm)
	return(GSVA::gsva(gsvaPar))
}

#' @export list.to.gmx
list.to.gmx <- function(list, wd, filename){
	n.obs <- sapply(list, length)
	seq.max <- seq_len(max(n.obs))
	mat <- sapply(list, "[", i = seq.max)
	write.table(mat, paste(wd, filename, sep = "/"), col.names = T, row.names = F, quote = F, sep = "\t")
}
