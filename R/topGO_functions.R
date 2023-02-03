#' @export
list2genes <- function(sig_list, background_list){
	background_list <- setdiff(background_list, sig_list)
	genes <- c(rep(1, length(sig_list)), rep(0, length(background_list)))
	names(genes) <- c(sig_list, background_list)
	return(genes)
}

#' @export
code2genes <- function(DE_df, code){
	x <- T
	for(i in 1:length(code)){
		x <- DE_df[,i] == code[i] & x
	}
	sig <- rownames(DE_df)[x]
	#sig <- rownames(DE_df)[DE_df[,1] == code[1] & DE_df[,2] == code[2] & DE_df[,3] == code[3]]
	set <- setdiff(rownames(DE_df), sig)
	genes <- c(rep(1, length(sig)), rep(0, length(set)))
	names(genes) <- c(sig, set)
	return(genes)
}

#' @export run.topGO
run.topGO <- function(geneList, method = "classic", ontology = "BP", size = 1000, tgd = NULL){
	## prepare data
	if(is.null(tgd)){
		tgd <- new("topGOdata", ontology = ontology, allGenes = factor(geneList), nodeSize=2,
		annot=annFUN.org, mapping="org.Hs.eg.db", ID = "alias")
	}
	## run tests
	if(method == "classic"){
		tst <- topGO::runTest(tgd, algorithm = "classic", statistic = "Fisher")
		tbl <- topGO::GenTable(tgd, classic = tst, orderBy = "classic", topNodes = size, numChar = 1000)
	} else if(method == "weight01"){
		tst <- topGO::runTest(tgd, algorithm = "weight01", statistic = "Fisher" )
		tbl <- topGO::GenTable(tgd, weight01 = tst, orderBy = "weight01", topNodes = size, numChar = 1000)
	} else if(method == "parentchild"){
		tst <- topGO::runTest(tgd, algorithm = "parentchild", statistic = "Fisher")
		tbl <- topGO::GenTable(tgd, parentchild = tst, orderBy = "parentchild", topNodes = size, numChar = 1000)
		} else {
		return("no matching method")
	}
	## create table
	tbl$pval <- topGO::score(tst)[tbl$GO.ID]
	return(list(tgd = tgd, tbl = tbl))
}

#' @export search.GO.tbl
search.GO.tbl <- function(tbl, searchList, p.value = 0.05, remove = F){
	if(remove){
		outSearch <- !grepl(paste(searchList, collapse="|"), tbl$Term)
		sig <- tbl[outSearch & tbl$pval < p.value,]
	} else {
		inSearch <- grepl(paste(searchList,collapse="|"), tbl$Term)
		sig <- tbl[inSearch & tbl$pval < p.value,]
	}
	return(sig)
}

#' @export GO2gene.signature
GO2gene.signature <- function(tgd, GOList, geneList){
	mygenes <- topGO::genesInTerm(tgd, GOList)
	temp <- unique(unlist(mygenes))
	final <- temp[temp %in% names(geneList)[geneList == 1]]
	return(final)
}

#' @export gene2GO.signature
gene2GO.signature <- function(tgd, GOList, geneList){
	temp <- topGO::inverseList(topGO::genesInTerm(tgd, GOList))
	final <- temp[geneList]
	return(final)
}

#' @export GO2DE.metrics
GO2DE.metrics <- function(tgd, tbl, geneList, DETable, absolute = F, metric = "log2FoldChange"){
	fold_changes <- vector(length = nrow(tbl))
	sig <- names(geneList)[geneList == 1]
	for(i in 1:nrow(tbl)){
		mygenes <- topGO::genesInTerm(tgd, tbl$GO.ID[i])
		temp <- unique(unlist(mygenes))
		final <- temp[temp %in% sig]
		if(metric == "signed logp"){
			vals <- -log10(DETable[["padj"]][rownames(DETable) %in% final]) * sign(DETable[["stat"]][rownames(DETable) %in% final])
		} else {
			vals <- DETable[[metric]][rownames(DETable) %in% final]
		}
		if(absolute){
			val <- mean(abs(vals))
		} else {
			val <- mean(vals)
		}
		fold_changes[i] <- val
	}
	return(fold_changes)
}