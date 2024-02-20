#' @export two.group.row.test
two.group.row.test <- function(data, labels, test = c("t", "w"), var_equal = FALSE, paired = FALSE, adjust_method = "fdr"){
	test <- match.arg(test)
	if(!is.factor(labels)) {
	 	message("labels are not factors, converting to factors")
	 	labels <- factor(labels)
	} 
	if(test == "t") ev <- paste("t.test(data[i,] ~ factor(labels), var.equal = ", var_equal, ", paired = ", paired,  ")", sep = "")
	else if(test == "w") ev <- paste("wilcox.test(data[i,] ~ factor(labels), paired = ", paired, ")", sep = "")
	df <- data.frame(matrix(nrow = nrow(data), ncol = 5))
	colnames(df) <- c("stat", "mean1", "mean2", "dm", "pvalue")
	rownames(df) <- rownames(data)
	inds <- which(biorunr::row.vars(data) > 0)
	for(i in inds){
		t_res <- eval(parse(text = ev))
		mean1 <- mean(data[i, labels == levels(labels)[1]])
		mean2 <- mean(data[i, labels == levels(labels)[2]])
		df[i,] <- c(t_res$statistic, mean1, mean2, mean1 - mean2, t_res$p.value)
	}
	df$padj <- p.adjust(df$pvalue, method = adjust_method)
	return(df)
}