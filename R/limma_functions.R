#' @export make.model.matrix
make.model.matrix <- function(labels){
	mm <- model.matrix(~ 0 + labels)
	colnames(mm) <- gsub("labels", "", colnames(mm))
	return(mm)
}

#' @export limma.fit.contrast
limma.fit.contrast <- function(data, mm, contrast = NULL){
	groups <- colnames(mm)
	if(is.null(contrast)) contrast <- paste(groups[1], "-", groups[2], sep = " ")
	contr <- do.call(limma::makeContrasts, list(contrast, levels = groups))
	fit <- limma::lmFit(data, mm)
	tmp <- limma::contrasts.fit(fit, contr)
	tmp <- limma::eBayes(tmp)
	return(tmp)
}

#' @export two.group.limma
two.group.limma <- function(data, labels, add_student = FALSE, col_name = "condition"){
	mm <- make.model.matrix(data[[col_name]])
	tmp <- limma.fit.contrast(t(data), mm)
	tt <- data.frame(limma::topTable(tmp, sort.by = "P", n = Inf))
	if(add_student) {
		tt$student.t <- as.numeric(tmp$coef/tmp$stdev.unscaled/tmp$sigma)[match(rownames(tt), rownames(tmp$coef))]
		tt$student.P.Value <- 2 * pt(-abs(tt$student.t), df = tmp$df.residual)
		tt$student.adj.P.Val <- p.adjust(tt$student.P.Value, method = "fdr") 
	}
	return(tt)
}
