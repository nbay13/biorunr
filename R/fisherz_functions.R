#from psych package, included here to avoid dependencies
#' @export fisherz
fisherz <- function(r) 0.5*log((1+r)/(1-r))  

#from psych package, included here to avoid dependencies
#' @export fisherz2r
fisherz2r <- function(z) (exp(2*z)-1)/(1+exp(2*z))

#from psych package, included here to avoid dependencies
#' @export r2t
r2t <- function(r, n) r * sqrt((n - 2)/(1 - r^2))

#' @export fisherz.avg.cor
fisherz.avg.cor <- function(cor_list, sizes = NULL, return_z = F){
	fisher_list <- lapply(cor_list, fisherz)
	if(!is.null(sizes)){
		if(length(sizes) != length(cor_list)) stop("Number of sizes does not equal number of correlation matrices in list")
		else {
			fisher_list <- lapply(1:length(fisher_list), function(i){
				fisher_list[[i]] * (sizes[i]/sum(sizes)) 
			})
			avg <- do.call("+", fisher_list)
		}
	} else avg <- do.call("+", fisher_list) / length(fisher_list)
	converted <- fisherz2r(avg)
	if(return_z) return(list(cor = converted, z = avg))
	else return(converted)
}

#' @export fisherz.pvalue
fisherz.pvalue <- function(fisherz_mat, n){
	return(2*pnorm(-abs(fisherz_mat * sqrt(n-3))))
}

#' @export fisherz.t.pvalue
fisherz.t.pvalue <- function(fisher_mat, n){
	t_mat <- r2t(fisher_mat, n)
	return(2*(pt(-abs(t_mat), n-2)))
}

#' @export fisherz.test
fisherz.test <- function(cor_list, sizes, weighted = F, method = "z"){
	if(weighted) res_list <- fisherz.avg.cor(cor_list, sizes = sizes, return_z = T)
	else res_list <- fisherz.avg.cor(cor_list, return_z = T)
	if(method == "z") res_list$p.value <- fisherz.pvalue(res_list$z, sum(sizes))
	else if(method == "t") res_list$p.value <- fisherz.t.pvalue(res_list$z, sum(sizes))
	else stop("Unknown test selection")
	return(res_list)
}
