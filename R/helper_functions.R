#' @export get.deltas
get.deltas <- function(df, id, group, comparison, vars, percent = F, return_all = F){
	filt_df <- df[df[,group] %in% comparison,]
	temp <- data.frame(df %>%
	dplyr::filter(get(group) %in% comparison) %>%
	dplyr::group_by(get(id)) %>%
	dplyr::summarise(check = all(comparison %in% get(group))))
	filt_df <- filt_df[filt_df[,id] %in% temp[temp[,2],1],]
	a <- filt_df[filt_df[,group] == comparison[1],vars, drop = F]
	b <- filt_df[filt_df[,group] == comparison[2],vars, drop = F]
	if(percent){
		delta <- (a-b) / (abs(b))
	} else {
		delta <- a - b
	}
	if(return_all){
		colnames(a) <- paste(comparison[1], vars, sep = ".")
		colnames(b) <- paste(comparison[2], vars, sep = ".")
		colnames(delta) <- paste("Delta", vars, sep = ".")
		delta <- data.frame(filt_df[filt_df[,group] == comparison[1],!colnames(filt_df) %in% c(group, vars), drop = F],a,b,delta)
	} else {
		colnames(delta) <- paste("Delta", vars, sep = ".")
		delta <- data.frame(filt_df[filt_df[,group] == comparison[1],!colnames(filt_df) %in% c(group, vars), drop = F],delta)
	}
	return(delta)
}
#' @export abs.max
abs.max <- function(x) max(abs(x))

#' @export lin.map
lin.map <- function(x) (x-min(x))/(max(x)-min(x))

#' @export qin.map
qin.map <- function(x, q = 0.99){(x-quantile(x,1-q))/(quantile(x, q)-quantile(x, 1-q))}
