#' @export factor.columns
factor.columns <- function(df, col_names = NULL, levels_list = NULL, chr_to_factor = FALSE){
	if(chr_to_factor){
		chr_col_names <- colnames(df)[which(sapply(df, is.character))]
		for(name in unique(c(chr_col_names, col_names))){
			if(name %in% names(levels_list)) df[[name]] <- factor(df[[name]], levels = levels_list[[name]])
			else df[[name]] <- factor(df[[name]])			
		}
	} else {
		if(any(!col_names %in% colnames(df))) message("You have provided column names not present in the data.frame")
		else if(any(!names(levels_list) %in% col_names))  message("You have provided levels for columns not marked for factor conversion")
		else {
			for(name in col_names){
				if(name %in% names(levels_list)) df[[name]] <- factor(df[[name]], levels = levels_list[[name]])
				else df[[name]] <- factor(df[[name]])
			}
		}	
	}
	return(df)
}

#' @export prop.table.2
prop.table.2 <- function(x, y, margin = 1, digits = 2) {
  if (is.null(digits)) return(prop.table(table(x, y), margin))
  else return(round(prop.table(table(x, y), margin), digits))
}

#' @export split.df
split.df <- function(df, labels, by_row = F){
	if(by_row){
		mats <-  split(data.frame(df), labels)	
	} else {
		temp <- split(data.frame(t(df)), labels)
		mats <- lapply(temp, t)			
	}
	return(mats)
}

#' @export match.features
match.features <- function(data_list, row = T, column = F){
	if(column){
		column_list <- lapply(data_list, colnames)
		shared_columns <- Reduce(intersect, column_list)
	}
	if(row){
		row_list <- lapply(data_list, rownames)
		shared_rows <- Reduce(intersect, row_list)
	} 
	out_list <- lapply(data_list, function(x){
		if(column) x <- x[,shared_columns]
		if(row) x <- x[shared_rows,]
		x
	})
	return(out_list)
}


#' @export melt.and.merge
melt.and.merge <- function(df_list, column_names = c("Var1", "Var2")){
	dims <- lapply(df_list, dim)
	if(length(unique(dims)) == 1){
		mel_list <- lapply(df_list, reshape2::melt)
		mel_list <- lapply(1:length(mel_list), function(i){
			colnames(mel_list[[i]])[ncol(mel_list[[i]])] <- paste(names(mel_list)[i], "value", sep = ".")
			mel_list[[i]]
		})
		out_df <- data.frame(
			setNames(mel_list[[1]][,1:2], column_names), 
			lapply(mel_list, function(x) {
				x[,ncol(x), drop = F]
			})
		)
		return(out_df)
	} else {
		cat("Datasets are not the same size and cannot be merged...\n")
		return(NULL)
	}
}

#' @export overlap.features
overlap.features <- function(dataset_list, row = T, col = F){
	overlap_list <- dataset_list
	if(row){
		row_features <- lapply(dataset_list, rownames)
		row_overlap <- Reduce(intersect, row_features)
		overlap_list <- lapply(overlap_list, function(x){
			x[match(row_overlap, rownames(x)),]
		})
	}
	if(col){
		col_features <- lapply(dataset_list, colnames)
		col_overlap <- Reduce(intersect, col_features)
		overlap_list <- lapply(overlap_list, function(x){
			x[,match(col_overlap, colnames(x))]
		})
	} 
	return(overlap_list)
}

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
		delta <- data.frame(filt_df[filt_df[,group] == comparison[1], group, drop = F],delta)
	}
	return(delta)
}

#' @export centered.breaks
centered.breaks <- function(pal, vals, mid = NULL){
	len <- length(pal)
	if(is.null(mid)) mid <- min(vals) + (max(vals) - min(vals))/2
	breaks <- c(
		seq(min(vals), mid, length.out=ceiling(len/2)+1),
		seq(mid, max(vals), length.out=ceiling(len/2)+1)[-1]
		)
	return(breaks)
}

#' @export abs.max
abs.max <- function(x) max(abs(x))

#' @export lin.map
lin.map <- function(x) (x-min(x))/(max(x)-min(x))

#' @export lin.scl
lin.scl <- function(x, lo, hi) (x-lo)/(hi-lo)

#' @export qin.map
qin.map <- function(x, q = 0.99) (x-quantile(x,1-q))/(quantile(x, q)-quantile(x, 1-q))

#' @export z.scale
z.scale <- function(x, mean, sd) (x - mean) / sd

#' @export get.conf.int.error
get.conf.int.error <- function(x, q = 0.95) sd(x) * qnorm(1 - (1-q)/2) / sqrt(length(x))

#' @export mean.z.score
mean.z.score <- function(exp_matrix, gene_list, sum = F, transpose = F){
	if(transpose){
		sig_matrix <- t(exp_matrix[,colnames(exp_matrix) %in% gene_list])
	} else {
		sig_matrix <- exp_matrix[rownames(exp_matrix) %in% gene_list,]
	}
	scores <- rowSums(scale(t(sig_matrix)))
	if(!sum) scores <- scores / ncol(exp_matrix)
	return(scores)
}

#' @export geom.mean
geom.mean <- function(x) exp(mean(log(x + 0.001))-0.001)

# https://davetang.org/muse/2014/07/07/quantile-normalisation-in-r/
#' @export quantile.normalization
quantile.normalization <- function(df){
	df_rank <- apply(df,2,rank,ties.method="min")
	df_sorted <- data.frame(apply(df, 2, sort))
	df_mean <- apply(df_sorted, 1, mean)
	df_final <- apply(df_rank, 2, index.to.mean, my_mean=df_mean)
	rownames(df_final) <- rownames(df)
	return(df_final)
}

#' @noRd index.to.mean
index.to.mean <- function(my_index, my_mean){
	return(my_mean[my_index])
}

#' @noRd .ls.objects
.ls.objects <- function (pos = 1, pattern, order.by,
                        decreasing=FALSE, head=FALSE, n=5) {
    napply <- function(names, fn) sapply(names, function(x)
                                         fn(get(x, pos = pos)))
    names <- ls(pos = pos, pattern = pattern)
    obj.class <- napply(names, function(x) as.character(class(x))[1])
    obj.mode <- napply(names, mode)
    obj.type <- ifelse(is.na(obj.class), obj.mode, obj.class)
    obj.prettysize <- napply(names, function(x) {
                           format(utils::object.size(x), units = "auto") })
    obj.size <- napply(names, object.size)
    obj.dim <- t(napply(names, function(x)
                        as.numeric(dim(x))[1:2]))
    vec <- is.na(obj.dim)[, 1] & (obj.type != "function")
    obj.dim[vec, 1] <- napply(names, length)[vec]
    out <- data.frame(obj.type, obj.size, obj.prettysize, obj.dim)
    names(out) <- c("Type", "Size", "PrettySize", "Length/Rows", "Columns")
    if (!missing(order.by))
        out <- out[order(out[[order.by]], decreasing=decreasing), ]
    if (head)
        out <- head(out, n)
    out
}

# shorthand
#' @export lsos
lsos <- function(..., n=10) .ls.objects(..., order.by="Size", decreasing=TRUE, head=TRUE, n=n)

#' @export corner
corner <- function(mat) mat[1:min(5,nrow(mat)),1:min(5,ncol(mat))]

#' @export ggplot.colors
ggplot.colors <- function(n=6, h=c(0, 360) +15){
  if ((diff(h)%%360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}
