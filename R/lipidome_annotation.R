# based on inverseList from BioCor
#' @export inverse.list
inverse.list <- function(x){
	stopifnot(length(names(x)) == length(x))
    stopifnot(all(sapply(x, function(x) {
        is.character(x) || is.na(x)
    })))
	values <- unlist(x, use.names = FALSE)
    names <- rep(names(x), lengths(x))
    split(names, values)
}

#' @export get.lipid.category
get.lipid.category <- function(species){
	category_list <- list(
		"Sterol" = c("CE"), 
		"Sphingolipid" = c("Cer", "LacCER", "HexCER", "LCER", "SM", "dhCer"), 
		"Glycerolipid" = c("DG", "TG", "DAG", "TAG"), 
		"Fatty.Acyl" = c("FA"),
		"Glycerophospholipid" = c("LPC", "LPE", "PC", "PE", "PG", "PI", "PS", "PA"), 
		"Ether" = c("PE.O", "PE.P")
	)
	inv_list <- inverse.list(category_list)
	return(unlist(sapply(species, function(x){
		if(x %in% names(inv_list)) inv_list[[x]]
		else NA
	})))
}

#' @export get.chain.group
get.chain.group <- function(lengths){
	group_list <- list(
		"SCFA" = as.character(c(1:4)), 
		"MCFA" = as.character(c(5:12)), 
		"LCFA" = as.character(c(13:21)), 
		"VLCFA" = as.character(c(22:32))
	)
	inv_list <- inverse.list(group_list)
	return(unlist(sapply(lengths, function(x){
		if(x %in% names(inv_list)) inv_list[[as.character(x)]]
		else NA
	})))
}

#' @export annotate.lipid.species
annotate.lipid.species <- function(input_names){
	# set up annotation vectors and lists
	two_chain <- c("PI", "PC", "PE", "PG", "PS", "PG", "DG", "DAG", "LacCER", "SM", "HexCER", "Cer", "dhCer", "PG", "PI", "PS")
	category_list <- list(
		"Sterol" = c("CE"), 
		"Sphingolipid" = c("Cer", "LacCER", "HexCER", "LCER", "SM", "dhCer"), 
		"Glycerolipid" = c("DG", "TG", "DAG", "TAG"), 
		"Fatty.Acyl" = c("FA"),
		"Glycerophospholipid" = c("LPC", "LPE", "PC", "PE", "PG", "PI", "PS", "PA"), 
		"Ether" = c("PE.O", "PE.P")
	)
	# clean up lipid name format
	lipid_names <- gsub(" |-|/|\\\\|:|;|_|~", "\\.", input_names)
	# split the details in the species names by periods
	temp <- strsplit(lipid_names, "\\.")
	# remove numbers from the class names
	class_name <- unlist(lapply(temp, function(x){
		return(gsub("[0-9]*","", x[[1]]))
	}))
	if(!all(class_name %in% unlist(category_list))){
		missing <- class_name[!class_name %in% unlist(category_list)]
		w.missing <- which(!class_name %in% unlist(category_list))
		prnt <- paste(w.missing, missing, sep = ": ", collapse = "\n")
		stop(paste("Unknown lipid classes...\n", prnt))
	} 
	
	# pre-define data.frame and fill with entries using for loop
	structure_anno <- data.frame(matrix(nrow = length(class_name), ncol = 5))
	colnames(structure_anno) <- c("Class","Total.Carbons", "Longest.Tail", "Total.DBs", "Saturation")
	structure_anno[,1] <- class_name
	# This script works based on specific class names. If new input data has different or new 
	# class names than the script or names must be edited.
	for(i in 1:length(class_name)){
		if(class_name[i] %in% c("TG", "TAG")){
			# For TAG class with three tails
			structure_anno[i,2] <- temp[[i]][2]
			structure_anno[i,3] <- gsub("[A-z]*", "", temp[[i]][4])
			structure_anno[i,4] <- temp[[i]][3]
			fa_dbs <- as.numeric(temp[[i]][5])
			if(fa_dbs > 1 | (as.numeric(structure_anno[i,4]) - fa_dbs) > 2){
				structure_anno[i,5] <- "PUFA"
			} else if(fa_dbs == 1 & (as.numeric(structure_anno[i,4]) - fa_dbs) < 2){
				structure_anno[i,5] <- "MUFA"
			} else if(fa_dbs == 0 & as.numeric(structure_anno[i,4]) == 1){
				structure_anno[i,5] <- "MUFA"
			}
			 else if(fa_dbs == 0 & as.numeric(structure_anno[i,4]) == 0){
				structure_anno[i,5] <- "SFA"
			} else {
				structure_anno[i,5] <- "Unknown"
			}
		} else if(class_name[i] == "PE" & temp[[i]][2] %in% c("P","O")){
			# For PE classes with two chains and extra character for PE-P and PE-O
			structure_anno[i,1] <- paste("PE", temp[[i]][2], sep = ".")
			structure_anno[i,2] <- sum(as.numeric(temp[[i]][3]), as.numeric(temp[[i]][5]))
			structure_anno[i,3] <- max(as.numeric(temp[[i]][3]), as.numeric(temp[[i]][5]))
			structure_anno[i,4] <- as.numeric(temp[[i]][4]) + as.numeric(temp[[i]][6])
			if(as.numeric(temp[[i]][4]) > 1 | as.numeric(temp[[i]][6]) > 1){
				structure_anno[i,5] <- "PUFA"
			} else if(as.numeric(temp[[i]][4]) == 1 | as.numeric(temp[[i]][6]) == 1){
				structure_anno[i,5] <- "MUFA"
			} else {
				structure_anno[i,5] <- "SFA"
			}
		} else if(class_name[i] %in% two_chain){
			# For other classes with two chains
			structure_anno[i,2] <- sum(as.numeric(gsub("d", "", temp[[i]][2])), as.numeric(temp[[i]][4]))
			structure_anno[i,3] <- max(as.numeric(gsub("d", "", temp[[i]][2])), as.numeric(temp[[i]][4]))
			structure_anno[i,4] <- as.numeric(temp[[i]][3]) + as.numeric(temp[[i]][5])
			if(structure_anno[i,1] == "Cer" & as.numeric(temp[[i]][3]) == 0) structure_anno[i,1] <- "dhCer"
			if(as.numeric(temp[[i]][3]) > 1 | as.numeric(temp[[i]][5]) > 1){
				structure_anno[i,5] <- "PUFA"
			} else if(as.numeric(temp[[i]][3]) == 1 | as.numeric(temp[[i]][5]) == 1){
				structure_anno[i,5] <- "MUFA"
			} else {
				structure_anno[i,5] <- "SFA"
			}
		} else {
			# For classes with one chain
			structure_anno[i,2] <- temp[[i]][2]
			structure_anno[i,3] <- temp[[i]][2]
			structure_anno[i,4] <- temp[[i]][3]
			if(as.numeric(temp[[i]][3]) > 1){
				structure_anno[i,5] <- "PUFA"
			} else if(as.numeric(temp[[i]][3]) == 1){
				structure_anno[i,5] <- "MUFA"
			} else {
				structure_anno[i,5] <- "SFA"
			}
		}
	}
	#rownames(structure_anno) <- lipid_names
	structure_anno$Species <- input_names
	structure_anno[,2:4] <- apply(structure_anno[,2:4], 2, as.numeric)
	structure_anno$Category <- get.lipid.category(structure_anno$Class)
	structure_anno$Chain <- get.chain.group(structure_anno$Longest.Tail)
	structure_anno[structure_anno$Class == "PA","Chain"] <- NA
	return(structure_anno[,c("Species", "Class", "Category", "Total.Carbons", "Longest.Tail", "Total.DBs", "Saturation", "Chain")])
}

#' @export get.tail.saturation
get.tail.saturation <- function(n_db){
	if(any(!is.numeric(n_db)) | any(n_db < 0) | any(!round(n_db) == n_db)) stop("Not positive whole numbers")
	res <- n_db
	res[n_db == 0] <- "SFA"
	res[n_db == 1] <- "MUFA"
	res[n_db > 1] <- "PUFA"
	return(res)
}

#' @export get.acyl.tails
get.acyl.tails <- function(input_names){
	# set up annotation vectors and lists
	two_chain <- c("PI", "PC", "PE", "PG", "PS", "PG", "DG", "DAG", "LacCER", "SM", "HexCER", "Cer", "dhCer", "PG", "PI", "PS")
	category_list <- list(
		"Sterol" = c("CE"), 
		"Sphingolipid" = c("Cer", "LacCER", "HexCER", "LCER", "SM", "dhCer"), 
		"Glycerolipid" = c("DG", "TG", "DAG", "TAG"), 
		"Fatty.Acyl" = c("FA"),
		"Glycerophospholipid" = c("LPC", "LPE", "PC", "PE", "PG", "PI", "PS", "PA"), 
		"Ether" = c("PE.O", "PE.P")
	)
	# clean up lipid name format
	lipid_names <- gsub(" |-|/|\\\\|:|;|_|~", "\\.", input_names)
	# split the details in the species names by periods
	temp <- strsplit(lipid_names, "\\.")
	# remove numbers from the class names
	class_name <- unlist(lapply(temp, function(x){
		return(gsub("[0-9]*","", x[[1]]))
	}))
	if(!all(class_name %in% unlist(category_list))){
		missing <- class_name[!class_name %in% unlist(category_list)]
		w.missing <- which(!class_name %in% unlist(category_list))
		prnt <- paste(w.missing, missing, sep = ": ", collapse = "\n")
		stop(paste("Unknown lipid classes...\n", prnt))
	} 
	
	# pre-define data.frame and fill with entries using for loop
	structure_anno <- data.frame(matrix(nrow = length(class_name), ncol = 3))
	extras <- list()
	colnames(structure_anno) <- c("Class","Total.Carbons", "Total.DBs")
	structure_anno[,1] <- class_name
	# This script works based on specific class names. If new input data has different or new 
	# class names than the script or names must be edited.
	for(i in 1:length(class_name)){
		if(class_name[i] %in% c("TG", "TAG")){
			# For TAG class with three tails
			structure_anno[i,2] <- gsub("[A-z]*", "", temp[[i]][4])
			fa_dbs <- as.numeric(temp[[i]][5])
			structure_anno[i,3] <- fa_dbs
		} else if(class_name[i] == "PE" & temp[[i]][2] %in% c("P","O")){
			# For PE classes with two chains and extra character for PE-P and PE-O
			structure_anno[i,1] <- paste("PE", temp[[i]][2], sep = ".")
			structure_anno[i,2] <- as.numeric(temp[[i]][3]) 
			extras[[as.character(i)]] <- c(paste("PE", temp[[i]][2], sep = "."), as.numeric(temp[[i]][5]), as.numeric(temp[[i]][6]))
			structure_anno[i,3] <- as.numeric(temp[[i]][4])
		} else if(class_name[i] %in% two_chain){
			# For other classes with two chains
			structure_anno[i,2] <- as.numeric(gsub("d", "", temp[[i]][2]))
			if(structure_anno[i,1] == "Cer" & as.numeric(temp[[i]][3]) == 0) structure_anno[i,1] <- "dhCer"
			extras[[as.character(i)]] <- c(structure_anno[i,1], as.numeric(temp[[i]][4]), as.numeric(temp[[i]][5]))
			structure_anno[i,3] <- as.numeric(temp[[i]][3])
		} else {
			# For classes with one chain
			structure_anno[i,2] <- temp[[i]][2]
			structure_anno[i,3] <- temp[[i]][3]
		}
	}
	#rownames(structure_anno) <- lipid_names
	extra_df <- data.frame(t(as.data.frame(extras)))
	colnames(extra_df) <- colnames(structure_anno)
	extra_df$Species <- input_names[as.numeric(names(extras))]
	structure_anno$Species <- input_names
	structure_anno <- data.frame(rbind(structure_anno, extra_df))
	structure_anno[,2:3] <- apply(structure_anno[,2:3], 2, as.numeric)
	structure_anno$Category <- get.lipid.category(structure_anno$Class)
	structure_anno$Chain <- get.chain.group(structure_anno$Total.Carbons)
	structure_anno$Saturation <- get.tail.saturation(structure_anno$Total.DBs)
	return(structure_anno[order(structure_anno$Species),c("Species", "Class", "Category", "Total.Carbons", "Chain", "Total.DBs", "Saturation")])
}

#' @export abundance.to.percent.total
abundance.to.percent.total <- function(in_filename, out_filename, directory, anno_column, id_var = "Short.ID", exclude_classes = c()){
	# read in file
	tbl <- read.table(paste0(directory, in_filename), sep = "\t", stringsAsFactors = F, header = F, fill = T)
	# split into annotation metadata and data
	lipids <- as.character(tbl[1,(anno_column+1):ncol(tbl)])
	lipid_anno <- annotate.lipid.species(lipids)
	anno <- tbl[-1,1:anno_column]
	data <- tbl[-1,(anno_column+1):ncol(tbl)]
	data <- apply(data, 2, as.numeric)
	data[is.na(data)] <- 0
	colnames(data) <- lipids
	colnames(anno) <- tbl[1,1:anno_column]
	# filter based on Class-level black list
	data <- data[,!lipid_anno$Class %in% exclude_classes]
	lipids <- lipids[!lipid_anno$Class %in% exclude_classes]
	lipid_anno <- lipid_anno[!lipid_anno$Class %in% exclude_classes,]
	# Calculate initial percent total values
	temp <- 
		data %>% 
	    t() %>% data.frame() %>%
	    dplyr::summarise(dplyr::across(colnames(.), ~./sum(.)*100)) %>%
	    magrittr::set_rownames(colnames(data)) %>%
	    plyr::rename(setNames(anno[[id_var]], colnames(.)), warn_duplicated = F)
	# Average replicates
	percent_total <- 
		temp %>% 
		t() %>% data.frame() %>% 
		dplyr::mutate(Sample = colnames(temp)) %>% 
		dplyr::group_by(Sample) %>% 
		dplyr::summarise(dplyr::across(setdiff(colnames(.), "Sample"), mean)) %>% 
		tibble::column_to_rownames('Sample') %>% 
		t() %>% data.frame() %>% 
		magrittr::set_rownames(colnames(data))
	# Check if averaging affected the sum, and re-calculate percent totals if necessary
	if(any(colSums(percent_total) != 100)){
		cat("Some columns no longer sum to 1 after averaging, re-calculating percent total...\n")
		percent_total <- 
			temp %>% 
			t() %>% data.frame() %>% 
			dplyr::mutate(Sample = colnames(temp)) %>% 
			dplyr::group_by(Sample) %>% 
			dplyr::summarise(dplyr::across(setdiff(colnames(.), "Sample"), mean)) %>% 
			tibble::column_to_rownames('Sample') %>% 
			t() %>% data.frame() %>%
			dplyr::summarise(dplyr::across(colnames(.), ~./sum(.)*100)) %>% 
			magrittr::set_rownames(colnames(data))
	}
	cat(paste("Saving data to: ", out_filename, "\nat: ", directory,"\n"))
	final_output <- data.frame(Species = lipids, lipid_anno[,colnames(lipid_anno) != "Species"], percent_total)
	rownames(final_output) <- final_output$Species
	write.table(final_output, paste0(directory, out_filename), sep = "\t", quote = F, row.names = F, col.names = T)
	return(final_output)
}

#' @export total.to.percent.class
total.to.percent.class <- function(lipid_mat, lipid_anno, out_filename, directory){
	combined_df <- data.frame(lipid_anno, lipid_mat)
	class_df <- 
		combined_df %>% 
		dplyr::group_by(Class) %>% 
		dplyr::summarise(Species = Species, dplyr::across(colnames(lipid_mat), list(~./sum(.)*100), .names = "{.col}")) %>% 
		dplyr::ungroup() %>% data.frame()
	rownames(class_df) <- class_df$Species
	class_df[is.na(class_df)] <- 0
	# reorder to match input, group_by will sort alphabetically
	class_df <- class_df[rownames(combined_df),]
	class_df <- data.frame(class_df, lipid_anno[,colnames(lipid_anno != "Species")])
	class_df <- class_df[,c(colnames(lipid_anno), colnames(lipid_mat))]
	cat(paste("Saving data to: ", out_filename, "\nat: ", directory,"\n"))
	write.table(class_df, paste0(directory, out_filename), sep = "\t", quote = F, row.names = F, col.names = T)
	return(class_df)
}

#' @export acyl.tail.to.bulk.species
acyl.tail.to.bulk.species <- function(lipid_mat, lipid_anno, out_filename, directory){
	lipid_mat[lipid_anno$Class %in% c("TG", "TAG"),] <- lipid_mat[lipid_anno$Class %in% c("TG", "TAG"),] / 3
	lipid_mat <- lipid_mat %>% dplyr::summarise(dplyr::across(colnames(lipid_mat), list(~./sum(.)*100), .names = "{.col}"))
	combined_df <- data.frame(lipid_anno, lipid_mat)
	bulk_df <- combined_df %>% dplyr::group_by(Class, Total.Carbons, Total.DBs) %>% dplyr::summarise(collapse = paste(Species, collapse = ","), dplyr::across(colnames(lipid_mat), sum)) %>% dplyr::ungroup() %>% dplyr::mutate(Species = paste(Class, Total.Carbons, Total.DBs, sep = ".")) %>% data.frame()
	rownames(bulk_df) <- bulk_df$Species
	bulk_df <- bulk_df[c("Species", "Class", "Total.Carbons", "Total.DBs", "collapse", colnames(lipid_mat))]
	cat(paste("Saving data to: ", out_filename, "\nat: ", directory,"\n"))
	write.table(bulk_df, paste0(directory, out_filename), sep = "\t", quote = F, row.names = F, col.names = T)
	return(bulk_df)
}

#' @export acyl.tail.to.bulk.species.by.tail
acyl.tail.to.bulk.species.by.tail <- function(lipid_mat, lipid_anno, out_filename, directory){
	two_chain_ambig <- c("PI", "PC", "PE", "PG", "PS", "PG", "DG", "DAG", "PG", "PI", "PS")
	lipid_mat[lipid_anno$Class %in% two_chain_ambig,] <- lipid_mat[lipid_anno$Class %in% two_chain_ambig,] / 2
	temp_mat <- lipid_mat[lipid_anno$Class %in% c("PE.O", "PE.P", "HexCER", "LacCER", "Cer", "dhCer", "SM"),][c(F,T),]
	lipid_mat <- lipid_mat[!lipid_anno$Class %in% c("PE.O", "PE.P", "HexCER", "LacCER", "Cer", "dhCer", "SM"),]
	temp_anno <- lipid_anno[lipid_anno$Class %in% c("PE.O", "PE.P", "HexCER", "LacCER", "Cer", "dhCer", "SM"),][c(F,T),]
	temp_anno$collapse <- temp_anno$Species
	lipid_anno <- lipid_anno[!lipid_anno$Class %in% c("PE.O", "PE.P", "HexCER", "LacCER", "Cer", "dhCer", "SM"),]
	combined_df <- data.frame(lipid_anno, lipid_mat)
	bulk_df <- combined_df %>% dplyr::group_by(Class, Total.Carbons, Total.DBs) %>% dplyr::summarise(collapse = paste(Species, collapse = ","), dplyr::across(colnames(lipid_mat), sum)) %>% dplyr::ungroup() %>% dplyr::mutate(Species = paste(Class, Total.Carbons, Total.DBs, sep = ".")) %>% data.frame()
	rownames(bulk_df) <- bulk_df$Species
	print(head(temp_anno))
	bulk_df <- bulk_df[c("Species", "Class", "Total.Carbons", "Total.DBs", "collapse", colnames(lipid_mat))]
	bulk_df <- data.frame(rbind(bulk_df, data.frame(temp_anno[,c("Species", "Class", "Total.Carbons", "Total.DBs", "collapse")], temp_mat)))
	bulk_df[,colnames(lipid_mat)] <- bulk_df[,colnames(lipid_mat)] %>% dplyr::summarise(dplyr::across(colnames(lipid_mat), list(~./sum(.)*100), .names = "{.col}"))
	cat(paste("Saving data to: ", out_filename, "\nat: ", directory,"\n"))
	write.table(bulk_df, paste0(directory, out_filename), sep = "\t", quote = F, row.names = F, col.names = T)
	return(bulk_df)
}

#' @export acyl.tail.to.bulk.class
acyl.tail.to.bulk.class <- function(lipid_mat, lipid_anno, out_filename, directory){
	lipid_mat[lipid_anno$Class %in% c("TG", "TAG"),] <- lipid_mat[lipid_anno$Class %in% c("TG", "TAG"),] / 3
	lipid_mat <- lipid_mat %>% dplyr::summarise(dplyr::across(colnames(lipid_mat), list(~./sum(.)*100), .names = "{.col}"))
	combined_df <- data.frame(lipid_anno, lipid_mat)
	bulk_df <- combined_df %>% dplyr::group_by(Class) %>% dplyr::summarise(collapse = paste(Species, collapse = ","), dplyr::across(colnames(lipid_mat), sum)) %>% data.frame()
	cat(paste("Saving data to: ", out_filename, "\nat: ", directory,"\n"))
	write.table(bulk_df, paste0(directory, out_filename), sep = "\t", quote = F, row.names = F, col.names = T)
	return(bulk_df)
}


