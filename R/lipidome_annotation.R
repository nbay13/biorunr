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

#' @export annotate.lipid.species
annotate.lipid.species <- function(lipid_names){
	# set up annotation vectors and lists
	two_chain <- c("PI", "PC", "PE", "PG", "PS", "PG", "DG", "LacCER", "SM", "HexCER", "Cer")
	category_list <- list(
		"Sterol" = c("CE"), 
		"Sphingolipid" = c("Cer", "LacCER", "HexCER", "LCER", "SM"), 
		"Glycerolipid" = c("DG", "TG"), 
		"Fatty.Acyl" = c("FA"),
		"Glycerophospholipid" = c("LPC", "LPE", "PC", "PE"), 
		"Ether" = c("PE.O", "PE.P")
	)
	inv_list <- inverse.list(category_list)
	# clean up lipid name format
	lipid_names <- gsub(" |-|/|\\\\|:|;|_|~", "\\.", lipid_names)
	# split the details in the species names by periods
	temp <- strsplit(lipid_names, "\\.")
	# remove numbers from the class names
	class_name <- unlist(lapply(temp, function(x){
		return(gsub("[0-9]*","", x[[1]]))
	}))
	# pre-define data.frame and fill with entries using for loop
	structure_anno <- data.frame(matrix(nrow = length(class_name), ncol = 4))
	colnames(structure_anno) <- c("Class","Longest.Carbon", "Total.DBs", "Saturation")
	structure_anno[,1] <- class_name
	# This script works based on specific class names. If new input data has different or new 
	# class names than the script or names must be edited.
	for(i in 1:length(class_name)){
		if(class_name[i] %in% c("TG", "TAG")){
			# For TAG class with three tails
			structure_anno[i,2] <- temp[[i]][2]
			structure_anno[i,3] <- temp[[i]][3]
			fa_dbs <- as.numeric(temp[[i]][5])
			if(fa_dbs > 1 | (as.numeric(structure_anno[i,3]) - fa_dbs) > 2){
				structure_anno[i,4] <- "PUFA"
			} else if(fa_dbs == 1 & (as.numeric(structure_anno[i,3]) - fa_dbs) < 2){
				structure_anno[i,4] <- "MUFA"
			} else if(fa_dbs == 0 & as.numeric(structure_anno[i,3]) == 1){
				structure_anno[i,4] <- "MUFA"
			}
			 else if(fa_dbs == 0 & as.numeric(structure_anno[i,3]) == 0){
				structure_anno[i,4] <- "SFA"
			} else {
				structure_anno[i,4] <- "Unknown"
			}
		} else if(class_name[i] == "PE" & temp[[i]][2] %in% c("P","O")){
			# For PE classes with two chains and extra character for PE-P and PE-O
			structure_anno[i,1] <- paste("PE", temp[[i]][2], sep = ".")
			structure_anno[i,2] <- max(as.numeric(temp[[i]][3]), as.numeric(temp[[i]][5]))
			structure_anno[i,3] <- as.numeric(temp[[i]][4]) + as.numeric(temp[[i]][6])
			if(as.numeric(temp[[i]][4]) > 1 | as.numeric(temp[[i]][6]) > 1){
				structure_anno[i,4] <- "PUFA"
			} else if(as.numeric(temp[[i]][4]) == 1 | as.numeric(temp[[i]][6]) == 1){
				structure_anno[i,4] <- "MUFA"
			} else {
				structure_anno[i,4] <- "SFA"
			}
		} else if(class_name[i] %in% two_chain){
			# For other classes with two chains
			structure_anno[i,2] <- max(as.numeric(gsub("d", "", temp[[i]][2])), as.numeric(temp[[i]][4]))
			structure_anno[i,3] <- as.numeric(temp[[i]][3]) + as.numeric(temp[[i]][5])
			if(as.numeric(temp[[i]][3]) > 1 | as.numeric(temp[[i]][5]) > 1){
				structure_anno[i,4] <- "PUFA"
			} else if(as.numeric(temp[[i]][3]) == 1 | as.numeric(temp[[i]][5]) == 1){
				structure_anno[i,4] <- "MUFA"
			} else {
				structure_anno[i,4] <- "SFA"
			}
		} else {
			# For classes with one chain
			structure_anno[i,2] <- temp[[i]][2]
			structure_anno[i,3] <- temp[[i]][3]
			if(as.numeric(temp[[i]][3]) > 1){
				structure_anno[i,4] <- "PUFA"
			} else if(as.numeric(temp[[i]][3]) == 1){
				structure_anno[i,4] <- "MUFA"
			} else {
				structure_anno[i,4] <- "SFA"
			}
		}
	}
	rownames(structure_anno) <- lipid_names
	structure_anno[,2:3] <- apply(structure_anno[,2:3], 2, as.numeric)
	structure_anno$Category <- unlist(apply(structure_anno, 1, function(x){
		if(x[1] %in% names(inv_list)) inv_list[[x[1]]]
		else NA
	}))
	return(structure_anno)
}
