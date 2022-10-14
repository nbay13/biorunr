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

#' @export list.to.gmx
list.to.gmx <- function(list, wd, filename){
	n.obs <- sapply(list, length)
	seq.max <- seq_len(max(n.obs))
	mat <- sapply(list, "[", i = seq.max)
	write.table(mat, paste(wd, filename, sep = "/"), col.names = T, row.names = F, quote = F, sep = "\t")
}


load.genesets.old <- function(){
	#cahoy_inds <- grep("CAHOY", names(MSigDB$C6_ONCOGENIC))
	#cahoy_genesets <- MSigDB$C6_ONCOGENIC[cahoy_inds]
	#names(cahoy_genesets)[names(cahoy_genesets) == "CAHOY_OLIGODENDROCUTIC"] <- "CAHOY_OLIGODENDROCYTIC"
	data(package = "GSVAdata", brainTxDbSets)
	cahoy_genesets <- brainTxDbSets[grep("up", names(brainTxDbSets))]
	names(cahoy_genesets) <- paste0("Cahoy_2008_", names(cahoy_genesets))
	names(cahoy_genesets) <- gsub("_up", "", names(cahoy_genesets))
	old_wd <- getwd()
	setwd("C:/Users/Nick/Desktop/GBM/gene_lists")

	garofano_list <- as.list(read.table("Garofano_2021_GBM_subtypes.txt", sep = "\t", header = T, stringsAsFactors = F))

	names(garofano_list) <- paste0("Garofano_2021_", names(garofano_list))

	bhaduri_df <- read.table("Bhaduri 2020 developmental cell types.txt", fill = T, sep = "\t",
	 stringsAsFactors = F, header = T)

	bhaduri_list <- split(bhaduri_df[,1], bhaduri_df[,2])

	dusart_df <- as.list(read.table("Dusart_2019_brain_cell_type_lists.txt", sep = "\t", fill = T, stringsAsFactors = F, header = T))
	dusart_hundred <- lapply(dusart_df, function(x){
		x[1:103]
	})
	names(dusart_hundred) <- paste0("Dusart_2020_", names(dusart_hundred))

	neftel_df <- read.table("Filtered_Neftel_2019_cellular_states_20220217.txt", sep = "\t", stringsAsFactors = F, header = T)
	neftel_df[neftel_df == "LPPR1"] <- "PLPPR1"
	neftel_df[neftel_df == "GPR56"] <- "ADGRG1"
	neftel_df[neftel_df == "LOC150568"] <- "LINC01102"
	neftel_df[neftel_df == "MLF1IP"] <- "CENPU"
	neftel_genesets <- list(Neftel_2019_MES2 = neftel_df[,1], Neftel_2019_MES1 = neftel_df[,2], Neftel_2019_AC = neftel_df[,3], Neftel_2019_OPC = neftel_df[,4], Neftel_2019_NPC1 = neftel_df[,5], Neftel_2019_NPC2 = neftel_df[,6], Neftel_2019_G1S = neftel_df[,7], Neftel_2019_G2M = neftel_df[,8])

	mcp_genes <- read.table("MCPcounter_genes_reformatted.txt", header = T, sep = "\t", stringsAsFactors = F)
	mcp_genesets <- split(as.character(mcp_genes[,2]), mcp_genes[,1])
	names(mcp_genesets) <- gsub("^", "MCP_", names(mcp_genesets))

	davoli_genes <- read.table("davoli_genes_reformatted.txt", header = T, sep = "\t")
	davoli_genesets <- split(as.character(davoli_genes[,2]), davoli_genes[,1])
	names(davoli_genesets) <- gsub("^", "Davoli_", names(davoli_genesets))

	setwd("C:/Users/Nick/Desktop/GBM/raw_data")
	mes_genes <- read.table("Verhaak_2017_MES.txt", sep = "\t", stringsAsFactors = F, header = F)[,1]
	pn_genes <- read.table("Verhaak_2017_PN.txt", sep = "\t", stringsAsFactors = F, header = F)[,1]
	cl_genes <- read.table("Verhaak_2017_CL.txt", sep = "\t", stringsAsFactors = F, header = F)[,1]
	pn_genes[pn_genes == "ZNF643"] <- "ZFP69B"
	cl_genes[cl_genes == "KIAA0494"] <- "EFCAB14"
	subtype_genesets <- list(Verhaak_2017_MES = mes_genes, Verhaak_2017_PN = pn_genes, Verhaak_2017_CL = cl_genes)

	files <- list.files()
	GO_110_files <- files[grep("110_GO", files)]
	GO_110_genes <- c()
	for(i in 1:length(GO_110_files)){
	  temp <- read.table(GO_110_files[i], header = F, stringsAsFactors = F, sep = "\t")
	  #print(length(setdiff(temp[,2], GO_110_genes)))
	  GO_110_genes <- unique(c(GO_110_genes, temp[,2]))
	}

	GO_011_files <- files[grep("011_GO", files)]
	GO_011_genes <- c()
	for(i in 1:length(GO_011_files)){
	  temp <- read.table(GO_011_files[i], header = F, stringsAsFactors = F, sep = "\t")
	  #print(length(setdiff(temp[,2], GO_011_genes)))
	  GO_011_genes <- unique(c(GO_011_genes, temp[,2]))
	}

	# uaing x instead of dash because R doesn't allow dash in variable name
	GO_0x1x1_files <- files[grep("0-1-1_GO", files)]
	GO_0x1x1_genes <- c()
	for(i in 1:length(GO_0x1x1_files)){
	  temp <- read.table(GO_0x1x1_files[i], header = F, stringsAsFactors = F, sep = "\t")
	  #print(length(setdiff(temp[,2], GO_0x1x1_genes)))
	  GO_0x1x1_genes <- unique(c(GO_0x1x1_genes, temp[,2]))
	}

	GO_genesets <- list(GO_immune = GO_110_genes, GO_neural = GO_011_genes, GO_metabolism = GO_0x1x1_genes)

	setwd("C:/Users/Nick/Desktop/GBM/gene_lists")

	GO_df <- read.table("UCLA Sample Type DE TME module genesets.txt", sep = "\t", stringsAsFactors = F, header = T)

	GO_genesets_v2 <- as.list(GO_df)
	GO_genesets_v2 <- lapply(GO_genesets_v2, na.omit)

	GO_df2 <- read.table("UCLA Sample Type DE only triplets TME module genesets.txt", sep = "\t", stringsAsFactors = F, header = T)

	GO_genesets_v3 <- as.list(GO_df2)
	GO_genesets_v3 <- lapply(GO_genesets_v3, na.omit)


	SREBF_df <- read.table("UCLA SREBF module genesets.txt", sep = "\t", stringsAsFactors = F, header = T)
	SREBF_genesets <- as.list(SREBF_df)
	SREBF_genesets <- lapply(SREBF_genesets, na.omit)

	SREBFv2_df <- read.table("UCLA SREBF module genesets v2.txt", sep = "\t", stringsAsFactors = F, header = T)
	SREBFv2_genesets <- as.list(SREBFv2_df)
	SREBFv2_genesets <- lapply(SREBFv2_genesets, na.omit)
	SREBFv2_genesets <- SREBFv2_genesets[names(SREBFv2_genesets) != "outgroup"]
	names(SREBFv2_genesets) <- paste0(names(SREBFv2_genesets), "_v2")

	lipid_df <- read.table("JMinami lipid genesets.txt", sep = "\t", stringsAsFactors = F, header = T)
	lipid_genesets <- as.list(lipid_df)
	lipid_genesets <- lapply(lipid_genesets, function(x)x[x != ""])

	estimate_df <- read.table("ESTIMATE_gene_signatures.txt", sep = "\t", stringsAsFactors = F, header = T)
	estimate_genesets <- list(ESTIMATE.Stromal = estimate_df[,1], ESTIMATE.Immune = estimate_df[,2])

	# filbin 2018 paper
	filbin_df <- read.table("Filbin_2018_H3K27M_differentiation_signatures.txt", sep = "\t", stringsAsFactors = F, header = T)
	# fixing aliases for compatibility with TOIL output
	filbin_df[,1][filbin_df[,1] == "MLF1IP"] <- "CENPU"
	filbin_df[,3][filbin_df[,3] == "KAL1"] <- "ANOS1"
	# venteicher 2017 paper
	venteicher_df <- read.table("Venteicher_2017_IDH_differentiation_signatures.txt", sep = "\t", stringsAsFactors = F, header = T)
	fixed_venteicher_df <- apply(venteicher_df, 2, function(x) {gsub("^\\s*", "", x)})
	# patel 2014 paper
	patel_genes <- read.table("patel_stem_signature.txt", header = T, sep = "\t", stringsAsFactors = F)[,1]
	patel_hypo_genes <- read.table("patel_hypoxia_signature.txt", header = F, sep = "\t", stringsAsFactors = F)[,1]
	patel_hypo_genes[patel_hypo_genes == "FAM115C"] <- "TCAF2"
	patel_hypo_genes[patel_hypo_genes == "C5orf62"] <- "SMIM3"

	liau_df <- read.table("Liau_2017_GBM_chromatin_remodeling_lists.gmx", sep = "\t", stringsAsFactors = F, header = T)
	suva_genes <- read.table("Suva_2014_stem_tfs.txt", sep = "\t", stringsAsFactors = F, header = F)[,1]
	zhang_df <- read.table("Zhang_2016_brain_signatures.txt", sep = "\t", stringsAsFactors = F, header = T)

	dar_df <- read.table("my_darmanis_gene_sigs.txt", sep = "\t", header = T, stringsAsFactors = F)
	dar_list <- as.list(dar_df)
	names(dar_list) <- paste0("Darmanis_2017_", names(dar_list))
	names(dar_list)[1] <- "Darmanis_2017_Astrocyte" 

	filbin_genesets <- list(Filbin_2018_cell_cycle = filbin_df[,1], Filbin_2018_OC = filbin_df[,2], Filbin_2018_AC = filbin_df[,3], Filbin_2018_OPC_shared = filbin_df[,4], Filbin_2018_OPC_variable = filbin_df[,5])
	venteicher_genesets <- list(Venteicher_2017_Oligo = fixed_venteicher_df[,1], Venteicher_2017_Astro = fixed_venteicher_df[,2], Venteicher_2017_Stem = fixed_venteicher_df[,3])
	patel_genesets <- list(Patel_2014_GSC = patel_genes, Patel_2014_hypoxia = patel_hypo_genes)
	liau_genesets <- list(Liau_2017_cell_cycle = liau_df[,4][liau_df[,4] != ""], Codega_2014_qNSC = liau_df[,6][liau_df[,6] != ""], Bobadilla_2015_qNSC = liau_df[,7][liau_df[,7] != ""], Martynoga_2013_qNSC = liau_df[,9][liau_df[,9] != ""])
	suva_geneset <- list(Suva_2014_stem = suva_genes)
	zhang_genesets <- list(Zhang_2016_Endothelial = zhang_df[zhang_df$Cell.type == "Endothelial",2], Zhang_2016_Astrocyte = zhang_df[zhang_df$Cell.type == "Astrocyte",2], 
		Zhang_2016_Oligodendrocyte = zhang_df[zhang_df$Cell.type == "Oligodendrocyte",2], Zhang_2016_Myeloid = zhang_df[zhang_df$Cell.type == "Myeloid",2], Zhang_2016_Neuron = zhang_df[zhang_df$Cell.type == "Neuron",2])

	# tme-d sig
	tme_signatures <- read.table("TME.D vs. TME.I DESeq2 lfcS signature.tsv", sep = "\t", header = T, stringsAsFactors = F)
	tme_v2_signatures <- read.table("TME.D vs. TME.I DESeq2 lfcS v2 signature.tsv", sep = "\t", header = T, stringsAsFactors = F)

	splt <- split(rownames(tme_signatures), tme_signatures[,1])
	splt_v2 <- split(rownames(tme_v2_signatures), tme_v2_signatures[,1])
	names(splt_v2) <- c("sTME.D", "sTME.I")


	ptme_signatures <- read.table("pTME.D vs. pTME.I DESeq2 lfcS signature.tsv", sep = "\t", header = T, stringsAsFactors = F)
	ptme_v2_signatures <- read.table("pTME.D vs. pTME.I DESeq2 lfcS v2 signature.tsv", sep = "\t", header = T, stringsAsFactors = F)

	psplt <- split(rownames(ptme_signatures), ptme_signatures[,1])
	psplt_v2 <- split(rownames(ptme_v2_signatures), ptme_v2_signatures[,1])
	names(psplt) <- c("pTME.D", "pTME.I")
	names(psplt_v2) <- c("psTME.D", "psTME.I")

	type_signatures <- as.list(read.table("Llaguno 2015 Type 1 Type 2 genesets.txt", sep = "\t", header = T, stringsAsFactors = F)[,1:8])
	type_table <- read.table("Llaguno 2015 Type 1 vs Type 2.txt", sep = "\t", header = T, stringsAsFactors = F)
	type_list <- list(Type1.mod = type_table$ID[type_table$log.Fold.Change < 0], Type2.mod = type_table$ID[type_table$log.Fold.Change > 0])


	creb_sigs <- as.list(read.table("CREB gene sets.txt", sep = "\t", header = T, stringsAsFactors = F))
	names(creb_sigs) <- c("Pardo_2017_Astro_CREB", "Pardo_2017_Neuro_CREB")

	connectivity <- as.list(read.table("Hai_2021_biorxiv_connectivity_signature.txt", sep = "\t", header = T, stringsAsFactors = F))
	names(connectivity) <- "Hai_2021_Connectivity"

	lgg_subtypes <- as.list(read.table("TCGA LGG 2015 subtype gene signatures.txt", sep = "\t", header = T, stringsAsFactors = F))

	scav_synth <- as.list(read.table("ScavSynthGenelist.txt", sep = "\t", header = T, stringsAsFactors = F))
	cell_types <- as.list(read.table("GBM scRNA Seuratv4 gbm.filtered.batch3 STACAS.integrated Cell.Type signatures.txt", sep = "\t", header = T, stringsAsFactors = F))

	names(cell_types) <- paste0("UCLA_", names(cell_types))

	diapause_sigs <- as.list(read.table("Rehman 2021 diapause signatures.txt", sep = "\t", header = T, stringsAsFactors = F))

	lipid_sigs <- as.list(read.table("GBM scRNA lipid gene clusters 20220718.txt", sep = "\t", header = T, stringsAsFactors = F))

	lipid_sigs <- lipid_sigs[!names(lipid_sigs) %in% c("Storage_lipid", "Structure_lipid")]

	curated_lipid_sigs <- as.list(read.table("20221004 GO storage and structure lipid signatures.txt", sep = "\t", header = T, stringsAsFactors = F))

	names(curated_lipid_sigs) <- paste(names(curated_lipid_sigs), "_lipid", sep = "")

	superset <- c(neftel_genesets, subtype_genesets, filbin_genesets, venteicher_genesets, patel_genesets, liau_genesets, suva_geneset, 
		zhang_genesets, cahoy_genesets, dusart_hundred, dar_list, bhaduri_list, GO_genesets, GO_genesets_v2, GO_genesets_v3, lipid_genesets, SREBF_genesets, SREBFv2_genesets, splt, splt_v2, psplt, psplt_v2, mcp_genesets, davoli_genesets, estimate_genesets, 
		type_signatures, type_list, garofano_list, creb_sigs, connectivity, lgg_subtypes, scav_synth, cell_types, diapause_sigs, lipid_sigs, curated_lipid_sigs)
	


	filt_sets <- lapply(superset, function(x){
		x[x != "" & !is.na(x)]
	})


	source("C:/Users/Nick/Desktop/GBM/scripts/FGSEA_functions.R")
	msig_df <- load.MSigDB("Homo sapiens")
	hallmarks <- get.MSigDB.genesets(msig_df, genesets = c("H$"))

	all_sets <- c(filt_sets, hallmarks)

	#list.to.gmx(all_sets, "C:/Users/Nick/Desktop/GBM/gene_lists/", "GBM master ssGSEA genesets new.txt")
	setwd(old_wd)
	return(filt_sets)
}


