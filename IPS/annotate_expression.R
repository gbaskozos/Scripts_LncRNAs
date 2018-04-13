library("org.Hs.eg.db")
library("data.table")
library(GenomicFeatures)
library("rtracklayer")
library(GenomicAlignments)
library(BiocParallel)
library(GenomicRanges)
library(gplots)
library(biomaRt)

options(stringsAsFactors=FALSE)
PATH <- "/home/george/Desktop/LncRNA_v2/IPS"

load(file = paste(PATH,"/toc_all.RData", sep=""))

load(file = paste(PATH,"/DE/norm_counts.RData", sep=""))

AD2 <- read.csv(file = paste(PATH,"/DE/res_NEURONS_vs_IPS_AD2.csv", sep=""), row.name=1)
AD4 <- read.csv(file = paste(PATH,"/DE/res_NEURONS_vs_IPS_AD4.csv", sep=""), row.name=1)
NHDF <- read.csv(file = paste(PATH,"/DE/res_NEURONS_vs_IPS_NHDF.csv", sep=""), row.name=1)

#################
# Known LncRNAs #
#################
#########################################
########### ENSEMBL LncRNAs #############
#########################################

nc_ensembl <- read.csv(file = paste(PATH,"/HS_genome/ens_lncs_annotation.csv", sep=""), row.name=1)

nc_ensembl <- nc_ensembl[-which(duplicated(nc_ensembl$ensembl_gene_id)),]

rownames(nc_ensembl) <- nc_ensembl$ensembl_gene_id

ens_lncs_names <- rownames(nc_ensembl)

cutoff <- 0
## Apply cut off to condition
toc_coff <- toc[apply(toc, 1, function(x) sum(x[grep("IPS", names(toc))] > cutoff) == length(grep("IPS", names(toc))) | sum(x[grep("NEURONS", names(toc))] > cutoff) == length(grep("NEURONS", names(toc)))), ]

expressed_nc_ensembl <- toc_coff[rownames(toc_coff) %in% as.character(nc_ensembl$ensembl_gene_id), ]

save(file = paste(PATH,"/DE/expressed_nc_ensembl.RData", sep=""), expressed_nc_ensembl)
write.csv(file = paste(PATH,"/DE/expressed_nc_ensembl_toc.csv", sep=""), expressed_nc_ensembl)

nc_rld <- rld_mat[rownames(rld) %in% rownames(expressed_nc_ensembl),] 

nc_ensembl_rld <- rld_mat[rownames(rld_mat) %in% rownames(expressed_nc_ensembl),] 

sig_nc_AD2 <- AD2[rownames(AD2) %in% rownames(expressed_nc_ensembl) & AD2$padj < 0.05 & !is.na(AD2$padj), ]
sig_nc_AD4 <- AD4[rownames(AD4) %in% rownames(expressed_nc_ensembl) & AD4$padj < 0.05 & !is.na(AD4$padj), ]
sig_nc_NHDF <- NHDF[rownames(NHDF) %in% rownames(expressed_nc_ensembl) & NHDF$padj < 0.05 & !is.na(NHDF$padj), ]

#Order by LFC, base mean and - padj
sig_nc_AD2 <- sig_nc_AD2[order(abs(sig_nc_AD2$log2FoldChange), sig_nc_AD2$baseMean, -sig_nc_AD2$padj),]
sig_nc_AD4 <- sig_nc_AD4[order(abs(sig_nc_AD4$log2FoldChange), sig_nc_AD4$baseMean, -sig_nc_AD4$padj),]
sig_nc_NHDF <- sig_NHDF[order(abs(sig_NHDF$log2FoldChange), sig_NHDF$baseMean, -sig_NHDF$padj),]

sig_nc_AD2$description <- nc_ensembl[rownames(sig_nc_AD2),]$description
sig_nc_AD4$description <- nc_ensembl[rownames(sig_nc_AD4),]$description
sig_nc_NHDF$description <- nc_ensembl[rownames(sig_NHDF),]$description

write.csv(file = paste(PATH,"/DE/sig_nc_AD2.csv", sep=""),sig_nc_AD2)
write.csv(file = paste(PATH,"/DE/sig_nc_AD4.csv", sep=""),sig_nc_AD4)
write.csv(file = paste(PATH,"/DE/sig_nc_NHDF.csv", sep=""),sig_NHDF)


pdf(file = paste(PATH,"/DE/expressed_lncs_heatmap.pdf", sep=""))
cluster <- heatmap.2(as.matrix(nc_ensembl_rld), col=redgreen(200), scale="row", trace="none", cexCol=1, cexRow=0.4, margins=c(5,7), key=FALSE, Colv=TRUE, Rowv=TRUE, main="ENSEMBL LncRNAs", dendrogram="column", symkey=TRUE, srtCol=42, adjCol=c(1,1), lmat=rbind(c(0,3), c(0,1), c(2,4)), lhei=c(2,7,1), lwid = c(1,7),labRow=nc_ensembl[rownames(nc_ensembl_rld),]$mgi_symbol) 
dev.off()

#########################################################################################################
antisense_ENSEMBL <- read.table( file = paste(PATH, "/antisense/lncs_antisense_ENSEMBL.bed", sep=""))

lincs_ENSEMBL <- read.table(file = paste(PATH, "/non_antisense/lncs_intergenic_closest_ENSEMBL.bed", sep=""))

#######################################################
antisense_mapping <- antisense_ENSEMBL[,c(4,16)]

antisense_mapping$antisense_symbol <- mapIds(org.Hs.eg.db, keys=as.character(antisense_mapping[,1]), column="SYMBOL", keytype="ENSEMBL", multiVals="first")

antisense_mapping$sense_symbol <- mapIds(org.Hs.eg.db, keys=as.character(antisense_mapping[,2]), column="SYMBOL", keytype="ENSEMBL", multiVals="first")

antisense_eset <- antisense_mapping

load(file = paste(PATH,"/pain_genes.RData", sep=""))

pain_genes <- unique(pain_genes[,c(2,1)])

pain_genes_symbols <- unique(as.character(pain_genes$HGNC.symbol))


AD2 <- read.csv(file = paste(PATH,"/DE/res_NEURONS_vs_IPS_AD2.csv", sep=""), row.name=1)
AD4 <- read.csv(file = paste(PATH,"/DE/res_NEURONS_vs_IPS_AD4.csv", sep=""), row.name=1)
NHDF <- read.csv(file = paste(PATH,"/DE/res_NEURONS_vs_IPS_NHDF.csv", sep=""), row.name=1)


ens_genes_AD2 <- AD2[grep("ENS", rownames(AD2)),]
ens_genes_AD4 <- AD4[grep("ENS", rownames(AD4)),]
ens_genes_NHDF <- NHDF[grep("ENS", rownames(NHDF)),]

######################################################################################################
# AD2
########
system("mkdir -p /home/george/Desktop/LncRNA_v2/IPS/antisense/ensembl_antisense/AD2")

pain <- ens_genes_AD2[rownames(ens_genes_AD2) %in% pain_genes$Gene.stable.ID.1,]
pain_sig <- subset(pain, pain$padj <.05 & !is.na(pain$padj))


sodium_channel <- ens_genes_AD2[grep("SCN", ens_genes_AD2$symbol),]
sodium_channel_sig <- subset(sodium_channel, sodium_channel$padj <.05 & !is.na(sodium_channel$padj))


potassium_channel <- ens_genes_AD2[grep("KCN", ens_genes_AD2$symbol),]
potassium_channel_sig <- subset(potassium_channel, potassium_channel$padj <.05 & !is.na(potassium_channel$padj))

chloride_channel <- ens_genes_AD2[grep("CLC", ens_genes_AD2$symbol),]
chloride_channel_sig <- subset(chloride_channel, chloride_channel$padj <.05 & !is.na(chloride_channel$padj))

calcium_channel <- ens_genes_AD2[grep("CAC", ens_genes_AD2$symbol),]
calcium_channel_sig <- subset(calcium_channel, calcium_channel$padj <.05 & !is.na(calcium_channel$padj))

trp_channel <- ens_genes_AD2[grep("TRP", ens_genes_AD2$symbol),]
trp_channel_sig <- subset(trp_channel, trp_channel$padj <.05 & !is.na(trp_channel$padj))

antisense_eset$lnc_lfc <- 0
antisense_eset$lnc_pvalue <- 0
antisense_eset$lnc_baseMean <- 0
antisense_eset$gene_lfc <- 0
antisense_eset$gene_pvalue <- 0



pb <- txtProgressBar(min = 0, max = nrow(antisense_eset), style = 3)
for (i in 1:nrow(antisense_eset)) {
	
	antisense_eset$gene_lfc[i] <- ens_genes_AD2[as.character(antisense_eset$V16[i]),2]
	antisense_eset$gene_pvalue[i] <- ens_genes_AD2[as.character(antisense_eset$V16[i]),6]
	
			
	antisense_eset$lnc_lfc[i] <- ens_genes_AD2[as.character(antisense_eset$V4[i]),2]
	antisense_eset$lnc_pvalue[i] <- ens_genes_AD2[as.character(antisense_eset$V4[i]),6]
	antisense_eset$lnc_baseMean[i] <- ens_genes_AD2[as.character(antisense_eset$V4[i]),1] 
	setTxtProgressBar(pb, i)
	
	}
close(pb)
cat("\n")


both_significant <- subset(antisense_eset, antisense_eset$gene_pvalue<.05 & antisense_eset$lnc_pvalue<.05)

anticorellated <- subset(antisense_eset, (antisense_eset$lnc_lfc * antisense_eset$gene_lfc) < 0)

anticorellated_sig <- subset(both_significant, (both_significant$lnc_lfc * both_significant$gene_lfc) < 0) 

pain_antisense <- subset(antisense_eset, antisense_eset$sense_symbol %in% pain_genes_symbols)

pain_anticorellated <- subset(pain_antisense, (pain_antisense$lnc_lfc * pain_antisense$gene_lfc) < 0)

sodium_channel_antisense <- subset(antisense_eset, antisense_eset$sense_symbol %in% sodium_channel$symbol)
sodium_channel_both_sig <- subset(sodium_channel_antisense, sodium_channel_antisense$gene_pvalue<.05 & sodium_channel_antisense$lnc_pvalue<.05)
write.csv(sodium_channel_antisense, file = paste(PATH, "/antisense/ensembl_antisense/AD2_sodium_channel_antisense.csv", sep=""))
write.csv(sodium_channel_both_sig, file = paste(PATH, "/antisense/ensembl_antisense/AD2_sodium_channel_antisense_both_sig.csv", sep=""))

potassium_channel_antisense <- subset(antisense_eset, antisense_eset$sense_symbol %in% potassium_channel$symbol)
potassium_channel_both_sig <- subset(potassium_channel_antisense, potassium_channel_antisense$gene_pvalue<.05 & potassium_channel_antisense$lnc_pvalue<.05)
write.csv(potassium_channel_antisense, file = paste(PATH, "/antisense/ensembl_antisense/AD2_potassium_channel_antisense.csv", sep=""))
write.csv(potassium_channel_both_sig, file = paste(PATH, "/antisense/ensembl_antisense/AD2_potassium_channel_antisense_both_sig.csv", sep=""))

chloride_channel_antisense <- subset(antisense_eset, antisense_eset$sense_symbol %in% chloride_channel$symbol)
chloride_channel_both_sig <- subset(chloride_channel_antisense, chloride_channel_antisense$gene_pvalue<.05 & chloride_channel_antisense$lnc_pvalue<.05)
write.csv(chloride_channel_antisense, file = paste(PATH, "/antisense/ensembl_antisense/AD2_chloride_channel_antisense.csv", sep=""))
write.csv(chloride_channel_both_sig, file = paste(PATH, "/antisense/ensembl_antisense/AD2_chloride_channel_antisense_both_sig.csv", sep=""))

calcium_channel_antisense <- subset(antisense_eset, antisense_eset$sense_symbol %in% calcium_channel$symbol)
calcium_channel_both_sig <- subset(calcium_channel_antisense, calcium_channel_antisense$gene_pvalue<.05 & calcium_channel_antisense$lnc_pvalue<.05)
write.csv(calcium_channel_antisense, file = paste(PATH, "/antisense/ensembl_antisense/AD2_calcium_channel_antisense.csv", sep=""))
write.csv(calcium_channel_both_sig, file = paste(PATH, "/antisense/ensembl_antisense/AD2_calcium_channel_antisense_both_sig.csv", sep=""))

trp_channel_antisense <- subset(antisense_eset, antisense_eset$sense_symbol %in% trp_channel$symbol)
trp_channel_both_sig <- subset(trp_channel_antisense, trp_channel_antisense$gene_pvalue<.05 & trp_channel_antisense$lnc_pvalue<.05)
write.csv(trp_channel_antisense, file = paste(PATH, "/antisense/ensembl_antisense/AD2_trp_channel_antisense.csv", sep=""))
write.csv(trp_channel_both_sig, file = paste(PATH, "/antisense/ensembl_antisense/AD2_trp_channel_antisense_both_sig.csv", sep=""))

write.csv(antisense_eset, file = paste(PATH, "/antisense/ensembl_antisense/AD2/antisense_all_expression.csv", sep=""))

write.csv(both_significant, file = paste(PATH, "/antisense/ensembl_antisense/AD2/antisense_both_sig.csv", sep=""))

write.csv(pain_antisense, file = paste(PATH, "/antisense/ensembl_antisense/AD2/pain_antisense.csv", sep=""))

write.csv(anticorellated, file = paste(PATH, "/antisense/ensembl_antisense/AD2/antisense_anticorellated.csv", sep=""))

write.csv(anticorellated_sig, file = paste(PATH, "/antisense/ensembl_antisense/AD2/antisense_anticorellated_sig.csv", sep=""))

write.csv(pain_anticorellated, file = paste(PATH, "/antisense/ensembl_antisense/AD2/pain_anticorellated.csv", sep=""))

########
# AD4
########
system("mkdir -p /home/george/Desktop/LncRNA_v2/IPS/antisense/ensembl_antisense/AD4")

pain <- ens_genes_AD4[rownames(ens_genes_AD4) %in% pain_genes$Gene.stable.ID.1,]
pain_sig <- subset(pain, pain$padj <.05 & !is.na(pain$padj))


sodium_channel <- ens_genes_AD4[grep("SCN", ens_genes_AD4$symbol),]
sodium_channel_sig <- subset(sodium_channel, sodium_channel$padj <.05 & !is.na(sodium_channel$padj))


potassium_channel <- ens_genes_AD4[grep("KCN", ens_genes_AD4$symbol),]
potassium_channel_sig <- subset(potassium_channel, potassium_channel$padj <.05 & !is.na(potassium_channel$padj))

chloride_channel <- ens_genes_AD4[grep("CLC", ens_genes_AD4$symbol),]
chloride_channel_sig <- subset(chloride_channel, chloride_channel$padj <.05 & !is.na(chloride_channel$padj))

calcium_channel <- ens_genes_AD4[grep("CAC", ens_genes_AD4$symbol),]
calcium_channel_sig <- subset(calcium_channel, calcium_channel$padj <.05 & !is.na(calcium_channel$padj))

trp_channel <- ens_genes_AD4[grep("TRP", ens_genes_AD4$symbol),]
trp_channel_sig <- subset(trp_channel, trp_channel$padj <.05 & !is.na(trp_channel$padj))

antisense_eset$lnc_lfc <- 0
antisense_eset$lnc_pvalue <- 0
antisense_eset$lnc_baseMean <- 0
antisense_eset$gene_lfc <- 0
antisense_eset$gene_pvalue <- 0



pb <- txtProgressBar(min = 0, max = nrow(antisense_eset), style = 3)
for (i in 1:nrow(antisense_eset)) {
	
	antisense_eset$gene_lfc[i] <- ens_genes_AD4[as.character(antisense_eset$V16[i]),2]
	antisense_eset$gene_pvalue[i] <- ens_genes_AD4[as.character(antisense_eset$V16[i]),6]
	
			
	antisense_eset$lnc_lfc[i] <- ens_genes_AD4[as.character(antisense_eset$V4[i]),2]
	antisense_eset$lnc_pvalue[i] <- ens_genes_AD4[as.character(antisense_eset$V4[i]),6]
	antisense_eset$lnc_baseMean[i] <- ens_genes_AD4[as.character(antisense_eset$V4[i]),1] 
	setTxtProgressBar(pb, i)
	
	}
close(pb)
cat("\n")


both_significant <- subset(antisense_eset, antisense_eset$gene_pvalue<.05 & antisense_eset$lnc_pvalue<.05)

anticorellated <- subset(antisense_eset, (antisense_eset$lnc_lfc * antisense_eset$gene_lfc) < 0)

anticorellated_sig <- subset(both_significant, (both_significant$lnc_lfc * both_significant$gene_lfc) < 0) 

pain_antisense <- subset(antisense_eset, antisense_eset$sense_symbol %in% pain_genes_symbols)

pain_anticorellated <- subset(pain_antisense, (pain_antisense$lnc_lfc * pain_antisense$gene_lfc) < 0)

sodium_channel_antisense <- subset(antisense_eset, antisense_eset$sense_symbol %in% sodium_channel$symbol)
sodium_channel_both_sig <- subset(sodium_channel_antisense, sodium_channel_antisense$gene_pvalue<.05 & sodium_channel_antisense$lnc_pvalue<.05)
write.csv(sodium_channel_antisense, file = paste(PATH, "/antisense/ensembl_antisense/AD4_sodium_channel_antisense.csv", sep=""))
write.csv(sodium_channel_both_sig, file = paste(PATH, "/antisense/ensembl_antisense/AD4_sodium_channel_antisense_both_sig.csv", sep=""))

potassium_channel_antisense <- subset(antisense_eset, antisense_eset$sense_symbol %in% potassium_channel$symbol)
potassium_channel_both_sig <- subset(potassium_channel_antisense, potassium_channel_antisense$gene_pvalue<.05 & potassium_channel_antisense$lnc_pvalue<.05)
write.csv(potassium_channel_antisense, file = paste(PATH, "/antisense/ensembl_antisense/AD4_potassium_channel_antisense.csv", sep=""))
write.csv(potassium_channel_both_sig, file = paste(PATH, "/antisense/ensembl_antisense/AD4_potassium_channel_antisense_both_sig.csv", sep=""))

chloride_channel_antisense <- subset(antisense_eset, antisense_eset$sense_symbol %in% chloride_channel$symbol)
chloride_channel_both_sig <- subset(chloride_channel_antisense, chloride_channel_antisense$gene_pvalue<.05 & chloride_channel_antisense$lnc_pvalue<.05)
write.csv(chloride_channel_antisense, file = paste(PATH, "/antisense/ensembl_antisense/AD4_chloride_channel_antisense.csv", sep=""))
write.csv(chloride_channel_both_sig, file = paste(PATH, "/antisense/ensembl_antisense/AD4_chloride_channel_antisense_both_sig.csv", sep=""))

calcium_channel_antisense <- subset(antisense_eset, antisense_eset$sense_symbol %in% calcium_channel$symbol)
calcium_channel_both_sig <- subset(calcium_channel_antisense, calcium_channel_antisense$gene_pvalue<.05 & calcium_channel_antisense$lnc_pvalue<.05)
write.csv(calcium_channel_antisense, file = paste(PATH, "/antisense/ensembl_antisense/AD4_calcium_channel_antisense.csv", sep=""))
write.csv(calcium_channel_both_sig, file = paste(PATH, "/antisense/ensembl_antisense/AD4_calcium_channel_antisense_both_sig.csv", sep=""))

trp_channel_antisense <- subset(antisense_eset, antisense_eset$sense_symbol %in% trp_channel$symbol)
trp_channel_both_sig <- subset(trp_channel_antisense, trp_channel_antisense$gene_pvalue<.05 & trp_channel_antisense$lnc_pvalue<.05)
write.csv(trp_channel_antisense, file = paste(PATH, "/antisense/ensembl_antisense/AD4_trp_channel_antisense.csv", sep=""))
write.csv(trp_channel_both_sig, file = paste(PATH, "/antisense/ensembl_antisense/AD4_trp_channel_antisense_both_sig.csv", sep=""))

write.csv(antisense_eset, file = paste(PATH, "/antisense/ensembl_antisense/AD4/antisense_all_expression.csv", sep=""))

write.csv(both_significant, file = paste(PATH, "/antisense/ensembl_antisense/AD4/antisense_both_sig.csv", sep=""))

write.csv(pain_antisense, file = paste(PATH, "/antisense/ensembl_antisense/AD4/pain_antisense.csv", sep=""))

write.csv(anticorellated, file = paste(PATH, "/antisense/ensembl_antisense/AD4/antisense_anticorellated.csv", sep=""))

write.csv(anticorellated_sig, file = paste(PATH, "/antisense/ensembl_antisense/AD4/antisense_anticorellated_sig.csv", sep=""))

write.csv(pain_anticorellated, file = paste(PATH, "/antisense/ensembl_antisense/AD4/pain_anticorellated.csv", sep=""))

####################
# MHDF #
####################
system("mkdir -p /home/george/Desktop/LncRNA_v2/IPS/antisense/ensembl_antisense/NHDF")

pain <- ens_genes_NHDF[rownames(ens_genes_NHDF) %in% pain_genes$Gene.stable.ID.1,]
pain_sig <- subset(pain, pain$padj <.05 & !is.na(pain$padj))


sodium_channel <- ens_genes_NHDF[grep("SCN", ens_genes_NHDF$symbol),]
sodium_channel_sig <- subset(sodium_channel, sodium_channel$padj <.05 & !is.na(sodium_channel$padj))


potassium_channel <- ens_genes_NHDF[grep("KCN", ens_genes_NHDF$symbol),]
potassium_channel_sig <- subset(potassium_channel, potassium_channel$padj <.05 & !is.na(potassium_channel$padj))

chloride_channel <- ens_genes_NHDF[grep("CLC", ens_genes_NHDF$symbol),]
chloride_channel_sig <- subset(chloride_channel, chloride_channel$padj <.05 & !is.na(chloride_channel$padj))

calcium_channel <- ens_genes_NHDF[grep("CAC", ens_genes_NHDF$symbol),]
calcium_channel_sig <- subset(calcium_channel, calcium_channel$padj <.05 & !is.na(calcium_channel$padj))

trp_channel <- ens_genes_NHDF[grep("TRP", ens_genes_NHDF$symbol),]
trp_channel_sig <- subset(trp_channel, trp_channel$padj <.05 & !is.na(trp_channel$padj))

antisense_eset$lnc_lfc <- 0
antisense_eset$lnc_pvalue <- 0
antisense_eset$lnc_baseMean <- 0
antisense_eset$gene_lfc <- 0
antisense_eset$gene_pvalue <- 0



pb <- txtProgressBar(min = 0, max = nrow(antisense_eset), style = 3)
for (i in 1:nrow(antisense_eset)) {
	
	antisense_eset$gene_lfc[i] <- ens_genes_NHDF[as.character(antisense_eset$V16[i]),2]
	antisense_eset$gene_pvalue[i] <- ens_genes_NHDF[as.character(antisense_eset$V16[i]),6]
	
			
	antisense_eset$lnc_lfc[i] <- ens_genes_NHDF[as.character(antisense_eset$V4[i]),2]
	antisense_eset$lnc_pvalue[i] <- ens_genes_NHDF[as.character(antisense_eset$V4[i]),6]
	antisense_eset$lnc_baseMean[i] <- ens_genes_NHDF[as.character(antisense_eset$V4[i]),1] 
	setTxtProgressBar(pb, i)
	
	}
close(pb)
cat("\n")


both_significant <- subset(antisense_eset, antisense_eset$gene_pvalue<.05 & antisense_eset$lnc_pvalue<.05)

anticorellated <- subset(antisense_eset, (antisense_eset$lnc_lfc * antisense_eset$gene_lfc) < 0)

anticorellated_sig <- subset(both_significant, (both_significant$lnc_lfc * both_significant$gene_lfc) < 0) 

pain_antisense <- subset(antisense_eset, antisense_eset$sense_symbol %in% pain_genes_symbols)

pain_anticorellated <- subset(pain_antisense, (pain_antisense$lnc_lfc * pain_antisense$gene_lfc) < 0)

sodium_channel_antisense <- subset(antisense_eset, antisense_eset$sense_symbol %in% sodium_channel$symbol)
sodium_channel_both_sig <- subset(sodium_channel_antisense, sodium_channel_antisense$gene_pvalue<.05 & sodium_channel_antisense$lnc_pvalue<.05)
write.csv(sodium_channel_antisense, file = paste(PATH, "/antisense/ensembl_antisense/NHDF_sodium_channel_antisense.csv", sep=""))
write.csv(sodium_channel_both_sig, file = paste(PATH, "/antisense/ensembl_antisense/NHDF_sodium_channel_antisense_both_sig.csv", sep=""))

potassium_channel_antisense <- subset(antisense_eset, antisense_eset$sense_symbol %in% potassium_channel$symbol)
potassium_channel_both_sig <- subset(potassium_channel_antisense, potassium_channel_antisense$gene_pvalue<.05 & potassium_channel_antisense$lnc_pvalue<.05)
write.csv(potassium_channel_antisense, file = paste(PATH, "/antisense/ensembl_antisense/NHDF_potassium_channel_antisense.csv", sep=""))
write.csv(potassium_channel_both_sig, file = paste(PATH, "/antisense/ensembl_antisense/NHDF_potassium_channel_antisense_both_sig.csv", sep=""))

chloride_channel_antisense <- subset(antisense_eset, antisense_eset$sense_symbol %in% chloride_channel$symbol)
chloride_channel_both_sig <- subset(chloride_channel_antisense, chloride_channel_antisense$gene_pvalue<.05 & chloride_channel_antisense$lnc_pvalue<.05)
write.csv(chloride_channel_antisense, file = paste(PATH, "/antisense/ensembl_antisense/NHDF_chloride_channel_antisense.csv", sep=""))
write.csv(chloride_channel_both_sig, file = paste(PATH, "/antisense/ensembl_antisense/NHDF_chloride_channel_antisense_both_sig.csv", sep=""))

calcium_channel_antisense <- subset(antisense_eset, antisense_eset$sense_symbol %in% calcium_channel$symbol)
calcium_channel_both_sig <- subset(calcium_channel_antisense, calcium_channel_antisense$gene_pvalue<.05 & calcium_channel_antisense$lnc_pvalue<.05)
write.csv(calcium_channel_antisense, file = paste(PATH, "/antisense/ensembl_antisense/NHDF_calcium_channel_antisense.csv", sep=""))
write.csv(calcium_channel_both_sig, file = paste(PATH, "/antisense/ensembl_antisense/NHDF_calcium_channel_antisense_both_sig.csv", sep=""))

trp_channel_antisense <- subset(antisense_eset, antisense_eset$sense_symbol %in% trp_channel$symbol)
trp_channel_both_sig <- subset(trp_channel_antisense, trp_channel_antisense$gene_pvalue<.05 & trp_channel_antisense$lnc_pvalue<.05)
write.csv(trp_channel_antisense, file = paste(PATH, "/antisense/ensembl_antisense/NHDF_trp_channel_antisense.csv", sep=""))
write.csv(trp_channel_both_sig, file = paste(PATH, "/antisense/ensembl_antisense/NHDF_trp_channel_antisense_both_sig.csv", sep=""))

write.csv(antisense_eset, file = paste(PATH, "/antisense/ensembl_antisense/NHDF/antisense_all_expression.csv", sep=""))

write.csv(both_significant, file = paste(PATH, "/antisense/ensembl_antisense/NHDF/antisense_both_sig.csv", sep=""))

write.csv(pain_antisense, file = paste(PATH, "/antisense/ensembl_antisense/NHDF/pain_antisense.csv", sep=""))

write.csv(anticorellated, file = paste(PATH, "/antisense/ensembl_antisense/NHDF/antisense_anticorellated.csv", sep=""))

write.csv(anticorellated_sig, file = paste(PATH, "/antisense/ensembl_antisense/NHDF/antisense_anticorellated_sig.csv", sep=""))

write.csv(pain_anticorellated, file = paste(PATH, "/antisense/ensembl_antisense/NHDF/pain_anticorellated.csv", sep=""))

#################
# Novel LncRNAs #
#################
lncRNAs_gtf <- read.table(file = paste(PATH,"/myRegs_cov/Novel_lncRNAs_IPS.gtf", sep=""), sep="\t")

load(paste(PATH,"/LncRNAs_toc_IPS.RData", sep=""))

toc_lncs <- LncRNAs_toc

lncs_names <- toc_lncs$name

name_mapping <- cbind(rownames(toc_lncs), toc_lncs$name)

rownames(name_mapping) <- toc_lncs$name

antisense <- read.table(file = paste(PATH,"/antisense/antisense_lncs_sorted.bed", sep=""),sep="\t", header=FALSE)

antisense_refseq <- read.table(file = paste(PATH,"/antisense/antisense_lncs_refseq_sorted.bed", sep=""),sep="\t", header=FALSE)

intronic <- read.table(file = paste(PATH,"/antisense/intronic.bed", sep=""),sep="\t", header=FALSE)

#######################################################
antisense_mapping <- antisense[,c(4,16)]

intronic_mapping <- intronic[, c(4,16)]

antisense_mapping <- unique(antisense_mapping)

antisense_mapping_refseq <- antisense_refseq[,c(4,16)]

antisense_mapping_refseq <- unique(antisense_mapping_refseq)

######################################################
antisense_mapping$symbol <- mapIds(org.Hs.eg.db, keys=as.character(antisense_mapping[,2]), column="SYMBOL", keytype="ENSEMBL", multiVals="first")

intronic_mapping$symbol <- mapIds(org.Hs.eg.db, keys=as.character(intronic_mapping[,2]), column="SYMBOL", keytype="ENSEMBL", multiVals="first")

intronic_eset <- intronic_mapping

antisense_eset <- antisense_mapping

antisense_names <- name_mapping[name_mapping[,2] %in% antisense_eset$V4, ]

rownames(antisense_names) <- antisense_names[,2]

load(file = paste(PATH,"/pain_genes.RData", sep=""))

pain_genes <- unique(pain_genes[,c(2,1)])

pain_genes_symbols <- unique(as.character(pain_genes$HGNC.symbol))


AD2 <- read.csv(file = paste(PATH,"/DE/res_NEURONS_vs_IPS_AD2.csv", sep=""), row.name=1)
AD4 <- read.csv(file = paste(PATH,"/DE/res_NEURONS_vs_IPS_AD4.csv", sep=""), row.name=1)
NHDF <- read.csv(file = paste(PATH,"/DE/res_NEURONS_vs_IPS_NHDF.csv", sep=""), row.name=1)


ens_genes_AD2 <- AD2[grep("ENS", rownames(AD2)),]
ens_genes_AD4 <- AD4[grep("ENS", rownames(AD4)),]
ens_genes_NHDF <- NHDF[grep("ENS", rownames(NHDF)),]

lncs_AD2 <- AD2[grep("LncRNA", rownames(AD2)),]
lncs_AD4 <- AD4[grep("LncRNA", rownames(AD4)),]
lncs_NHDF <- NHDF[grep("LncRNA", rownames(NHDF)),]

################
# SNI vs Sham AD2 #
################
system("mkdir /home/george/Desktop/LncRNA_v2/IPS/antisense/AD2")

pain <- ens_genes_AD2[rownames(ens_genes_AD2) %in% pain_genes$Gene.stable.ID.1,]
pain_sig <- subset(pain, pain$padj <.05 & !is.na(pain$padj))
write.csv(pain, file = paste(PATH, "/AD2_pain.csv", sep=""))
write.csv(pain_sig, file = paste(PATH, "/AD2_pain_sig.csv", sep=""))

sodium_channel <- ens_genes_AD2[grep("SCN", ens_genes_AD2$symbol),]
sodium_channel_sig <- subset(sodium_channel, sodium_channel$padj <.05 & !is.na(sodium_channel$padj))
write.csv(pain, file = paste(PATH, "/AD2_sodium_channel.csv", sep=""))
write.csv(pain_sig, file = paste(PATH, "/AD2_sodium_channel_sig.csv", sep=""))

potassium_channel <- ens_genes_AD2[grep("KCN", ens_genes_AD2$symbol),]
potassium_channel_sig <- subset(potassium_channel, potassium_channel$padj <.05 & !is.na(potassium_channel$padj))
write.csv(potassium_channel, file = paste(PATH, "/AD2_potassium_channel.csv", sep=""))
write.csv(potassium_channel_sig, file = paste(PATH, "/AD2_potassium_channel_sig.csv", sep=""))

chloride_channel <- ens_genes_AD2[grep("CLC", ens_genes_AD2$symbol),]
chloride_channel_sig <- subset(chloride_channel, chloride_channel$padj <.05 & !is.na(chloride_channel$padj))
write.csv(chloride_channel, file = paste(PATH, "/AD2_chloride_channel.csv", sep=""))
write.csv(chloride_channel_sig, file = paste(PATH, "/AD2_chloride_channel_sig.csv", sep=""))

calcium_channel <- ens_genes_AD2[grep("CAC", ens_genes_AD2$symbol),]
calcium_channel_sig <- subset(calcium_channel, calcium_channel$padj <.05 & !is.na(calcium_channel$padj))
write.csv(calcium_channel, file = paste(PATH, "/AD2_calcium_channel.csv", sep=""))
write.csv(calcium_channel_sig, file = paste(PATH, "/AD2_calcium_channel_sig.csv", sep=""))

trp_channel <- ens_genes_AD2[grep("TRP", ens_genes_AD2$symbol),]
trp_channel_sig <- subset(trp_channel, trp_channel$padj <.05 & !is.na(trp_channel$padj))
write.csv(trp_channel, file = paste(PATH, "/AD2_trp_channel.csv", sep=""))
write.csv(trp_channel_sig, file = paste(PATH, "/AD2_trp_channel_sig.csv", sep=""))

antisense_eset$lnc_lfc <- 0
antisense_eset$lnc_pvalue <- 0
antisense_eset$lnc_baseMean <- 0
antisense_eset$gene_lfc <- 0
antisense_eset$gene_pvalue <- 0



pb <- txtProgressBar(min = 0, max = nrow(antisense_eset), style = 3)
for (i in 1:nrow(antisense_eset)) {
	if (grepl("^N", as.character(antisense_eset$V16[i]))) { 
	if (length(grep(as.character(antisense_eset$symbol[i]), ens_genes_AD2$symbol)) >0 ) {  	
	antisense_eset$gene_lfc[i] <- ens_genes_AD2[grep(as.character(antisense_eset$symbol[i]), ens_genes_AD2$symbol),2]
	antisense_eset$gene_pvalue[i] <- ens_genes_AD2[grep(as.character(antisense_eset$symbol[i]), ens_genes_AD2$symbol),6]}
	}
	else {	
	antisense_eset$gene_lfc[i] <- ens_genes_AD2[as.character(antisense_eset$V16[i]),2]
	antisense_eset$gene_pvalue[i] <- ens_genes_AD2[as.character(antisense_eset$V16[i]),6]
	}
	antisense_eset$ID[i] <- antisense_names[as.character(antisense_eset$V4)[i],1]
		
	antisense_eset$lnc_lfc[i] <- lncs_AD2[as.character(antisense_eset$ID[i]),2]
	antisense_eset$lnc_pvalue[i] <- lncs_AD2[as.character(antisense_eset$ID[i]),6]
	antisense_eset$lnc_baseMean[i] <- lncs_AD2[as.character(antisense_eset$ID[i]),1] 
	setTxtProgressBar(pb, i)
	
	}
close(pb)
cat("\n")


both_significant <- subset(antisense_eset, antisense_eset$gene_pvalue<.05 & antisense_eset$lnc_pvalue<.05)

anticorellated <- subset(antisense_eset, (antisense_eset$lnc_lfc * antisense_eset$gene_lfc) < 0)

anticorellated_sig <- subset(both_significant, (both_significant$lnc_lfc * both_significant$gene_lfc) < 0) 

pain_antisense <- subset(antisense_eset, antisense_eset$symbol %in% pain_genes_symbols)

pain_anticorellated <- subset(pain_antisense, (pain_antisense$lnc_lfc * pain_antisense$gene_lfc) < 0)


sodium_channel_antisense <- subset(antisense_eset, antisense_eset$symbol %in% sodium_channel$symbol)
sodium_channel_both_sig <- subset(sodium_channel_antisense, sodium_channel_antisense$gene_pvalue<.05 & sodium_channel_antisense$lnc_pvalue<.05)
write.csv(sodium_channel_antisense, file = paste(PATH, "/AD2_sodium_channel_antisense.csv", sep=""))
write.csv(sodium_channel_both_sig, file = paste(PATH, "/AD2_sodium_channel_antisense_both_sig.csv", sep=""))

potassium_channel_antisense <- subset(antisense_eset, antisense_eset$symbol %in% potassium_channel$symbol)
potassium_channel_both_sig <- subset(potassium_channel_antisense, potassium_channel_antisense$gene_pvalue<.05 & potassium_channel_antisense$lnc_pvalue<.05)
write.csv(potassium_channel_antisense, file = paste(PATH, "/AD2_potassium_channel_antisense.csv", sep=""))
write.csv(potassium_channel_both_sig, file = paste(PATH, "/AD2_potassium_channel_antisense_both_sig.csv", sep=""))

chloride_channel_antisense <- subset(antisense_eset, antisense_eset$symbol %in% chloride_channel$symbol)
chloride_channel_both_sig <- subset(chloride_channel_antisense, chloride_channel_antisense$gene_pvalue<.05 & chloride_channel_antisense$lnc_pvalue<.05)
write.csv(chloride_channel_antisense, file = paste(PATH, "/AD2_chloride_channel_antisense.csv", sep=""))
write.csv(chloride_channel_both_sig, file = paste(PATH, "/AD2_chloride_channel_antisense_both_sig.csv", sep=""))

calcium_channel_antisense <- subset(antisense_eset, antisense_eset$symbol %in% calcium_channel$symbol)
calcium_channel_both_sig <- subset(calcium_channel_antisense, calcium_channel_antisense$gene_pvalue<.05 & calcium_channel_antisense$lnc_pvalue<.05)
write.csv(calcium_channel_antisense, file = paste(PATH, "/AD2_calcium_channel_antisense.csv", sep=""))
write.csv(calcium_channel_both_sig, file = paste(PATH, "/AD2_calcium_channel_antisense_both_sig.csv", sep=""))

trp_channel_antisense <- subset(antisense_eset, antisense_eset$symbol %in% trp_channel$symbol)
trp_channel_both_sig <- subset(trp_channel_antisense, trp_channel_antisense$gene_pvalue<.05 & trp_channel_antisense$lnc_pvalue<.05)
write.csv(trp_channel_antisense, file = paste(PATH, "/AD2_trp_channel_antisense.csv", sep=""))
write.csv(trp_channel_both_sig, file = paste(PATH, "/AD2_trp_channel_antisense_both_sig.csv", sep=""))

write.csv(antisense_eset, file = paste(PATH, "/antisense/AD2/antisense_all_expression.csv", sep=""))

write.csv(both_significant, file = paste(PATH, "/antisense/AD2/antisense_both_sig.csv", sep=""))

write.csv(pain_antisense, file = paste(PATH, "/antisense/AD2/pain_antisense.csv", sep=""))

write.csv(anticorellated, file = paste(PATH, "/antisense/AD2/antisense_anticorellated.csv", sep=""))

write.csv(anticorellated_sig, file = paste(PATH, "/antisense/AD2/antisense_anticorellated_sig.csv", sep=""))

write.csv(pain_anticorellated, file = paste(PATH, "/antisense/AD2/pain_anticorellated.csv", sep=""))

#################################################
load(paste(PATH,"/DE/norm_counts.RData", sep=""))

rld_mat_lncs <-  rld_mat[grep("LncRNA", rownames(rld_mat)),]
rld_mat <-  rld_mat[grep("ENS", rownames(rld_mat)),]

colnames(rld_mat) <- paste("genes_", colnames(rld_mat), sep="")
colnames(rld_mat_lncs) <- paste("lncs_", colnames(rld_mat_lncs), sep="")

antisense_counts <- antisense_eset[as.character(antisense_eset$V16) %in% rownames(rld_mat),] 

antisense_counts <- antisense_counts[,c(1:3,9)]

genes_counts <- rld_mat[as.character(antisense_counts$V16),]

lncs_counts <- rld_mat_lncs[as.character(antisense_counts$ID),]

antisense_counts <- cbind(antisense_counts, lncs_counts, genes_counts)

for (i in 1:nrow(antisense_counts)) {
antisense_counts$cor[i] <- cor.test(x=as.numeric(antisense_counts[i,c(5:16)]), y=as.numeric(antisense_counts[i,c(17:28)]))$estimate  

antisense_counts$cor_pvalue[i] <- cor.test(x=as.numeric(antisense_counts[i,c(5:16)]), y=as.numeric(antisense_counts[i,c(17:28)]))$p.value 
}

antisense_counts$cor_pvalue <- p.adjust(antisense_counts$cor_pvalue, method = "BH")

save(file = paste(PATH, "/antisense/antisense_counts.RData", sep=""), antisense_counts)


#################################################
################
# AD4 #
################
system("mkdir /home/george/Desktop/LncRNA_v2/IPS/antisense/AD4")

pain <- ens_genes_AD4[rownames(ens_genes_AD4) %in% pain_genes$Gene.stable.ID.1,]
pain_sig <- subset(pain, pain$padj <.05 & !is.na(pain$padj))
write.csv(pain, file = paste(PATH, "/AD4_pain.csv", sep=""))
write.csv(pain_sig, file = paste(PATH, "/AD4_pain_sig.csv", sep=""))

sodium_channel <- ens_genes_AD4[grep("SCN", ens_genes_AD4$symbol),]
sodium_channel_sig <- subset(sodium_channel, sodium_channel$padj <.05 & !is.na(sodium_channel$padj))
write.csv(pain, file = paste(PATH, "/AD4_sodium_channel.csv", sep=""))
write.csv(pain_sig, file = paste(PATH, "/AD4_sodium_channel_sig.csv", sep=""))

potassium_channel <- ens_genes_AD4[grep("KCN", ens_genes_AD4$symbol),]
potassium_channel_sig <- subset(potassium_channel, potassium_channel$padj <.05 & !is.na(potassium_channel$padj))
write.csv(potassium_channel, file = paste(PATH, "/AD4_potassium_channel.csv", sep=""))
write.csv(potassium_channel_sig, file = paste(PATH, "/AD4_potassium_channel_sig.csv", sep=""))

chloride_channel <- ens_genes_AD4[grep("CLC", ens_genes_AD4$symbol),]
chloride_channel_sig <- subset(chloride_channel, chloride_channel$padj <.05 & !is.na(chloride_channel$padj))
write.csv(chloride_channel, file = paste(PATH, "/AD4_chloride_channel.csv", sep=""))
write.csv(chloride_channel_sig, file = paste(PATH, "/AD4_chloride_channel_sig.csv", sep=""))

calcium_channel <- ens_genes_AD4[grep("CAC", ens_genes_AD4$symbol),]
calcium_channel_sig <- subset(calcium_channel, calcium_channel$padj <.05 & !is.na(calcium_channel$padj))
write.csv(calcium_channel, file = paste(PATH, "/AD4_calcium_channel.csv", sep=""))
write.csv(calcium_channel_sig, file = paste(PATH, "/AD4_calcium_channel_sig.csv", sep=""))

trp_channel <- ens_genes_AD4[grep("TRP", ens_genes_AD4$symbol),]
trp_channel_sig <- subset(trp_channel, trp_channel$padj <.05 & !is.na(trp_channel$padj))
write.csv(trp_channel, file = paste(PATH, "/AD4_trp_channel.csv", sep=""))
write.csv(trp_channel_sig, file = paste(PATH, "/AD4_trp_channel_sig.csv", sep=""))

antisense_eset$lnc_lfc <- 0
antisense_eset$lnc_pvalue <- 0
antisense_eset$lnc_baseMean <- 0
antisense_eset$gene_lfc <- 0
antisense_eset$gene_pvalue <- 0



pb <- txtProgressBar(min = 0, max = nrow(antisense_eset), style = 3)
for (i in 1:nrow(antisense_eset)) {
	if (grepl("^N", as.character(antisense_eset$V16[i]))) { 
	if (length(grep(as.character(antisense_eset$symbol[i]), ens_genes_AD4$symbol)) >0 ) {  	
	antisense_eset$gene_lfc[i] <- ens_genes_AD4[grep(as.character(antisense_eset$symbol[i]), ens_genes_AD4$symbol),2]
	antisense_eset$gene_pvalue[i] <- ens_genes_AD4[grep(as.character(antisense_eset$symbol[i]), ens_genes_AD4$symbol),6]}
	}
	else {	
	antisense_eset$gene_lfc[i] <- ens_genes_AD4[as.character(antisense_eset$V16[i]),2]
	antisense_eset$gene_pvalue[i] <- ens_genes_AD4[as.character(antisense_eset$V16[i]),6]
	}
	antisense_eset$ID[i] <- antisense_names[as.character(antisense_eset$V4)[i],1]
		
	antisense_eset$lnc_lfc[i] <- lncs_AD4[as.character(antisense_eset$ID[i]),2]
	antisense_eset$lnc_pvalue[i] <- lncs_AD4[as.character(antisense_eset$ID[i]),6]
	antisense_eset$lnc_baseMean[i] <- lncs_AD4[as.character(antisense_eset$ID[i]),1] 
	setTxtProgressBar(pb, i)
	
	}
close(pb)
cat("\n")


both_significant <- subset(antisense_eset, antisense_eset$gene_pvalue<.05 & antisense_eset$lnc_pvalue<.05)

anticorellated <- subset(antisense_eset, (antisense_eset$lnc_lfc * antisense_eset$gene_lfc) < 0)

anticorellated_sig <- subset(both_significant, (both_significant$lnc_lfc * both_significant$gene_lfc) < 0) 

pain_antisense <- subset(antisense_eset, antisense_eset$symbol %in% pain_genes_symbols)

pain_anticorellated <- subset(pain_antisense, (pain_antisense$lnc_lfc * pain_antisense$gene_lfc) < 0)


sodium_channel_antisense <- subset(antisense_eset, antisense_eset$symbol %in% sodium_channel$symbol)
sodium_channel_both_sig <- subset(sodium_channel_antisense, sodium_channel_antisense$gene_pvalue<.05 & sodium_channel_antisense$lnc_pvalue<.05)
write.csv(sodium_channel_antisense, file = paste(PATH, "/AD4_sodium_channel_antisense.csv", sep=""))
write.csv(sodium_channel_both_sig, file = paste(PATH, "/AD4_sodium_channel_antisense_both_sig.csv", sep=""))

potassium_channel_antisense <- subset(antisense_eset, antisense_eset$symbol %in% potassium_channel$symbol)
potassium_channel_both_sig <- subset(potassium_channel_antisense, potassium_channel_antisense$gene_pvalue<.05 & potassium_channel_antisense$lnc_pvalue<.05)
write.csv(potassium_channel_antisense, file = paste(PATH, "/AD4_potassium_channel_antisense.csv", sep=""))
write.csv(potassium_channel_both_sig, file = paste(PATH, "/AD4_potassium_channel_antisense_both_sig.csv", sep=""))

chloride_channel_antisense <- subset(antisense_eset, antisense_eset$symbol %in% chloride_channel$symbol)
chloride_channel_both_sig <- subset(chloride_channel_antisense, chloride_channel_antisense$gene_pvalue<.05 & chloride_channel_antisense$lnc_pvalue<.05)
write.csv(chloride_channel_antisense, file = paste(PATH, "/AD4_chloride_channel_antisense.csv", sep=""))
write.csv(chloride_channel_both_sig, file = paste(PATH, "/AD4_chloride_channel_antisense_both_sig.csv", sep=""))

calcium_channel_antisense <- subset(antisense_eset, antisense_eset$symbol %in% calcium_channel$symbol)
calcium_channel_both_sig <- subset(calcium_channel_antisense, calcium_channel_antisense$gene_pvalue<.05 & calcium_channel_antisense$lnc_pvalue<.05)
write.csv(calcium_channel_antisense, file = paste(PATH, "/AD4_calcium_channel_antisense.csv", sep=""))
write.csv(calcium_channel_both_sig, file = paste(PATH, "/AD4_calcium_channel_antisense_both_sig.csv", sep=""))

trp_channel_antisense <- subset(antisense_eset, antisense_eset$symbol %in% trp_channel$symbol)
trp_channel_both_sig <- subset(trp_channel_antisense, trp_channel_antisense$gene_pvalue<.05 & trp_channel_antisense$lnc_pvalue<.05)
write.csv(trp_channel_antisense, file = paste(PATH, "/AD4_trp_channel_antisense.csv", sep=""))
write.csv(trp_channel_both_sig, file = paste(PATH, "/AD4_trp_channel_antisense_both_sig.csv", sep=""))

write.csv(antisense_eset, file = paste(PATH, "/antisense/AD4/antisense_all_expression.csv", sep=""))

write.csv(both_significant, file = paste(PATH, "/antisense/AD4/antisense_both_sig.csv", sep=""))

write.csv(pain_antisense, file = paste(PATH, "/antisense/AD4/pain_antisense.csv", sep=""))

write.csv(anticorellated, file = paste(PATH, "/antisense/AD4/antisense_anticorellated.csv", sep=""))

write.csv(anticorellated_sig, file = paste(PATH, "/antisense/AD4/antisense_anticorellated_sig.csv", sep=""))

write.csv(pain_anticorellated, file = paste(PATH, "/antisense/AD4/pain_anticorellated.csv", sep=""))

#################################################

################
# NHDF #
################
system("mkdir /home/george/Desktop/LncRNA_v2/IPS/antisense/NHDF")

pain <- ens_genes_NHDF[rownames(ens_genes_NHDF) %in% pain_genes$Gene.stable.ID.1,]
pain_sig <- subset(pain, pain$padj <.05 & !is.na(pain$padj))
write.csv(pain, file = paste(PATH, "/NHDF_pain.csv", sep=""))
write.csv(pain_sig, file = paste(PATH, "/NHDF_pain_sig.csv", sep=""))

sodium_channel <- ens_genes_NHDF[grep("SCN", ens_genes_NHDF$symbol),]
sodium_channel_sig <- subset(sodium_channel, sodium_channel$padj <.05 & !is.na(sodium_channel$padj))
write.csv(pain, file = paste(PATH, "/NHDF_sodium_channel.csv", sep=""))
write.csv(pain_sig, file = paste(PATH, "/NHDF_sodium_channel_sig.csv", sep=""))

potassium_channel <- ens_genes_NHDF[grep("KCN", ens_genes_NHDF$symbol),]
potassium_channel_sig <- subset(potassium_channel, potassium_channel$padj <.05 & !is.na(potassium_channel$padj))
write.csv(potassium_channel, file = paste(PATH, "/NHDF_potassium_channel.csv", sep=""))
write.csv(potassium_channel_sig, file = paste(PATH, "/NHDF_potassium_channel_sig.csv", sep=""))

chloride_channel <- ens_genes_NHDF[grep("CLC", ens_genes_NHDF$symbol),]
chloride_channel_sig <- subset(chloride_channel, chloride_channel$padj <.05 & !is.na(chloride_channel$padj))
write.csv(chloride_channel, file = paste(PATH, "/NHDF_chloride_channel.csv", sep=""))
write.csv(chloride_channel_sig, file = paste(PATH, "/NHDF_chloride_channel_sig.csv", sep=""))

calcium_channel <- ens_genes_NHDF[grep("CAC", ens_genes_NHDF$symbol),]
calcium_channel_sig <- subset(calcium_channel, calcium_channel$padj <.05 & !is.na(calcium_channel$padj))
write.csv(calcium_channel, file = paste(PATH, "/NHDF_calcium_channel.csv", sep=""))
write.csv(calcium_channel_sig, file = paste(PATH, "/NHDF_calcium_channel_sig.csv", sep=""))

trp_channel <- ens_genes_NHDF[grep("TRP", ens_genes_NHDF$symbol),]
trp_channel_sig <- subset(trp_channel, trp_channel$padj <.05 & !is.na(trp_channel$padj))
write.csv(trp_channel, file = paste(PATH, "/NHDF_trp_channel.csv", sep=""))
write.csv(trp_channel_sig, file = paste(PATH, "/NHDF_trp_channel_sig.csv", sep=""))

antisense_eset$lnc_lfc <- 0
antisense_eset$lnc_pvalue <- 0
antisense_eset$lnc_baseMean <- 0
antisense_eset$gene_lfc <- 0
antisense_eset$gene_pvalue <- 0



pb <- txtProgressBar(min = 0, max = nrow(antisense_eset), style = 3)
for (i in 1:nrow(antisense_eset)) {
	if (grepl("^N", as.character(antisense_eset$V16[i]))) { 
	if (length(grep(as.character(antisense_eset$symbol[i]), ens_genes_NHDF$symbol)) >0 ) {  	
	antisense_eset$gene_lfc[i] <- ens_genes_NHDF[grep(as.character(antisense_eset$symbol[i]), ens_genes_NHDF$symbol),2]
	antisense_eset$gene_pvalue[i] <- ens_genes_NHDF[grep(as.character(antisense_eset$symbol[i]), ens_genes_NHDF$symbol),6]}
	}
	else {	
	antisense_eset$gene_lfc[i] <- ens_genes_NHDF[as.character(antisense_eset$V16[i]),2]
	antisense_eset$gene_pvalue[i] <- ens_genes_NHDF[as.character(antisense_eset$V16[i]),6]
	}
	antisense_eset$ID[i] <- antisense_names[as.character(antisense_eset$V4)[i],1]
		
	antisense_eset$lnc_lfc[i] <- lncs_NHDF[as.character(antisense_eset$ID[i]),2]
	antisense_eset$lnc_pvalue[i] <- lncs_NHDF[as.character(antisense_eset$ID[i]),6]
	antisense_eset$lnc_baseMean[i] <- lncs_NHDF[as.character(antisense_eset$ID[i]),1] 
	setTxtProgressBar(pb, i)
	
	}
close(pb)
cat("\n")


both_significant <- subset(antisense_eset, antisense_eset$gene_pvalue<.05 & antisense_eset$lnc_pvalue<.05)

anticorellated <- subset(antisense_eset, (antisense_eset$lnc_lfc * antisense_eset$gene_lfc) < 0)

anticorellated_sig <- subset(both_significant, (both_significant$lnc_lfc * both_significant$gene_lfc) < 0) 

pain_antisense <- subset(antisense_eset, antisense_eset$symbol %in% pain_genes_symbols)

pain_anticorellated <- subset(pain_antisense, (pain_antisense$lnc_lfc * pain_antisense$gene_lfc) < 0)


sodium_channel_antisense <- subset(antisense_eset, antisense_eset$symbol %in% sodium_channel$symbol)
sodium_channel_both_sig <- subset(sodium_channel_antisense, sodium_channel_antisense$gene_pvalue<.05 & sodium_channel_antisense$lnc_pvalue<.05)
write.csv(sodium_channel_antisense, file = paste(PATH, "/NHDF_sodium_channel_antisense.csv", sep=""))
write.csv(sodium_channel_both_sig, file = paste(PATH, "/NHDF_sodium_channel_antisense_both_sig.csv", sep=""))

potassium_channel_antisense <- subset(antisense_eset, antisense_eset$symbol %in% potassium_channel$symbol)
potassium_channel_both_sig <- subset(potassium_channel_antisense, potassium_channel_antisense$gene_pvalue<.05 & potassium_channel_antisense$lnc_pvalue<.05)
write.csv(potassium_channel_antisense, file = paste(PATH, "/NHDF_potassium_channel_antisense.csv", sep=""))
write.csv(potassium_channel_both_sig, file = paste(PATH, "/NHDF_potassium_channel_antisense_both_sig.csv", sep=""))

chloride_channel_antisense <- subset(antisense_eset, antisense_eset$symbol %in% chloride_channel$symbol)
chloride_channel_both_sig <- subset(chloride_channel_antisense, chloride_channel_antisense$gene_pvalue<.05 & chloride_channel_antisense$lnc_pvalue<.05)
write.csv(chloride_channel_antisense, file = paste(PATH, "/NHDF_chloride_channel_antisense.csv", sep=""))
write.csv(chloride_channel_both_sig, file = paste(PATH, "/NHDF_chloride_channel_antisense_both_sig.csv", sep=""))

calcium_channel_antisense <- subset(antisense_eset, antisense_eset$symbol %in% calcium_channel$symbol)
calcium_channel_both_sig <- subset(calcium_channel_antisense, calcium_channel_antisense$gene_pvalue<.05 & calcium_channel_antisense$lnc_pvalue<.05)
write.csv(calcium_channel_antisense, file = paste(PATH, "/NHDF_calcium_channel_antisense.csv", sep=""))
write.csv(calcium_channel_both_sig, file = paste(PATH, "/NHDF_calcium_channel_antisense_both_sig.csv", sep=""))

trp_channel_antisense <- subset(antisense_eset, antisense_eset$symbol %in% trp_channel$symbol)
trp_channel_both_sig <- subset(trp_channel_antisense, trp_channel_antisense$gene_pvalue<.05 & trp_channel_antisense$lnc_pvalue<.05)
write.csv(trp_channel_antisense, file = paste(PATH, "/NHDF_trp_channel_antisense.csv", sep=""))
write.csv(trp_channel_both_sig, file = paste(PATH, "/NHDF_trp_channel_antisense_both_sig.csv", sep=""))

write.csv(antisense_eset, file = paste(PATH, "/antisense/NHDF/antisense_all_expression.csv", sep=""))

write.csv(both_significant, file = paste(PATH, "/antisense/NHDF/antisense_both_sig.csv", sep=""))

write.csv(pain_antisense, file = paste(PATH, "/antisense/NHDF/pain_antisense.csv", sep=""))

write.csv(anticorellated, file = paste(PATH, "/antisense/NHDF/antisense_anticorellated.csv", sep=""))

write.csv(anticorellated_sig, file = paste(PATH, "/antisense/NHDF/antisense_anticorellated_sig.csv", sep=""))

write.csv(pain_anticorellated, file = paste(PATH, "/antisense/NHDF/pain_anticorellated.csv", sep=""))












##################################
########### Intronic #############
##################################
################
# AD2 #
################
system("mkdir /home/george/Desktop/LncRNA_v2/IPS/intronic")
system("mkdir /home/george/Desktop/LncRNA_v2/IPS/intronic/AD2")

intronic_names <- name_mapping[name_mapping[,2] %in% intronic_eset$V4, ]

rownames(intronic_names) <- intronic_names[,2]


pain <- ens_genes_AD2[rownames(ens_genes_AD2) %in% pain_genes$Gene.stable.ID.1,]
pain_sig <- subset(pain, pain$padj <.05 & !is.na(pain$padj))

intronic_eset$lnc_lfc <- 0
intronic_eset$lnc_pvalue <- 0
intronic_eset$lnc_baseMean <- 0
intronic_eset$gene_lfc <- 0
intronic_eset$gene_pvalue <- 0



pb <- txtProgressBar(min = 0, max = nrow(intronic_eset), style = 3)
for (i in 1:nrow(intronic_eset)) {
	if (grepl("^N", as.character(intronic_eset$V16[i]))) { 
	if (length(grep(as.character(intronic_eset$symbol[i]), ens_genes_AD2$symbol)) >0 ) {  	
	intronic_eset$gene_lfc[i] <- ens_genes_AD2[grep(as.character(intronic_eset$symbol[i]), ens_genes_AD2$symbol),2]
	intronic_eset$gene_pvalue[i] <- ens_genes_AD2[grep(as.character(intronic_eset$symbol[i]), ens_genes_AD2$symbol),6]}
	}
	else {	
	intronic_eset$gene_lfc[i] <- ens_genes_AD2[as.character(intronic_eset$V16[i]),2]
	intronic_eset$gene_pvalue[i] <- ens_genes_AD2[as.character(intronic_eset$V16[i]),6]
	}
	intronic_eset$ID[i] <- intronic_names[as.character(intronic_eset$V4)[i],1]
		
	intronic_eset$lnc_lfc[i] <- lncs_AD2[as.character(intronic_eset$ID[i]),2]
	intronic_eset$lnc_pvalue[i] <- lncs_AD2[as.character(intronic_eset$ID[i]),6]
	intronic_eset$lnc_baseMean[i] <- lncs_AD2[as.character(intronic_eset$ID[i]),1] 
	setTxtProgressBar(pb, i)
	
	}
close(pb)
cat("\n")


both_significant <- subset(intronic_eset, intronic_eset$gene_pvalue<.05 & intronic_eset$lnc_pvalue<.05)

anticorellated <- subset(intronic_eset, (intronic_eset$lnc_lfc * intronic_eset$gene_lfc) < 0)

anticorellated_sig <- subset(both_significant, (both_significant$lnc_lfc * both_significant$gene_lfc) < 0) 

pain_intronic <- subset(intronic_eset, intronic_eset$symbol %in% pain_genes_symbols)

pain_anticorellated <- subset(pain_intronic, (pain_intronic$lnc_lfc * pain_intronic$gene_lfc) < 0)


write.csv(intronic_eset, file = paste(PATH, "/intronic/AD2/intronic_all_expression.csv", sep=""))

write.csv(both_significant, file = paste(PATH, "/intronic/AD2/intronic_both_sig.csv", sep=""))

write.csv(pain_intronic, file = paste(PATH, "/intronic/AD2/pain_intronic.csv", sep=""))

write.csv(anticorellated, file = paste(PATH, "/intronic/AD2/intronic_anticorellated.csv", sep=""))

write.csv(anticorellated_sig, file = paste(PATH, "/intronic/AD2/intronic_anticorellated_sig.csv", sep=""))

write.csv(pain_anticorellated, file = paste(PATH, "/intronic/AD2/pain_anticorellated.csv", sep=""))

#################################################
load(paste(PATH,"/DE/norm_counts.RData", sep=""))

rld_mat_lncs <-  rld_mat[grep("LncRNA", rownames(rld_mat)),]
rld_mat <-  rld_mat[grep("ENS", rownames(rld_mat)),]

colnames(rld_mat) <- paste("genes_", colnames(rld_mat), sep="")
colnames(rld_mat_lncs) <- paste("lncs_", colnames(rld_mat_lncs), sep="")

intronic_counts <- intronic_eset[as.character(intronic_eset$V16) %in% rownames(rld_mat),] 

intronic_counts <- intronic_counts[,c(1:3,9)]

genes_counts <- rld_mat[as.character(intronic_counts$V16),]

lncs_counts <- rld_mat_lncs[as.character(intronic_counts$ID),]

intronic_counts <- cbind(intronic_counts, lncs_counts, genes_counts)

for (i in 1:nrow(intronic_counts)) {
intronic_counts$cor[i] <- cor.test(x=as.numeric(intronic_counts[i,c(5:16)]), y=as.numeric(intronic_counts[i,c(17:28)]))$estimate  

intronic_counts$cor_pvalue[i] <- cor.test(x=as.numeric(intronic_counts[i,c(5:16)]), y=as.numeric(intronic_counts[i,c(17:28)]))$p.value 
}

intronic_counts$cor_pvalue <- p.adjust(intronic_counts$cor_pvalue, method = "BH")

save(file = paste(PATH, "/intronic/intronic_counts.RData", sep=""), intronic_counts)


#################################################
################
# AD4 #
################
system("mkdir /home/george/Desktop/LncRNA_v2/IPS/intronic/AD4")

pain <- ens_genes_AD4[rownames(ens_genes_AD4) %in% pain_genes$Gene.stable.ID.1,]
pain_sig <- subset(pain, pain$padj <.05 & !is.na(pain$padj))


intronic_eset$lnc_lfc <- 0
intronic_eset$lnc_pvalue <- 0
intronic_eset$lnc_baseMean <- 0
intronic_eset$gene_lfc <- 0
intronic_eset$gene_pvalue <- 0



pb <- txtProgressBar(min = 0, max = nrow(intronic_eset), style = 3)
for (i in 1:nrow(intronic_eset)) {
	if (grepl("^N", as.character(intronic_eset$V16[i]))) { 
	if (length(grep(as.character(intronic_eset$symbol[i]), ens_genes_AD4$symbol)) >0 ) {  	
	intronic_eset$gene_lfc[i] <- ens_genes_AD4[grep(as.character(intronic_eset$symbol[i]), ens_genes_AD4$symbol),2]
	intronic_eset$gene_pvalue[i] <- ens_genes_AD4[grep(as.character(intronic_eset$symbol[i]), ens_genes_AD4$symbol),6]}
	}
	else {	
	intronic_eset$gene_lfc[i] <- ens_genes_AD4[as.character(intronic_eset$V16[i]),2]
	intronic_eset$gene_pvalue[i] <- ens_genes_AD4[as.character(intronic_eset$V16[i]),6]
	}
	intronic_eset$ID[i] <- intronic_names[as.character(intronic_eset$V4)[i],1]
		
	intronic_eset$lnc_lfc[i] <- lncs_AD4[as.character(intronic_eset$ID[i]),2]
	intronic_eset$lnc_pvalue[i] <- lncs_AD4[as.character(intronic_eset$ID[i]),6]
	intronic_eset$lnc_baseMean[i] <- lncs_AD4[as.character(intronic_eset$ID[i]),1] 
	setTxtProgressBar(pb, i)
	
	}
close(pb)
cat("\n")


both_significant <- subset(intronic_eset, intronic_eset$gene_pvalue<.05 & intronic_eset$lnc_pvalue<.05)

anticorellated <- subset(intronic_eset, (intronic_eset$lnc_lfc * intronic_eset$gene_lfc) < 0)

anticorellated_sig <- subset(both_significant, (both_significant$lnc_lfc * both_significant$gene_lfc) < 0) 

pain_intronic <- subset(intronic_eset, intronic_eset$symbol %in% pain_genes_symbols)

pain_anticorellated <- subset(pain_intronic, (pain_intronic$lnc_lfc * pain_intronic$gene_lfc) < 0)


write.csv(intronic_eset, file = paste(PATH, "/intronic/AD4/intronic_all_expression.csv", sep=""))

write.csv(both_significant, file = paste(PATH, "/intronic/AD4/intronic_both_sig.csv", sep=""))

write.csv(pain_intronic, file = paste(PATH, "/intronic/AD4/pain_intronic.csv", sep=""))

write.csv(anticorellated, file = paste(PATH, "/intronic/AD4/intronic_anticorellated.csv", sep=""))

write.csv(anticorellated_sig, file = paste(PATH, "/intronic/AD4/intronic_anticorellated_sig.csv", sep=""))

write.csv(pain_anticorellated, file = paste(PATH, "/intronic/AD4/pain_anticorellated.csv", sep=""))

#################################################

################
# NHDF #
################
system("mkdir /home/george/Desktop/LncRNA_v2/IPS/intronic/NHDF")

pain <- ens_genes_NHDF[rownames(ens_genes_NHDF) %in% pain_genes$Gene.stable.ID.1,]
pain_sig <- subset(pain, pain$padj <.05 & !is.na(pain$padj))

intronic_eset$lnc_lfc <- 0
intronic_eset$lnc_pvalue <- 0
intronic_eset$lnc_baseMean <- 0
intronic_eset$gene_lfc <- 0
intronic_eset$gene_pvalue <- 0



pb <- txtProgressBar(min = 0, max = nrow(intronic_eset), style = 3)
for (i in 1:nrow(intronic_eset)) {
	if (grepl("^N", as.character(intronic_eset$V16[i]))) { 
	if (length(grep(as.character(intronic_eset$symbol[i]), ens_genes_NHDF$symbol)) >0 ) {  	
	intronic_eset$gene_lfc[i] <- ens_genes_NHDF[grep(as.character(intronic_eset$symbol[i]), ens_genes_NHDF$symbol),2]
	intronic_eset$gene_pvalue[i] <- ens_genes_NHDF[grep(as.character(intronic_eset$symbol[i]), ens_genes_NHDF$symbol),6]}
	}
	else {	
	intronic_eset$gene_lfc[i] <- ens_genes_NHDF[as.character(intronic_eset$V16[i]),2]
	intronic_eset$gene_pvalue[i] <- ens_genes_NHDF[as.character(intronic_eset$V16[i]),6]
	}
	intronic_eset$ID[i] <- intronic_names[as.character(intronic_eset$V4)[i],1]
		
	intronic_eset$lnc_lfc[i] <- lncs_NHDF[as.character(intronic_eset$ID[i]),2]
	intronic_eset$lnc_pvalue[i] <- lncs_NHDF[as.character(intronic_eset$ID[i]),6]
	intronic_eset$lnc_baseMean[i] <- lncs_NHDF[as.character(intronic_eset$ID[i]),1] 
	setTxtProgressBar(pb, i)
	
	}
close(pb)
cat("\n")


both_significant <- subset(intronic_eset, intronic_eset$gene_pvalue<.05 & intronic_eset$lnc_pvalue<.05)

anticorellated <- subset(intronic_eset, (intronic_eset$lnc_lfc * intronic_eset$gene_lfc) < 0)

anticorellated_sig <- subset(both_significant, (both_significant$lnc_lfc * both_significant$gene_lfc) < 0) 

pain_intronic <- subset(intronic_eset, intronic_eset$symbol %in% pain_genes_symbols)

pain_anticorellated <- subset(pain_intronic, (pain_intronic$lnc_lfc * pain_intronic$gene_lfc) < 0)


write.csv(intronic_eset, file = paste(PATH, "/intronic/NHDF/intronic_all_expression.csv", sep=""))

write.csv(both_significant, file = paste(PATH, "/intronic/NHDF/intronic_both_sig.csv", sep=""))

write.csv(pain_intronic, file = paste(PATH, "/intronic/NHDF/pain_intronic.csv", sep=""))

write.csv(anticorellated, file = paste(PATH, "/intronic/NHDF/intronic_anticorellated.csv", sep=""))

write.csv(anticorellated_sig, file = paste(PATH, "/intronic/NHDF/intronic_anticorellated_sig.csv", sep=""))

write.csv(pain_anticorellated, file = paste(PATH, "/intronic/NHDF/pain_anticorellated.csv", sep=""))




