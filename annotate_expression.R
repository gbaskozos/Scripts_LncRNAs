library("org.Mm.eg.db")
library("data.table")
library(GenomicFeatures)
library("rtracklayer")
library(GenomicAlignments)
library(BiocParallel)
library(GenomicRanges)
library(gplots)
library(biomaRt)

options(stringsAsFactors=FALSE)
PATH <- "/home/george/Desktop/LncRNA_v2/Mouse"

load(file = paste(PATH,"/toc_all.RData", sep=""))

load(file = paste(PATH,"/DE/norm_counts.RData", sep=""))

balbc <- read.csv(file = paste(PATH,"/DE/res_sni_vs_sham_BALB.c.csv", sep=""), row.name=1)
b10.d2 <- read.csv(file = paste(PATH,"/DE/res_sni_vs_sham_B10.D2.csv", sep=""), row.name=1)
diff_sni_b10d2_vs_balbc <- read.csv(file = paste(PATH,"/DE/res_sni_diff_B10.D2_vs_BALB.c.csv", sep=""), row.name=1)

#################
# Known LncRNAs #
#################
#########################################
########### ENSEMBL LncRNAs #############
#########################################

nc_ensembl <- read.csv(file = paste(PATH,"/mm10_genome/ens_lncs_annotation.csv", sep=""), row.name=1)

rownames(nc_ensembl) <- nc_ensembl$ensembl_gene_id

ens_lncs_names <- rownames(nc_ensembl)

cutoff <- 0
## Apply cut off to strain
toc_coff <- toc[apply(toc, 1, function(x) sum(x[grep("BALB.c", names(toc))] > cutoff) == length(grep("BALB.c", names(toc))) | sum(x[grep("B10.D2", names(toc))] > cutoff) == length(grep("B10.D2", names(toc)))), ]

expressed_nc_ensembl <- toc_coff[rownames(toc_coff) %in% as.character(nc_ensembl$ensembl_gene_id), ]

save(file = paste(PATH,"/DE/expressed_nc_ensembl.RData", sep=""), expressed_nc_ensembl)
write.csv(file = paste(PATH,"/DE/expressed_nc_ensembl_toc.csv", sep=""), expressed_nc_ensembl)

nc_rld <- rld_mat[rownames(rld) %in% rownames(expressed_nc_ensembl),] 

nc_ensembl_rld <- rld_mat[rownames(rld_mat) %in% rownames(expressed_nc_ensembl),] 

sig_nc_BALB.c <- balbc[rownames(balbc) %in% rownames(expressed_nc_ensembl) & balbc$padj < 0.05 & !is.na(balbc$padj), ]
sig_nc_B10.D2 <- b10.d2[rownames(b10.d2) %in% rownames(expressed_nc_ensembl) & b10.d2$padj < 0.05 & !is.na(b10.d2$padj), ]
sig_diff_sni_b10d2_vs_balbc <- diff_sni_b10d2_vs_balbc[rownames(diff_sni_b10d2_vs_balbc) %in% rownames(expressed_nc_ensembl) & diff_sni_b10d2_vs_balbc$padj < 0.05 & !is.na(diff_sni_b10d2_vs_balbc$padj), ]

#Order by LFC, base mean and - padj
sig_nc_BALB.c <- sig_nc_BALB.c[order(abs(sig_nc_BALB.c$log2FoldChange), sig_nc_BALB.c$baseMean, -sig_nc_BALB.c$padj),]
sig_nc_B10.D2 <- sig_nc_B10.D2[order(abs(sig_nc_B10.D2$log2FoldChange), sig_nc_B10.D2$baseMean, -sig_nc_B10.D2$padj),]
sig_diff_sni_b10d2_vs_balbc <- sig_diff_sni_b10d2_vs_balbc[order(abs(sig_diff_sni_b10d2_vs_balbc$log2FoldChange), sig_diff_sni_b10d2_vs_balbc$baseMean, -sig_diff_sni_b10d2_vs_balbc$padj),]

sig_nc_BALB.c$description <- nc_ensembl[rownames(sig_nc_BALB.c),]$description
sig_nc_B10.D2$description <- nc_ensembl[rownames(sig_nc_B10.D2),]$description
sig_diff_sni_b10d2_vs_balbc$description <- nc_ensembl[rownames(sig_diff_sni_b10d2_vs_balbc),]$description

write.csv(file = paste(PATH,"/DE/sig_nc_BALB.c.csv", sep=""),sig_nc_BALB.c)
write.csv(file = paste(PATH,"/DE/sig_nc_B10.D2.csv", sep=""),sig_nc_B10.D2)
write.csv(file = paste(PATH,"/DE/sig_diff_sni_b10d2_vs_balbc.csv", sep=""),sig_diff_sni_b10d2_vs_balbc)


pdf(file = paste(PATH,"/DE/expressed_lncs_heatmap.pdf", sep=""))
cluster <- heatmap.2(as.matrix(nc_ensembl_rld), col=redgreen(200), scale="row", trace="none", cexCol=1, cexRow=0.4, margins=c(5,7), key=FALSE, Colv=TRUE, Rowv=TRUE, main="ENSEMBL LncRNAs", dendrogram="column", symkey=TRUE, srtCol=42, adjCol=c(1,1), lmat=rbind(c(0,3), c(0,1), c(2,4)), lhei=c(2,7,1), lwid = c(1,7),labRow=nc_ensembl[rownames(nc_ensembl_rld),]$mgi_symbol) 
dev.off()

#########################################################################################################
antisense_ENSEMBL <- read.table( file = paste(PATH, "/antisense/lncs_antisense_ENSEMBL.bed", sep=""))

lincs_ENSEMBL <- read.table(file = paste(PATH, "/non_antisense/lncs_intergenic_closest_ENSEMBL.bed", sep=""))

#######################################################
antisense_mapping <- antisense_ENSEMBL[,c(4,16)]

antisense_mapping$antisense_symbol <- mapIds(org.Mm.eg.db, keys=as.character(antisense_mapping[,1]), column="SYMBOL", keytype="ENSEMBL", multiVals="first")

antisense_mapping$sense_symbol <- mapIds(org.Mm.eg.db, keys=as.character(antisense_mapping[,2]), column="SYMBOL", keytype="ENSEMBL", multiVals="first")

antisense_eset <- antisense_mapping

pain_genes <- read.csv(file = paste(PATH,"/../pain_genes_homologs.csv", sep=""))

pain_genes <- unique(pain_genes[,c(2:3)])

pain_genes_symbols <- unique(as.character(pain_genes$MGI.symbol))


BALB.c <- read.csv(file = paste(PATH,"/DE/res_sni_vs_sham_BALB.c.csv", sep=""), row.name=1)
B10.D2 <- read.csv(file = paste(PATH,"/DE/res_sni_vs_sham_B10.D2.csv", sep=""), row.name=1)
Diff_B10.D2_vs_BALB.c <- read.csv(file = paste(PATH,"/DE/res_sni_diff_B10.D2_vs_BALB.c.csv", sep=""), row.name=1)


ens_genes_sham_vs_sni_BALB.c <- BALB.c[grep("ENS", rownames(BALB.c)),]
ens_genes_sham_vs_sni_B10.D2 <- B10.D2[grep("ENS", rownames(B10.D2)),]
ens_genes_sni_diff_B10.D2_vs_BALB.c <- Diff_B10.D2_vs_BALB.c[grep("ENS", rownames(Diff_B10.D2_vs_BALB.c)),]

######################################################################################################
# BALB/c
########
system("mkdir -p /home/george/Desktop/LncRNA_v2/Mouse/antisense/ensembl_antisense/sni_vs_sham_balb.c")

pain <- ens_genes_sham_vs_sni_BALB.c[rownames(ens_genes_sham_vs_sni_BALB.c) %in% pain_genes$Gene.ID,]
pain_sig <- subset(pain, pain$padj <.05 & !is.na(pain$padj))


sodium_channel <- ens_genes_sham_vs_sni_BALB.c[grep("Scn", ens_genes_sham_vs_sni_BALB.c$symbol),]
sodium_channel_sig <- subset(sodium_channel, sodium_channel$padj <.05 & !is.na(sodium_channel$padj))


potassium_channel <- ens_genes_sham_vs_sni_BALB.c[grep("Kcn", ens_genes_sham_vs_sni_BALB.c$symbol),]
potassium_channel_sig <- subset(potassium_channel, potassium_channel$padj <.05 & !is.na(potassium_channel$padj))

chloride_channel <- ens_genes_sham_vs_sni_BALB.c[grep("Clc", ens_genes_sham_vs_sni_BALB.c$symbol),]
chloride_channel_sig <- subset(chloride_channel, chloride_channel$padj <.05 & !is.na(chloride_channel$padj))

calcium_channel <- ens_genes_sham_vs_sni_BALB.c[grep("Cac", ens_genes_sham_vs_sni_BALB.c$symbol),]
calcium_channel_sig <- subset(calcium_channel, calcium_channel$padj <.05 & !is.na(calcium_channel$padj))

trp_channel <- ens_genes_sham_vs_sni_BALB.c[grep("Trp", ens_genes_sham_vs_sni_BALB.c$symbol),]
trp_channel_sig <- subset(trp_channel, trp_channel$padj <.05 & !is.na(trp_channel$padj))

antisense_eset$lnc_lfc <- 0
antisense_eset$lnc_pvalue <- 0
antisense_eset$lnc_baseMean <- 0
antisense_eset$gene_lfc <- 0
antisense_eset$gene_pvalue <- 0



pb <- txtProgressBar(min = 0, max = nrow(antisense_eset), style = 3)
for (i in 1:nrow(antisense_eset)) {
	
	antisense_eset$gene_lfc[i] <- ens_genes_sham_vs_sni_BALB.c[as.character(antisense_eset$V16[i]),2]
	antisense_eset$gene_pvalue[i] <- ens_genes_sham_vs_sni_BALB.c[as.character(antisense_eset$V16[i]),6]
	
			
	antisense_eset$lnc_lfc[i] <- ens_genes_sham_vs_sni_BALB.c[as.character(antisense_eset$V4[i]),2]
	antisense_eset$lnc_pvalue[i] <- ens_genes_sham_vs_sni_BALB.c[as.character(antisense_eset$V4[i]),6]
	antisense_eset$lnc_baseMean[i] <- ens_genes_sham_vs_sni_BALB.c[as.character(antisense_eset$V4[i]),1] 
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
write.csv(sodium_channel_antisense, file = paste(PATH, "/antisense/ensembl_antisense/balbc_sodium_channel_antisense.csv", sep=""))
write.csv(sodium_channel_both_sig, file = paste(PATH, "/antisense/ensembl_antisense/balbc_sodium_channel_antisense_both_sig.csv", sep=""))

potassium_channel_antisense <- subset(antisense_eset, antisense_eset$sense_symbol %in% potassium_channel$symbol)
potassium_channel_both_sig <- subset(potassium_channel_antisense, potassium_channel_antisense$gene_pvalue<.05 & potassium_channel_antisense$lnc_pvalue<.05)
write.csv(potassium_channel_antisense, file = paste(PATH, "/antisense/ensembl_antisense/balbc_potassium_channel_antisense.csv", sep=""))
write.csv(potassium_channel_both_sig, file = paste(PATH, "/antisense/ensembl_antisense/balbc_potassium_channel_antisense_both_sig.csv", sep=""))

chloride_channel_antisense <- subset(antisense_eset, antisense_eset$sense_symbol %in% chloride_channel$symbol)
chloride_channel_both_sig <- subset(chloride_channel_antisense, chloride_channel_antisense$gene_pvalue<.05 & chloride_channel_antisense$lnc_pvalue<.05)
write.csv(chloride_channel_antisense, file = paste(PATH, "/antisense/ensembl_antisense/balbc_chloride_channel_antisense.csv", sep=""))
write.csv(chloride_channel_both_sig, file = paste(PATH, "/antisense/ensembl_antisense/balbc_chloride_channel_antisense_both_sig.csv", sep=""))

calcium_channel_antisense <- subset(antisense_eset, antisense_eset$sense_symbol %in% calcium_channel$symbol)
calcium_channel_both_sig <- subset(calcium_channel_antisense, calcium_channel_antisense$gene_pvalue<.05 & calcium_channel_antisense$lnc_pvalue<.05)
write.csv(calcium_channel_antisense, file = paste(PATH, "/antisense/ensembl_antisense/balbc_calcium_channel_antisense.csv", sep=""))
write.csv(calcium_channel_both_sig, file = paste(PATH, "/antisense/ensembl_antisense/balbc_calcium_channel_antisense_both_sig.csv", sep=""))

trp_channel_antisense <- subset(antisense_eset, antisense_eset$sense_symbol %in% trp_channel$symbol)
trp_channel_both_sig <- subset(trp_channel_antisense, trp_channel_antisense$gene_pvalue<.05 & trp_channel_antisense$lnc_pvalue<.05)
write.csv(trp_channel_antisense, file = paste(PATH, "/antisense/ensembl_antisense/balbc_trp_channel_antisense.csv", sep=""))
write.csv(trp_channel_both_sig, file = paste(PATH, "/antisense/ensembl_antisense/balbc_trp_channel_antisense_both_sig.csv", sep=""))

write.csv(antisense_eset, file = paste(PATH, "/antisense/ensembl_antisense/sni_vs_sham_balb.c/antisense_all_expression.csv", sep=""))

write.csv(both_significant, file = paste(PATH, "/antisense/ensembl_antisense/sni_vs_sham_balb.c/antisense_both_sig.csv", sep=""))

write.csv(pain_antisense, file = paste(PATH, "/antisense/ensembl_antisense/sni_vs_sham_balb.c/pain_antisense.csv", sep=""))

write.csv(anticorellated, file = paste(PATH, "/antisense/ensembl_antisense/sni_vs_sham_balb.c/antisense_anticorellated.csv", sep=""))

write.csv(anticorellated_sig, file = paste(PATH, "/antisense/ensembl_antisense/sni_vs_sham_balb.c/antisense_anticorellated_sig.csv", sep=""))

write.csv(pain_anticorellated, file = paste(PATH, "/antisense/ensembl_antisense/sni_vs_sham_balb.c/pain_anticorellated.csv", sep=""))

########
# B10.D2
########
system("mkdir -p /home/george/Desktop/LncRNA_v2/Mouse/antisense/ensembl_antisense/sni_vs_sham_b10.d2")

pain <- ens_genes_sham_vs_sni_B10.D2[rownames(ens_genes_sham_vs_sni_B10.D2) %in% pain_genes$Gene.ID,]
pain_sig <- subset(pain, pain$padj <.05 & !is.na(pain$padj))


sodium_channel <- ens_genes_sham_vs_sni_B10.D2[grep("Scn", ens_genes_sham_vs_sni_B10.D2$symbol),]
sodium_channel_sig <- subset(sodium_channel, sodium_channel$padj <.05 & !is.na(sodium_channel$padj))


potassium_channel <- ens_genes_sham_vs_sni_B10.D2[grep("Kcn", ens_genes_sham_vs_sni_B10.D2$symbol),]
potassium_channel_sig <- subset(potassium_channel, potassium_channel$padj <.05 & !is.na(potassium_channel$padj))

chloride_channel <- ens_genes_sham_vs_sni_B10.D2[grep("Clc", ens_genes_sham_vs_sni_B10.D2$symbol),]
chloride_channel_sig <- subset(chloride_channel, chloride_channel$padj <.05 & !is.na(chloride_channel$padj))

calcium_channel <- ens_genes_sham_vs_sni_B10.D2[grep("Cac", ens_genes_sham_vs_sni_B10.D2$symbol),]
calcium_channel_sig <- subset(calcium_channel, calcium_channel$padj <.05 & !is.na(calcium_channel$padj))

trp_channel <- ens_genes_sham_vs_sni_B10.D2[grep("Trp", ens_genes_sham_vs_sni_B10.D2$symbol),]
trp_channel_sig <- subset(trp_channel, trp_channel$padj <.05 & !is.na(trp_channel$padj))

antisense_eset$lnc_lfc <- 0
antisense_eset$lnc_pvalue <- 0
antisense_eset$lnc_baseMean <- 0
antisense_eset$gene_lfc <- 0
antisense_eset$gene_pvalue <- 0



pb <- txtProgressBar(min = 0, max = nrow(antisense_eset), style = 3)
for (i in 1:nrow(antisense_eset)) {
	
	antisense_eset$gene_lfc[i] <- ens_genes_sham_vs_sni_B10.D2[as.character(antisense_eset$V16[i]),2]
	antisense_eset$gene_pvalue[i] <- ens_genes_sham_vs_sni_B10.D2[as.character(antisense_eset$V16[i]),6]
	
			
	antisense_eset$lnc_lfc[i] <- ens_genes_sham_vs_sni_B10.D2[as.character(antisense_eset$V4[i]),2]
	antisense_eset$lnc_pvalue[i] <- ens_genes_sham_vs_sni_B10.D2[as.character(antisense_eset$V4[i]),6]
	antisense_eset$lnc_baseMean[i] <- ens_genes_sham_vs_sni_B10.D2[as.character(antisense_eset$V4[i]),1] 
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
write.csv(sodium_channel_antisense, file = paste(PATH, "/antisense/ensembl_antisense/b10.d2_sodium_channel_antisense.csv", sep=""))
write.csv(sodium_channel_both_sig, file = paste(PATH, "/antisense/ensembl_antisense/b10.d2_sodium_channel_antisense_both_sig.csv", sep=""))

potassium_channel_antisense <- subset(antisense_eset, antisense_eset$sense_symbol %in% potassium_channel$symbol)
potassium_channel_both_sig <- subset(potassium_channel_antisense, potassium_channel_antisense$gene_pvalue<.05 & potassium_channel_antisense$lnc_pvalue<.05)
write.csv(potassium_channel_antisense, file = paste(PATH, "/antisense/ensembl_antisense/b10.d2_potassium_channel_antisense.csv", sep=""))
write.csv(potassium_channel_both_sig, file = paste(PATH, "/antisense/ensembl_antisense/b10.d2_potassium_channel_antisense_both_sig.csv", sep=""))

chloride_channel_antisense <- subset(antisense_eset, antisense_eset$sense_symbol %in% chloride_channel$symbol)
chloride_channel_both_sig <- subset(chloride_channel_antisense, chloride_channel_antisense$gene_pvalue<.05 & chloride_channel_antisense$lnc_pvalue<.05)
write.csv(chloride_channel_antisense, file = paste(PATH, "/antisense/ensembl_antisense/b10.d2_chloride_channel_antisense.csv", sep=""))
write.csv(chloride_channel_both_sig, file = paste(PATH, "/antisense/ensembl_antisense/b10.d2_chloride_channel_antisense_both_sig.csv", sep=""))

calcium_channel_antisense <- subset(antisense_eset, antisense_eset$sense_symbol %in% calcium_channel$symbol)
calcium_channel_both_sig <- subset(calcium_channel_antisense, calcium_channel_antisense$gene_pvalue<.05 & calcium_channel_antisense$lnc_pvalue<.05)
write.csv(calcium_channel_antisense, file = paste(PATH, "/antisense/ensembl_antisense/b10.d2_calcium_channel_antisense.csv", sep=""))
write.csv(calcium_channel_both_sig, file = paste(PATH, "/antisense/ensembl_antisense/b10.d2_calcium_channel_antisense_both_sig.csv", sep=""))

trp_channel_antisense <- subset(antisense_eset, antisense_eset$sense_symbol %in% trp_channel$symbol)
trp_channel_both_sig <- subset(trp_channel_antisense, trp_channel_antisense$gene_pvalue<.05 & trp_channel_antisense$lnc_pvalue<.05)
write.csv(trp_channel_antisense, file = paste(PATH, "/antisense/ensembl_antisense/b10.d2_trp_channel_antisense.csv", sep=""))
write.csv(trp_channel_both_sig, file = paste(PATH, "/antisense/ensembl_antisense/b10.d2_trp_channel_antisense_both_sig.csv", sep=""))

write.csv(antisense_eset, file = paste(PATH, "/antisense/ensembl_antisense/sni_vs_sham_b10.d2/antisense_all_expression.csv", sep=""))

write.csv(both_significant, file = paste(PATH, "/antisense/ensembl_antisense/sni_vs_sham_b10.d2/antisense_both_sig.csv", sep=""))

write.csv(pain_antisense, file = paste(PATH, "/antisense/ensembl_antisense/sni_vs_sham_b10.d2/pain_antisense.csv", sep=""))

write.csv(anticorellated, file = paste(PATH, "/antisense/ensembl_antisense/sni_vs_sham_b10.d2/antisense_anticorellated.csv", sep=""))

write.csv(anticorellated_sig, file = paste(PATH, "/antisense/ensembl_antisense/sni_vs_sham_b10.d2/antisense_anticorellated_sig.csv", sep=""))

write.csv(pain_anticorellated, file = paste(PATH, "/antisense/ensembl_antisense/sni_vs_sham_b10.d2/pain_anticorellated.csv", sep=""))

####################
# B10.D2 vs BALB/c #
####################
system("mkdir -p /home/george/Desktop/LncRNA_v2/Mouse/antisense/ensembl_antisense/diff_B10.D2_vs_BALB.c")

pain <- ens_genes_sni_diff_B10.D2_vs_BALB.c[rownames(ens_genes_sni_diff_B10.D2_vs_BALB.c) %in% pain_genes$Gene.ID,]
pain_sig <- subset(pain, pain$padj <.05 & !is.na(pain$padj))


sodium_channel <- ens_genes_sni_diff_B10.D2_vs_BALB.c[grep("Scn", ens_genes_sni_diff_B10.D2_vs_BALB.c$symbol),]
sodium_channel_sig <- subset(sodium_channel, sodium_channel$padj <.05 & !is.na(sodium_channel$padj))


potassium_channel <- ens_genes_sni_diff_B10.D2_vs_BALB.c[grep("Kcn", ens_genes_sni_diff_B10.D2_vs_BALB.c$symbol),]
potassium_channel_sig <- subset(potassium_channel, potassium_channel$padj <.05 & !is.na(potassium_channel$padj))

chloride_channel <- ens_genes_sni_diff_B10.D2_vs_BALB.c[grep("Clc", ens_genes_sni_diff_B10.D2_vs_BALB.c$symbol),]
chloride_channel_sig <- subset(chloride_channel, chloride_channel$padj <.05 & !is.na(chloride_channel$padj))

calcium_channel <- ens_genes_sni_diff_B10.D2_vs_BALB.c[grep("Cac", ens_genes_sni_diff_B10.D2_vs_BALB.c$symbol),]
calcium_channel_sig <- subset(calcium_channel, calcium_channel$padj <.05 & !is.na(calcium_channel$padj))

trp_channel <- ens_genes_sni_diff_B10.D2_vs_BALB.c[grep("Trp", ens_genes_sni_diff_B10.D2_vs_BALB.c$symbol),]
trp_channel_sig <- subset(trp_channel, trp_channel$padj <.05 & !is.na(trp_channel$padj))

antisense_eset$lnc_lfc <- 0
antisense_eset$lnc_pvalue <- 0
antisense_eset$lnc_baseMean <- 0
antisense_eset$gene_lfc <- 0
antisense_eset$gene_pvalue <- 0



pb <- txtProgressBar(min = 0, max = nrow(antisense_eset), style = 3)
for (i in 1:nrow(antisense_eset)) {
	
	antisense_eset$gene_lfc[i] <- ens_genes_sni_diff_B10.D2_vs_BALB.c[as.character(antisense_eset$V16[i]),2]
	antisense_eset$gene_pvalue[i] <- ens_genes_sni_diff_B10.D2_vs_BALB.c[as.character(antisense_eset$V16[i]),6]
	
			
	antisense_eset$lnc_lfc[i] <- ens_genes_sni_diff_B10.D2_vs_BALB.c[as.character(antisense_eset$V4[i]),2]
	antisense_eset$lnc_pvalue[i] <- ens_genes_sni_diff_B10.D2_vs_BALB.c[as.character(antisense_eset$V4[i]),6]
	antisense_eset$lnc_baseMean[i] <- ens_genes_sni_diff_B10.D2_vs_BALB.c[as.character(antisense_eset$V4[i]),1] 
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
write.csv(sodium_channel_antisense, file = paste(PATH, "/antisense/ensembl_antisense/diff_B10.D2_vs_BALB.c_sodium_channel_antisense.csv", sep=""))
write.csv(sodium_channel_both_sig, file = paste(PATH, "/antisense/ensembl_antisense/diff_B10.D2_vs_BALB.c_sodium_channel_antisense_both_sig.csv", sep=""))

potassium_channel_antisense <- subset(antisense_eset, antisense_eset$sense_symbol %in% potassium_channel$symbol)
potassium_channel_both_sig <- subset(potassium_channel_antisense, potassium_channel_antisense$gene_pvalue<.05 & potassium_channel_antisense$lnc_pvalue<.05)
write.csv(potassium_channel_antisense, file = paste(PATH, "/antisense/ensembl_antisense/diff_B10.D2_vs_BALB.c_potassium_channel_antisense.csv", sep=""))
write.csv(potassium_channel_both_sig, file = paste(PATH, "/antisense/ensembl_antisense/diff_B10.D2_vs_BALB.c_potassium_channel_antisense_both_sig.csv", sep=""))

chloride_channel_antisense <- subset(antisense_eset, antisense_eset$sense_symbol %in% chloride_channel$symbol)
chloride_channel_both_sig <- subset(chloride_channel_antisense, chloride_channel_antisense$gene_pvalue<.05 & chloride_channel_antisense$lnc_pvalue<.05)
write.csv(chloride_channel_antisense, file = paste(PATH, "/antisense/ensembl_antisense/diff_B10.D2_vs_BALB.c_chloride_channel_antisense.csv", sep=""))
write.csv(chloride_channel_both_sig, file = paste(PATH, "/antisense/ensembl_antisense/diff_B10.D2_vs_BALB.c_chloride_channel_antisense_both_sig.csv", sep=""))

calcium_channel_antisense <- subset(antisense_eset, antisense_eset$sense_symbol %in% calcium_channel$symbol)
calcium_channel_both_sig <- subset(calcium_channel_antisense, calcium_channel_antisense$gene_pvalue<.05 & calcium_channel_antisense$lnc_pvalue<.05)
write.csv(calcium_channel_antisense, file = paste(PATH, "/antisense/ensembl_antisense/diff_B10.D2_vs_BALB.c_calcium_channel_antisense.csv", sep=""))
write.csv(calcium_channel_both_sig, file = paste(PATH, "/antisense/ensembl_antisense/diff_B10.D2_vs_BALB.c_calcium_channel_antisense_both_sig.csv", sep=""))

trp_channel_antisense <- subset(antisense_eset, antisense_eset$sense_symbol %in% trp_channel$symbol)
trp_channel_both_sig <- subset(trp_channel_antisense, trp_channel_antisense$gene_pvalue<.05 & trp_channel_antisense$lnc_pvalue<.05)
write.csv(trp_channel_antisense, file = paste(PATH, "/antisense/ensembl_antisense/diff_B10.D2_vs_BALB.c_trp_channel_antisense.csv", sep=""))
write.csv(trp_channel_both_sig, file = paste(PATH, "/antisense/ensembl_antisense/diff_B10.D2_vs_BALB.c_trp_channel_antisense_both_sig.csv", sep=""))

write.csv(antisense_eset, file = paste(PATH, "/antisense/ensembl_antisense/diff_B10.D2_vs_BALB.c/antisense_all_expression.csv", sep=""))

write.csv(both_significant, file = paste(PATH, "/antisense/ensembl_antisense/diff_B10.D2_vs_BALB.c/antisense_both_sig.csv", sep=""))

write.csv(pain_antisense, file = paste(PATH, "/antisense/ensembl_antisense/diff_B10.D2_vs_BALB.c/pain_antisense.csv", sep=""))

write.csv(anticorellated, file = paste(PATH, "/antisense/ensembl_antisense/diff_B10.D2_vs_BALB.c/antisense_anticorellated.csv", sep=""))

write.csv(anticorellated_sig, file = paste(PATH, "/antisense/ensembl_antisense/diff_B10.D2_vs_BALB.c/antisense_anticorellated_sig.csv", sep=""))

write.csv(pain_anticorellated, file = paste(PATH, "/antisense/ensembl_antisense/diff_B10.D2_vs_BALB.c/pain_anticorellated.csv", sep=""))

#################
# Novel LncRNAs #
#################
lncRNAs_gtf <- read.table(file = paste(PATH,"/myRegs_cov/Novel_lncRNAs_mouse_DRG.gtf", sep=""), sep="\t")

toc_lncs <- read.csv(file = paste(PATH,"/LncRNAs_toc_mouse_DRG.csv", sep=""), row.name=1)

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
antisense_mapping$symbol <- mapIds(org.Mm.eg.db, keys=as.character(antisense_mapping[,2]), column="SYMBOL", keytype="ENSEMBL", multiVals="first")

antisense_mapping_refseq$symbol <- mapIds(org.Mm.eg.db, keys=as.character(antisense_mapping_refseq[,2]), column="SYMBOL", keytype="REFSEQ", multiVals="first")

intronic_mapping$symbol <- mapIds(org.Mm.eg.db, keys=as.character(intronic_mapping[,2]), column="SYMBOL", keytype="ENSEMBL", multiVals="first")

intronic_eset <- intronic_mapping

antisense_eset <- antisense_mapping

antisense_eset_only_refseq <- antisense_mapping_refseq[!antisense_mapping_refseq$V4 %in% antisense_mapping$V4, ] 

antisense_eset <- rbind(antisense_eset, antisense_eset_only_refseq)

antisense_names <- name_mapping[name_mapping[,2] %in% antisense_eset$V4, ]

rownames(antisense_names) <- antisense_names[,2]

pain_genes <- read.csv(file = paste(PATH,"/../pain_genes_homologs.csv", sep=""))

pain_genes <- unique(pain_genes[,c(2:3)])

pain_genes_symbols <- unique(as.character(pain_genes$MGI.symbol))


BALB.c <- read.csv(file = paste(PATH,"/DE/res_sni_vs_sham_BALB.c.csv", sep=""), row.name=1)
B10.D2 <- read.csv(file = paste(PATH,"/DE/res_sni_vs_sham_B10.D2.csv", sep=""), row.name=1)
Diff_B10.D2_vs_BALB.c <- read.csv(file = paste(PATH,"/DE/res_sni_diff_B10.D2_vs_BALB.c.csv", sep=""), row.name=1)


ens_genes_sham_vs_sni_BALB.c <- BALB.c[grep("ENS", rownames(BALB.c)),]
ens_genes_sham_vs_sni_B10.D2 <- B10.D2[grep("ENS", rownames(B10.D2)),]
ens_genes_sni_diff_B10.D2_vs_BALB.c <- Diff_B10.D2_vs_BALB.c[grep("ENS", rownames(Diff_B10.D2_vs_BALB.c)),]

lncs_sham_vs_sni_BALB.c <- BALB.c[grep("LncRNA", rownames(BALB.c)),]
lncs_sham_vs_sni_B10.D2 <- B10.D2[grep("LncRNA", rownames(B10.D2)),]
lncs_sni_diff_B10.D2_vs_BALB.c <- Diff_B10.D2_vs_BALB.c[grep("LncRNA", rownames(Diff_B10.D2_vs_BALB.c)),]

################
# SNI vs Sham BALB.c #
################
system("mkdir /home/george/Desktop/LncRNA_v2/Mouse/antisense/sni_vs_sham_balb.c")

pain <- ens_genes_sham_vs_sni_BALB.c[rownames(ens_genes_sham_vs_sni_BALB.c) %in% pain_genes$Gene.ID,]
pain_sig <- subset(pain, pain$padj <.05 & !is.na(pain$padj))
write.csv(pain, file = paste(PATH, "/balbc_pain.csv", sep=""))
write.csv(pain_sig, file = paste(PATH, "/balbc_pain_sig.csv", sep=""))

sodium_channel <- ens_genes_sham_vs_sni_BALB.c[grep("Scn", ens_genes_sham_vs_sni_BALB.c$symbol),]
sodium_channel_sig <- subset(sodium_channel, sodium_channel$padj <.05 & !is.na(sodium_channel$padj))
write.csv(pain, file = paste(PATH, "/balbc_sodium_channel.csv", sep=""))
write.csv(pain_sig, file = paste(PATH, "/balbc_sodium_channel_sig.csv", sep=""))

potassium_channel <- ens_genes_sham_vs_sni_BALB.c[grep("Kcn", ens_genes_sham_vs_sni_BALB.c$symbol),]
potassium_channel_sig <- subset(potassium_channel, potassium_channel$padj <.05 & !is.na(potassium_channel$padj))
write.csv(potassium_channel, file = paste(PATH, "/balbc_potassium_channel.csv", sep=""))
write.csv(potassium_channel_sig, file = paste(PATH, "/balbc_potassium_channel_sig.csv", sep=""))

chloride_channel <- ens_genes_sham_vs_sni_BALB.c[grep("Clc", ens_genes_sham_vs_sni_BALB.c$symbol),]
chloride_channel_sig <- subset(chloride_channel, chloride_channel$padj <.05 & !is.na(chloride_channel$padj))
write.csv(chloride_channel, file = paste(PATH, "/balbc_chloride_channel.csv", sep=""))
write.csv(chloride_channel_sig, file = paste(PATH, "/balbc_chloride_channel_sig.csv", sep=""))

calcium_channel <- ens_genes_sham_vs_sni_BALB.c[grep("Cac", ens_genes_sham_vs_sni_BALB.c$symbol),]
calcium_channel_sig <- subset(calcium_channel, calcium_channel$padj <.05 & !is.na(calcium_channel$padj))
write.csv(calcium_channel, file = paste(PATH, "/balbc_calcium_channel.csv", sep=""))
write.csv(calcium_channel_sig, file = paste(PATH, "/balbc_calcium_channel_sig.csv", sep=""))

trp_channel <- ens_genes_sham_vs_sni_BALB.c[grep("Trp", ens_genes_sham_vs_sni_BALB.c$symbol),]
trp_channel_sig <- subset(trp_channel, trp_channel$padj <.05 & !is.na(trp_channel$padj))
write.csv(trp_channel, file = paste(PATH, "/balbc_trp_channel.csv", sep=""))
write.csv(trp_channel_sig, file = paste(PATH, "/balbc_trp_channel_sig.csv", sep=""))

antisense_eset$lnc_lfc <- 0
antisense_eset$lnc_pvalue <- 0
antisense_eset$lnc_baseMean <- 0
antisense_eset$gene_lfc <- 0
antisense_eset$gene_pvalue <- 0



pb <- txtProgressBar(min = 0, max = nrow(antisense_eset), style = 3)
for (i in 1:nrow(antisense_eset)) {
	if (grepl("^N", as.character(antisense_eset$V16[i]))) { 
	if (length(grep(as.character(antisense_eset$symbol[i]), ens_genes_sham_vs_sni_BALB.c$symbol)) >0 ) {  	
	antisense_eset$gene_lfc[i] <- ens_genes_sham_vs_sni_BALB.c[grep(as.character(antisense_eset$symbol[i]), ens_genes_sham_vs_sni_BALB.c$symbol),2]
	antisense_eset$gene_pvalue[i] <- ens_genes_sham_vs_sni_BALB.c[grep(as.character(antisense_eset$symbol[i]), ens_genes_sham_vs_sni_BALB.c$symbol),6]}
	}
	else {	
	antisense_eset$gene_lfc[i] <- ens_genes_sham_vs_sni_BALB.c[as.character(antisense_eset$V16[i]),2]
	antisense_eset$gene_pvalue[i] <- ens_genes_sham_vs_sni_BALB.c[as.character(antisense_eset$V16[i]),6]
	}
	antisense_eset$ID[i] <- antisense_names[as.character(antisense_eset$V4)[i],1]
		
	antisense_eset$lnc_lfc[i] <- lncs_sham_vs_sni_BALB.c[as.character(antisense_eset$ID[i]),2]
	antisense_eset$lnc_pvalue[i] <- lncs_sham_vs_sni_BALB.c[as.character(antisense_eset$ID[i]),6]
	antisense_eset$lnc_baseMean[i] <- lncs_sham_vs_sni_BALB.c[as.character(antisense_eset$ID[i]),1] 
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
write.csv(sodium_channel_antisense, file = paste(PATH, "/balbc_sodium_channel_antisense.csv", sep=""))
write.csv(sodium_channel_both_sig, file = paste(PATH, "/balbc_sodium_channel_antisense_both_sig.csv", sep=""))

potassium_channel_antisense <- subset(antisense_eset, antisense_eset$symbol %in% potassium_channel$symbol)
potassium_channel_both_sig <- subset(potassium_channel_antisense, potassium_channel_antisense$gene_pvalue<.05 & potassium_channel_antisense$lnc_pvalue<.05)
write.csv(potassium_channel_antisense, file = paste(PATH, "/balbc_potassium_channel_antisense.csv", sep=""))
write.csv(potassium_channel_both_sig, file = paste(PATH, "/balbc_potassium_channel_antisense_both_sig.csv", sep=""))

chloride_channel_antisense <- subset(antisense_eset, antisense_eset$symbol %in% chloride_channel$symbol)
chloride_channel_both_sig <- subset(chloride_channel_antisense, chloride_channel_antisense$gene_pvalue<.05 & chloride_channel_antisense$lnc_pvalue<.05)
write.csv(chloride_channel_antisense, file = paste(PATH, "/balbc_chloride_channel_antisense.csv", sep=""))
write.csv(chloride_channel_both_sig, file = paste(PATH, "/balbc_chloride_channel_antisense_both_sig.csv", sep=""))

calcium_channel_antisense <- subset(antisense_eset, antisense_eset$symbol %in% calcium_channel$symbol)
calcium_channel_both_sig <- subset(calcium_channel_antisense, calcium_channel_antisense$gene_pvalue<.05 & calcium_channel_antisense$lnc_pvalue<.05)
write.csv(calcium_channel_antisense, file = paste(PATH, "/balbc_calcium_channel_antisense.csv", sep=""))
write.csv(calcium_channel_both_sig, file = paste(PATH, "/balbc_calcium_channel_antisense_both_sig.csv", sep=""))

trp_channel_antisense <- subset(antisense_eset, antisense_eset$symbol %in% trp_channel$symbol)
trp_channel_both_sig <- subset(trp_channel_antisense, trp_channel_antisense$gene_pvalue<.05 & trp_channel_antisense$lnc_pvalue<.05)
write.csv(trp_channel_antisense, file = paste(PATH, "/balbc_trp_channel_antisense.csv", sep=""))
write.csv(trp_channel_both_sig, file = paste(PATH, "/balbc_trp_channel_antisense_both_sig.csv", sep=""))

write.csv(antisense_eset, file = paste(PATH, "/antisense/sni_vs_sham_balb.c/antisense_all_expression.csv", sep=""))

write.csv(both_significant, file = paste(PATH, "/antisense/sni_vs_sham_balb.c/antisense_both_sig.csv", sep=""))

write.csv(pain_antisense, file = paste(PATH, "/antisense/sni_vs_sham_balb.c/pain_antisense.csv", sep=""))

write.csv(anticorellated, file = paste(PATH, "/antisense/sni_vs_sham_balb.c/antisense_anticorellated.csv", sep=""))

write.csv(anticorellated_sig, file = paste(PATH, "/antisense/sni_vs_sham_balb.c/antisense_anticorellated_sig.csv", sep=""))

write.csv(pain_anticorellated, file = paste(PATH, "/antisense/sni_vs_sham_balb.c/pain_anticorellated.csv", sep=""))

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
antisense_counts$cor[i] <- cor.test(x=as.numeric(antisense_counts[i,c(5:24)]), y=as.numeric(antisense_counts[i,c(25:44)]))$estimate  

antisense_counts$cor_pvalue[i] <- cor.test(x=as.numeric(antisense_counts[i,c(5:24)]), y=as.numeric(antisense_counts[i,c(25:44)]))$p.value 
}

antisense_counts$cor_pvalue <- p.adjust(antisense_counts$cor_pvalue, method = "BH")

save(file = paste(PATH, "/antisense/antisense_counts.RData", sep=""), antisense_counts)


#################################################
################
# SNI vs Sham B10.D2 #
################
system("mkdir /home/george/Desktop/LncRNA_v2/Mouse/antisense/sni_vs_sham_b10.d2")

pain <- ens_genes_sham_vs_sni_B10.D2[rownames(ens_genes_sham_vs_sni_B10.D2) %in% pain_genes$Gene.ID,]
pain_sig <- subset(pain, pain$padj <.05 & !is.na(pain$padj))
write.csv(pain, file = paste(PATH, "/b10.d2_pain.csv", sep=""))
write.csv(pain_sig, file = paste(PATH, "/b10.d2_pain_sig.csv", sep=""))

sodium_channel <- ens_genes_sham_vs_sni_B10.D2[grep("Scn", ens_genes_sham_vs_sni_B10.D2$symbol),]
sodium_channel_sig <- subset(sodium_channel, sodium_channel$padj <.05 & !is.na(sodium_channel$padj))
write.csv(pain, file = paste(PATH, "/b10.d2_sodium_channel.csv", sep=""))
write.csv(pain_sig, file = paste(PATH, "/b10.d2_sodium_channel_sig.csv", sep=""))

potassium_channel <- ens_genes_sham_vs_sni_B10.D2[grep("Kcn", ens_genes_sham_vs_sni_B10.D2$symbol),]
potassium_channel_sig <- subset(potassium_channel, potassium_channel$padj <.05 & !is.na(potassium_channel$padj))
write.csv(potassium_channel, file = paste(PATH, "/b10.d2_potassium_channel.csv", sep=""))
write.csv(potassium_channel_sig, file = paste(PATH, "/b10.d2_potassium_channel_sig.csv", sep=""))

chloride_channel <- ens_genes_sham_vs_sni_B10.D2[grep("Clc", ens_genes_sham_vs_sni_B10.D2$symbol),]
chloride_channel_sig <- subset(chloride_channel, chloride_channel$padj <.05 & !is.na(chloride_channel$padj))
write.csv(chloride_channel, file = paste(PATH, "/b10.d2_chloride_channel.csv", sep=""))
write.csv(chloride_channel_sig, file = paste(PATH, "/b10.d2_chloride_channel_sig.csv", sep=""))

calcium_channel <- ens_genes_sham_vs_sni_B10.D2[grep("Cac", ens_genes_sham_vs_sni_B10.D2$symbol),]
calcium_channel_sig <- subset(calcium_channel, calcium_channel$padj <.05 & !is.na(calcium_channel$padj))
write.csv(calcium_channel, file = paste(PATH, "/b10.d2_calcium_channel.csv", sep=""))
write.csv(calcium_channel_sig, file = paste(PATH, "/b10.d2_calcium_channel_sig.csv", sep=""))

trp_channel <- ens_genes_sham_vs_sni_B10.D2[grep("Trp", ens_genes_sham_vs_sni_B10.D2$symbol),]
trp_channel_sig <- subset(trp_channel, trp_channel$padj <.05 & !is.na(trp_channel$padj))
write.csv(trp_channel, file = paste(PATH, "/b10.d2_trp_channel.csv", sep=""))
write.csv(trp_channel_sig, file = paste(PATH, "/b10.d2_trp_channel_sig.csv", sep=""))

antisense_eset$lnc_lfc <- 0
antisense_eset$lnc_pvalue <- 0
antisense_eset$lnc_baseMean <- 0
antisense_eset$gene_lfc <- 0
antisense_eset$gene_pvalue <- 0



pb <- txtProgressBar(min = 0, max = nrow(antisense_eset), style = 3)
for (i in 1:nrow(antisense_eset)) {
	if (grepl("^N", as.character(antisense_eset$V16[i]))) { 
	if (length(grep(as.character(antisense_eset$symbol[i]), ens_genes_sham_vs_sni_B10.D2$symbol)) >0 ) {  	
	antisense_eset$gene_lfc[i] <- ens_genes_sham_vs_sni_B10.D2[grep(as.character(antisense_eset$symbol[i]), ens_genes_sham_vs_sni_B10.D2$symbol),2]
	antisense_eset$gene_pvalue[i] <- ens_genes_sham_vs_sni_B10.D2[grep(as.character(antisense_eset$symbol[i]), ens_genes_sham_vs_sni_B10.D2$symbol),6]}
	}
	else {	
	antisense_eset$gene_lfc[i] <- ens_genes_sham_vs_sni_B10.D2[as.character(antisense_eset$V16[i]),2]
	antisense_eset$gene_pvalue[i] <- ens_genes_sham_vs_sni_B10.D2[as.character(antisense_eset$V16[i]),6]
	}
	antisense_eset$ID[i] <- antisense_names[as.character(antisense_eset$V4)[i],1]
		
	antisense_eset$lnc_lfc[i] <- lncs_sham_vs_sni_B10.D2[as.character(antisense_eset$ID[i]),2]
	antisense_eset$lnc_pvalue[i] <- lncs_sham_vs_sni_B10.D2[as.character(antisense_eset$ID[i]),6]
	antisense_eset$lnc_baseMean[i] <- lncs_sham_vs_sni_B10.D2[as.character(antisense_eset$ID[i]),1] 
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
write.csv(sodium_channel_antisense, file = paste(PATH, "/b10.d2_sodium_channel_antisense.csv", sep=""))
write.csv(sodium_channel_both_sig, file = paste(PATH, "/b10.d2_sodium_channel_antisense_both_sig.csv", sep=""))

potassium_channel_antisense <- subset(antisense_eset, antisense_eset$symbol %in% potassium_channel$symbol)
potassium_channel_both_sig <- subset(potassium_channel_antisense, potassium_channel_antisense$gene_pvalue<.05 & potassium_channel_antisense$lnc_pvalue<.05)
write.csv(potassium_channel_antisense, file = paste(PATH, "/b10.d2_potassium_channel_antisense.csv", sep=""))
write.csv(potassium_channel_both_sig, file = paste(PATH, "/b10.d2_potassium_channel_antisense_both_sig.csv", sep=""))

chloride_channel_antisense <- subset(antisense_eset, antisense_eset$symbol %in% chloride_channel$symbol)
chloride_channel_both_sig <- subset(chloride_channel_antisense, chloride_channel_antisense$gene_pvalue<.05 & chloride_channel_antisense$lnc_pvalue<.05)
write.csv(chloride_channel_antisense, file = paste(PATH, "/b10.d2_chloride_channel_antisense.csv", sep=""))
write.csv(chloride_channel_both_sig, file = paste(PATH, "/b10.d2_chloride_channel_antisense_both_sig.csv", sep=""))

calcium_channel_antisense <- subset(antisense_eset, antisense_eset$symbol %in% calcium_channel$symbol)
calcium_channel_both_sig <- subset(calcium_channel_antisense, calcium_channel_antisense$gene_pvalue<.05 & calcium_channel_antisense$lnc_pvalue<.05)
write.csv(calcium_channel_antisense, file = paste(PATH, "/b10.d2_calcium_channel_antisense.csv", sep=""))
write.csv(calcium_channel_both_sig, file = paste(PATH, "/b10.d2_calcium_channel_antisense_both_sig.csv", sep=""))

trp_channel_antisense <- subset(antisense_eset, antisense_eset$symbol %in% trp_channel$symbol)
trp_channel_both_sig <- subset(trp_channel_antisense, trp_channel_antisense$gene_pvalue<.05 & trp_channel_antisense$lnc_pvalue<.05)
write.csv(trp_channel_antisense, file = paste(PATH, "/b10.d2_trp_channel_antisense.csv", sep=""))
write.csv(trp_channel_both_sig, file = paste(PATH, "/b10.d2_trp_channel_antisense_both_sig.csv", sep=""))

write.csv(antisense_eset, file = paste(PATH, "/antisense/sni_vs_sham_b10.d2/antisense_all_expression.csv", sep=""))

write.csv(both_significant, file = paste(PATH, "/antisense/sni_vs_sham_b10.d2/antisense_both_sig.csv", sep=""))

write.csv(pain_antisense, file = paste(PATH, "/antisense/sni_vs_sham_b10.d2/pain_antisense.csv", sep=""))

write.csv(anticorellated, file = paste(PATH, "/antisense/sni_vs_sham_b10.d2/antisense_anticorellated.csv", sep=""))

write.csv(anticorellated_sig, file = paste(PATH, "/antisense/sni_vs_sham_b10.d2/antisense_anticorellated_sig.csv", sep=""))

write.csv(pain_anticorellated, file = paste(PATH, "/antisense/sni_vs_sham_b10.d2/pain_anticorellated.csv", sep=""))

#################################################

################
# dif SNI vs Sham B10.D2 vs Balb.c #
################
system("mkdir /home/george/Desktop/LncRNA_v2/Mouse/antisense/diff_B10.D2_vs_BALB.c")

pain <- ens_genes_sni_diff_B10.D2_vs_BALB.c[rownames(ens_genes_sni_diff_B10.D2_vs_BALB.c) %in% pain_genes$Gene.ID,]
pain_sig <- subset(pain, pain$padj <.05 & !is.na(pain$padj))
write.csv(pain, file = paste(PATH, "/diff_B10.D2_vs_BALB.c_pain.csv", sep=""))
write.csv(pain_sig, file = paste(PATH, "/diff_B10.D2_vs_BALB.c_pain_sig.csv", sep=""))

sodium_channel <- ens_genes_sni_diff_B10.D2_vs_BALB.c[grep("Scn", ens_genes_sni_diff_B10.D2_vs_BALB.c$symbol),]
sodium_channel_sig <- subset(sodium_channel, sodium_channel$padj <.05 & !is.na(sodium_channel$padj))
write.csv(pain, file = paste(PATH, "/diff_B10.D2_vs_BALB.c_sodium_channel.csv", sep=""))
write.csv(pain_sig, file = paste(PATH, "/diff_B10.D2_vs_BALB.c_sodium_channel_sig.csv", sep=""))

potassium_channel <- ens_genes_sni_diff_B10.D2_vs_BALB.c[grep("Kcn", ens_genes_sni_diff_B10.D2_vs_BALB.c$symbol),]
potassium_channel_sig <- subset(potassium_channel, potassium_channel$padj <.05 & !is.na(potassium_channel$padj))
write.csv(potassium_channel, file = paste(PATH, "/diff_B10.D2_vs_BALB.c_potassium_channel.csv", sep=""))
write.csv(potassium_channel_sig, file = paste(PATH, "/diff_B10.D2_vs_BALB.c_potassium_channel_sig.csv", sep=""))

chloride_channel <- ens_genes_sni_diff_B10.D2_vs_BALB.c[grep("Clc", ens_genes_sni_diff_B10.D2_vs_BALB.c$symbol),]
chloride_channel_sig <- subset(chloride_channel, chloride_channel$padj <.05 & !is.na(chloride_channel$padj))
write.csv(chloride_channel, file = paste(PATH, "/diff_B10.D2_vs_BALB.c_chloride_channel.csv", sep=""))
write.csv(chloride_channel_sig, file = paste(PATH, "/diff_B10.D2_vs_BALB.c_chloride_channel_sig.csv", sep=""))

calcium_channel <- ens_genes_sni_diff_B10.D2_vs_BALB.c[grep("Cac", ens_genes_sni_diff_B10.D2_vs_BALB.c$symbol),]
calcium_channel_sig <- subset(calcium_channel, calcium_channel$padj <.05 & !is.na(calcium_channel$padj))
write.csv(calcium_channel, file = paste(PATH, "/diff_B10.D2_vs_BALB.c_calcium_channel.csv", sep=""))
write.csv(calcium_channel_sig, file = paste(PATH, "/diff_B10.D2_vs_BALB.c_calcium_channel_sig.csv", sep=""))

trp_channel <- ens_genes_sni_diff_B10.D2_vs_BALB.c[grep("Trp", ens_genes_sni_diff_B10.D2_vs_BALB.c$symbol),]
trp_channel_sig <- subset(trp_channel, trp_channel$padj <.05 & !is.na(trp_channel$padj))
write.csv(trp_channel, file = paste(PATH, "/diff_B10.D2_vs_BALB.c_trp_channel.csv", sep=""))
write.csv(trp_channel_sig, file = paste(PATH, "/diff_B10.D2_vs_BALB.c_trp_channel_sig.csv", sep=""))

antisense_eset$lnc_lfc <- 0
antisense_eset$lnc_pvalue <- 0
antisense_eset$lnc_baseMean <- 0
antisense_eset$gene_lfc <- 0
antisense_eset$gene_pvalue <- 0



pb <- txtProgressBar(min = 0, max = nrow(antisense_eset), style = 3)
for (i in 1:nrow(antisense_eset)) {
	if (grepl("^N", as.character(antisense_eset$V16[i]))) { 
	if (length(grep(as.character(antisense_eset$symbol[i]), ens_genes_sni_diff_B10.D2_vs_BALB.c$symbol)) >0 ) {  	
	antisense_eset$gene_lfc[i] <- ens_genes_sni_diff_B10.D2_vs_BALB.c[grep(as.character(antisense_eset$symbol[i]), ens_genes_sni_diff_B10.D2_vs_BALB.c$symbol),2]
	antisense_eset$gene_pvalue[i] <- ens_genes_sni_diff_B10.D2_vs_BALB.c[grep(as.character(antisense_eset$symbol[i]), ens_genes_sni_diff_B10.D2_vs_BALB.c$symbol),6]}
	}
	else {	
	antisense_eset$gene_lfc[i] <- ens_genes_sni_diff_B10.D2_vs_BALB.c[as.character(antisense_eset$V16[i]),2]
	antisense_eset$gene_pvalue[i] <- ens_genes_sni_diff_B10.D2_vs_BALB.c[as.character(antisense_eset$V16[i]),6]
	}
	antisense_eset$ID[i] <- antisense_names[as.character(antisense_eset$V4)[i],1]
		
	antisense_eset$lnc_lfc[i] <- lncs_sni_diff_B10.D2_vs_BALB.c[as.character(antisense_eset$ID[i]),2]
	antisense_eset$lnc_pvalue[i] <- lncs_sni_diff_B10.D2_vs_BALB.c[as.character(antisense_eset$ID[i]),6]
	antisense_eset$lnc_baseMean[i] <- lncs_sni_diff_B10.D2_vs_BALB.c[as.character(antisense_eset$ID[i]),1] 
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
write.csv(sodium_channel_antisense, file = paste(PATH, "/diff_B10.D2_vs_BALB.c_sodium_channel_antisense.csv", sep=""))
write.csv(sodium_channel_both_sig, file = paste(PATH, "/diff_B10.D2_vs_BALB.c_sodium_channel_antisense_both_sig.csv", sep=""))

potassium_channel_antisense <- subset(antisense_eset, antisense_eset$symbol %in% potassium_channel$symbol)
potassium_channel_both_sig <- subset(potassium_channel_antisense, potassium_channel_antisense$gene_pvalue<.05 & potassium_channel_antisense$lnc_pvalue<.05)
write.csv(potassium_channel_antisense, file = paste(PATH, "/diff_B10.D2_vs_BALB.c_potassium_channel_antisense.csv", sep=""))
write.csv(potassium_channel_both_sig, file = paste(PATH, "/diff_B10.D2_vs_BALB.c_potassium_channel_antisense_both_sig.csv", sep=""))

chloride_channel_antisense <- subset(antisense_eset, antisense_eset$symbol %in% chloride_channel$symbol)
chloride_channel_both_sig <- subset(chloride_channel_antisense, chloride_channel_antisense$gene_pvalue<.05 & chloride_channel_antisense$lnc_pvalue<.05)
write.csv(chloride_channel_antisense, file = paste(PATH, "/diff_B10.D2_vs_BALB.c_chloride_channel_antisense.csv", sep=""))
write.csv(chloride_channel_both_sig, file = paste(PATH, "/diff_B10.D2_vs_BALB.c_chloride_channel_antisense_both_sig.csv", sep=""))

calcium_channel_antisense <- subset(antisense_eset, antisense_eset$symbol %in% calcium_channel$symbol)
calcium_channel_both_sig <- subset(calcium_channel_antisense, calcium_channel_antisense$gene_pvalue<.05 & calcium_channel_antisense$lnc_pvalue<.05)
write.csv(calcium_channel_antisense, file = paste(PATH, "/diff_B10.D2_vs_BALB.c_calcium_channel_antisense.csv", sep=""))
write.csv(calcium_channel_both_sig, file = paste(PATH, "/diff_B10.D2_vs_BALB.c_calcium_channel_antisense_both_sig.csv", sep=""))

trp_channel_antisense <- subset(antisense_eset, antisense_eset$symbol %in% trp_channel$symbol)
trp_channel_both_sig <- subset(trp_channel_antisense, trp_channel_antisense$gene_pvalue<.05 & trp_channel_antisense$lnc_pvalue<.05)
write.csv(trp_channel_antisense, file = paste(PATH, "/diff_B10.D2_vs_BALB.c_trp_channel_antisense.csv", sep=""))
write.csv(trp_channel_both_sig, file = paste(PATH, "/diff_B10.D2_vs_BALB.c_trp_channel_antisense_both_sig.csv", sep=""))

write.csv(antisense_eset, file = paste(PATH, "/antisense/diff_B10.D2_vs_BALB.c/antisense_all_expression.csv", sep=""))

write.csv(both_significant, file = paste(PATH, "/antisense/diff_B10.D2_vs_BALB.c/antisense_both_sig.csv", sep=""))

write.csv(pain_antisense, file = paste(PATH, "/antisense/diff_B10.D2_vs_BALB.c/pain_antisense.csv", sep=""))

write.csv(anticorellated, file = paste(PATH, "/antisense/diff_B10.D2_vs_BALB.c/antisense_anticorellated.csv", sep=""))

write.csv(anticorellated_sig, file = paste(PATH, "/antisense/diff_B10.D2_vs_BALB.c/antisense_anticorellated_sig.csv", sep=""))

write.csv(pain_anticorellated, file = paste(PATH, "/antisense/diff_B10.D2_vs_BALB.c/pain_anticorellated.csv", sep=""))












##################################
########### Intronic #############
##################################
################
# SNI vs Sham BALB.c #
################
system("mkdir /home/george/Desktop/LncRNA_v2/Mouse/intronic")
system("mkdir /home/george/Desktop/LncRNA_v2/Mouse/intronic/sni_vs_sham_balb.c")

intronic_names <- name_mapping[name_mapping[,2] %in% intronic_eset$V4, ]

rownames(intronic_names) <- intronic_names[,2]


pain <- ens_genes_sham_vs_sni_BALB.c[rownames(ens_genes_sham_vs_sni_BALB.c) %in% pain_genes$Gene.ID,]
pain_sig <- subset(pain, pain$padj <.05 & !is.na(pain$padj))

intronic_eset$lnc_lfc <- 0
intronic_eset$lnc_pvalue <- 0
intronic_eset$lnc_baseMean <- 0
intronic_eset$gene_lfc <- 0
intronic_eset$gene_pvalue <- 0



pb <- txtProgressBar(min = 0, max = nrow(intronic_eset), style = 3)
for (i in 1:nrow(intronic_eset)) {
	if (grepl("^N", as.character(intronic_eset$V16[i]))) { 
	if (length(grep(as.character(intronic_eset$symbol[i]), ens_genes_sham_vs_sni_BALB.c$symbol)) >0 ) {  	
	intronic_eset$gene_lfc[i] <- ens_genes_sham_vs_sni_BALB.c[grep(as.character(intronic_eset$symbol[i]), ens_genes_sham_vs_sni_BALB.c$symbol),2]
	intronic_eset$gene_pvalue[i] <- ens_genes_sham_vs_sni_BALB.c[grep(as.character(intronic_eset$symbol[i]), ens_genes_sham_vs_sni_BALB.c$symbol),6]}
	}
	else {	
	intronic_eset$gene_lfc[i] <- ens_genes_sham_vs_sni_BALB.c[as.character(intronic_eset$V16[i]),2]
	intronic_eset$gene_pvalue[i] <- ens_genes_sham_vs_sni_BALB.c[as.character(intronic_eset$V16[i]),6]
	}
	intronic_eset$ID[i] <- intronic_names[as.character(intronic_eset$V4)[i],1]
		
	intronic_eset$lnc_lfc[i] <- lncs_sham_vs_sni_BALB.c[as.character(intronic_eset$ID[i]),2]
	intronic_eset$lnc_pvalue[i] <- lncs_sham_vs_sni_BALB.c[as.character(intronic_eset$ID[i]),6]
	intronic_eset$lnc_baseMean[i] <- lncs_sham_vs_sni_BALB.c[as.character(intronic_eset$ID[i]),1] 
	setTxtProgressBar(pb, i)
	
	}
close(pb)
cat("\n")


both_significant <- subset(intronic_eset, intronic_eset$gene_pvalue<.05 & intronic_eset$lnc_pvalue<.05)

anticorellated <- subset(intronic_eset, (intronic_eset$lnc_lfc * intronic_eset$gene_lfc) < 0)

anticorellated_sig <- subset(both_significant, (both_significant$lnc_lfc * both_significant$gene_lfc) < 0) 

pain_intronic <- subset(intronic_eset, intronic_eset$symbol %in% pain_genes_symbols)

pain_anticorellated <- subset(pain_intronic, (pain_intronic$lnc_lfc * pain_intronic$gene_lfc) < 0)


write.csv(intronic_eset, file = paste(PATH, "/intronic/sni_vs_sham_balb.c/intronic_all_expression.csv", sep=""))

write.csv(both_significant, file = paste(PATH, "/intronic/sni_vs_sham_balb.c/intronic_both_sig.csv", sep=""))

write.csv(pain_intronic, file = paste(PATH, "/intronic/sni_vs_sham_balb.c/pain_intronic.csv", sep=""))

write.csv(anticorellated, file = paste(PATH, "/intronic/sni_vs_sham_balb.c/intronic_anticorellated.csv", sep=""))

write.csv(anticorellated_sig, file = paste(PATH, "/intronic/sni_vs_sham_balb.c/intronic_anticorellated_sig.csv", sep=""))

write.csv(pain_anticorellated, file = paste(PATH, "/intronic/sni_vs_sham_balb.c/pain_anticorellated.csv", sep=""))

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
intronic_counts$cor[i] <- cor.test(x=as.numeric(intronic_counts[i,c(5:24)]), y=as.numeric(intronic_counts[i,c(25:44)]))$estimate  

intronic_counts$cor_pvalue[i] <- cor.test(x=as.numeric(intronic_counts[i,c(5:24)]), y=as.numeric(intronic_counts[i,c(25:44)]))$p.value 
}

intronic_counts$cor_pvalue <- p.adjust(intronic_counts$cor_pvalue, method = "BH")

save(file = paste(PATH, "/intronic/intronic_counts.RData", sep=""), intronic_counts)


#################################################
################
# SNI vs Sham B10.D2 #
################
system("mkdir /home/george/Desktop/LncRNA_v2/Mouse/intronic/sni_vs_sham_b10.d2")

pain <- ens_genes_sham_vs_sni_B10.D2[rownames(ens_genes_sham_vs_sni_B10.D2) %in% pain_genes$Gene.ID,]
pain_sig <- subset(pain, pain$padj <.05 & !is.na(pain$padj))


intronic_eset$lnc_lfc <- 0
intronic_eset$lnc_pvalue <- 0
intronic_eset$lnc_baseMean <- 0
intronic_eset$gene_lfc <- 0
intronic_eset$gene_pvalue <- 0



pb <- txtProgressBar(min = 0, max = nrow(intronic_eset), style = 3)
for (i in 1:nrow(intronic_eset)) {
	if (grepl("^N", as.character(intronic_eset$V16[i]))) { 
	if (length(grep(as.character(intronic_eset$symbol[i]), ens_genes_sham_vs_sni_B10.D2$symbol)) >0 ) {  	
	intronic_eset$gene_lfc[i] <- ens_genes_sham_vs_sni_B10.D2[grep(as.character(intronic_eset$symbol[i]), ens_genes_sham_vs_sni_B10.D2$symbol),2]
	intronic_eset$gene_pvalue[i] <- ens_genes_sham_vs_sni_B10.D2[grep(as.character(intronic_eset$symbol[i]), ens_genes_sham_vs_sni_B10.D2$symbol),6]}
	}
	else {	
	intronic_eset$gene_lfc[i] <- ens_genes_sham_vs_sni_B10.D2[as.character(intronic_eset$V16[i]),2]
	intronic_eset$gene_pvalue[i] <- ens_genes_sham_vs_sni_B10.D2[as.character(intronic_eset$V16[i]),6]
	}
	intronic_eset$ID[i] <- intronic_names[as.character(intronic_eset$V4)[i],1]
		
	intronic_eset$lnc_lfc[i] <- lncs_sham_vs_sni_B10.D2[as.character(intronic_eset$ID[i]),2]
	intronic_eset$lnc_pvalue[i] <- lncs_sham_vs_sni_B10.D2[as.character(intronic_eset$ID[i]),6]
	intronic_eset$lnc_baseMean[i] <- lncs_sham_vs_sni_B10.D2[as.character(intronic_eset$ID[i]),1] 
	setTxtProgressBar(pb, i)
	
	}
close(pb)
cat("\n")


both_significant <- subset(intronic_eset, intronic_eset$gene_pvalue<.05 & intronic_eset$lnc_pvalue<.05)

anticorellated <- subset(intronic_eset, (intronic_eset$lnc_lfc * intronic_eset$gene_lfc) < 0)

anticorellated_sig <- subset(both_significant, (both_significant$lnc_lfc * both_significant$gene_lfc) < 0) 

pain_intronic <- subset(intronic_eset, intronic_eset$symbol %in% pain_genes_symbols)

pain_anticorellated <- subset(pain_intronic, (pain_intronic$lnc_lfc * pain_intronic$gene_lfc) < 0)


write.csv(intronic_eset, file = paste(PATH, "/intronic/sni_vs_sham_b10.d2/intronic_all_expression.csv", sep=""))

write.csv(both_significant, file = paste(PATH, "/intronic/sni_vs_sham_b10.d2/intronic_both_sig.csv", sep=""))

write.csv(pain_intronic, file = paste(PATH, "/intronic/sni_vs_sham_b10.d2/pain_intronic.csv", sep=""))

write.csv(anticorellated, file = paste(PATH, "/intronic/sni_vs_sham_b10.d2/intronic_anticorellated.csv", sep=""))

write.csv(anticorellated_sig, file = paste(PATH, "/intronic/sni_vs_sham_b10.d2/intronic_anticorellated_sig.csv", sep=""))

write.csv(pain_anticorellated, file = paste(PATH, "/intronic/sni_vs_sham_b10.d2/pain_anticorellated.csv", sep=""))

#################################################

################
# dif SNI vs Sham B10.D2 vs Balb.c #
################
system("mkdir /home/george/Desktop/LncRNA_v2/Mouse/intronic/diff_B10.D2_vs_BALB.c")

pain <- ens_genes_sni_diff_B10.D2_vs_BALB.c[rownames(ens_genes_sni_diff_B10.D2_vs_BALB.c) %in% pain_genes$Gene.ID,]
pain_sig <- subset(pain, pain$padj <.05 & !is.na(pain$padj))

intronic_eset$lnc_lfc <- 0
intronic_eset$lnc_pvalue <- 0
intronic_eset$lnc_baseMean <- 0
intronic_eset$gene_lfc <- 0
intronic_eset$gene_pvalue <- 0



pb <- txtProgressBar(min = 0, max = nrow(intronic_eset), style = 3)
for (i in 1:nrow(intronic_eset)) {
	if (grepl("^N", as.character(intronic_eset$V16[i]))) { 
	if (length(grep(as.character(intronic_eset$symbol[i]), ens_genes_sni_diff_B10.D2_vs_BALB.c$symbol)) >0 ) {  	
	intronic_eset$gene_lfc[i] <- ens_genes_sni_diff_B10.D2_vs_BALB.c[grep(as.character(intronic_eset$symbol[i]), ens_genes_sni_diff_B10.D2_vs_BALB.c$symbol),2]
	intronic_eset$gene_pvalue[i] <- ens_genes_sni_diff_B10.D2_vs_BALB.c[grep(as.character(intronic_eset$symbol[i]), ens_genes_sni_diff_B10.D2_vs_BALB.c$symbol),6]}
	}
	else {	
	intronic_eset$gene_lfc[i] <- ens_genes_sni_diff_B10.D2_vs_BALB.c[as.character(intronic_eset$V16[i]),2]
	intronic_eset$gene_pvalue[i] <- ens_genes_sni_diff_B10.D2_vs_BALB.c[as.character(intronic_eset$V16[i]),6]
	}
	intronic_eset$ID[i] <- intronic_names[as.character(intronic_eset$V4)[i],1]
		
	intronic_eset$lnc_lfc[i] <- lncs_sni_diff_B10.D2_vs_BALB.c[as.character(intronic_eset$ID[i]),2]
	intronic_eset$lnc_pvalue[i] <- lncs_sni_diff_B10.D2_vs_BALB.c[as.character(intronic_eset$ID[i]),6]
	intronic_eset$lnc_baseMean[i] <- lncs_sni_diff_B10.D2_vs_BALB.c[as.character(intronic_eset$ID[i]),1] 
	setTxtProgressBar(pb, i)
	
	}
close(pb)
cat("\n")


both_significant <- subset(intronic_eset, intronic_eset$gene_pvalue<.05 & intronic_eset$lnc_pvalue<.05)

anticorellated <- subset(intronic_eset, (intronic_eset$lnc_lfc * intronic_eset$gene_lfc) < 0)

anticorellated_sig <- subset(both_significant, (both_significant$lnc_lfc * both_significant$gene_lfc) < 0) 

pain_intronic <- subset(intronic_eset, intronic_eset$symbol %in% pain_genes_symbols)

pain_anticorellated <- subset(pain_intronic, (pain_intronic$lnc_lfc * pain_intronic$gene_lfc) < 0)


write.csv(intronic_eset, file = paste(PATH, "/intronic/diff_B10.D2_vs_BALB.c/intronic_all_expression.csv", sep=""))

write.csv(both_significant, file = paste(PATH, "/intronic/diff_B10.D2_vs_BALB.c/intronic_both_sig.csv", sep=""))

write.csv(pain_intronic, file = paste(PATH, "/intronic/diff_B10.D2_vs_BALB.c/pain_intronic.csv", sep=""))

write.csv(anticorellated, file = paste(PATH, "/intronic/diff_B10.D2_vs_BALB.c/intronic_anticorellated.csv", sep=""))

write.csv(anticorellated_sig, file = paste(PATH, "/intronic/diff_B10.D2_vs_BALB.c/intronic_anticorellated_sig.csv", sep=""))

write.csv(pain_anticorellated, file = paste(PATH, "/intronic/diff_B10.D2_vs_BALB.c/pain_anticorellated.csv", sep=""))




