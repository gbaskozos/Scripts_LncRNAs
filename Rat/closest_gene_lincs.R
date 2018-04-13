library("AnnotationDbi")
library("org.Rn.eg.db")
library("Hmisc")

options(stringsAsFactors=FALSE)
PATH <- "/home/george/Desktop/LncRNA_v2/Rat"

load(file = paste(PATH,"/toc_all.RData", sep=""))

load(file = paste(PATH,"/DE/norm_counts.RData", sep=""))

res_snt_vs_sham <- read.csv(file = paste(PATH,"/DE/res_snt_vs_sham.csv", sep=""), row.name=1)

expressed_nc_ensembl <- read.csv(file = paste(PATH,"/DE/expressed_nc_ensembl_toc.csv", sep=""), row.name=1)

lincs <- read.table(paste(PATH, "/non_antisense/lncs_intergenic_closest.bed", sep=""), sep="\t")

lincs_ens <- read.table(paste(PATH, "/non_antisense/lncs_intergenic_closest_ENSEMBL.bed", sep=""), sep="\t")

pain_genes <- read.csv(file = paste(PATH,"/../pain_genes_homologs.csv", sep=""))

pain_genes <- unique(pain_genes[,c(5,4)])

pain_genes_symbols <- unique(as.character(pain_genes$RGD.symbol))


##############################
# ENSEMBL anotated Linc RNAs #
##############################
lincs_closest_ens <- lincs_ens[,c(4,16,25)]

lincs_closest_ens <- lincs_closest_ens[lincs_closest_ens$V4 %in% rownames(expressed_nc_ensembl), ] 

rld_mat_lincs_ens <- rld_mat[rownames(rld) %in% as.character(lincs_closest_ens$V4),] 

rld_mat_genes <- rld_mat[rownames(rld) %in% as.character(lincs_closest_ens$V16),]

names(lincs_closest_ens) <- c("lincRNA", "ENSEMBL_gene_id", "Distance")

lincs_closest_ens$gene_symbol <- mapIds(org.Rn.eg.db, keys=as.character(lincs_closest_ens$ENSEMBL_gene_id), column="SYMBOL", keytype="ENSEMBL", multiVals="first")

lincs_closest_ens$LincRNA_symbol <- mapIds(org.Rn.eg.db, keys=as.character(lincs_closest_ens$lincRNA), column="SYMBOL", keytype="ENSEMBL", multiVals="first")

lincs_closest_ens <- lincs_closest_ens[as.character(lincs_closest_ens$ENSEMBL_gene_id) %in% rownames(rld_mat),]

lincs_closest_ens$baseMean <- res_snt_vs_sham[as.character(lincs_closest_ens$lincRNA),]$baseMean

lincs_closest_ens$lnc_lfc_res_snt_vs_sham <- res_snt_vs_sham[as.character(lincs_closest_ens$lincRNA),]$log2FoldChange

lincs_closest_ens$lnc_pvalue_res_snt_vs_sham <- res_snt_vs_sham[as.character(lincs_closest_ens$lincRNA),]$padj

lincs_closest_ens$gene_lfc_res_snt_vs_sham <- res_snt_vs_sham[as.character(lincs_closest_ens$ENSEMBL_gene_id),]$log2FoldChange

lincs_closest_ens$gene_pvalue_res_snt_vs_sham <- res_snt_vs_sham[as.character(lincs_closest_ens$ENSEMBL_gene_id),]$padj


colnames(rld_mat_genes) <- paste("genes_", colnames(rld_mat_genes), sep="")
colnames(rld_mat_lincs_ens) <- paste("lncs_", colnames(rld_mat_lincs_ens), sep="")

genes_counts <- matrix(nrow=nrow(lincs_closest_ens), ncol = ncol(rld_mat_genes), data = 0)

lincs_counts <- matrix(nrow=nrow(lincs_closest_ens), ncol = ncol(rld_mat_lincs_ens), data = 0)

#rld_mat_genes[lincs_closest_ens$ENSEMBL_gene_id,]


for (i in 1:nrow(lincs_closest_ens)) {

genes_counts[i,] <- rld_mat_genes[as.character(lincs_closest_ens$ENSEMBL_gene_id[i]),]

lincs_counts[i,] <- rld_mat_lincs_ens[as.character(lincs_closest_ens$lincRNA[i]),]
}

lincs_closest_ens <- cbind(lincs_closest_ens, lincs_counts, genes_counts)

names(lincs_closest_ens) <- c(names(lincs_closest_ens)[1:10], colnames(rld_mat_lincs_ens), colnames(rld_mat_genes))

for (i in 1:nrow(lincs_closest_ens)) {
lincs_closest_ens$cor[i] <- cor.test(x=as.numeric(lincs_closest_ens[i,c(11:18)]), y=as.numeric(lincs_closest_ens[i,c(19:26)]))$estimate  

lincs_closest_ens$cor_pvalue[i] <- cor.test(x=as.numeric(lincs_closest_ens[i,c(11:18)]), y=as.numeric(lincs_closest_ens[i,c(19:26)]))$p.value 
}

lincs_closest_ens$cor_pvalue <- p.adjust(lincs_closest_ens$cor_pvalue, method = "BH")

lincs_pain <- lincs_closest_ens[as.character(lincs_closest_ens$ENSEMBL_gene_id) %in% as.character(pain_genes$Gene.ID),]

sig_true <- lincs_closest_ens[lincs_closest_ens$cor_pvalue < 0.05 & !is.na(lincs_closest_ens$cor_pvalue),]

sig_true_high <- lincs_closest_ens[lincs_closest_ens$cor_pvalue < 0.05 & !is.na(lincs_closest_ens$cor_pvalue) & abs(lincs_closest_ens$cor) > 0.8,]


##############
res_snt_vs_sham_sig_ens <- lincs_closest_ens[lincs_closest_ens$lnc_pvalue_res_snt_vs_sham < 0.05 & !is.na(lincs_closest_ens$lnc_pvalue_res_snt_vs_sham),]


pdf(paste (PATH, "/Distance_ensembl_lincRNAs.pdf", sep= ""))
boxplot(abs(lincs_closest_ens$Distance), abs(sig_true$Distance), abs(res_snt_vs_sham_sig_ens$Distance), pars=list(ylim=c(0,1000000)), names=c("All LincRNAs", "Sig. correlated", "DE SNT vs SHAM"), main = "Distance between ENSEMBL genes and LincRNAs" , cex.axis=0.5, las=1)
axis(side=2, at=c(floor(quantile(abs(lincs_closest_ens$Distance))[2]),floor(median(abs(lincs_closest_ens$Distance)))), las=1, cex.axis=0.6, padj=-0.5)

dev.off()

sink(paste (PATH, "/Distance_ENSEMBL_lincs_sig_test.txt", sep= ""))
ks.test(sig_true_high$Distance, abs(lincs_closest_ens$Distance))
sink()

pdf(paste (PATH, "/Distance_ensembl_lincRNAs_cor.pdf", sep= ""))
plot(abs(lincs_closest_ens$Distance), lincs_closest_ens$cor, main="Scatterplot of correlation vs distane", 
  	xlab="Distance (absolute)", ylab="Pearson's correlation coefficient", pch=19)
abline(lm(lincs_closest_ens$cor~abs(lincs_closest_ens$Distance)), col="red") # regression line (y~x) 
legend("topright", title="Correlation", c(paste("R= ",cor.test(abs(lincs_closest_ens$Distance), lincs_closest_ens$cor)$estimate, sep= ""), paste("P.value= ", cor.test(abs(lincs_closest_ens$Distance), lincs_closest_ens$cor)$p.value, sep="")))
dev.off()

lincs_pain_cor <- lincs_pain[lincs_pain$cor_pvalue < 0.05 & !is.na(lincs_pain$cor_pvalue) & abs(lincs_pain$cor) > 0.8,]

res_snt_vs_sham_sig_pain <- lincs_pain[lincs_pain$lnc_pvalue_res_snt_vs_sham < 0.05 & !is.na(lincs_pain$lnc_pvalue_res_snt_vs_sham) & lincs_pain$gene_pvalue_res_snt_vs_sham < 0.05,]


write.csv(file = paste(PATH, "/DE/lincs_ENSEMBL_pain.csv", sep=""), lincs_pain)

write.csv(file = paste(PATH, "/DE/lincs_ENSEMBL_closest.csv", sep=""), lincs_closest_ens)

write.csv(file = paste(PATH, "/DE/lincs_ENSEMBL_pain_cor.csv", sep=""), lincs_pain_cor)

write.csv(file = paste(PATH, "/DE/lincs_ENSEMBL_res_snt_vs_sham_sig.csv", sep=""), res_snt_vs_sham_sig_ens[order(-abs(res_snt_vs_sham_sig_ens$lnc_lfc_res_snt_vs_sham), -res_snt_vs_sham_sig_ens$baseMean),])

write.csv(file = paste(PATH, "/DE/lincs_ENSEMBL_cor.csv", sep=""), sig_true_high)


###################
# Novel Linc RNAs #
###################
lincs_closest <- lincs[,c(4,16,25)]

toc_lncs <- read.csv(file = paste(PATH,"/LncRNAs_toc_rat_DRG.csv", sep=""), row.name=1)

lncs_names <- toc_lncs$name

name_mapping <- cbind(rownames(toc_lncs), as.character(toc_lncs$name))

rownames(name_mapping) <- as.character(toc_lncs$name)

lincs_closest$lincRNA <- name_mapping[as.character(lincs_closest$V4),1]

lincs_closest <- data.frame(lincRNA = lincs_closest$lincRNA, ENSEMBL_gene_id = lincs_closest$V16, Distance = lincs_closest$V25, gene_symbol = mapIds(org.Rn.eg.db, keys=as.character(lincs_closest$V16), column="SYMBOL", keytype="ENSEMBL", multiVals="first"), LincRNA_symbol = lincs_closest$V4)

lincs_closest <- lincs_closest[as.character(lincs_closest$lincRNA) %in% rownames(rld_mat),]

rld_mat_lincs <- rld_mat[rownames(rld) %in% as.character(lincs_closest$lincRNA),] 

rld_mat_genes <- rld_mat[rownames(rld) %in% as.character(lincs_closest$ENSEMBL_gene_id),]

lincs_closest$baseMean <- res_snt_vs_sham[as.character(lincs_closest$lincRNA),]$baseMean

lincs_closest$lnc_lfc_res_snt_vs_sham <- res_snt_vs_sham[as.character(lincs_closest$lincRNA),]$log2FoldChange

lincs_closest$lnc_pvalue_res_snt_vs_sham <- res_snt_vs_sham[as.character(lincs_closest$lincRNA),]$padj

lincs_closest$gene_lfc_res_snt_vs_sham <- res_snt_vs_sham[as.character(lincs_closest$ENSEMBL_gene_id),]$log2FoldChange

lincs_closest$gene_pvalue_res_snt_vs_sham <- res_snt_vs_sham[as.character(lincs_closest$ENSEMBL_gene_id),]$padj


colnames(rld_mat_genes) <- paste("genes_", colnames(rld_mat_genes), sep="")
colnames(rld_mat_lincs) <- paste("lncs_", colnames(rld_mat_lincs), sep="")

genes_counts <- matrix(nrow=nrow(lincs_closest), ncol = ncol(rld_mat_genes), data = 0)

lincs_counts <- matrix(nrow=nrow(lincs_closest), ncol = ncol(rld_mat_lincs), data = 0)


for (i in 1:nrow(lincs_closest)) {

if (as.character(lincs_closest$ENSEMBL_gene_id[i] %in% rownames(rld_mat_genes))) {
genes_counts[i,] <- rld_mat_genes[as.character(lincs_closest$ENSEMBL_gene_id[i]),]
}
lincs_counts[i,] <- rld_mat_lincs[as.character(lincs_closest$lincRNA[i]),]
}

lincs_closest <- cbind(lincs_closest, lincs_counts, genes_counts)

names(lincs_closest) <- c(names(lincs_closest)[1:10], colnames(rld_mat_lincs), colnames(rld_mat_genes))




for (i in 1:nrow(lincs_closest)) {
lincs_closest$cor[i] <- cor.test(x=as.numeric(lincs_closest[i,c(11:18)]), y=as.numeric(lincs_closest[i,c(19:26)]))$estimate  

lincs_closest$cor_pvalue[i] <- cor.test(x=as.numeric(lincs_closest[i,c(11:18)]), y=as.numeric(lincs_closest[i,c(19:26)]))$p.value 
}

lincs_closest$cor_pvalue <- p.adjust(lincs_closest$cor_pvalue, method = "BH")

lincs_pain <- lincs_closest[as.character(lincs_closest$ENSEMBL_gene_id) %in% as.character(pain_genes$Gene.ID),]


sig_true <- lincs_closest[lincs_closest$cor_pvalue < 0.05 & !is.na(lincs_closest$cor_pvalue),]

sig_true_high <- lincs_closest[lincs_closest$cor_pvalue < 0.05 & !is.na(lincs_closest$cor_pvalue) & abs(lincs_closest$cor) > 0.8,]


##############
res_snt_vs_sham_sig <- lincs_closest[lincs_closest$lnc_pvalue_res_snt_vs_sham < 0.05 & !is.na(lincs_closest$lnc_pvalue_res_snt_vs_sham),]


pdf(paste (PATH, "/Distance_novel_lincRNAs.pdf", sep= ""))
boxplot(abs(lincs_closest$Distance), abs(sig_true$Distance), abs(res_snt_vs_sham_sig$Distance),  pars=list(ylim=c(0,1000000)), names=c("All LincRNAs", "Sig. correlated", "DE SNT vs SHAM"), main = "Distance between ENSEMBL genes and LincRNAs" , cex.axis=0.5, las=1)
axis(side=2, at=c(floor(quantile(abs(lincs_closest$Distance))[2]),floor(median(abs(lincs_closest$Distance)))), las=1, cex.axis=0.6, padj=-0.5)

dev.off()

sink(paste (PATH, "/Distance_novel_lincs_sig_test.txt", sep= ""))
ks.test(sig_true_high$Distance, abs(lincs_closest$Distance))
sink()

pdf(paste (PATH, "/Distance_novel_lincRNAs_cor.pdf", sep= ""))
plot(abs(lincs_closest$Distance), lincs_closest$cor, main="Scatterplot of correlation vs distane", 
  	xlab="Distance (absolute)", ylab="Pearson's correlation coefficient", pch=19)
abline(lm(lincs_closest$cor~abs(lincs_closest$Distance)), col="red") # regression line (y~x) 
legend("topright", title="Correlation", c(paste("R= ",cor.test(abs(lincs_closest$Distance), lincs_closest$cor)$estimate, sep= ""), paste("P.value= ", cor.test(abs(lincs_closest$Distance), lincs_closest$cor)$p.value, sep="")))
dev.off()

lincs_pain_cor <- lincs_pain[lincs_pain$cor_pvalue < 0.05 & !is.na(lincs_pain$cor_pvalue) & abs(lincs_pain$cor) > 0.8,]

res_snt_vs_sham_sig_pain <- lincs_pain[lincs_pain$lnc_pvalue_res_snt_vs_sham < 0.05 & !is.na(lincs_pain$lnc_pvalue_res_snt_vs_sham) & lincs_pain$gene_pvalue_res_snt_vs_sham < 0.05,]

write.csv(file = paste(PATH, "/DE/lincs_pain.csv", sep=""), lincs_pain)

write.csv(file = paste(PATH, "/DE/lincs_closest.csv", sep=""), lincs_closest)

write.csv(file = paste(PATH, "/DE/lincs_pain_cor.csv", sep=""), lincs_pain_cor)

write.csv(file = paste(PATH, "/DE/lincs_res_snt_vs_sham_sig.csv", sep=""), res_snt_vs_sham_sig)

write.csv(file = paste(PATH, "/DE/lincs_cor.csv", sep=""), sig_true_high)

