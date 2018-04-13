library("GenomicRanges")
library("DESeq2")
library("gplots")
library("ggplot2")
library("RColorBrewer")
library(genefilter)
library("vsn")
library("BiocParallel")
library("AnnotationDbi")
library("org.Hs.eg.db")
library("sva")
library(Tmisc)
library(calibrate)
library(lattice)
library(plotrix)
library(WGCNA)

n_cores <- detectCores() - 1

multicoreParam <- MulticoreParam(workers = n_cores, progressbar = TRUE)
options(mc.cores=n_cores)
register(multicoreParam)

allowWGCNAThreads(n_cores)

options(stringsAsFactors=FALSE)

PATH <- "/home/george/Desktop/LncRNA_v2/IPS"

load(paste(PATH,"/LncRNAs_toc_IPS.RData", sep=""))

load(paste(PATH,"/ens_genes_toc_all.RData", sep=""))

toc_genes <- toc

colData <- data.frame(condition = c(rep("IPS",6), rep("NEURONS",6)), cell_line = rep(c(rep("AD2",2), rep("AD4",2), rep("NHDF",2)),2))

colData$condition <- as.factor(colData$condition)

colData$cell_line <- as.factor(colData$cell_line)

name_mapping <- data.frame(Lnc_id = rownames(LncRNAs_toc), Lnc_name = LncRNAs_toc$name)

toc_lncs <- LncRNAs_toc[, -13]

toc <- rbind(toc_genes, toc_lncs)

save(toc, file = paste(PATH,"/toc_all.RData", sep=""))

########
#DESeq2#
########

load(file = paste(PATH,"/toc_all.RData", sep=""))

dds <- DESeqDataSetFromMatrix(countData = toc, colData = colData, design = ~ cell_line*condition)

dds$condition <- relevel(dds$condition,"IPS")
dds$cell_line <- relevel(dds$cell_line,"AD2")

dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)

norm_counts <- counts(dds, normalized=TRUE)
colnames(norm_counts) <- names(toc)

####################
# Sample distances #
####################
rld <- rlog(dds, blind=FALSE)
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
rld_mat <- assay(rld)
colnames(rld_mat) <- names(toc)

save(rld_mat, rld, norm_counts, file = paste(PATH,"/DE/norm_counts.RData", sep=""))

hmcol <- colorRampPalette(brewer.pal(9,"GnBu"))(100)

distsRL <- dist(t(rld_mat))
mat <- as.matrix(distsRL)
rownames(mat) <- colnames(mat) <- names(toc)
hc <- hclust(distsRL)

pdf(file = paste(PATH,"/DE/heatmap_all_genes.pdf", sep=""))
heatmap.2(mat, col = rev(hmcol), Rowv=as.dendrogram(hc), trace="none", symm=TRUE, margin=c(13,13))
dev.off()


pdf(file = paste(PATH,"/DE/pca_all_genes.pdf", sep=""))
d = plotPCA(rld, intgroup=c("condition", "cell_line"), returnData=TRUE, ntop = 10000)
percentVar <- round(100 * attr(d, "percentVar"))
ggplot(d, aes(x=PC1,y=PC2, color=condition, shape=cell_line, label=names(toc))) + geom_point(size=10) + xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance")) + ggtitle("PCA Neurons vs IPS") + theme(axis.text.y = element_text(size= 17, face="bold"), axis.title.y = element_text(size=17, face="bold"), axis.title.x = element_text(size=17, face="bold"), axis.text.x = element_text(size= 14, angle=45, hjust=1, face="bold"), legend.title=NULL, legend.text=element_text(size=17, face="bold"), legend.key.size=unit(0.5, "cm"), plot.title=element_text(size=25, face="bold", hjust=0.5))
dev.off()


#####################
dds <- nbinomWaldTest(dds)


resultsNames(dds)


# the condition effect for AD2 (the main effect)
res_NEURONS_vs_IPS_AD2 <- results(dds, name = "condition_NEURONS_vs_IPS")


# the condition effect for AD4
# this is, by definition, the main effect *plus* the interaction term
# (the extra condition effect in AD4 compared to AD2).
res_NEURONS_vs_IPS_AD4 <- results(dds, list( c("condition_NEURONS_vs_IPS","cell_lineAD4.conditionNEURONS")))


res_NEURONS_vs_IPS_NHDF <- results(dds, list( c("condition_NEURONS_vs_IPS","cell_lineNHDF.conditionNEURONS")))


# the interaction term, answering: is the condition effect *different* across strains?
diff_AD4_vs_AD2 <- results(dds, name="cell_lineAD4.conditionNEURONS")

diff_NHDF_vs_AD2 <- results(dds, name="cell_lineNHDF.conditionNEURONS")

# cell line
AD4_vs_AD2 <- results(dds, name="cell_line_AD4_vs_AD2")

NHDF_vs_AD2 <- results(dds, name="cell_line_NHDF_vs_AD2")


rownames(name_mapping) <- name_mapping$Lnc_id

res_NEURONS_vs_IPS_AD2$symbol <- mapIds(org.Hs.eg.db, keys=row.names(res_NEURONS_vs_IPS_AD2), column="SYMBOL", keytype="ENSEMBL", multiVals="first")
res_NEURONS_vs_IPS_AD4$symbol <- mapIds(org.Hs.eg.db, keys=row.names(res_NEURONS_vs_IPS_AD4), column="SYMBOL", keytype="ENSEMBL", multiVals="first")
res_NEURONS_vs_IPS_NHDF$symbol <- mapIds(org.Hs.eg.db, keys=row.names(res_NEURONS_vs_IPS_NHDF), column="SYMBOL", keytype="ENSEMBL", multiVals="first")

diff_AD4_vs_AD2$symbol <- mapIds(org.Hs.eg.db, keys=row.names(diff_AD4_vs_AD2), column="SYMBOL", keytype="ENSEMBL", multiVals="first")
diff_NHDF_vs_AD2$symbol <- mapIds(org.Hs.eg.db, keys=row.names(diff_NHDF_vs_AD2), column="SYMBOL", keytype="ENSEMBL", multiVals="first")

AD4_vs_AD2$symbol <- mapIds(org.Hs.eg.db, keys=row.names(AD4_vs_AD2), column="SYMBOL", keytype="ENSEMBL", multiVals="first")
NHDF_vs_AD2$symbol <- mapIds(org.Hs.eg.db, keys=row.names(NHDF_vs_AD2), column="SYMBOL", keytype="ENSEMBL", multiVals="first")

res_NEURONS_vs_IPS_AD2[grep("LncRNA", rownames(res_NEURONS_vs_IPS_AD2)),]$symbol <- name_mapping[as.character(rownames(res_NEURONS_vs_IPS_AD2[grep("LncRNA", rownames(res_NEURONS_vs_IPS_AD2)),])),]$Lnc_name
res_NEURONS_vs_IPS_AD4[grep("LncRNA", rownames(res_NEURONS_vs_IPS_AD4)),]$symbol <- name_mapping[as.character(rownames(res_NEURONS_vs_IPS_AD4[grep("LncRNA", rownames(res_NEURONS_vs_IPS_AD4)),])),]$Lnc_name
res_NEURONS_vs_IPS_NHDF[grep("LncRNA", rownames(res_NEURONS_vs_IPS_NHDF)),]$symbol <- name_mapping[as.character(rownames(res_NEURONS_vs_IPS_NHDF[grep("LncRNA", rownames(res_NEURONS_vs_IPS_NHDF)),])),]$Lnc_name


sig_AD2 <- subset(res_NEURONS_vs_IPS_AD2, padj<.05)

sig_AD4 <- subset(res_NEURONS_vs_IPS_AD4, padj<.05)

sig_NHDF <- subset(res_NEURONS_vs_IPS_NHDF, padj<.05)

hsig_AD2 <- subset(res_NEURONS_vs_IPS_AD2, padj<.05 & abs(log2FoldChange)>1)

hsig_AD4 <- subset(res_NEURONS_vs_IPS_AD4, padj<.05 & abs(log2FoldChange)>1)

hsig_NHDF <- subset(res_NEURONS_vs_IPS_NHDF, padj<.05 & abs(log2FoldChange)>1)


write.csv(file = paste(PATH,"/DE/res_NEURONS_vs_IPS_AD2.csv", sep=""),res_NEURONS_vs_IPS_AD2[order(-abs(res_NEURONS_vs_IPS_AD2$log2FoldChange), res_NEURONS_vs_IPS_AD2$padj),])
write.csv(file = paste(PATH,"/DE/res_NEURONS_vs_IPS_AD4.csv", sep=""),res_NEURONS_vs_IPS_AD4[order(-abs(res_NEURONS_vs_IPS_AD4$log2FoldChange), res_NEURONS_vs_IPS_AD4$padj),])
write.csv(file = paste(PATH,"/DE/res_NEURONS_vs_IPS_NHDF.csv", sep=""),res_NEURONS_vs_IPS_NHDF[order(-abs(res_NEURONS_vs_IPS_NHDF$log2FoldChange), res_NEURONS_vs_IPS_NHDF$padj),])

write.csv(file = paste(PATH,"/DE/sig_AD2.csv", sep=""), sig_AD2[order(-abs(sig_AD2$log2FoldChange), sig_AD2$padj),])
write.csv(file = paste(PATH,"/DE/sig_AD4.csv", sep=""), sig_AD4[order(-abs(sig_AD4$log2FoldChange), sig_AD4$padj),])
write.csv(file = paste(PATH,"/DE/sig_NHDF.csv", sep=""), sig_NHDF[order(-abs(sig_NHDF$log2FoldChange), sig_NHDF$padj),])

write.csv(file = paste(PATH,"/DE/hsig_AD2.csv", sep=""), hsig_AD2[order(-abs(hsig_AD2$log2FoldChange), hsig_AD2$padj),])
write.csv(file = paste(PATH,"/DE/hsig_AD4.csv", sep=""), hsig_AD4[order(-abs(hsig_AD4$log2FoldChange), hsig_AD4$padj),])
write.csv(file = paste(PATH,"/DE/hsig_NHDF.csv", sep=""), hsig_NHDF[order(-abs(hsig_NHDF$log2FoldChange), hsig_NHDF$padj),])



