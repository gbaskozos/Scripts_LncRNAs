library("GenomicRanges")
library("DESeq2")
library("gplots")
library("ggplot2")
library("RColorBrewer")
library(genefilter)
library("vsn")
library("BiocParallel")
library("AnnotationDbi")
library("org.Rn.eg.db")
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

PATH <- "/home/george/Desktop/LncRNA_v2/Rat"

load(paste(PATH,"/LncRNAs_toc_rat_DRG.RData", sep=""))

load(paste(PATH,"/ens_genes_toc_all.RData", sep=""))

toc_genes <- toc

colData <- data.frame(condition = c(rep("SHAM",4), rep("D21_SNT",4)))

colData$condition <- as.factor(colData$condition)

name_mapping <- data.frame(Lnc_id = rownames(LncRNAs_toc), Lnc_name = LncRNAs_toc$name)

toc_lncs <- LncRNAs_toc[, -9]

toc <- rbind(toc_genes, toc_lncs)

save(toc, file = paste(PATH,"/toc_all.RData", sep=""))

########
#DESeq2#
########

load(file = paste(PATH,"/toc_all.RData", sep=""))

dds <- DESeqDataSetFromMatrix(countData = toc, colData = colData, design = ~ condition)

dds$condition <- relevel(dds$condition,"SHAM")

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
rownames(mat) <- colnames(mat) <- paste(names(toc), colData$condition, sep="")
hc <- hclust(distsRL)

pdf(file = paste(PATH,"/DE/heatmap_all_genes.pdf", sep=""))
heatmap.2(mat, col = rev(hmcol), Rowv=as.dendrogram(hc), trace="none", symm=TRUE, margin=c(13,13), cexRow=0.8, cexCol=0.8)
dev.off()


pdf(file = paste(PATH,"/DE/pca_all_genes.pdf", sep=""))
d = plotPCA(rld, intgroup=c("condition"), returnData=TRUE, ntop = 10000)
percentVar <- round(100 * attr(d, "percentVar"))
ggplot(d, aes(x=PC1,y=PC2, color=condition, label=names(toc))) + geom_point(size=10) + xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance")) + ggtitle("PCA SNT vs Sham rat") + theme(axis.text.y = element_text(size= 17, face="bold"), axis.title.y = element_text(size=17, face="bold"), axis.title.x = element_text(size=17, face="bold"), axis.text.x = element_text(size= 14, angle=45, hjust=1, face="bold"), legend.title=NULL, legend.text=element_text(size=17, face="bold"), legend.key.size=unit(0.5, "cm"), plot.title=element_text(size=25, face="bold", hjust=0.5))
dev.off()


#####################
dds <- nbinomWaldTest(dds)


resultsNames(dds)


# the condition effect
res_snt_vs_sham <- results(dds, name = "condition_D21_SNT_vs_SHAM")

rownames(name_mapping) <- name_mapping$Lnc_id

res_snt_vs_sham$symbol <- mapIds(org.Rn.eg.db, keys=row.names(res_snt_vs_sham), column="SYMBOL", keytype="ENSEMBL", multiVals="first")

res_snt_vs_sham[grep("LncRNA", rownames(res_snt_vs_sham)),]$symbol <- name_mapping[as.character(rownames(res_snt_vs_sham[grep("LncRNA", rownames(res_snt_vs_sham)),])),]$Lnc_name

sig_snt_vs_sham <- subset(res_snt_vs_sham, padj<.05)

hsig_snt_vs_sham <- subset(res_snt_vs_sham, padj<.05 & abs(log2FoldChange)>1)


write.csv(file = paste(PATH,"/DE/res_snt_vs_sham.csv", sep=""),res_snt_vs_sham[order(-abs(res_snt_vs_sham$log2FoldChange), res_snt_vs_sham$padj),])

write.csv(file = paste(PATH,"/DE/sig_snt_vs_sham.csv", sep=""), sig_snt_vs_sham[order(-abs(sig_snt_vs_sham$log2FoldChange), sig_snt_vs_sham$padj),])

write.csv(file = paste(PATH,"/DE/highly_snt_vs_sham.csv", sep=""), hsig_snt_vs_sham[order(-abs(hsig_snt_vs_sham$log2FoldChange), hsig_snt_vs_sham$padj),])

###########
## WGCNA ##
###########
## Not enough samples (n=8) to do WGCNA 

