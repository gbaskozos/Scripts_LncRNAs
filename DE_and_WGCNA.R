library("GenomicRanges")
library("DESeq2")
library("gplots")
library("ggplot2")
library("RColorBrewer")
library(genefilter)
library("vsn")
library("BiocParallel")
library("AnnotationDbi")
library("org.Mm.eg.db")
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

PATH <- "/home/george/Desktop/LncRNA_v2/Mouse"

load(paste(PATH,"/LncRNAs_toc_mouse_DRG.RData", sep=""))

load(paste(PATH,"/ens_genes_toc_all.RData", sep=""))

toc_genes <- toc

colData <- data.frame(strain = c(rep("BALB.c",5), rep("B10.D2",3), "BALB.c", rep("B10.D2",2), rep("BALB.c",3), "B10.D2", "BALB.c", rep("B10.D2",4)), sex = c("M", "M", rep("F",3), rep("M",2), "F", "M", "F", "M", "M", rep("F",3), "M", rep("F",2), "M", "M"), condition = c(rep("SHAM", 11), rep("SNI",9)))

colData$strain <- as.factor(colData$strain)

colData$sex <- as.factor(colData$sex)

colData$condition <- as.factor(colData$condition)

colData$group <- factor(paste(colData$strain, colData$condition, sep="_"))

name_mapping <- data.frame(Lnc_id = rownames(LncRNAs_toc), Lnc_name = LncRNAs_toc$name)

toc_lncs <- LncRNAs_toc[, -21]

toc <- rbind(toc_genes, toc_lncs)

save(toc, file = paste(PATH,"/toc_all.RData", sep=""))

########
#DESeq2#
########

load(file = paste(PATH,"/toc_all.RData", sep=""))

dds <- DESeqDataSetFromMatrix(countData = toc, colData = colData, design = ~ sex + strain*condition)

dds$condition <- relevel(dds$condition,"SHAM")
dds$strain <- relevel(dds$strain,"BALB.c")

dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)

norm_counts <- counts(dds, normalized=TRUE)
colnames(norm_counts) <- names(toc)

####################
# Sample distances #
####################
rld <- rlog(dds, blind=TRUE)
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
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
d = plotPCA(rld, intgroup=c("condition", "strain"), returnData=TRUE, ntop = 10000)
percentVar <- round(100 * attr(d, "percentVar"))
ggplot(d, aes(x=PC1,y=PC2, color=condition, shape=strain, label=names(toc))) + geom_point(size=10) + xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance")) + ggtitle("PCA SNI vs Sham mouse") + theme(axis.text.y = element_text(size= 17, face="bold"), axis.title.y = element_text(size=17, face="bold"), axis.title.x = element_text(size=17, face="bold"), axis.text.x = element_text(size= 14, angle=45, hjust=1, face="bold"), legend.title=NULL, legend.text=element_text(size=17, face="bold"), legend.key.size=unit(0.5, "cm"), plot.title=element_text(size=25, face="bold", hjust=0.5))
dev.off()


#####################
dds <- nbinomWaldTest(dds)


resultsNames(dds)


# the condition effect for BALB.c (the main effect)
res_sham_vs_sni_BALB.c <- results(dds, name = "condition_SNI_vs_SHAM")


# the condition effect for genotype B10.D2
# this is, by definition, the main effect *plus* the interaction term
# (the extra condition effect in B10.D2 compared to Balb.c).
res_sham_vs_sni_B10.D2 <- results(dds, list( c("condition_SNI_vs_SHAM","strainB10.D2.conditionSNI")))

# the interaction term, answering: is the condition effect *different* across strains?
res_sni_diff_B10.D2_vs_BALB.c <- results(dds, name="strainB10.D2.conditionSNI")

# the sex effect for genotype BALB.c
res_sex <- results(dds, contrast=c("sex", "F", "M"))

# the strain effect in sham condition
res_strain <- results(dds, name="strain_B10.D2_vs_BALB.c")



rownames(name_mapping) <- name_mapping$Lnc_id

res_sham_vs_sni_B10.D2$symbol <- mapIds(org.Mm.eg.db, keys=row.names(res_sham_vs_sni_B10.D2), column="SYMBOL", keytype="ENSEMBL", multiVals="first")
res_sham_vs_sni_BALB.c$symbol <- mapIds(org.Mm.eg.db, keys=row.names(res_sham_vs_sni_BALB.c), column="SYMBOL", keytype="ENSEMBL", multiVals="first")
res_sni_diff_B10.D2_vs_BALB.c$symbol <- mapIds(org.Mm.eg.db, keys=row.names(res_sni_diff_B10.D2_vs_BALB.c), column="SYMBOL", keytype="ENSEMBL", multiVals="first")
res_sex$symbol <- mapIds(org.Mm.eg.db, keys=row.names(res_sex), column="SYMBOL", keytype="ENSEMBL", multiVals="first")
res_strain$symbol <- mapIds(org.Mm.eg.db, keys=row.names(res_strain), column="SYMBOL", keytype="ENSEMBL", multiVals="first")


res_sham_vs_sni_B10.D2[grep("LncRNA", rownames(res_sham_vs_sni_B10.D2)),]$symbol <- name_mapping[as.character(rownames(res_sham_vs_sni_B10.D2[grep("LncRNA", rownames(res_sham_vs_sni_B10.D2)),])),]$Lnc_name
res_sham_vs_sni_BALB.c[grep("LncRNA", rownames(res_sham_vs_sni_BALB.c)),]$symbol <- name_mapping[as.character(rownames(res_sham_vs_sni_BALB.c[grep("LncRNA", rownames(res_sham_vs_sni_BALB.c)),])),]$Lnc_name
res_sni_diff_B10.D2_vs_BALB.c[grep("LncRNA", rownames(res_sni_diff_B10.D2_vs_BALB.c)),]$symbol <- name_mapping[as.character(rownames(res_sni_diff_B10.D2_vs_BALB.c[grep("LncRNA", rownames(res_sni_diff_B10.D2_vs_BALB.c)),])),]$Lnc_name
res_sex[grep("LncRNA", rownames(res_sex)),]$symbol <- name_mapping[as.character(rownames(res_sex[grep("LncRNA", rownames(res_sex)),])),]$Lnc_name
res_strain[grep("LncRNA", rownames(res_strain)),]$symbol <- name_mapping[as.character(rownames(res_strain[grep("LncRNA", rownames(res_strain)),])),]$Lnc_name

sig_b10d2 <- subset(res_sham_vs_sni_B10.D2, padj<.05)

sig_balbc <- subset(res_sham_vs_sni_BALB.c, padj<.05)

hsig_b10d2 <- subset(res_sham_vs_sni_B10.D2, padj<.05 & abs(log2FoldChange)>1)

hsig_balbc <- subset(res_sham_vs_sni_BALB.c, padj<.05 & abs(log2FoldChange)>1)

write.csv(file = paste(PATH,"/DE/res_sni_vs_sham_B10.D2.csv", sep=""),res_sham_vs_sni_B10.D2[order(-abs(res_sham_vs_sni_B10.D2$log2FoldChange), res_sham_vs_sni_B10.D2$padj),])
write.csv(file = paste(PATH,"/DE/res_sni_vs_sham_BALB.c.csv", sep=""),res_sham_vs_sni_BALB.c[order(-abs(res_sham_vs_sni_BALB.c$log2FoldChange), res_sham_vs_sni_BALB.c$padj),])
write.csv(file = paste(PATH,"/DE/res_sni_diff_B10.D2_vs_BALB.c.csv", sep=""),res_sni_diff_B10.D2_vs_BALB.c[order(-abs(res_sni_diff_B10.D2_vs_BALB.c$log2FoldChange), res_sni_diff_B10.D2_vs_BALB.c$padj),])
write.csv(file = paste(PATH,"/DE/res_sex.csv", sep=""),res_sex[order(abs(res_sex$log2FoldChange)),])
write.csv(file = paste(PATH,"/DE/res_strain.csv", sep=""),res_strain[order(abs(res_strain$log2FoldChange)),])

write.csv(file = paste(PATH,"/DE/sig_sni_vs_sham_B10.D2.csv", sep=""), sig_b10d2[order(-abs(sig_b10d2$log2FoldChange), sig_b10d2$padj),])
write.csv(file = paste(PATH,"/DE/sig_sni_vs_sham_BALB.c.csv", sep=""), sig_balbc[order(-abs(sig_balbc$log2FoldChange), sig_balbc$padj),])

write.csv(file = paste(PATH,"/DE/highly_sig_sni_vs_sham_B10.D2.csv", sep=""), hsig_b10d2[order(-abs(hsig_b10d2$log2FoldChange), hsig_b10d2$padj),])
write.csv(file = paste(PATH,"/DE/highly_sig_sni_vs_sham_BALB.c.csv", sep=""), hsig_balbc[order(-abs(hsig_balbc$log2FoldChange), hsig_balbc$padj),])

###########
## WGCNA ##
###########
ffun=filterfun(pOverA(p = 0.25, A = 10))
filt=genefilter(toc,ffun)
toc_coff=toc[filt,]

dds_coff <- DESeqDataSetFromMatrix(countData = toc_coff, colData = colData, design = ~ sex + strain*condition)

dds_coff$condition <- relevel(dds_coff$condition,"SHAM")
dds_coff$strain <- relevel(dds_coff$strain,"BALB.c")

dds_coff <- estimateSizeFactors(dds_coff)
dds_coff <- estimateDispersions(dds_coff)

vsd_coff <- varianceStabilizingTransformation(dds_coff, blind=FALSE)
vsd_mat_coff <- assay(vsd_coff)
colnames(vsd_mat_coff) <- names(toc_coff)

save(vsd_mat_coff, vsd_coff, file = paste(PATH,"/Network/norm_counts_coff.RData", sep=""))

data_matrix <- vsd_mat_coff

WGCNA_matrix <- t(data_matrix[order(apply(data_matrix,1,mad), decreasing=T)[1:(nrow(data_matrix)/4)],])

s = abs(bicor(WGCNA_matrix))

powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(WGCNA_matrix, powerVector = powers, verbose = 5)

pdf(file = paste(PATH,"/Network/pick_soft_threshold_b.pdf", sep=""))
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab='Soft Threshold (power)',ylab='Scale Free Topology Model Fit,signed R^2',
     type='n', main = paste('Scale independence'));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=1,col='red'); abline(h=0.90,col='red')
dev.off()

gc()

beta <- 5
a <- s^beta
w <- 1-a

#create gene tree by average linkage hierarchical clustering 
geneTree = hclust(as.dist(w), method = 'average')

#module identification using dynamic tree cut algorithm
modules = cutreeDynamic(dendro = geneTree, distM = w, deepSplit = 2, pamRespectsDendro = FALSE, cutHeight = 0.995, minClusterSize = 30)

#assign module colours
module.colours = labels2colors(modules)


library(ape)
#calculate eigengenes
MEList = moduleEigengenes(WGCNA_matrix, colors = module.colours)
MEs = MEList$eigengenes

#calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);

#cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = 'average')

MEDissThres = 0.2
# Plot the result
pdf(file = paste(PATH,"/Network/unmerged_modules_clustering.pdf", sep=""))
plot(METree, main = "Clustering of module eigengenes",
xlab = "", sub = "")
abline(h=MEDissThres, col = "red")
dev.off()

save(MEs, MEList, METree, modules, module.colours, geneTree, file = paste(PATH, "/Network/WGCNA.RData", sep=""))

#Merging modules
merge = mergeCloseModules(WGCNA_matrix, module.colours, cutHeight = MEDissThres, verbose = 3)
merged.colours = merge$colors;
mergedMEs = merge$newMEs;

pdf(file = paste(PATH, "/Network/Gene_dendro_merged_modules.pdf", sep=""))
plotDendroAndColors(geneTree, cbind(module.colours, merged.colours),
c("Dynamic Tree Cut", "Merged dynamic"),
dendroLabels = FALSE, hang = 0.03,
addGuide = TRUE, guideHang = 0.05)
dev.off()

module.colours = merged.colours
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(100));
moduleLabels = match(module.colours, colorOrder)-1;
MEs = mergedMEs;

module_annotation <- data.frame(moduleLabels = unique(moduleLabels), module.colours = unique(module.colours))

save(MEs, moduleLabels, module.colours, geneTree, module_annotation, file = paste(PATH, "/Network/merged_modules.RData", sep=""))

#reference genes = all ENSEMBL top 1/4 MAD genes 
ref_genes <- grep("ENS", colnames(WGCNA_matrix), value=TRUE)
ref_matrix <- WGCNA_matrix[,grep("ENS", colnames(WGCNA_matrix))]
ref_labels <- moduleLabels[grep("ENS", colnames(WGCNA_matrix))]
ref_colours <- module.colours[grep("ENS", colnames(WGCNA_matrix))]

#create data frame for GO analysis
library(org.Mm.eg.db)
GO = toTable(org.Mm.egGO); ENSEMBL = toTable(org.Mm.egENSEMBL)
GO_data_frame = data.frame(GO$go_id, GO$Evidence, ENSEMBL$ensembl_id[match(GO$gene_id,ENSEMBL$gene_id)])

GO_data_frame <- na.omit(GO_data_frame)

#create GOAllFrame object
library(AnnotationDbi)
GO_ALLFrame = GOAllFrame(GOFrame(GO_data_frame, organism = 'Mus musculus'))

#create gene set
library(GSEABase)
gsc <- GeneSetCollection(GO_ALLFrame, setType = GOCollection())

#perform GO enrichment analysis and save results to list - this make take several minutes
library(GOstats)
GSEAGO = vector('list',length(unique(ref_labels)))
for(i in length(unique(ref_labels))){
  GSEAGO[[i]] = summary(hyperGTest(GSEAGOHyperGParams(name = 'Mus musculus GO', 
              geneSetCollection = gsc, geneIds = colnames(ref_matrix)[ref_labels==unique(ref_labels)[[i]]], 
              universeGeneIds = ref_genes, ontology = 'BP', pvalueCutoff = 0.05, 
              conditional = FALSE, testDirection = 'over')))
  print(i)
print(paste("Module:", unique(ref_labels)[[i]]))
}

ref_module_annotation <- data.frame(moduleLabels = unique(ref_labels), module.colours = unique(ref_colours))

cutoff_size = 100

GO_module_name = rep(NA,length(unique(unique(ref_labels))))
for (i in 1:length(unique(unique(ref_labels)))){
  GO_module_name[i] = 
    GSEAGO[[i]][GSEAGO[[i]]$Size<cutoff_size,
    ][which(GSEAGO[[i]][GSEAGO[[i]]$Size<cutoff_size,]$Count==max(GSEAGO[[i]][GSEAGO[[i]]$
    Size<cutoff_size,]$Count)),7]
}

module_annotation <- ref_module_annotation

module_annotation$GO_module  <- GO_module_name

save(gsc, GO_module_name, GSEAGO, module_annotation, file = paste(PATH,"/Network/GO_modules.pdf", sep=""))

#LncRNAs in top 1/4 MAD genes 
LncRNA_matrix <- WGCNA_matrix[, grep("LncRNA", colnames(WGCNA_matrix))]
LncRNA_labels <- moduleLabels[grep("LncRNA", colnames(WGCNA_matrix))]
LncRNA_colours <- module.colours[grep("LncRNA", colnames(WGCNA_matrix))]


module_annotation <- module_annotation[order(module_annotation$moduleLabels), ]

rownames(module_annotation) <- module_annotation$moduleLabels

write.csv(module_annotation, file = paste(PATH,"/Network/module_annotation.csv", sep=""))

GO_lncRNAs <- module_annotation[as.character(LncRNA_labels),]$GO_module

LncRNAs_GO_module <- data.frame(LncRNA = colnames(LncRNA_matrix), GO_module = GO_lncRNAs, module.label = LncRNA_labels, module.colour = LncRNA_colours)

LncRNAs_GO_module$module.membership_weight <- NA

MM <- abs(bicor(LncRNA_matrix, MEs))

for(i in 1:nrow(LncRNAs_GO_module)) {

LncRNAs_GO_module$module.membership_weight[i] <- MM[LncRNAs_GO_module$LncRNA[i], paste0("ME",LncRNAs_GO_module$module.colour[i])]

}

write.csv(LncRNAs_GO_module, file = paste(PATH,"/Network/LncRNAs_GO_module.csv", sep=""))

##DO THIS FOR ENSEMBL LNCS AS WELL

