library(biomaRt)
library("org.Mm.eg.db")
library("org.Rn.eg.db")
library("org.Hs.eg.db")

system("mkdir /home/george/Desktop/LncRNA_v2/homology")


pain_genes <- read.csv("/home/george/Desktop/LncRNA_v2/pain_genes_homologs.csv")

mouse <- useMart("ensembl",dataset="mmusculus_gene_ensembl")
rat <- useMart("ensembl",dataset="rnorvegicus_gene_ensembl")
human <- useMart("ensembl",dataset="hsapiens_gene_ensembl")

sni_vs_sham_BALB.c <- read.csv(file = "/home/george/Desktop/LncRNA_v2/Mouse/DE/res_sni_vs_sham_BALB.c.csv", row.name=1)
sni_vs_sham_B10.D2 <- read.csv(file = "/home/george/Desktop/LncRNA_v2/Mouse/DE/res_sni_vs_sham_B10.D2.csv", row.name=1)
sni_diff_B10.D2_vs_BALB.c <- read.csv(file = "/home/george/Desktop/LncRNA_v2/Mouse/DE/res_sni_diff_B10.D2_vs_BALB.c.csv", row.name=1)
ens_genes_sni_vs_sham_BALB.c <- sni_vs_sham_BALB.c[grep("ENS", rownames(sni_vs_sham_BALB.c)), ]
ens_genes_sni_vs_sham_B10.D2 <- sni_vs_sham_B10.D2[grep("ENS", rownames(sni_vs_sham_B10.D2)), ]
mouse_lncs_balbc <- sni_vs_sham_BALB.c[grep("LncRNA", rownames(sni_vs_sham_BALB.c)), ]
mouse_lncs_b10d2 <- sni_vs_sham_B10.D2[grep("LncRNA", rownames(sni_vs_sham_B10.D2)), ]

snt_vs_sham_rat <- read.csv("/home/george/Desktop/LncRNA_v2/Rat/DE/res_snt_vs_sham.csv", row.names=1)
rat_genes <- snt_vs_sham_rat[grep("ENS", rownames(snt_vs_sham_rat)), ]
rat_lncs <- snt_vs_sham_rat[grep("LncRNA", rownames(snt_vs_sham_rat)), ]

neurons_vs_ips <- read.csv("/home/george/Desktop/LncRNA_v2/IPS/DE/res_NEURONS_vs_IPS_AD2.csv", row.names=1)
IPS_genes <- neurons_vs_ips[grep("ENS", rownames(neurons_vs_ips)), ]
IPS_lncs <- neurons_vs_ips[grep("LncRNA", rownames(neurons_vs_ips)), ]

sig_de_balbc <- ens_genes_sni_vs_sham_BALB.c[ens_genes_sni_vs_sham_BALB.c$padj < 0.05 & !is.na(ens_genes_sni_vs_sham_BALB.c$padj),]

sig_de_b10.d2 <- ens_genes_sni_vs_sham_B10.D2[ens_genes_sni_vs_sham_B10.D2$padj < 0.05 & !is.na(ens_genes_sni_vs_sham_B10.D2$padj),]

common <- sig_de_balbc[rownames(sig_de_balbc) %in% rownames(sig_de_b10.d2),]

only_balbc <- sig_de_balbc[!(rownames(sig_de_balbc) %in% rownames(sig_de_b10.d2)),]  

only_b10.d2 <- sig_de_b10.d2[!(rownames(sig_de_b10.d2) %in% rownames(sig_de_balbc)),]


write.csv(file = "/home/george/Desktop/LncRNA_v2/homology/common_mouse_de.csv", common)
write.csv(file = "/home/george/Desktop/LncRNA_v2/homology/only_balbc_de.csv", only_balbc)
write.csv(file = "/home/george/Desktop/LncRNA_v2/homology/only_b10d2_de.csv", only_b10.d2)




###############
# Mouse - Rat #
###############
mm_ortholog_rat_up <- getLDS(attributes = "ensembl_gene_id", filters = "ensembl_gene_id", values = as.character(rownames(rat_genes[rat_genes$padj < 0.05 & rat_genes$log2FoldChange > 0 & !is.na(rat_genes$padj),])), mart = rat, attributesL = c('mgi_symbol', 'ensembl_gene_id'), martL = mouse)

mm_ortholog_rat_down <- getLDS(attributes = "ensembl_gene_id", filters = "ensembl_gene_id", values = as.character(rownames(rat_genes[rat_genes$padj < 0.05 & rat_genes$log2FoldChange < 0 & !is.na(rat_genes$padj),])), mart = rat, attributesL = c('mgi_symbol', 'ensembl_gene_id'), martL = mouse)  

rat_mm_orthologs <- rbind(mm_ortholog_rat_up, mm_ortholog_rat_down)
write.csv(file = "/home/george/Desktop/LncRNA_v2/homology/DE_rat_mm_orthologs.csv", rat_mm_orthologs)

BALB.c_rat_up <- ens_genes_sni_vs_sham_BALB.c[rownames(ens_genes_sni_vs_sham_BALB.c) %in% as.character(mm_ortholog_rat_up$Gene.stable.ID.1) & ens_genes_sni_vs_sham_BALB.c$padj < 0.05 & ens_genes_sni_vs_sham_BALB.c$log2FoldChange > 0 & !is.na(ens_genes_sni_vs_sham_BALB.c$padj),]
BALB.c_rat_down <- ens_genes_sni_vs_sham_BALB.c[rownames(ens_genes_sni_vs_sham_BALB.c) %in% as.character(mm_ortholog_rat_down$Gene.stable.ID.1) & ens_genes_sni_vs_sham_BALB.c$padj < 0.05 & ens_genes_sni_vs_sham_BALB.c$log2FoldChange < 0 & !is.na(ens_genes_sni_vs_sham_BALB.c$padj),]

BALB.c_rat_orthologs <- rbind(BALB.c_rat_up, BALB.c_rat_down)

write.csv(BALB.c_rat_orthologs, file = "/home/george/Desktop/LncRNA_v2/homology/DE_BALB.c_rat_orthologs.csv")

B10.D2_rat_up <- ens_genes_sni_vs_sham_B10.D2[rownames(ens_genes_sni_vs_sham_B10.D2) %in% as.character(mm_ortholog_rat_up$Gene.stable.ID.1) & ens_genes_sni_vs_sham_B10.D2$padj < 0.05 & ens_genes_sni_vs_sham_B10.D2$log2FoldChange > 0 & !is.na(ens_genes_sni_vs_sham_B10.D2$padj),]
B10.D2_rat_down <- ens_genes_sni_vs_sham_B10.D2[rownames(ens_genes_sni_vs_sham_B10.D2) %in% as.character(mm_ortholog_rat_down$Gene.stable.ID.1) & ens_genes_sni_vs_sham_B10.D2$padj < 0.05 & ens_genes_sni_vs_sham_B10.D2$log2FoldChange < 0 & !is.na(ens_genes_sni_vs_sham_B10.D2$padj),]

B10.D2_rat_orthologs <- rbind(B10.D2_rat_up, B10.D2_rat_down)

write.csv(B10.D2_rat_orthologs, file = "/home/george/Desktop/LncRNA_v2/homology/DE_B10.D2_rat_orthologs.csv")

###############
# Mouse - IPS #
###############
mm_ortholog_IPS_up <- getLDS(attributes = "ensembl_gene_id", filters = "ensembl_gene_id", values = as.character(rownames(IPS_genes[IPS_genes$padj < 0.05 & IPS_genes$log2FoldChange > 0 & !is.na(IPS_genes$padj),])), mart = human, attributesL = c('mgi_symbol', 'ensembl_gene_id'), martL = mouse)

mm_ortholog_IPS_down <- getLDS(attributes = "ensembl_gene_id", filters = "ensembl_gene_id", values = as.character(rownames(IPS_genes[IPS_genes$padj < 0.05 & IPS_genes$log2FoldChange < 0 & !is.na(IPS_genes$padj),])), mart = human, attributesL = c('mgi_symbol', 'ensembl_gene_id'), martL = mouse)  

IPS_mm_orthologs <- rbind(mm_ortholog_IPS_up, mm_ortholog_IPS_down)
write.csv(file = "/home/george/Desktop/LncRNA_v2/homology/DE_IPS_mm_orthologs.csv", IPS_mm_orthologs)

BALB.c_IPS_up <- ens_genes_sni_vs_sham_BALB.c[rownames(ens_genes_sni_vs_sham_BALB.c) %in% as.character(mm_ortholog_IPS_up$Gene.stable.ID.1) & ens_genes_sni_vs_sham_BALB.c$padj < 0.05 & ens_genes_sni_vs_sham_BALB.c$log2FoldChange > 0 & !is.na(ens_genes_sni_vs_sham_BALB.c$padj),]
BALB.c_IPS_down <- ens_genes_sni_vs_sham_BALB.c[rownames(ens_genes_sni_vs_sham_BALB.c) %in% as.character(mm_ortholog_IPS_down$Gene.stable.ID.1) & ens_genes_sni_vs_sham_BALB.c$padj < 0.05 & ens_genes_sni_vs_sham_BALB.c$log2FoldChange < 0 & !is.na(ens_genes_sni_vs_sham_BALB.c$padj),]

BALB.c_IPS_orthologs <- rbind(BALB.c_IPS_up, BALB.c_IPS_down)

write.csv(BALB.c_IPS_orthologs, file = "/home/george/Desktop/LncRNA_v2/homology/DE_BALB.c_IPS_orthologs.csv")

B10.D2_IPS_up <- ens_genes_sni_vs_sham_B10.D2[rownames(ens_genes_sni_vs_sham_B10.D2) %in% as.character(mm_ortholog_IPS_up$Gene.stable.ID.1) & ens_genes_sni_vs_sham_B10.D2$padj < 0.05 & ens_genes_sni_vs_sham_B10.D2$log2FoldChange > 0 & !is.na(ens_genes_sni_vs_sham_B10.D2$padj),]
B10.D2_IPS_down <- ens_genes_sni_vs_sham_B10.D2[rownames(ens_genes_sni_vs_sham_B10.D2) %in% as.character(mm_ortholog_IPS_down$Gene.stable.ID.1) & ens_genes_sni_vs_sham_B10.D2$padj < 0.05 & ens_genes_sni_vs_sham_B10.D2$log2FoldChange < 0 & !is.na(ens_genes_sni_vs_sham_B10.D2$padj),]

B10.D2_IPS_orthologs <- rbind(B10.D2_IPS_up, B10.D2_IPS_down)

write.csv(B10.D2_IPS_orthologs, file = "/home/george/Desktop/LncRNA_v2/homology/DE_B10.D2_IPS_orthologs.csv")

##################
##Novel LncsRNAs##
##################
#####################################
#Antisense of the orthologous genes #
#####################################

###############
# Mouse - Rat #
###############
rat_antisense <- read.csv("/home/george/Desktop/LncRNA_v2/Rat/antisense/snt_vs_sham/antisense_all_expression.csv", row.names=1)

balbc_antisense <- read.csv("/home/george/Desktop/LncRNA_v2/Mouse/antisense/sni_vs_sham_balb.c/antisense_all_expression.csv", row.names=1)

b10d2_antisense <- read.csv("/home/george/Desktop/LncRNA_v2/Mouse/antisense/sni_vs_sham_b10.d2/antisense_all_expression.csv", row.names=1)

mm_ortholog_rat_antisense <- getLDS(attributes = "ensembl_gene_id", filters = "ensembl_gene_id", values = as.character(rat_antisense$V16), mart = rat, attributesL = c('mgi_symbol', 'ensembl_gene_id'), martL = mouse)

balbc_antisense_rat_orthologs <- balbc_antisense[as.character(balbc_antisense$V16) %in% as.character(mm_ortholog_rat_antisense$Gene.stable.ID.1), ]

b10d2_antisense_rat_orthologs <- b10d2_antisense[as.character(b10d2_antisense$V16) %in% as.character(mm_ortholog_rat_antisense$Gene.stable.ID.1), ]

rat_antisense_mm_orthologs <- rat_antisense[as.character(rat_antisense$V16) %in% as.character(mm_ortholog_rat_antisense$Gene.stable.ID), ]

write.csv(balbc_antisense_rat_orthologs, file = "/home/george/Desktop/LncRNA_v2/homology/balbc_antisense_rat_orthologs.csv")
write.csv(b10d2_antisense_rat_orthologs, file = "/home/george/Desktop/LncRNA_v2/homology/b10d2_antisense_rat_orthologs.csv")
write.csv(rat_antisense_mm_orthologs, file = "/home/george/Desktop/LncRNA_v2/homology/rat_antisense_mm_orthologs.csv")

###############
# Mouse - IPS #
###############
IPS_antisense <- read.csv("/home/george/Desktop/LncRNA_v2/IPS/antisense/AD2/antisense_all_expression.csv", row.names=1)

mm_ortholog_IPS_antisense <- getLDS(attributes = "ensembl_gene_id", filters = "ensembl_gene_id", values = as.character(IPS_antisense$V16), mart = human, attributesL = c('mgi_symbol', 'ensembl_gene_id'), martL = mouse)

balbc_antisense_IPS_orthologs <- balbc_antisense[as.character(balbc_antisense$V16) %in% as.character(mm_ortholog_IPS_antisense$Gene.stable.ID.1), ]

b10d2_antisense_IPS_orthologs <- b10d2_antisense[as.character(b10d2_antisense$V16) %in% as.character(mm_ortholog_IPS_antisense$Gene.stable.ID.1), ]

IPS_antisense_mm_orthologs <- IPS_antisense[as.character(IPS_antisense$V16) %in% as.character(mm_ortholog_IPS_antisense$Gene.stable.ID), ]

write.csv(balbc_antisense_IPS_orthologs, file = "/home/george/Desktop/LncRNA_v2/homology/balbc_antisense_IPS_orthologs.csv")
write.csv(b10d2_antisense_IPS_orthologs, file = "/home/george/Desktop/LncRNA_v2/homology/b10d2_antisense_IPS_orthologs.csv")
write.csv(IPS_antisense_mm_orthologs, file = "/home/george/Desktop/LncRNA_v2/homology/IPS_antisense_mm_orthologs.csv")

######################################################################################
#LncRNAs found in syntenic conserved blocks between Mm10, Hg38 and Rn5 lifted to Rn6 #
######################################################################################
syntenic_mouse <- read.table("/home/george/Desktop/LncRNA_v2/Mouse_syntenic_lncRNAs.bed", sep="\t")
syntenic_rat <- read.table("/home/george/Desktop/LncRNA_v2/Rat_syntenic_lncRNAs.bed", sep="\t")
syntenic_IPS <- read.table("/home/george/Desktop/LncRNA_v2/Hs_syntenic_lncRNAs.bed", sep="\t")

syntenic_mouse <- syntenic_mouse[,c(4,16)]
syntenic_rat <- syntenic_rat[,c(4,16)]
syntenic_IPS <- syntenic_IPS[,c(4,16)]

write.csv(syntenic_mouse, file = "/home/george/Desktop/LncRNA_v2/homology/syntenic_mouse.csv")
write.csv(syntenic_rat, file = "/home/george/Desktop/LncRNA_v2/homology/syntenic_rat.csv")
write.csv(syntenic_IPS, file = "/home/george/Desktop/LncRNA_v2/homology/syntenic_IPS.csv")

common_blocks_mouse_rat <- intersect(syntenic_mouse$V16, syntenic_rat$V16)

syntenic_mouse_mm10_rn6 <- syntenic_mouse[syntenic_mouse$V16 %in% common_blocks_mouse_rat, ]
syntenic_rat_mm10_rn6 <- syntenic_rat[syntenic_rat$V16 %in% common_blocks_mouse_rat, ]

write.csv(syntenic_mouse_mm10_rn6, file = "/home/george/Desktop/LncRNA_v2/homology/syntenic_mouse_mm10_rn6.csv")
write.csv(syntenic_rat_mm10_rn6, file = "/home/george/Desktop/LncRNA_v2/homology/syntenic_rat_mm10_rn6.csv")

common_blocks_mouse_hs <- intersect(syntenic_mouse$V16, syntenic_IPS$V16)

syntenic_mouse_mm10_hg38 <- syntenic_mouse[syntenic_mouse$V16 %in% common_blocks_mouse_hs, ]
syntenic_IPS_mm10_hg38 <- syntenic_IPS[syntenic_IPS$V16 %in% common_blocks_mouse_hs, ]

write.csv(syntenic_mouse_mm10_hg38, file = "/home/george/Desktop/LncRNA_v2/homology/syntenic_mouse_mm10_hg38.csv")
write.csv(syntenic_IPS_mm10_hg38, file = "/home/george/Desktop/LncRNA_v2/homology/syntenic_IPS_mm10_hg38.csv")





