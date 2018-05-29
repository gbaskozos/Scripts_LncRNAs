library("org.Mm.eg.db")
library("org.Rn.eg.db")
library("org.Hs.eg.db")
library(Hmisc)
library("ggplot2")
library("gplots")
library(Tmisc)
library(calibrate)
library(lattice)
library(plotrix)
library(parallel)
library("RColorBrewer")
library("DESeq2")
library(grid)

################
# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
################
n_cores =3 

options(stringsAsFactors=FALSE)

PATH_mouse <- "/home/george/Desktop/LncRNA_v2/Mouse"
PATH_rat <- "/home/george/Desktop/LncRNA_v2/Rat"
PATH_ips <- "/home/george/Desktop/LncRNA_v2/IPS"

PATH_out <- "/home/george/Desktop/LncRNA_v2/Figures"

##################
#### Figure 2 ####
##################
antisense_mouse <- read.csv(paste0(PATH_mouse, "/antisense/sni_vs_sham_balb.c/antisense_all_expression.csv"))

intronic_mouse <- read.csv(paste0(PATH_mouse, "/intronic/sni_vs_sham_balb.c/intronic_all_expression.csv"))

lincs_mouse <- read.csv(paste0(PATH_mouse, "/DE/lincs_closest.csv"))

antisense_rat <- read.csv(paste0(PATH_rat, "/antisense/snt_vs_sham/antisense_all_expression.csv"))

intronic_rat <- read.csv(paste0(PATH_rat, "/intronic/snt_vs_sham/intronic_all_expression.csv"))

lincs_rat <- read.csv(paste0(PATH_rat, "/DE/lincs_closest.csv"))

mouse_lncs_IDs <- unique(c(antisense_mouse$ID, intronic_mouse$ID, lincs_mouse$lincRNA))
mouse_lncs_names <- unique(c(antisense_mouse$V4, intronic_mouse$V4, lincs_mouse$LincRNA_symbol))

rat_lncs_IDs <- unique(c(antisense_rat$ID, intronic_rat$ID, lincs_rat$lincRNA))
rat_lncs_names <- unique(c(antisense_rat$V4, intronic_rat$V4, lincs_rat$LincRNA_symbol))

load(paste0(PATH_mouse, "/myRegs_cov/rRegs_complete_coff.RData"))
mouse_lncs <- rRegs_complete_coff_nc[mouse_lncs_names[mouse_lncs_names %in% names(rRegs_complete_coff_nc)]]
nr_of_exons_mouse <- unlist(mclapply(mouse_lncs, function(x) length(x), mc.cores=n_cores, mc.cleanup=TRUE, mc.allow.recursive=TRUE))

load(paste0(PATH_rat, "/myRegs_cov/rRegs_complete_coff.RData"))
rat_lncs <- rRegs_complete_coff_nc[rat_lncs_names[rat_lncs_names %in% names(rRegs_complete_coff_nc)]]
nr_of_exons_rat <- unlist(mclapply(rat_lncs, function(x) length(x), mc.cores=n_cores, mc.cleanup=TRUE, mc.allow.recursive=TRUE))


pdf(file = paste(PATH_out,"/Figure_2.pdf", sep=""), width = 9, height = 9)

par(oma = c(0, 0, 2, 0), mar = c(3, 3, 2, 2), mgp = c(2, 1, 0), xpd = NA)
layout(matrix(c(1,2,3,4,5,5), 3, 2, byrow = TRUE))

#Proportions mouse
slices <- c(length(unique(antisense_mouse$ID)), length(unique(lincs_mouse$lincRNA)), length(unique(intronic_mouse$ID))) 
lbls <- c("Antisense", "Intergenic", "Intronic")

pct <- round(slices/sum(slices)*100)
lbls <- paste(lbls, pct) # add percents to labels 
lbls <- paste(lbls,"%",sep="") # ad % to labels 
pie(slices,labels = lbls, col=rainbow(length(lbls)), main="Mouse DRG", cex=1.1, cex.main=1.2, radius = 0.8 )

#Proportions rat
slices <- c(length(unique(antisense_rat$ID)), length(unique(lincs_rat$lincRNA)), length(unique(intronic_rat$ID))) 
lbls <- c("Antisense", "Intergenic", "Intronic")
pct <- round(slices/sum(slices)*100)
lbls <- paste(lbls, pct) # add percents to labels 
lbls <- paste(lbls,"%",sep="") # ad % to labels 
pie(slices,labels = lbls, col=rainbow(length(lbls)), main="Rat DRG", cex=1.1, cex.main=1.2, radius = 0.8 )


#Nr of exons mouse
plot(table(nr_of_exons_mouse),type="l",col="red", xaxt="n", xlab="Nr of exons", ylab="Novel LncRNAs", main = "Nr of exons in Mouse", cex.main=1.2, cex.lab=1.1)
axis(1, at=seq(1,max(nr_of_exons_mouse),1),labels=NULL, col.axis="black", las=0, cex.axis=1.1, cex.lab=1.1)
legend("topright", title=NULL, legend = c(paste("Median nr of exons: ", median(nr_of_exons_mouse), sep=""), c(paste("Mode: ", names(sort(-table(nr_of_exons_mouse)))[1], sep=""
))), cex=1.1, bg="white", xjust=0)

#Nr of exons rat
plot(table(nr_of_exons_rat),type="l",col="red", xaxt="n", xlab="Nr of exons", ylab="Novel LncRNAs", main = "Nr of exons in Rat", cex.main=1.2, cex.lab=1.1)
axis(1, at=seq(1,max(nr_of_exons_rat),1),labels=NULL, col.axis="black", las=0, cex.axis=1.1, cex.lab=1.1)
legend("topright", title=NULL, legend = c(paste("Median nr of exons: ", median(nr_of_exons_rat), sep=""), c(paste("Mode: ", names(sort(-table(nr_of_exons_rat)))[1], sep=""
))), cex=1.1, bg="white", xjust=0)



#TSS
tss_lncs_distance <- read.table(file = paste0(PATH_mouse, "/tss/lncs_tss_distance.bed"), sep="\t", header=FALSE)

dsBase.iqr <- tss_lncs_distance
 
# Create a variable/vector/collection of the column names you want to remove outliers on.
vars <- c("V22")

Outliers <- c()
 
for(i in vars){
 
  # Get the Min/Max values
  max <- quantile(dsBase.iqr[,i],0.75, na.rm=TRUE) + (IQR(dsBase.iqr[,i], na.rm=TRUE) * 1.5 )
  min <- quantile(dsBase.iqr[,i],0.25, na.rm=TRUE) - (IQR(dsBase.iqr[,i], na.rm=TRUE) * 1.5 )
  
  # Get the id's using which
  idx <- which(dsBase.iqr[,i] < min | dsBase.iqr[,i] > max)
  
  # Output the number of outliers in each variable
  print(paste(i, length(idx), sep=''))
  
  # Append the outliers list
  Outliers <- c(Outliers, idx) 
}
 
Outliers <- sort(Outliers)
dsBase.iqr <- dsBase.iqr[-Outliers,]

d_lncs_no <- density(dsBase.iqr$V22)
plot(d_lncs_no, main="Distance between novel LncRNAs and Transcription Start Sites in Mouse DRG", xlab="distance from TSS", cex.main = 1.2, cex.lab=1.1, cex.axis=1.1)
polygon(d_lncs_no, col="red", border="blue")
legend("topleft", title="Distances > 1.5xIQR were removed", legend=c(paste("Mean: ", round(mean(dsBase.iqr$V22), digits=2), sep="")), cex=1.1, bg="white")

title("Genomic context of novel LncRNAs", outer=TRUE, cex=1.4)

dev.off()

#Data for figure 2
tss_mouse_fig2 <- data.frame(LncRNA_name = tss_lncs_distance$V4, distance_TSS = tss_lncs_distance$V22)
rownames(tss_mouse_fig2) <- tss_mouse_fig2$LncRNA_name

mouse_fig2 <- data.frame(LncRNA_ID = c(unique(antisense_mouse$ID), unique(intronic_mouse$ID), unique(lincs_mouse$lincRNA)), LncRNA_name = c(unique(antisense_mouse$V4), unique(intronic_mouse$V4), unique(lincs_mouse$LincRNA_symbol)),  genomic_context = c(rep("Antisense", length(unique(antisense_mouse$ID))), rep("Intronic", length(unique(intronic_mouse$ID))), rep("Intergenic", length(unique(lincs_mouse$lincRNA)))), nr_of_exons = nr_of_exons_mouse[c(unique(antisense_mouse$V4), unique(intronic_mouse$V4), unique(lincs_mouse$LincRNA_symbol))], distance_form_TSS = tss_mouse_fig2[c(unique(antisense_mouse$V4), unique(intronic_mouse$V4), unique(lincs_mouse$LincRNA_symbol)),]$distance_TSS)

write.csv(file = paste(PATH_out,"/Figure_2_mouse.csv", sep=""), mouse_fig2)

rat_fig2 <- data.frame(LncRNA_ID = c(unique(antisense_rat$ID), unique(intronic_rat$ID), unique(lincs_rat$lincRNA)), LncRNA_name = c(unique(antisense_rat$V4), unique(intronic_rat$V4), unique(lincs_rat$LincRNA_symbol)),  genomic_context = c(rep("Antisense", length(unique(antisense_rat$ID))), rep("Intronic", length(unique(intronic_rat$ID))), rep("Intergenic", length(unique(lincs_rat$lincRNA)))), nr_of_exons = nr_of_exons_rat[c(unique(antisense_rat$V4), unique(intronic_rat$V4), unique(lincs_rat$LincRNA_symbol))])

write.csv(file = paste(PATH_out,"/Figure_2_rat.csv", sep=""), rat_fig2)

#######################################
############## Figure 3 ###############
#######################################


cutoff <- 0
load(paste(PATH_mouse,"/toc_all.RData", sep=""))
toc_coff_mouse <- toc[apply(toc, 1, function(x) sum(x[grep("BALB.c", names(toc))] > cutoff) == length(grep("BALB.c", names(toc))) | sum(x[grep("B10.D2", names(toc))] > cutoff) == length(grep("B10.D2", names(toc)))), ]
nc_ensembl <- read.csv(file = paste(PATH_mouse,"/mm10_genome/ens_lncs_annotation.csv", sep=""), row.name=1)

toc_mouse_genes <- toc_coff_mouse[grepl("ENS", rownames(toc_coff_mouse)) & !rownames(toc_coff_mouse) %in% nc_ensembl$ensembl_gene_id,]
toc_mouse_lncRNAs <- toc_coff_mouse[grepl("LncRNA", rownames(toc_coff_mouse)) | rownames(toc_coff_mouse) %in% nc_ensembl$ensembl_gene_id,]

load(paste(PATH_rat,"/toc_all.RData", sep=""))
toc_coff_rat <- toc[apply(toc, 1, function(x) sum(x[c(1:4)] > cutoff) == 4 | sum(x[c(5:8)] > cutoff) == 4), ]
nc_ensembl <- read.csv(file = paste(PATH_rat,"/rn6_genome/ensembl_nc.csv", sep=""), row.name=1)

toc_rat_genes <- toc_coff_rat[grepl("ENS", rownames(toc_coff_rat)) & !rownames(toc_coff_rat) %in% nc_ensembl$ensembl_gene_id,]
toc_rat_lncRNAs <- toc_coff_rat[grepl("LncRNA", rownames(toc_coff_rat)) | rownames(toc_coff_rat) %in% nc_ensembl$ensembl_gene_id,]


med_mouse_genes <- apply(toc_mouse_genes,2, median)
med_rat_genes <- apply(toc_rat_genes,2, median)


med_mouse_lncs <- apply(toc_mouse_lncRNAs,2, median)
med_rat_lncs <- apply(toc_rat_lncRNAs,2, median)

n <- length(c(med_mouse_genes, med_rat_genes))

sample_stats <- data.frame(sample = rep(c(names(med_mouse_genes), names(med_rat_genes)),2), cat = c(rep("Genes",n), rep("LncRNAs",n)), counts = c(med_mouse_genes, med_rat_genes, med_mouse_lncs, med_rat_lncs))

sample_stats$group <- rep( c(rep("Mouse",20), rep("Rat", 8)),2)

aggdata <- aggregate(sample_stats$counts, by=list(sample_stats$group, sample_stats$cat), FUN=mean, na.rm=TRUE)

aggdata_sd <- aggregate(sample_stats$counts, by=list(sample_stats$group, sample_stats$cat), FUN=std.error, na.rm=TRUE)

aggdata <- cbind(aggdata, aggdata_sd$x)

data_sem <- transform(aggdata, lower=x-aggdata_sd$x, upper=x+aggdata_sd$x)

names(data_sem) <- c("name", "Legend", "median", "sd", "lower", "upper")

data_sem$cat <- c("Genes","Genes","LncRNAs", "LncRNAs")

wilcox.test(med_mouse_genes, med_mouse_lncs)$p.value < 0.001
wilcox.test(med_rat_genes, med_rat_lncs)$p.value < 0.001

p1 <- ggplot(data_sem, aes(x = name, y = median, fill = Legend)) +  geom_bar(stat="identity", position=position_dodge(), width=1) + geom_errorbar(aes(ymax=upper, ymin=lower), position=position_dodge(0.9), data=data_sem, width = 0.5, size = 0.15) + xlab("Group average") + ylab("Median normalised read-pairs + SEM") + ggtitle("LncRNAs (Novel + annotated) vs p.coding genes") + theme(axis.text.y = element_text(size= 12, face="bold"), axis.title.y = element_text(size=15, face="bold"), axis.line.y = element_line(color="black", size = 0.3), axis.title.x = element_blank(), axis.text.x = element_text(size= 12, angle=45, hjust=1, face="bold"), legend.title=NULL, legend.text=element_text(size=15, face="bold"), legend.key.size=unit(0.5, "cm"), plot.title=element_text(size=15, face="bold", hjust = 0.5)) + scale_y_continuous(limits=c(0, max(data_sem$upper)), breaks = seq(from = 0, to = max(data_sem$upper), by=20) ) 

load(paste0(PATH_mouse, "/tissue_specific/lncs_neuron_types_counts.RData"))

tissue_specificity <- read.csv(paste0(PATH_mouse, "/tissue_specific/tissue_specificity.csv"), row.name=1)

tissue_specific <- read.csv(paste0(PATH_mouse, "/tissue_specific/tissue_specific.csv"), row.name=1)

lnc_ti <- rownames(tissue_specificity)[grep("LncRNA", rownames(tissue_specificity))]

d <- plotPCA(vsd[rownames(vsd) %in% lnc_ti], intgroup="neuron_type", ntop = 500, returnData=TRUE)
percentVar <- round(100 * attr(d, "percentVar"))
p2 <- ggplot(d, aes(x=PC1,y=PC2, color=neuron_type, label=names(toc))) + geom_point(size=4) + geom_text(size=1, check_overlap=TRUE, vjust = 0, nudge_y = 0.5) + xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance")) 

vsd_mat <- assay(vsd[rownames(vsd) %in% lnc_ti])

rv <- rowVars(vsd_mat)
top_500 <- vsd_mat[order(-rv)[1:500], ]

ti_data <- read.csv(file = paste0(PATH_mouse,"/tissue_specific/ti_data.csv"), row.name=1)

p3 <- ggplot(ti_data, aes(x = Gene_biotype, y = percentage, fill = Specificity)) + geom_bar(stat = "identity") + geom_text(data=ti_data, aes(x = Gene_biotype, y = 100- (0.5*percentage), label = paste0(round(percentage, digits = 2),"%")), size=4) + ggtitle("Neuron-subtype specificity") + theme(axis.text.y = element_text(size= 12, face="bold"), axis.title.y = element_text(size=15, face="bold"), axis.line.y = element_line(color="black", size = 0.3), axis.title.x = element_blank(), axis.text.x = element_text(size= 12, angle=45, hjust=1, face="bold"), legend.title=NULL, legend.text=element_text(size=10, face="bold"), legend.key.size=unit(0.5, "cm"), plot.title=element_text(size=15, face="bold", hjust = 0.5)) + scale_y_continuous(limits=c(0, 103), breaks = seq(from = 0, to = 100, by=5) ) 

tissue_specific <- read.csv(file = paste0(PATH_mouse, "/tissue_specific/tissue_specific_lncRNAs.csv"), row.name=1)

tissue_specific$t.index <- factor(tissue_specific$t.index, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10") )

tau_freq <- data.frame(Neuron_subtype = c("1. MHN", "2. MHN (MI, IS)", "3. C-LTMR", "4. MHN (IS)", "5. MHN (IS)", "6. MHN", "7. MHN (NS)", "8. MR", "9.MHN", "10. MR"), frequency = as.matrix(table(tissue_specific$t.index)))

tau_freq$Neuron_subtype <- factor(tau_freq$Neuron_subtype, levels = c("1. MHN", "2. MHN (MI, IS)", "3. C-LTMR", "4. MHN (IS)", "5. MHN (IS)", "6. MHN", "7. MHN (NS)", "8. MR", "9.MHN", "10. MR") )

p4 <- ggplot(data=tau_freq, aes(x=Neuron_subtype, y=frequency)) + geom_bar(stat="identity", fill="steelblue") + geom_text(aes(label=frequency), vjust=1, hjust=1.2, color="white", size=4, fontface='bold') + ggtitle("Subtype-specific novel LncRNAs") + coord_flip() + theme(axis.text.y = element_text(size= 12, face="bold"), axis.title.y = element_blank(), axis.title.x = element_text(size=15, face="bold"), axis.text.x = element_text(size= 12, angle=45, hjust=1, face="bold"), legend.title=NULL, legend.text=element_text(size=15, face="bold"), legend.key.size=unit(0.5, "cm"), plot.title=element_text(size=15, face="bold", hjust = 0.5))

pdf(file = paste(PATH_out,"/Figure_3_test.pdf", sep=""), width = 12, height = 12)
print(multiplot(p1, p3, p2, p4, cols=2))
dev.off()

# Data for figure 3

write.csv(file = paste0(PATH_out, "/Figure_3_genes_lncs.csv"), sample_stats) 
write.csv(file = paste0(PATH_out, "/Figure_3_top_500_pca.csv"), top_500)
write.csv(file = paste0(PATH_out, "/Figure_3_tissue_specific.csv"), tissue_specific)
write.csv(file = paste0(PATH_out, "/Figure_3_tissue_index.csv"), ti_data)



#######################################
############## Figure 4 ###############
#######################################
antisense_ips <- read.csv(paste0(PATH_ips, "/antisense/AD2/antisense_all_expression.csv"))

intronic_ips <- read.csv(paste0(PATH_ips, "/intronic/AD2/intronic_all_expression.csv"))

lincs_ips <- read.csv(paste0(PATH_ips, "/DE/lincs_closest.csv"))

ips_lncs_IDs <- unique(c(antisense_ips$ID, intronic_ips$ID, lincs_ips$lincRNA))
ips_lncs_names <- unique(c(antisense_ips$V4, intronic_ips$V4, lincs_ips$LincRNA_symbol))

load(paste0(PATH_ips, "/myRegs_cov/rRegs_complete_coff.RData"))
ips_lncs <- rRegs_complete_coff_nc[ips_lncs_names[ips_lncs_names %in% names(rRegs_complete_coff_nc)]]
nr_of_exons_ips <- unlist(mclapply(ips_lncs, function(x) length(x), mc.cores=n_cores, mc.cleanup=TRUE, mc.allow.recursive=TRUE))

slices <- c(length(unique(antisense_ips$ID)), length(unique(lincs_ips$lincRNA)), length(unique(intronic_ips$ID))) 
lbls <- c("Antisense", "Intergenic", "Intronic")
pct <- round(slices/sum(slices)*100)
lbls <- paste(lbls, pct) # add percents to labels 
lbls <- paste(lbls,"%",sep="") # ad % to labels 

load(paste(PATH_ips,"/toc_all.RData", sep=""))
toc_coff_ips <- toc[apply(toc, 1, function(x) sum(x[grep("IPS", names(toc))] > cutoff) == length(grep("IPS", names(toc))) | sum(x[grep("NEURONS", names(toc))] > cutoff) == length(grep("NEURONS", names(toc)))), ]
nc_ensembl <- read.csv(file = paste(PATH_ips,"/HS_genome/ensembl_nc.csv", sep=""), row.name=1)

toc_ips_genes <- toc_coff_ips[grepl("ENS", rownames(toc_coff_ips)) & !rownames(toc_coff_ips) %in% nc_ensembl$ensembl_gene_id,]
toc_ips_lncRNAs <- toc_coff_ips[grepl("LncRNA", rownames(toc_coff_ips)) | rownames(toc_coff_ips) %in% nc_ensembl$ensembl_gene_id,]
med_ips_genes <- apply(toc_ips_genes,2, median)

med_ips_lncs <- apply(toc_ips_lncRNAs,2, median)

n <- length(med_ips_lncs)

sample_stats <- data.frame(sample = rep(names(med_ips_genes),2), cat = c(rep("Genes",n), rep("LncRNAs",n)), counts = c(med_ips_genes, med_ips_lncs))

sample_stats$group <- rep( c(rep("HS.IPSc", 12)),2)

aggdata <- aggregate(sample_stats$counts, by=list(sample_stats$group, sample_stats$cat), FUN=mean, na.rm=TRUE)

aggdata_sd <- aggregate(sample_stats$counts, by=list(sample_stats$group, sample_stats$cat), FUN=std.error, na.rm=TRUE)

aggdata <- cbind(aggdata, aggdata_sd$x)

data_sem <- transform(aggdata, lower=x-aggdata_sd$x, upper=x+aggdata_sd$x)

names(data_sem) <- c("name", "Legend", "median", "sd", "lower", "upper")

data_sem$cat <- c("Genes","LncRNAs")

wilcox.test(med_ips_genes, med_ips_lncs)$p.value < 0.001

p1 <- ggplot(data_sem, aes(x = name, y = median, fill = Legend)) +  geom_bar(stat="identity", position=position_dodge(), width=1) + geom_errorbar(aes(ymax=upper, ymin=lower), position=position_dodge(0.9), data=data_sem, width = 0.5, size = 0.15) + xlab("Group average") + ylab("Median normalised read-pairs + SEM") + ggtitle("LncRNAs (Novel + annotated) vs ENSEMBL genes") + theme(axis.text.y = element_text(size= 12, face="bold"), axis.title.y = element_text(size=15, face="bold"), axis.line.y = element_line(color="black", size = 0.3), axis.title.x = element_blank(), axis.text.x = element_text(size= 12, angle=45, hjust=1, face="bold"), legend.title=NULL, legend.text=element_text(size=15, face="bold"), legend.key.size=unit(0.5, "cm"), plot.title=element_text(size=15, face="bold", hjust = 0.5)) + scale_y_continuous(limits=c(0, max(data_sem$upper)), breaks = seq(from = 0, to = max(data_sem$upper), by=20) ) 

pdf(file = paste(PATH_out,"/Figure_4_genes_lncs.pdf", sep=""), width = 6, height = 6)
print(p1)
dev.off()


pdf(file = paste(PATH_out,"/Figure_4_upper_pannel.pdf", sep=""), width = 9, height = 3)

par(oma = c(0, 0, 2, 0), mar = c(3, 3, 2, 2), mgp = c(2, 1, 0), xpd = NA)
layout(matrix(c(1,2), 1, 2, byrow = TRUE))

pie(slices,labels = lbls, col=rainbow(length(lbls)), main="IPSc derived Neurons", cex=1.1, cex.main=1.2, radius = 0.8 )

#Nr of exons ips
plot(table(nr_of_exons_ips),type="l",col="red", xaxt="n", xlab="Nr of exons", ylab="Novel LncRNAs", main = "Nr of exons in IPSc", cex.main=1.2, cex.lab=1.1)
axis(1, at=seq(1,max(nr_of_exons_ips),1),labels=NULL, col.axis="black", las=0, cex.axis=1.1, cex.lab=1.1)
legend("topright", title=NULL, legend = c(paste("Median nr of exons: ", median(nr_of_exons_ips), sep=""), c(paste("Mode: ", names(sort(-table(nr_of_exons_ips)))[1], sep=""
))), cex=1.1, bg="white", xjust=0)

dev.off()


antisense <- read.csv(paste0(PATH_ips, "/antisense/AD2/antisense_all_expression.csv"), row.name=1)
antisense_ens <- read.csv(paste0(PATH_ips, "/antisense/ensembl_antisense/AD2/antisense_all_expression.csv"), row.name=1)
load(file=paste0(PATH_ips, "/DE/norm_counts.RData"))
nc_ensembl <- read.csv(file = paste0(PATH_ips, "/HS_genome/ens_lncs_annotation.csv"), row.name=1)
load(paste0(PATH_ips, "/pain_genes.RData"))
lincs_closest <- read.csv(file = paste0(PATH_ips, "/DE/lincs_closest.csv"), row.name=1)
lincs_closest_ens <- read.csv(file = paste0(PATH_ips, "/DE/lincs_ENSEMBL_closest.csv"), row.name=1)

both_significant <- subset(antisense, antisense$gene_pvalue<.05 & antisense$lnc_pvalue<.05)

lnc_significant <- subset(antisense, antisense$lnc_pvalue<.05)

both_significant_ens <- subset(antisense_ens, antisense_ens$gene_pvalue<.05 & antisense_ens$lnc_pvalue<.05)

lnc_significant_ens <- subset(antisense_ens, antisense_ens$lnc_pvalue<.05)

res_NEURONS_vs_IPS_AD2 <- read.csv(file = paste0(PATH_ips, "/DE/res_NEURONS_vs_IPS_AD2.csv"), row.name=1)

pain <- res_NEURONS_vs_IPS_AD2[rownames(res_NEURONS_vs_IPS_AD2) %in% pain_genes$Gene.stable.ID.1,]
sodium_channel <- res_NEURONS_vs_IPS_AD2[grep("SCN", res_NEURONS_vs_IPS_AD2$symbol),]
potassium_channel <- res_NEURONS_vs_IPS_AD2[grep("KCN", res_NEURONS_vs_IPS_AD2$symbol),]
chloride_channel <- res_NEURONS_vs_IPS_AD2[grep("CLC", res_NEURONS_vs_IPS_AD2$symbol),]
calcium_channel <- res_NEURONS_vs_IPS_AD2[grep("CAC", res_NEURONS_vs_IPS_AD2$symbol),]
trp_channel <- res_NEURONS_vs_IPS_AD2[grep("TRP", res_NEURONS_vs_IPS_AD2$symbol),]

pain_antisense <- subset(antisense, antisense$V16 %in% rownames(pain))
sodium_channel_antisense <- subset(antisense, antisense$symbol %in% sodium_channel$symbol)
potassium_channel_antisense <- subset(antisense, antisense$symbol %in% potassium_channel$symbol)
chloride_channel_antisense <- subset(antisense, antisense$symbol %in% chloride_channel$symbol)
calcium_channel_antisense <- subset(antisense, antisense$symbol %in% calcium_channel$symbol)
trp_channel_antisense <- subset(antisense, antisense$symbol %in% trp_channel$symbol)

pain_antisense_ens <- subset(antisense_ens, antisense_ens$V16 %in% rownames(pain))
sodium_channel_antisense_ens <- subset(antisense_ens, antisense_ens$sense_symbol %in% sodium_channel$symbol)
potassium_channel_antisense_ens <- subset(antisense_ens, antisense_ens$sense_symbol %in% potassium_channel$symbol)
chloride_channel_antisense_ens <- subset(antisense_ens, antisense_ens$sense_symbol %in% chloride_channel$symbol)
calcium_channel_antisense_ens <- subset(antisense_ens, antisense_ens$sense_symbol %in% calcium_channel$symbol)
trp_channel_antisense_ens <- subset(antisense_ens, antisense_ens$sense_symbol %in% trp_channel$symbol)

pdf(file = paste0(PATH_out,"/Figure_4_antisense_B.pdf"), width = 12, height = 6)
#dev.new(width=5, height=5)
par(oma = c(0, 0, 1, 0), mar = c(8, 3, 2, 2), mgp = c(2, 1, 0))
layout(matrix(c(1,2), 1, 2, byrow = TRUE))

with(antisense, plot(lnc_lfc, gene_lfc, pch=20, col="grey", cex=0.9, main="Genes and novel antisense LncRNAs", xlab="Novel lncRNA log2 fold change Neurons vs IPSc", ylab="Sense gene log2 fold change Neurons vs IPSc", cex.lab=0.7))
legend("topleft", legend = c("Both DE (p.value < 0.05)", "Only LncRNA DE", "Pain gene", "Sodium channel", "Potassium channel", "Chloride channel", "Calcium channel", "Trp channel"), pch=c(20,20,7, 21,22,23,24,25), col=c("green1","red", "blue3", "blue3", "blue3", "blue3", "blue3","blue3"  ), cex=0.7)

abline(h = 0, lty=3)
abline(v = 0, lty=3)

with(subset(antisense, antisense$ID %in% lnc_significant$ID),points(lnc_lfc, gene_lfc, pch=20, col="red",cex=0.9))

with(subset(antisense, antisense$ID %in% both_significant$ID),points(lnc_lfc, gene_lfc, pch=20, col="green1", cex=0.9))

with(subset(antisense, antisense$ID %in% pain_antisense$ID),points(lnc_lfc, gene_lfc, pch=7, col="blue3", cex=0.9))
with(subset(antisense, antisense$ID %in% sodium_channel_antisense$ID),points(lnc_lfc, gene_lfc, pch=21, col="blue3", cex=0.9))
with(subset(antisense, antisense$ID %in% potassium_channel_antisense$ID),points(lnc_lfc, gene_lfc, pch=22, col="blue3", cex=0.9))
with(subset(antisense, antisense$ID %in% chloride_channel_antisense$ID),points(lnc_lfc, gene_lfc, pch=23, col="blue3", cex=0.9))
with(subset(antisense, antisense$ID %in% calcium_channel_antisense$ID),points(lnc_lfc, gene_lfc, pch=24, col="blue3", cex=0.9))
with(subset(antisense, antisense$ID %in% trp_channel_antisense$ID),points(lnc_lfc, gene_lfc, pch=25, col="blue3", cex=0.9))

with(antisense_ens, plot(lnc_lfc, gene_lfc, pch=20, col="grey", cex=0.9, main="Genes and ENSEMBL antisense LncRNAs", xlab="Novel lncRNA log2 fold change Neurons vs IPSc", ylab="Sense gene log2 fold change Neurons vs IPSc", cex.lab=0.7))
legend("topleft", legend = c("Both DE (p.value < 0.05)", "Only LncRNA DE", "Pain gene", "Sodium channel", "Potassium channel", "Chloride channel", "Calcium channel", "Trp channel"), pch=c(20,20,7, 21,22,23,24,25), col=c("green1","red", "blue3", "blue3", "blue3", "blue3", "blue3","blue3"  ), cex=0.7)

abline(h = 0, lty=3)
abline(v = 0, lty=3)

with(subset(antisense_ens, antisense_ens$V4 %in% lnc_significant_ens$V4),points(lnc_lfc, gene_lfc, pch=20, col="red",cex=0.9))

with(subset(antisense_ens, antisense_ens$V16 %in% both_significant_ens$V16),points(lnc_lfc, gene_lfc, pch=20, col="green1", cex=0.9))

with(subset(antisense_ens, antisense_ens$V4 %in% pain_antisense_ens$V4),points(lnc_lfc, gene_lfc, pch=7, col="blue3", cex=0.9))
with(subset(antisense_ens, antisense_ens$V4 %in% sodium_channel_antisense_ens$V4),points(lnc_lfc, gene_lfc, pch=21, col="blue3", cex=0.9))
with(subset(antisense_ens, antisense_ens$V4 %in% potassium_channel_antisense_ens$V4),points(lnc_lfc, gene_lfc, pch=22, col="blue3", cex=0.9))
with(subset(antisense_ens, antisense_ens$V4 %in% chloride_channel_antisense_ens$V4),points(lnc_lfc, gene_lfc, pch=23, col="blue3", cex=0.9))
with(subset(antisense_ens, antisense_ens$V4 %in% calcium_channel_antisense_ens$V4),points(lnc_lfc, gene_lfc, pch=24, col="blue3", cex=0.9))
with(subset(antisense_ens, antisense_ens$V4 %in% trp_channel_antisense_ens$V4),points(lnc_lfc, gene_lfc, pch=25, col="blue3", cex=0.9))

dev.off()

sig_intergenic <- res_NEURONS_vs_IPS_AD2[rownames(res_NEURONS_vs_IPS_AD2) %in% lincs_closest$lincRNA & res_NEURONS_vs_IPS_AD2$padj < 0.05 & !is.na(res_NEURONS_vs_IPS_AD2$padj),]
sig_antisense <- res_NEURONS_vs_IPS_AD2[rownames(res_NEURONS_vs_IPS_AD2) %in% antisense$ID & res_NEURONS_vs_IPS_AD2$padj < 0.05 & !is.na(res_NEURONS_vs_IPS_AD2$padj),]

sig_intergenic_ens <- res_NEURONS_vs_IPS_AD2[rownames(res_NEURONS_vs_IPS_AD2) %in% lincs_closest_ens$lincRNA & res_NEURONS_vs_IPS_AD2$padj < 0.05 & !is.na(res_NEURONS_vs_IPS_AD2$padj),]
sig_antisense_ens <- res_NEURONS_vs_IPS_AD2[rownames(res_NEURONS_vs_IPS_AD2) %in% antisense_ens$V4 & res_NEURONS_vs_IPS_AD2$padj < 0.05 & !is.na(res_NEURONS_vs_IPS_AD2$padj),]

sig_intergenic <- lincs_closest[lincs_closest$lincRNA %in% rownames(sig_intergenic), ]
sig_antisense <- antisense[antisense$ID %in% rownames(sig_antisense), ]

sig_intergenic_ens <- lincs_closest_ens[lincs_closest_ens$lincRNA %in% rownames(sig_intergenic_ens), ]
sig_antisense_ens <- antisense_ens[antisense_ens$V4 %in% rownames(sig_antisense_ens), ]

sig_intergenic <- data.frame(LncRNA = sig_intergenic$lincRNA, LncRNA_symbol=sig_intergenic$LincRNA_symbol, closest_gene= sig_intergenic$ENSEMBL_gene_id, sig_intergenic$gene_symbol)
sig_intergenic_ens <- data.frame(LncRNA = sig_intergenic_ens$lincRNA, LncRNA_symbol=sig_intergenic_ens$LincRNA_symbol, closest_gene= sig_intergenic_ens$ENSEMBL_gene_id, sig_intergenic_ens$gene_symbol)


sig_annot <- unique(c(sig_intergenic$LncRNA, sig_intergenic_ens$LncRNA))


lmat <- rbind(c(0,4,5), c(0,2,1), c(0,3,0))
lhei <- c(2,8,2)
lwid <- c(2,8,1)



 
pdf(file = paste(PATH_out,"/Figure4_intergenic.pdf", sep=""))

heatmap.2(as.matrix(rld_mat[rownames(rld_mat) %in% sig_annot,]), col=redgreen(200), scale="row", trace="none", cexCol=1.2, cexRow=0.5, margins=c(1,1), key=FALSE, Colv=TRUE, Rowv=TRUE, main="Intergenic LncRNAs DE in Neurons vs IPSc", dendrogram="column", symkey=TRUE, srtCol=42, adjCol=c(1,1), lmat=lmat, lhei=lhei, lwid = lwid, labCol=c("AD2_IPS", "AD2_IPS", "AD4_IPS", "AD4_IPS", "NHDF_IPS", "NHDF_IPS", "AD2_NEURONS", "AD2_NEURONS", "AD4_NEURONS", "AD4_NEURONS", "NHDF_NEURONS", "NHDF_NEURONS")
, labRow=NA, RowSideColors=c(rep("green",length(unique(sig_intergenic$LncRNA))), rep("blue",length(unique(sig_intergenic_ens$LncRNA)))))

legend(x = 0.85, y = -0.1, xpd=TRUE, legend = c("Novel", "Annotated"), col = c("green", "blue"), lty= 1, lwd = 5, cex=0.9, bty = "n")

dev.off()

haglr_qpcr <- read.csv(file = paste0(PATH_ips,"/haglr_qpcr.csv"), header=FALSE)
haglr_rna_seq <- norm_counts["ENSG00000224189",]

haglr <- data.frame(HAGLR = haglr_qpcr$V2)

neurons_vs_ips <- data.frame(ID = rep("HAGLR",2) , condition =  c("IPSc", "Neurons"), expression = c(mean(haglr[1:4,1]), mean(haglr[5:8,1])), sem = c(std.error(haglr[1:4,1]), std.error(haglr[5:8,1])))
neurons_vs_ips <- transform(neurons_vs_ips, lower=expression-sem, upper=expression+sem)

neurons_vs_ips_rna_seq <- data.frame(ID = rep("HAGLR",2) , condition =  c("IPSc", "Neurons"), expression = c(mean(haglr_rna_seq[1:6]), mean(haglr_rna_seq[7:12])), sem = c(std.error(haglr_rna_seq[1:6]), std.error(haglr_rna_seq[7:12])))
neurons_vs_ips_rna_seq <- transform(neurons_vs_ips_rna_seq, lower=expression-sem, upper=expression+sem)

neurons_vs_ips$assay <- "qPCR"
neurons_vs_ips_rna_seq$assay <- "RNA-seq"

validation <- cbind(neurons_vs_ips, neurons_vs_ips_rna_seq)

p1 <- ggplot(neurons_vs_ips, aes(x = ID, y = expression, fill = condition)) +  geom_bar(stat="identity", position=position_dodge(), width=1) + geom_errorbar(aes(ymax=upper, ymin=lower), position=position_dodge(0.9), data=neurons_vs_ips, width = 0.5, size = 0.15) + xlab("Group average") + ylab("Relative expression + SEM") + ggtitle("qPCR") + theme(axis.text.y = element_text(size= 12, face="bold"), axis.title.y = element_text(size=15, face="bold"), axis.line.y = element_line(color="black", size = 0.3), axis.title.x = element_blank(), axis.text.x = element_text(size= 12, angle=45, hjust=1, face="bold"), legend.title=NULL, legend.text=element_text(size=15, face="bold"), legend.key.size=unit(0.5, "cm"), plot.title=element_text(size=15, face="bold", hjust = 0.5))

p2 <- ggplot(neurons_vs_ips_rna_seq, aes(x = ID, y = expression, fill = condition)) +  geom_bar(stat="identity", position=position_dodge(), width=1) + geom_errorbar(aes(ymax=upper, ymin=lower), position=position_dodge(0.9), data=neurons_vs_ips_rna_seq, width = 0.5, size = 0.15) + xlab("Group average") + ylab("Relative expression + SEM") + ggtitle("RNA-seq") + theme(axis.text.y = element_text(size= 12, face="bold"), axis.title.y = element_text(size=15, face="bold"), axis.line.y = element_line(color="black", size = 0.3), axis.title.x = element_blank(), axis.text.x = element_text(size= 12, angle=45, hjust=1, face="bold"), legend.title=NULL, legend.text=element_text(size=15, face="bold"), legend.key.size=unit(0.5, "cm"), plot.title=element_text(size=15, face="bold", hjust = 0.5))

pdf(file = paste(PATH_out,"/Figure_haglr.pdf", sep=""))
grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 2, heights = unit(c(0.5, 5), "null"))))

grid.text("Relative expression IPSc vs Neurons", vp = viewport(layout.pos.row = 1, layout.pos.col = 1:2), gp=gpar(fontsize=20, col="black", fontface="bold"))
print(p1, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))         
print(p2, vp = viewport(layout.pos.row = 2, layout.pos.col = 2))
dev.off()


#Data for figure 4
ips_fig4 <- data.frame(LncRNA_ID = c(unique(antisense_ips$ID), unique(intronic_ips$ID), unique(lincs_ips$lincRNA)), LncRNA_name = c(unique(antisense_ips$V4), unique(intronic_ips$V4), unique(lincs_ips$LincRNA_symbol)),  genomic_context = c(rep("Antisense", length(unique(antisense_ips$ID))), rep("Intronic", length(unique(intronic_ips$ID))), rep("Intergenic", length(unique(lincs_ips$lincRNA)))), nr_of_exons = nr_of_exons_rat[c(unique(antisense_ips$V4), unique(intronic_ips$V4), unique(lincs_ips$LincRNA_symbol))])

write.csv(file = paste(PATH_out,"/Figure_4_ips.csv", sep=""), ips_fig4)

write.csv(file = paste0(PATH_out, "/Figure_4_antisense_novel.csv"), antisense) 
write.csv(file = paste0(PATH_out, "/Figure_4_antisense_ENSEMBL.csv"), antisense_ens) 
write.csv(file = paste0(PATH_out, "/Figure_4_intergenic.csv"), rld_mat[rownames(rld_mat) %in% sig_annot,]) 


write.csv(file = paste0(PATH_out, "/Figure_haglr.csv"), validation) 


#######################################
############## Figure 5 ###############
#######################################
#A: Von Frey, B: PCA mouse, C: Volcano plot LncRNAs, D: Network modules
load(paste0(PATH_mouse,"/DE/norm_counts.RData"))
mouse_rld <- rld

load(paste0(PATH_rat,"/DE/norm_counts.RData"))
rat_rld <- rld

nc_ensembl_mouse <- read.csv(file = paste(PATH_mouse,"/mm10_genome/ens_lncs_annotation.csv", sep=""), row.name=1)

nc_ensembl_rat <- read.csv(file = paste(PATH_rat,"/rn6_genome/ensembl_nc.csv", sep=""), row.name=1)

d = plotPCA(mouse_rld, intgroup=c("condition", "strain"), returnData=TRUE, ntop = 10000)
percentVar <- round(100 * attr(d, "percentVar"))
p1 <- ggplot(d, aes(x=PC1,y=PC2, color=condition, shape=strain, label=colnames(mouse_rld))) + geom_point(size=8) + xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance")) + ggtitle("PCA SNI vs Sham Mouse DRG") + theme(axis.text.y = element_text(size= 17, face="bold"), axis.title.y = element_text(size=17, face="bold"), axis.title.x = element_text(size=17, face="bold"), axis.text.x = element_text(size= 14, angle=45, hjust=1, face="bold"), legend.title=NULL, legend.text=element_text(size=17, face="bold"), legend.key.size=unit(0.5, "cm"), plot.title=element_text(size=20, face="bold", hjust=0.5))

d = plotPCA(rat_rld, intgroup=c("condition"), returnData=TRUE, ntop = 5000)
percentVar <- round(100 * attr(d, "percentVar"))
p2 <- ggplot(d, aes(x=PC1,y=PC2, color=condition, label=colnames(rat_rld))) + geom_point(size=8) + xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance")) + ggtitle("PCA SNT vs Sham Rat DRG") + theme(axis.text.y = element_text(size= 17, face="bold"), axis.title.y = element_text(size=17, face="bold"), axis.title.x = element_text(size=17, face="bold"), axis.text.x = element_text(size= 14, angle=45, hjust=1, face="bold"), legend.title=NULL, legend.text=element_text(size=17, face="bold"), legend.key.size=unit(0.5, "cm"), plot.title=element_text(size=20, face="bold", hjust=0.5))

pdf(file = paste(PATH_out,"/Figure_5_PCA.pdf", sep=""), width = 12, height = 6)
print(multiplot(p1, p2, cols=2))
dev.off()

balbc <- read.csv(paste0(PATH_mouse, "/DE/res_sni_vs_sham_BALB.c.csv"), row.name=1)
balbc_fdr <- balbc[!is.na(balbc$padj),]
b10d2 <- read.csv(paste0(PATH_mouse, "/DE/res_sni_vs_sham_B10.D2.csv"), row.name=1)
b10d2_fdr <- b10d2[!is.na(b10d2$padj),]


pdf(file = paste0(PATH_out, "/Figure5_volcano_C.pdf"), height=4, width=12)
par(oma = c(0, 0, 2, 0), mar = c(2, 2, 2, 2), mgp = c(2, 1, 0), xpd = FALSE)
layout(matrix(c(1,2,3), 1, 3, byrow = TRUE))

with(balbc_fdr, plot(log2FoldChange, -log10(padj), pch=20,col="grey",main="Volcano plot LncRNAs SNI vs Sham in BALB.c Mouse", xlim=c(-4,4), ylim=c(0,30), xlab="Log 2 Fold Change", ylab="-log10(Adjusted P.value)"))
with(subset(balbc_fdr, padj < 0.05), points(log2FoldChange, -log10(padj), pch=20, col="lightgreen"))
with(subset(balbc_fdr, grepl("LncRNA", rownames(balbc_fdr)) & padj < 0.05), points(log2FoldChange, -log10(padj), pch=20, col="lightcoral"))
with(subset(balbc_fdr, padj < 0.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(padj), pch=20, col="green", cex=1.2))
with(subset(balbc_fdr, grepl("LncRNA", rownames(balbc_fdr)) & padj < 0.05 & abs(log2FoldChange)>1),points(log2FoldChange, -log10(padj), pch=20, col="red", cex=1.2))
legend("topleft", c("Genes with Adj.Pvalue < 0.05","& LFC > 1", "Novel LncRNAs with Adj.Pvalue < 0.05", "& LFC > 1"), pch=c(20,20,20), col=c("lightgreen","green", "lightcoral", "red"))

with(b10d2_fdr, plot(log2FoldChange, -log10(padj), pch=20,col="grey",main="Volcano plot LncRNAs SNI vs Sham in B10.D2 Mouse", xlim=c(-4,4), ylim=c(0,30), xlab="Log 2 Fold Change", ylab="-log10(Adjusted P.value)"))
with(subset(b10d2_fdr, padj < 0.05), points(log2FoldChange, -log10(padj), pch=20, col="lightgreen"))
with(subset(b10d2_fdr, grepl("LncRNA", rownames(b10d2_fdr)) & padj < 0.05), points(log2FoldChange, -log10(padj), pch=20, col="lightcoral"))
with(subset(b10d2_fdr, padj < 0.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(padj), pch=20, col="green", cex=1.2))
with(subset(b10d2_fdr, grepl("LncRNA", rownames(b10d2_fdr)) & padj < 0.05 & abs(log2FoldChange)>1),points(log2FoldChange, -log10(padj), pch=20, col="red", cex=1.2))
legend("topleft", c("Genes with Adj.Pvalue < 0.05","& LFC > 1", "Novel LncRNAs with Adj.Pvalue < 0.05", "& LFC > 1"), pch=c(20,20,20), col=c("lightgreen","green", "lightcoral", "red"))

rat <- read.csv(paste0(PATH_rat, "/DE/res_snt_vs_sham.csv"), row.name=1)
rat_fdr <- rat[!is.na(rat$padj),]

with(rat_fdr, plot(log2FoldChange, -log10(padj), pch=20,col="grey",main="Volcano plot LncRNAs SNT vs Sham in Rat DRG", xlim=c(-5,5), ylim=c(0,30), xlab="Log 2 Fold Change", ylab="-log10(Adjusted P.value)"))
with(subset(rat_fdr, padj < 0.05), points(log2FoldChange, -log10(padj), pch=20, col="lightgreen"))
with(subset(rat_fdr, grepl("LncRNA", rownames(rat_fdr)) & padj < 0.05), points(log2FoldChange, -log10(padj), pch=20, col="lightcoral"))
with(subset(rat_fdr, padj < 0.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(padj), pch=20, col="green", cex=1.2))
with(subset(rat_fdr, padj < 0.05 & abs(log2FoldChange)>1 & grepl("LncRNA", rownames(rat_fdr))),points(log2FoldChange, -log10(padj), pch=20, col="red", cex=1.2))
legend("topleft", c("Genes with Adj.Pvalue < 0.05","& LFC > 1", "Novel LncRNAs with Adj.Pvalue < 0.05", "& LFC > 1"), pch=c(20,20,20), col=c("lightgreen","green", "lightcoral", "red"))
dev.off()

#modules
load(paste0(PATH_mouse,"/Network/LncRNAs_MM.RData"))
LncRNA_modules <- read.csv(paste0(PATH_mouse,"/Network/LncRNAs_GO_module.csv"), row.name=1)

antisense_mouse <- read.csv(paste0(PATH_mouse, "/antisense/sni_vs_sham_balb.c/antisense_all_expression.csv"))

intronic_mouse <- read.csv(paste0(PATH_mouse, "/intronic/sni_vs_sham_balb.c/intronic_all_expression.csv"))

lincs_mouse <- read.csv(paste0(PATH_mouse, "/DE/lincs_closest.csv"))

antisense_MM <- MM[rownames(MM) %in% antisense_mouse$ID, ]
intronic_MM <- MM[rownames(MM) %in% intronic_mouse$ID, ]
intergenic_MM <- MM[rownames(MM) %in% lincs_mouse$lincRNA, ]

MM_annot <- rbind(antisense_MM, intronic_MM, intergenic_MM)

module_annotation[match(gsub("ME","", colnames(MM_annot)), module_annotation$module.colours),]$GO_module

module_annotation[module_annotation$GO_module=="proteasome-mediated ubiquitin-dependent protein catabolic process",]$GO_module <- "p.m. u.d protein catabolic process"

module_annotation[module_annotation$GO_module=="ribonucleoprotein complex subunit organization",]$GO_module <- "ribonucleoprotein comp. sub. organization"

lmat <- rbind(c(0,4,5), c(3,2,0), c(0,1,0))
lhei <- c(2,6,0.3)
lwid <- c(2,9,7)

hmcol <- colorRampPalette(brewer.pal(9,"OrRd"))(100)

pdf(file = paste0(PATH_out, "/figure5_MM_D.pdf"), heigh=7, width=7)

heatmap.2(as.matrix(t(MM_annot)), col=hmcol, scale="column", trace="none", cexRow=1.2, margins=c(1,1), key=TRUE, keysize=0.5, key.par = list(cex=0.5), Colv=TRUE, Rowv=TRUE, main="LncRNAs GO module membership", cex.main = 0.4, dendrogram="both", symkey=TRUE, srtCol=42, adjCol=c(1,1), labRow=module_annotation[match(gsub("ME","", colnames(MM_annot)), module_annotation$module.colours),]$GO_module, labCol=NA, ColSideColors=c(rep("blue3",nrow(antisense_MM)), rep("maroon3",nrow(intronic_MM)), rep("darkgoldenrod1",nrow(intergenic_MM))), lmat=lmat, lhei=lhei, lwid=lwid)


legend(x = 0.6, y = -0.12, xpd=TRUE, legend = c("Antisense", "Intergenic", "Intronic"), col = c("blue3", "maroon3", "darkgoldenrod1"), lty= 1, lwd = 5, cex=0.5, bty = "n")


dev.off()


module_freq <- as.data.frame(table(LncRNA_modules$GO_module))
names(module_freq) <- c("Module_GO", "frequency")

p2 <- ggplot(data=module_freq, aes(x=Module_GO, y=frequency)) + geom_bar(stat="identity", fill="steelblue") + geom_text(aes(label=frequency), vjust=0.5, hjust=0, color="black", size=2, fontface='bold') + ggtitle("GO module annotation for novel LncRNAs") + coord_flip() + theme(axis.text.y = element_text(size= 7, angle=25), axis.title.y = element_blank(), axis.title.x = element_text(size=15, face="bold"), axis.text.x = element_text(size= 10, angle=45, hjust=1, face="bold"), legend.title=NULL, legend.text=element_text(size=15, face="bold"), legend.key.size=unit(0.5, "cm"), plot.title=element_text(size=15, face="bold", hjust = 1.5))

pdf(file = paste0(PATH_out, "/figure5_Modules_E.pdf"), heigh=6, width=6)
print(p2)
dev.off()

#Data figure 5
rv <- rowVars(assay(mouse_rld))
top_10000 <- assay(mouse_rld)[order(-rv)[1:10000], ]
write.csv(file = paste0(PATH_out,"/Figure5_mousePCA.csv"), top_10000)
rv <- rowVars(assay(rat_rld))
top_5000 <- assay(rat_rld)[order(-rv)[1:5000], ]
write.csv(file = paste0(PATH_out,"/Figure5_ratPCA.csv"), top_5000)
write.csv(file = paste0(PATH_out, "/Figure5_MM.csv"), MM_annot)
write.csv(file = paste0(PATH_out, "/Figure5_modules.csv"), LncRNA_modules)


#######################################
############## Figure 6 ###############
#######################################

antisense_rat <- read.csv(paste0(PATH_rat, "/antisense/snt_vs_sham/antisense_all_expression.csv"), row.name=1)
antisense_ens_rat <- read.csv(paste0(PATH_rat, "/antisense/ensembl_antisense/snt_vs_sham/antisense_all_expression.csv"), row.name=1)

lincs_closest_rat <- read.csv(file = paste0(PATH_rat, "/DE/lincs_closest.csv"), row.name=1)
lincs_closest_ens_rat <- read.csv(file = paste0(PATH_rat, "/DE/lincs_ENSEMBL_closest.csv"), row.name=1)

pain_genes <- read.csv(paste0(PATH_mouse,"/../pain_genes_homologs.csv"), row.name=1)

load(file=paste0(PATH_rat, "/DE/norm_counts.RData"))
rat_rld <- rld_mat

both_significant_rat <- subset(antisense_rat, antisense_rat$gene_pvalue<.05 & antisense_rat$lnc_pvalue<.05)

lnc_significant_rat <- subset(antisense_rat, antisense_rat$lnc_pvalue<.05)

both_significant_ens_rat <- subset(antisense_ens_rat, antisense_ens_rat$gene_pvalue<.05 & antisense_ens_rat$lnc_pvalue<.05)

lnc_significant_ens_rat <- subset(antisense_ens_rat, antisense_ens_rat$lnc_pvalue<.05)

rat <- read.csv(file = paste0(PATH_rat, "/DE/res_snt_vs_sham.csv"), row.name=1)

pain_rat <- rat[rownames(rat) %in% pain_genes$Gene.ID.1,]
sodium_channel_rat <- rat[grep("Scn", rat$symbol),]
potassium_channel_rat <- rat[grep("Kcn", rat$symbol),]
chloride_channel_rat <- rat[grep("Clc", rat$symbol),]
calcium_channel_rat <- rat[grep("Cac", rat$symbol),]
trp_channel_rat <- rat[grep("Trp", rat$symbol),]

pain_antisense_rat <- subset(antisense_rat, antisense_rat$V16 %in% rownames(pain_rat))
sodium_channel_antisense_rat <- subset(antisense_rat, antisense_rat$symbol %in% sodium_channel_rat$symbol)
potassium_channel_antisense_rat <- subset(antisense_rat, antisense_rat$symbol %in% potassium_channel_rat$symbol)
chloride_channel_antisense_rat <- subset(antisense_rat, antisense_rat$symbol %in% chloride_channel_rat$symbol)
calcium_channel_antisense_rat <- subset(antisense_rat, antisense_rat$symbol %in% calcium_channel_rat$symbol)
trp_channel_antisense_rat <- subset(antisense_rat, antisense_rat$symbol %in% trp_channel_rat$symbol)

pain_antisense_ens_rat <- subset(antisense_ens_rat, antisense_ens_rat$V16 %in% rownames(pain_rat))
sodium_channel_antisense_ens_rat <- subset(antisense_ens_rat, antisense_ens_rat$sense_symbol %in% sodium_channel_rat$symbol)
potassium_channel_antisense_ens_rat <- subset(antisense_ens_rat, antisense_ens_rat$sense_symbol %in% potassium_channel_rat$symbol)
chloride_channel_antisense_ens_rat <- subset(antisense_ens_rat, antisense_ens_rat$sense_symbol %in% chloride_channel_rat$symbol)
calcium_channel_antisense_ens_rat <- subset(antisense_ens_rat, antisense_ens_rat$sense_symbol %in% calcium_channel_rat$symbol)
trp_channel_antisense_ens_rat <- subset(antisense_ens_rat, antisense_ens_rat$sense_symbol %in% trp_channel_rat$symbol)

###############################################################################

pdf(file = paste(PATH_out,"/Figure_6_antisense.pdf", sep=""), width = 6, height = 6)

par(oma = c(0, 0, 1, 0), mar = c(8, 3, 2, 2), mgp = c(2, 1, 0))
#layout(matrix(c(1,2), 1, 2, byrow = TRUE))

with(antisense_rat, plot(lnc_lfc, gene_lfc, pch=20, col="grey", cex=0.9, main="Novel antisense LncRNAs (Rat)", xlab="Novel lncRNA log2 fold change SNT vs Sham", ylab="Sense gene log2 fold change SNT vs Sham", cex.lab=0.7))
legend("topleft", legend = c("Both DE (p.value < 0.05)", "Only LncRNA DE", "Pain gene", "Sodium channel", "Potassium channel", "Chloride channel", "Calcium channel", "Trp channel"), pch=c(20,20,7, 21,22,23,24,25), col=c("green1","red", "blue3", "blue3", "blue3", "blue3", "blue3","blue3"  ), cex=0.7)

abline(h = 0, lty=3)
abline(v = 0, lty=3)

with(subset(antisense_rat, antisense_rat$ID %in% lnc_significant_rat$ID),points(lnc_lfc, gene_lfc, pch=20, col="red",cex=0.9))

with(subset(antisense_rat, antisense_rat$ID %in% both_significant_rat$ID),points(lnc_lfc, gene_lfc, pch=20, col="green1", cex=0.9))

with(subset(antisense_rat, antisense_rat$ID %in% pain_antisense_rat$ID),points(lnc_lfc, gene_lfc, pch=7, col="blue3", cex=0.9))
with(subset(antisense_rat, antisense_rat$ID %in% sodium_channel_antisense_rat$ID),points(lnc_lfc, gene_lfc, pch=21, col="blue3", cex=0.9))
with(subset(antisense_rat, antisense_rat$ID %in% potassium_channel_antisense_rat$ID),points(lnc_lfc, gene_lfc, pch=22, col="blue3", cex=0.9))
with(subset(antisense_rat, antisense_rat$ID %in% chloride_channel_antisense_rat$ID),points(lnc_lfc, gene_lfc, pch=23, col="blue3", cex=0.9))
with(subset(antisense_rat, antisense_rat$ID %in% calcium_channel_antisense_rat$ID),points(lnc_lfc, gene_lfc, pch=24, col="blue3", cex=0.9))
with(subset(antisense_rat, antisense_rat$ID %in% trp_channel_antisense_rat$ID),points(lnc_lfc, gene_lfc, pch=25, col="blue3", cex=0.9))

dev.off()

sig_intergenic <- rat[rownames(rat) %in% lincs_closest_rat$lincRNA & rat$padj < 0.05 & !is.na(rat$padj),]

sig_intergenic_ens <- rat[rownames(rat) %in% lincs_closest_ens_rat$lincRNA & rat$padj < 0.05 & !is.na(rat$padj),]

sig_intergenic <- lincs_closest_rat[lincs_closest_rat$lincRNA %in% rownames(sig_intergenic), ]
sig_intergenic_ens <- lincs_closest_ens_rat[lincs_closest_ens_rat$lincRNA %in% rownames(sig_intergenic_ens), ]

sig_intergenic <- data.frame(LncRNA = sig_intergenic$lincRNA, LncRNA_symbol=sig_intergenic$LincRNA_symbol, closest_gene= sig_intergenic$ENSEMBL_gene_id, sig_intergenic$gene_symbol)
sig_intergenic_ens <- data.frame(LncRNA = sig_intergenic_ens$lincRNA, LncRNA_symbol=sig_intergenic_ens$LincRNA_symbol, closest_gene= sig_intergenic_ens$ENSEMBL_gene_id, sig_intergenic_ens$gene_symbol)


sig_annot <- unique(c(sig_intergenic$LncRNA, sig_intergenic_ens$LncRNA))


lmat <- rbind(c(0,4,5), c(0,2,1), c(0,3,0))
lhei <- c(2,8,2)
lwid <- c(0.3,8,1)

 
pdf(file = paste(PATH_out,"/Figure6_intergenic.pdf", sep=""), width = 6, height = 6)

heatmap.2(as.matrix(rat_rld[rownames(rat_rld) %in% sig_annot,]), col=redgreen(200), scale="row", trace="none", cexCol=1.4, cexRow=0.5, cex.main=1, margins=c(1,1), key=FALSE, Colv=TRUE, Rowv=TRUE, main="Intergenic LncRNAs DE in SNT vs Sham (Rat DRG)", dendrogram="column", symkey=TRUE, srtCol=42, adjCol=c(1,1), lmat=lmat, lhei=lhei, lwid = lwid, labCol=c(rep("Sham",4),rep("SNT",4)), labRow=NA, RowSideColors=c(rep("green",length(unique(sig_intergenic$LncRNA))), rep("blue",length(unique(sig_intergenic_ens$LncRNA)))))

legend(x = 0.8, y = -0.1, xpd=TRUE, legend = c("Novel", "Annotated"), col = c("green", "blue"), lty= 1, lwd = 5, cex=0.9, bty = "n")

dev.off()

#Data for figure 6
write.csv(file = paste0(PATH_out,"/Figure6_antisense.csv"), antisense_rat)
write.csv(file = paste0(PATH_out, "/Figure6_intergenic.csv"), rat_rld[rownames(rat_rld) %in% sig_annot,])


#######################################
############## Figure 7 ###############
#######################################
#Mouse Pain model vs Sham
pain_genes <- read.csv(paste0(PATH_mouse,"/../pain_genes_homologs.csv"), row.name=1)

antisense <- read.csv(paste0(PATH_mouse, "/antisense/sni_vs_sham_balb.c/antisense_all_expression.csv"), row.name=1)
antisense_ens <- read.csv(paste0(PATH_mouse, "/antisense/ensembl_antisense/sni_vs_sham_balb.c/antisense_all_expression.csv"), row.name=1)

antisense_b10d2 <- read.csv(paste0(PATH_mouse, "/antisense/sni_vs_sham_b10.d2/antisense_all_expression.csv"), row.name=1)
antisense_ens_b10d2 <- read.csv(paste0(PATH_mouse, "/antisense/ensembl_antisense/sni_vs_sham_b10.d2/antisense_all_expression.csv"), row.name=1)

lincs_closest <- read.csv(file = paste0(PATH_mouse, "/DE/lincs_closest.csv"), row.name=1)
lincs_closest_ens <- read.csv(file = paste0(PATH_mouse, "/DE/lincs_ENSEMBL_closest.csv"), row.name=1)

load(file=paste0(PATH_mouse, "/DE/norm_counts.RData"))
mouse_rld <- rld_mat

both_significant <- subset(antisense, antisense$gene_pvalue<.05 & antisense$lnc_pvalue<.05)

lnc_significant <- subset(antisense, antisense$lnc_pvalue<.05)

both_significant_ens <- subset(antisense_ens, antisense_ens$gene_pvalue<.05 & antisense_ens$lnc_pvalue<.05)

lnc_significant_ens <- subset(antisense_ens, antisense_ens$lnc_pvalue<.05)

balbc <- read.csv(file = paste0(PATH_mouse, "/DE/res_sni_vs_sham_BALB.c.csv"), row.name=1)

pain <- balbc[rownames(balbc) %in% pain_genes$Gene.ID,]
sodium_channel <- balbc[grep("Scn", balbc$symbol),]
potassium_channel <- balbc[grep("Kcn", balbc$symbol),]
chloride_channel <- balbc[grep("Clc", balbc$symbol),]
calcium_channel <- balbc[grep("Cac", balbc$symbol),]
trp_channel <- balbc[grep("Trp", balbc$symbol),]

pain_antisense <- subset(antisense, antisense$V16 %in% rownames(pain))
sodium_channel_antisense <- subset(antisense, antisense$symbol %in% sodium_channel$symbol)
potassium_channel_antisense <- subset(antisense, antisense$symbol %in% potassium_channel$symbol)
chloride_channel_antisense <- subset(antisense, antisense$symbol %in% chloride_channel$symbol)
calcium_channel_antisense <- subset(antisense, antisense$symbol %in% calcium_channel$symbol)
trp_channel_antisense <- subset(antisense, antisense$symbol %in% trp_channel$symbol)

pain_antisense_ens <- subset(antisense_ens, antisense_ens$V16 %in% rownames(pain))
sodium_channel_antisense_ens <- subset(antisense_ens, antisense_ens$sense_symbol %in% sodium_channel$symbol)
potassium_channel_antisense_ens <- subset(antisense_ens, antisense_ens$sense_symbol %in% potassium_channel$symbol)
chloride_channel_antisense_ens <- subset(antisense_ens, antisense_ens$sense_symbol %in% chloride_channel$symbol)
calcium_channel_antisense_ens <- subset(antisense_ens, antisense_ens$sense_symbol %in% calcium_channel$symbol)
trp_channel_antisense_ens <- subset(antisense_ens, antisense_ens$sense_symbol %in% trp_channel$symbol)

###################################################################################################################
both_significant_b10d2 <- subset(antisense_b10d2, antisense_b10d2$gene_pvalue<.05 & antisense_b10d2$lnc_pvalue<.05)

lnc_significant_b10d2 <- subset(antisense_b10d2, antisense_b10d2$lnc_pvalue<.05)

both_significant_ens_b10d2 <- subset(antisense_ens_b10d2, antisense_ens_b10d2$gene_pvalue<.05 & antisense_ens_b10d2$lnc_pvalue<.05)

lnc_significant_ens_b10d2 <- subset(antisense_ens_b10d2, antisense_ens_b10d2$lnc_pvalue<.05)

b10d2 <- read.csv(file = paste0(PATH_mouse, "/DE/res_sni_vs_sham_B10.D2.csv"), row.name=1)

pain_b10d2 <- b10d2[rownames(b10d2) %in% pain_genes$Gene.ID,]
sodium_channel_b10d2 <- b10d2[grep("Scn", b10d2$symbol),]
potassium_channel_b10d2 <- b10d2[grep("Kcn", b10d2$symbol),]
chloride_channel_b10d2 <- b10d2[grep("Clc", b10d2$symbol),]
calcium_channel_b10d2 <- b10d2[grep("Cac", b10d2$symbol),]
trp_channel_b10d2 <- b10d2[grep("Trp", b10d2$symbol),]

pain_antisense_b10d2 <- subset(antisense_b10d2, antisense_b10d2$V16 %in% rownames(pain_b10d2))
sodium_channel_antisense_b10d2 <- subset(antisense_b10d2, antisense_b10d2$symbol %in% sodium_channel_b10d2$symbol)
potassium_channel_antisense_b10d2 <- subset(antisense_b10d2, antisense_b10d2$symbol %in% potassium_channel_b10d2$symbol)
chloride_channel_antisense_b10d2 <- subset(antisense_b10d2, antisense_b10d2$symbol %in% chloride_channel_b10d2$symbol)
calcium_channel_antisense_b10d2 <- subset(antisense_b10d2, antisense_b10d2$symbol %in% calcium_channel_b10d2$symbol)
trp_channel_antisense_b10d2 <- subset(antisense_b10d2, antisense_b10d2$symbol %in% trp_channel_b10d2$symbol)

pain_antisense_ens_b10d2 <- subset(antisense_ens_b10d2, antisense_ens_b10d2$V16 %in% rownames(pain_b10d2))
sodium_channel_antisense_ens_b10d2 <- subset(antisense_ens_b10d2, antisense_ens_b10d2$sense_symbol %in% sodium_channel_b10d2$symbol)
potassium_channel_antisense_ens_b10d2 <- subset(antisense_ens_b10d2, antisense_ens_b10d2$sense_symbol %in% potassium_channel_b10d2$symbol)
chloride_channel_antisense_ens_b10d2 <- subset(antisense_ens_b10d2, antisense_ens_b10d2$sense_symbol %in% chloride_channel_b10d2$symbol)
calcium_channel_antisense_ens_b10d2 <- subset(antisense_ens_b10d2, antisense_ens_b10d2$sense_symbol %in% calcium_channel_b10d2$symbol)
trp_channel_antisense_ens_b10d2 <- subset(antisense_ens_b10d2, antisense_ens_b10d2$sense_symbol %in% trp_channel_b10d2$symbol)

#########################################################################################################################




pdf(file = paste0(PATH_out,"/Figure_7_antisense.pdf"), width = 12, height = 12)
#dev.new(width=5, height=5)
par(oma = c(0, 0, 1, 0), mar = c(8, 3, 2, 2), mgp = c(2, 1, 0))
layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE))


with(antisense_ens, plot(lnc_lfc, gene_lfc, pch=20, col="grey", cex=1, main="ENSEMBL antisense LncRNAs (BALB/c mouse)", xlab="Novel lncRNA log2 fold change SNI vs Sham", ylab="Sense gene log2 fold change SNI vs Sham", cex.lab=1))
legend("topleft", legend = c("Both DE (p.value < 0.05)", "Only LncRNA DE", "Pain gene", "Sodium channel", "Potassium channel", "Chloride channel", "Calcium channel", "Trp channel"), pch=c(20,20,7, 21,22,23,24,25), col=c("green1","red", "blue3", "blue3", "blue3", "blue3", "blue3","blue3"  ), cex=0.7)

abline(h = 0, lty=3)
abline(v = 0, lty=3)

with(subset(antisense_ens, antisense_ens$V4 %in% lnc_significant_ens$V4),points(lnc_lfc, gene_lfc, pch=20, col="red",cex=1))

with(subset(antisense_ens, antisense_ens$V16 %in% both_significant_ens$V16),points(lnc_lfc, gene_lfc, pch=20, col="green1", cex=1))

with(subset(antisense_ens, antisense_ens$V4 %in% pain_antisense_ens$V4),points(lnc_lfc, gene_lfc, pch=7, col="blue3", cex=1))
with(subset(antisense_ens, antisense_ens$V4 %in% sodium_channel_antisense_ens$V4),points(lnc_lfc, gene_lfc, pch=21, col="blue3", cex=1))
with(subset(antisense_ens, antisense_ens$V4 %in% potassium_channel_antisense_ens$V4),points(lnc_lfc, gene_lfc, pch=22, col="blue3", cex=1))
with(subset(antisense_ens, antisense_ens$V4 %in% chloride_channel_antisense_ens$V4),points(lnc_lfc, gene_lfc, pch=23, col="blue3", cex=1))
with(subset(antisense_ens, antisense_ens$V4 %in% calcium_channel_antisense_ens$V4),points(lnc_lfc, gene_lfc, pch=24, col="blue3", cex=1))
with(subset(antisense_ens, antisense_ens$V4 %in% trp_channel_antisense_ens$V4),points(lnc_lfc, gene_lfc, pch=25, col="blue3", cex=1))

with(antisense_ens_b10d2, plot(lnc_lfc, gene_lfc, pch=20, col="grey", cex=1, main="ENSEMBL antisense LncRNAs (B10.D2 mouse)", xlab="Novel lncRNA log2 fold change SNI vs Sham", ylab="Sense gene log2 fold change SNI vs Sham", cex.lab=1))
legend("topleft", legend = c("Both DE (p.value < 0.05)", "Only LncRNA DE", "Pain gene", "Sodium channel", "Potassium channel", "Chloride channel", "Calcium channel", "Trp channel"), pch=c(20,20,7, 21,22,23,24,25), col=c("green1","red", "blue3", "blue3", "blue3", "blue3", "blue3","blue3"  ), cex=0.7)

abline(h = 0, lty=3)
abline(v = 0, lty=3)

with(subset(antisense_ens_b10d2, antisense_ens_b10d2$V4 %in% lnc_significant_ens_b10d2$V4),points(lnc_lfc, gene_lfc, pch=20, col="red",cex=1))

with(subset(antisense_ens_b10d2, antisense_ens_b10d2$V16 %in% both_significant_ens_b10d2$V16),points(lnc_lfc, gene_lfc, pch=20, col="green1", cex=1))

with(subset(antisense_ens_b10d2, antisense_ens_b10d2$V4 %in% pain_antisense_ens_b10d2$V4),points(lnc_lfc, gene_lfc, pch=7, col="blue3", cex=1))
with(subset(antisense_ens_b10d2, antisense_ens_b10d2$V4 %in% sodium_channel_antisense_ens_b10d2$V4),points(lnc_lfc, gene_lfc, pch=21, col="blue3", cex=1))
with(subset(antisense_ens_b10d2, antisense_ens_b10d2$V4 %in% potassium_channel_antisense_ens_b10d2$V4),points(lnc_lfc, gene_lfc, pch=22, col="blue3", cex=1))
with(subset(antisense_ens_b10d2, antisense_ens_b10d2$V4 %in% chloride_channel_antisense_ens_b10d2$V4),points(lnc_lfc, gene_lfc, pch=23, col="blue3", cex=1))
with(subset(antisense_ens_b10d2, antisense_ens_b10d2$V4 %in% calcium_channel_antisense_ens_b10d2$V4),points(lnc_lfc, gene_lfc, pch=24, col="blue3", cex=1))
with(subset(antisense_ens_b10d2, antisense_ens_b10d2$V4 %in% trp_channel_antisense_ens_b10d2$V4),points(lnc_lfc, gene_lfc, pch=25, col="blue3", cex=1))

##################################################
with(antisense, plot(lnc_lfc, gene_lfc, pch=20, col="grey", cex=1, main="Novel antisense LncRNAs (BALB/c mouse)", xlab="Novel lncRNA log2 fold change SNI vs Sham", ylab="Sense gene log2 fold change SNI vs Sham", cex.lab=1))
legend("topleft", legend = c("Both DE (p.value < 0.05)", "Only LncRNA DE", "Pain gene", "Sodium channel", "Potassium channel", "Chloride channel", "Calcium channel", "Trp channel"), pch=c(20,20,7, 21,22,23,24,25), col=c("green1","red", "blue3", "blue3", "blue3", "blue3", "blue3","blue3"  ), cex=0.7)

abline(h = 0, lty=3)
abline(v = 0, lty=3)

with(subset(antisense, antisense$ID %in% lnc_significant$ID),points(lnc_lfc, gene_lfc, pch=20, col="red",cex=1))

with(subset(antisense, antisense$ID %in% both_significant$ID),points(lnc_lfc, gene_lfc, pch=20, col="green1", cex=1))

with(subset(antisense, antisense$ID %in% pain_antisense$ID),points(lnc_lfc, gene_lfc, pch=7, col="blue3", cex=1))
with(subset(antisense, antisense$ID %in% sodium_channel_antisense$ID),points(lnc_lfc, gene_lfc, pch=21, col="blue3", cex=1))
with(subset(antisense, antisense$ID %in% potassium_channel_antisense$ID),points(lnc_lfc, gene_lfc, pch=22, col="blue3", cex=1))
with(subset(antisense, antisense$ID %in% chloride_channel_antisense$ID),points(lnc_lfc, gene_lfc, pch=23, col="blue3", cex=1))
with(subset(antisense, antisense$ID %in% calcium_channel_antisense$ID),points(lnc_lfc, gene_lfc, pch=24, col="blue3", cex=1))
with(subset(antisense, antisense$ID %in% trp_channel_antisense$ID),points(lnc_lfc, gene_lfc, pch=25, col="blue3", cex=1))


with(antisense_b10d2, plot(lnc_lfc, gene_lfc, pch=20, col="grey", cex=1, main="Novel antisense LncRNAs (B10.D2 mouse)", xlab="Novel lncRNA log2 fold change SNI vs Sham", ylab="Sense gene log2 fold change SNI vs Sham", cex.lab=1))
legend("topleft", legend = c("Both DE (p.value < 0.05)", "Only LncRNA DE", "Pain gene", "Sodium channel", "Potassium channel", "Chloride channel", "Calcium channel", "Trp channel"), pch=c(20,20,7, 21,22,23,24,25), col=c("green1","red", "blue3", "blue3", "blue3", "blue3", "blue3","blue3"  ), cex=0.7)

abline(h = 0, lty=3)
abline(v = 0, lty=3)

with(subset(antisense_b10d2, antisense_b10d2$ID %in% lnc_significant_b10d2$ID),points(lnc_lfc, gene_lfc, pch=20, col="red",cex=1))

with(subset(antisense_b10d2, antisense_b10d2$ID %in% both_significant_b10d2$ID),points(lnc_lfc, gene_lfc, pch=20, col="green1", cex=1))

with(subset(antisense_b10d2, antisense_b10d2$ID %in% pain_antisense_b10d2$ID),points(lnc_lfc, gene_lfc, pch=7, col="blue3", cex=1))
with(subset(antisense_b10d2, antisense_b10d2$ID %in% sodium_channel_antisense_b10d2$ID),points(lnc_lfc, gene_lfc, pch=21, col="blue3", cex=1))
with(subset(antisense_b10d2, antisense_b10d2$ID %in% potassium_channel_antisense_b10d2$ID),points(lnc_lfc, gene_lfc, pch=22, col="blue3", cex=1))
with(subset(antisense_b10d2, antisense_b10d2$ID %in% chloride_channel_antisense_b10d2$ID),points(lnc_lfc, gene_lfc, pch=23, col="blue3", cex=1))
with(subset(antisense_b10d2, antisense_b10d2$ID %in% calcium_channel_antisense_b10d2$ID),points(lnc_lfc, gene_lfc, pch=24, col="blue3", cex=1))
with(subset(antisense_b10d2, antisense_b10d2$ID %in% trp_channel_antisense_b10d2$ID),points(lnc_lfc, gene_lfc, pch=25, col="blue3", cex=1))

dev.off()




intergenic <- balbc[rownames(balbc) %in% lincs_closest$lincRNA,]


intergenic_ens <- balbc[rownames(balbc) %in% lincs_closest_ens$lincRNA,]


intergenic <- lincs_closest[lincs_closest$lincRNA %in% rownames(intergenic), ]


intergenic_ens <- lincs_closest_ens[lincs_closest_ens$lincRNA %in% rownames(intergenic_ens), ]


intergenic <- data.frame(LncRNA = intergenic$lincRNA, LncRNA_symbol=intergenic$LincRNA_symbol, closest_gene= intergenic$ENSEMBL_gene_id, symbol=intergenic$gene_symbol)
intergenic_ens <- data.frame(LncRNA = intergenic_ens$lincRNA, LncRNA_symbol=intergenic_ens$LincRNA_symbol, closest_gene= intergenic_ens$ENSEMBL_gene_id, symbol=intergenic_ens$gene_symbol)


intergenic_all <- unique(c(intergenic$LncRNA, intergenic_ens$LncRNA))


lmat <- rbind(c(0,4,5), c(0,2,1), c(0,3,0))
lhei <- c(4,8,1.5)
lwid <- c(2,8,1)

pdf(file = paste(PATH_out,"/Figure7_intergenic_balbc.pdf", sep=""))

heatmap.2(as.matrix(rld_mat[rownames(rld_mat) %in% intergenic_all,]), col=redgreen(200), scale="row", trace="none", cexCol=1.2, cexRow=0.5, margins=c(2,2), key=FALSE, Colv=TRUE, Rowv=TRUE, main="Intergenic LncRNAs in Mouse DRG", dendrogram="column", symkey=TRUE, srtCol=42, adjCol=c(1,1), lmat=lmat, lhei=lhei, lwid = lwid, labCol = substr(colnames(rld_mat), 10, nchar(colnames(rld_mat))), labRow=NA, RowSideColors=c(rep("green",length(unique(intergenic$LncRNA))), rep("blue",length(unique(intergenic_ens$LncRNA)))))

legend(x = 0.85, y = -0.1, xpd=TRUE, legend = c("Novel", "Annotated"), col = c("green", "blue"), lty= 1, lwd = 5, cex=0.9, bty = "n")

dev.off()

lncs_all <- c(intergenic_all, antisense_ens$V4, antisense$V4)
eset <- t(rld_mat[rownames(rld_mat) %in% lncs_all,])

distsRL <- dist(eset)
hc <- hclust(distsRL)
labels(hc) <- substr(labels(hc), 10, nchar(labels(hc)))

library(dendextend)
dend <- as.dendrogram(hc)

# Color the branches based on the clusters:
dend <- color_branches(dend, k=4, col = c("magenta", "blue", "magenta", "blue"))
 
dend <- hang.dendrogram(dend,hang_height=2)
dend <- set(dend, "labels_cex", 0.75)
#labels(dend) <- substr(colnames(rld_mat), 10, nchar(colnames(rld_mat)))
pdf(file = paste0(PATH_out, "/Figure7_clustering.pdf"))
par(mar = c(2,2,2,4), cex=1.2)
plot(dend, 
     main = "Samples clustering based on LncRNAs", 
     horiz =  TRUE,  nodePar = list(cex = .01), cex.main = 0.95)
dev.off()




#Data for figure 7
write.csv(file = paste0(PATH_out,"/Figure7_antisense_balbc.csv"), antisense)
write.csv(file = paste0(PATH_out,"/Figure7_antisense_b10d2.csv"), antisense_b10d2)
write.csv(file = paste0(PATH_out,"/Figure7_antisense_balbc_ensembl.csv"), antisense_ens)
write.csv(file = paste0(PATH_out,"/Figure7_antisense_b10d2_ensembl.csv"), antisense_ens_b10d2)
write.csv(file = paste0(PATH_out, "/Figure7_intergenic.csv"), rld_mat[rownames(rld_mat) %in% intergenic_all,])
write.csv(file = paste0(PATH_out, "/Figure7_clustering.csv"), rld_mat[rownames(rld_mat) %in% lncs_all,])


#######################################
############## Figure 9 ###############
#######################################
load(file=paste0(PATH_mouse, "/DE/norm_counts.RData"))
qpcr <- read.csv(file = paste0(PATH_mouse,"/qpcr_sni_vs_sham.csv"), row.name=1)
rna_seq <- norm_counts[rownames(qpcr),]

qpcr_balbc <- qpcr[,grep("BALB.c", names(qpcr))]

qpcr_b10d2 <- qpcr[,grep("B10.D2", names(qpcr))]

rna_seq_balbc <- rna_seq[,grep("BALB.c", colnames(rna_seq))]
rna_seq_balbc <- sweep( rna_seq_balbc, 1, apply(rna_seq_balbc[,c(1:6)], 1 ,mean), `/`) 

rna_seq_b10d2 <- rna_seq[,grep("B10.D2", colnames(rna_seq))]
rna_seq_b10d2 <- sweep( rna_seq_b10d2, 1, apply(rna_seq_b10d2[,c(1:5)], 1 ,mean), `/`) 

qpcr_balbc <- qpcr_balbc[, colnames(qpcr_balbc) %in% colnames(rna_seq_balbc)]

qpcr_b10d2 <- qpcr_b10d2[, colnames(qpcr_b10d2) %in% colnames(rna_seq_b10d2)]

qpcr_balbc <- sweep(qpcr_balbc, 1, apply(qpcr_balbc[,c(1:6)], 1 ,mean), `/`) 

qpcr_b10d2 <- sweep(qpcr_b10d2, 1, apply(qpcr_b10d2[,c(1:5)], 1 ,mean), `/`) 

colData_qpcr_balbc <- data.frame(condition = c(rep("SHAM",6), rep("SNI",4)), sex = c(rep("M",2), rep("F",3), rep("M",2), rep("F",2), "M"))


colData_qpcr_b10d2 <- data.frame(condition = c(rep("SHAM",5), rep("SNI",5)), sex = c(rep("M",2), rep("F",2),"M", rep("F",3), rep("M",2)))




p.value_balbc <- NULL

for(i in 1:nrow(qpcr)) {
p.value_balbc[i] <- anova(lm(as.numeric(qpcr_balbc[i,]) ~ colData_qpcr_balbc$sex + colData_qpcr_balbc$condition))[[5]][2]
}

p.value_b10d2 <- NULL

for(i in 1:nrow(qpcr)) {
p.value_b10d2[i] <- anova(lm(as.numeric(qpcr_b10d2[i,]) ~ colData_qpcr_b10d2$sex + colData_qpcr_b10d2$condition))[[5]][2]
}

qpcr_balbc$p.value <- p.value_balbc

qpcr_b10d2$p.value <- p.value_b10d2


qpcr_balbc <- data.frame(ID = rep(rownames(qpcr_balbc),2) , condition =  c(rep("SHAM",length(rownames(qpcr_balbc))), rep("SNI",length(rownames(qpcr_balbc)))), expression = c(apply(qpcr_balbc[,c(1:6)],1, mean, na.rm = TRUE), apply(qpcr_balbc[,c(7:10)],1, mean, na.rm = TRUE)), sem = c(apply(qpcr_balbc[,c(1:6)],1, std.error, na.rm = TRUE), apply(qpcr_balbc[,c(7:10)],1, std.error, na.rm = TRUE)))

qpcr_balbc <- transform(qpcr_balbc, lower=expression-sem, upper=expression+sem)

rna_seq_balbc <- data.frame(ID = rep(rownames(rna_seq_balbc),2) , condition =  c(rep("SHAM",length(rownames(rna_seq_balbc))), rep("SNI",length(rownames(rna_seq_balbc)))), expression = c(apply(rna_seq_balbc[,c(1:6)],1, mean, na.rm = TRUE), apply(rna_seq_balbc[,c(7:10)],1, mean, na.rm = TRUE)), sem = c(apply(rna_seq_balbc[,c(1:6)],1, std.error, na.rm = TRUE), apply(rna_seq_balbc[,c(7:10)],1, std.error, na.rm = TRUE)))

rna_seq_balbc <- transform(rna_seq_balbc, lower=expression-sem, upper=expression+sem)

qpcr_b10d2 <- data.frame(ID = rep(rownames(qpcr_b10d2),2) , condition =  c(rep("SHAM",length(rownames(qpcr_b10d2))), rep("SNI",length(rownames(qpcr_b10d2)))), expression = c(apply(qpcr_b10d2[,c(1:5)],1, mean, na.rm = TRUE), apply(qpcr_b10d2[,c(6:10)],1, mean, na.rm = TRUE)), sem = c(apply(qpcr_b10d2[,c(1:5)],1, std.error, na.rm = TRUE), apply(qpcr_b10d2[,c(6:10)],1, std.error, na.rm = TRUE)))

qpcr_b10d2 <- transform(qpcr_b10d2, lower=expression-sem, upper=expression+sem)

rna_seq_b10d2 <- data.frame(ID = rep(rownames(rna_seq_b10d2),2) , condition =  c(rep("SHAM",length(rownames(rna_seq_b10d2))), rep("SNI",length(rownames(rna_seq_b10d2)))), expression = c(apply(rna_seq_b10d2[,c(1:5)],1, mean, na.rm = TRUE), apply(rna_seq_b10d2[,c(6:10)],1, mean, na.rm = TRUE)), sem = c(apply(rna_seq_b10d2[,c(1:5)],1, std.error, na.rm = TRUE), apply(rna_seq_b10d2[,c(6:10)],1, std.error, na.rm = TRUE)))

rna_seq_b10d2 <- transform(rna_seq_b10d2, lower=expression-sem, upper=expression+sem)


p1 <- ggplot(qpcr_balbc, aes(x = ID, y = expression, fill = condition)) +  geom_bar(stat="identity", position=position_dodge(), width=1) + geom_errorbar(aes(ymax=upper, ymin=lower), position=position_dodge(0.9), data=qpcr_balbc, width = 0.5, size = 0.15) + xlab("Group average") + ylab("Relative expression (qPCR) + SEM") + ggtitle("Relative expression SNI vs sham (BALB/c mouse)") + theme(axis.text.y = element_text(size= 12, face="bold"), axis.title.y = element_text(size=15, face="bold"), axis.line.y = element_line(color="black", size = 0.3), axis.title.x = element_blank(), axis.text.x = element_text(size= 12, angle=45, hjust=1, face="bold"), legend.title=NULL, legend.text=element_text(size=15, face="bold"), legend.key.size=unit(0.5, "cm"), plot.title=element_text(size=15, face="bold", hjust = 0.5)) + scale_y_continuous(breaks=seq(0, 3, 0.25))

p2 <- ggplot(rna_seq_balbc, aes(x = ID, y = expression, fill = condition)) +  geom_bar(stat="identity", position=position_dodge(), width=1) + geom_errorbar(aes(ymax=upper, ymin=lower), position=position_dodge(0.9), data=rna_seq_balbc, width = 0.5, size = 0.15) + xlab("Group average") + ylab("Relative expression (RNA-seq) + SEM") + ggtitle("Relative expression SNI vs sham (BALB/c mouse)") + theme(axis.text.y = element_text(size= 12, face="bold"), axis.title.y = element_text(size=15, face="bold"), axis.line.y = element_line(color="black", size = 0.3), axis.title.x = element_blank(), axis.text.x = element_text(size= 12, angle=45, hjust=1, face="bold"), legend.title=NULL, legend.text=element_text(size=15, face="bold"), legend.key.size=unit(0.5, "cm"), plot.title=element_text(size=15, face="bold", hjust = 0.5)) + scale_y_continuous(breaks=seq(0, 3, 0.25))


p3 <- ggplot(qpcr_b10d2, aes(x = ID, y = expression, fill = condition)) +  geom_bar(stat="identity", position=position_dodge(), width=1) + geom_errorbar(aes(ymax=upper, ymin=lower), position=position_dodge(0.9), data=qpcr_b10d2, width = 0.5, size = 0.15) + xlab("Group average") + ylab("Relative expression (qPCR) + SEM") + ggtitle("Relative expression SNI vs sham (B10.D2 mouse)") + theme(axis.text.y = element_text(size= 12, face="bold"), axis.title.y = element_text(size=15, face="bold"), axis.line.y = element_line(color="black", size = 0.3), axis.title.x = element_blank(), axis.text.x = element_text(size= 12, angle=45, hjust=1, face="bold"), legend.title=NULL, legend.text=element_text(size=15, face="bold"), legend.key.size=unit(0.5, "cm"), plot.title=element_text(size=15, face="bold", hjust = 0.5)) + scale_y_continuous(breaks=seq(0, 3, 0.25))

p4 <- ggplot(rna_seq_b10d2, aes(x = ID, y = expression, fill = condition)) +  geom_bar(stat="identity", position=position_dodge(), width=1) + geom_errorbar(aes(ymax=upper, ymin=lower), position=position_dodge(0.9), data=rna_seq_b10d2, width = 0.5, size = 0.15) + xlab("Group average") + ylab("Relative expression (RNA-seq) + SEM") + ggtitle("Relative expression SNI vs sham (B10.D2 mouse)") + theme(axis.text.y = element_text(size= 12, face="bold"), axis.title.y = element_text(size=15, face="bold"), axis.line.y = element_line(color="black", size = 0.3), axis.title.x = element_blank(), axis.text.x = element_text(size= 12, angle=45, hjust=1, face="bold"), legend.title=NULL, legend.text=element_text(size=15, face="bold"), legend.key.size=unit(0.5, "cm"), plot.title=element_text(size=15, face="bold", hjust = 0.5)) + scale_y_continuous(breaks=seq(0, 3, 0.25))


pdf(file = paste(PATH_out,"/Figure_9_new.pdf", sep=""), width = 12, height=12)
grid.newpage()
pushViewport(viewport(layout = grid.layout(3, 2, heights = unit(c(0.5, 5,5), "null"))))

grid.text("Relative expression SNI vs Sham", vp = viewport(layout.pos.row = 1, layout.pos.col = 1:2), gp=gpar(fontsize=18, col="black", fontface="bold"))
print(p2, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))         
print(p1, vp = viewport(layout.pos.row = 2, layout.pos.col = 2))
print(p4, vp = viewport(layout.pos.row = 3, layout.pos.col = 1))
print(p3, vp = viewport(layout.pos.row = 3, layout.pos.col = 2))

dev.off()

qpcr$p.value_balbc <- p.value_balbc
qpcr$p.value_b10d2 <- p.value_b10d2


#Data for figure 9
write.csv(file = paste0(PATH_out, "/Figure_9_qpcr.csv"), qpcr)
write.csv(file = paste0(PATH_out, "/Figure_9_rnaseq.csv"), rna_seq)

###############
# HAGLR mouse #
###############
load(file=paste0(PATH_mouse, "/DE/norm_counts.RData"))

rna_seq <- norm_counts["ENSMUSG00000075277",]

rna_seq_balbc <- rna_seq[grep("BALB.c", names(rna_seq))]
rna_seq_balbc <- rna_seq_balbc / mean(rna_seq_balbc[c(1:6)])

rna_seq_b10d2 <- rna_seq[grep("B10.D2", names(rna_seq))]
rna_seq_b10d2 <- rna_seq_b10d2 / mean(rna_seq_b10d2[c(1:5)]) 


rna_seq_balbc <- data.frame(ID = rep("Haglr",2) , condition =  c("SHAM", "SNI"), expression = c(mean(rna_seq_balbc[c(1:6)]),  mean(rna_seq_balbc[c(7:10)])), sem = c(std.error(rna_seq_balbc[c(1:6)]),  std.error(rna_seq_balbc[c(7:10)]))) 

rna_seq_balbc <- transform(rna_seq_balbc, lower=expression-sem, upper=expression+sem)

rna_seq_b10d2 <- data.frame(ID = rep("Haglr",2) , condition =  c("SHAM", "SNI"), expression = c(mean(rna_seq_b10d2[c(1:5)]),  mean(rna_seq_b10d2[c(6:10)])), sem = c(std.error(rna_seq_b10d2[c(1:5)]),  std.error(rna_seq_b10d2[c(6:10)]))) 

rna_seq_b10d2 <- transform(rna_seq_b10d2, lower=expression-sem, upper=expression+sem)


p1 <- ggplot(rna_seq_balbc, aes(x = ID, y = expression, fill = condition)) +  geom_bar(stat="identity", position=position_dodge(), width=1) + geom_errorbar(aes(ymax=upper, ymin=lower), position=position_dodge(0.9), data=rna_seq_balbc, width = 0.5, size = 0.15) + xlab("Group average") + ylab("Relative expression (RNA-seq) + SEM") + ggtitle("Relative expression SNI vs sham (BALB/c mouse)") + theme(axis.text.y = element_text(size= 12, face="bold"), axis.title.y = element_text(size=15, face="bold"), axis.line.y = element_line(color="black", size = 0.3), axis.title.x = element_blank(), axis.text.x = element_text(size= 12, angle=45, hjust=1, face="bold"), legend.title=NULL, legend.text=element_text(size=15, face="bold"), legend.key.size=unit(0.5, "cm"), plot.title=element_text(size=15, face="bold", hjust = 0.5)) + scale_y_continuous(breaks=seq(0, 3, 0.25))

p2 <- ggplot(rna_seq_b10d2, aes(x = ID, y = expression, fill = condition)) +  geom_bar(stat="identity", position=position_dodge(), width=1) + geom_errorbar(aes(ymax=upper, ymin=lower), position=position_dodge(0.9), data=rna_seq_b10d2, width = 0.5, size = 0.15) + xlab("Group average") + ylab("Relative expression (RNA-seq) + SEM") + ggtitle("Relative expression SNI vs sham (B10.D2 mouse)") + theme(axis.text.y = element_text(size= 12, face="bold"), axis.title.y = element_text(size=15, face="bold"), axis.line.y = element_line(color="black", size = 0.3), axis.title.x = element_blank(), axis.text.x = element_text(size= 12, angle=45, hjust=1, face="bold"), legend.title=NULL, legend.text=element_text(size=15, face="bold"), legend.key.size=unit(0.5, "cm"), plot.title=element_text(size=15, face="bold", hjust = 0.5)) + scale_y_continuous(breaks=seq(0, 3, 0.25))

pdf(file = paste(PATH_out,"/supplemental_mouse_haglr.pdf", sep=""), width = 12, height=12)
grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 2, heights = unit(c(0.5, 5), "null"))))

grid.text("Relative expression SNI vs Sham", vp = viewport(layout.pos.row = 1, layout.pos.col = 1:2), gp=gpar(fontsize=18, col="black", fontface="bold"))
print(p1, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
print(p2, vp = viewport(layout.pos.row = 2, layout.pos.col = 2))
dev.off()

