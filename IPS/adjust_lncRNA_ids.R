library(GenomicFeatures)
library("rtracklayer")
library(biovizBase)
library(GenomicAlignments)
library(GenomicRanges)
library("AnnotationDbi")
library("org.Hs.eg.db")
library(gtools)

options(stringsAsFactors=FALSE)

PATH_results <- "/home/george/Desktop/LncRNA_v2/IPS"

load(paste(PATH_results,"/myRegs_cov/rRegs_complete_coff.RData", sep=""))

antisense <- read.table(paste(PATH_results,"/antisense/antisense_lncs_sorted.bed", sep=""),sep="\t", header=FALSE)

antisense_refseq <- read.table(paste(PATH_results,"/antisense/antisense_lncs_refseq_sorted.bed", sep=""),sep="\t", header=FALSE)

intronic <- read.table(paste(PATH_results,"/antisense/intronic.bed", sep=""),sep="\t", header=FALSE)

lincs <- read.table(paste(PATH_results,"/non_antisense/lncs_intergenic_closest.bed", sep=""),sep="\t", header=FALSE)

lnc_names <- unique(c(antisense$V4, antisense_refseq$V4, intronic$V4, lincs$V4))

rRegs_complete_coff_nc <- rRegs_complete_coff_nc[names(rRegs_complete_coff_nc) %in% lnc_names]

export(asGFF(rRegs_complete_coff_nc), paste(PATH_results,"/myRegs_cov/rRegs_complete_coff_nc.gtf", sep=""), format="gtf")

lncRNA_gtf <- import.gff(paste(PATH_results,"/myRegs_cov/rRegs_complete_coff_nc.gtf", sep=""), format="gtf")

name_mapping_coff_nc <- name_mapping_coff[name_mapping_coff$Name %in% names(rRegs_complete_coff_nc),]

name_mapping <- name_mapping_coff_nc[mixedorder(name_mapping_coff_nc$ID),]

lncRNAs <- read.table(paste(PATH_results,"/myRegs_cov/rRegs_complete_coff_nc.gtf", sep=""), sep="\t")

old_ids <- lncRNA_gtf$ID[1:length(rRegs_complete_coff_nc)]

new_ids <- name_mapping$ID

new_ids <- gsub("mRNA", "LncRNA", new_ids)

dict <- data.frame(old_ids = old_ids, new_ids=new_ids, row.names=old_ids)

to_translate <- lncRNAs$V9

temp <- strsplit(to_translate, split=";", fixed = TRUE)

old_names <- as.data.frame(matrix(unlist(temp), ncol=2, byrow=TRUE))

to_translate <- c(old_names[1:length(rRegs_complete_coff_nc),1],  tail(old_names[,2], - length(rRegs_complete_coff_nc)))

to_translate <- gsub("ID ", "", to_translate, fixed=TRUE)

to_translate <- gsub(" Parent ", "", to_translate, fixed=TRUE)

to_translate[1:length(rRegs_complete_coff_nc)] <- new_ids

pb <- txtProgressBar(min = 0, max = length(to_translate) - length(rRegs_complete_coff_nc), style = 3)
for (i in length(rRegs_complete_coff_nc):length(to_translate)) {

to_translate[i] <- dict[to_translate[[i]],]$new_ids
setTxtProgressBar(pb, i)
}

new_names_lncRNAs <- paste("ID ", to_translate[1:length(rRegs_complete_coff_nc)],";", old_names[1:length(rRegs_complete_coff_nc),]$V2,";", sep="")

new_names_exons <- paste(tail(old_names[,1], - length(rRegs_complete_coff_nc)),"; Parent ", tail(to_translate, - length(rRegs_complete_coff_nc)),";", sep="")

new_names <- c(new_names_lncRNAs, new_names_exons)

lncRNAs$V9 <- new_names

lncRNAs_bed <-  read.table(paste(PATH_results,"/HS_genome/rRegs_complete_coff_nc.bed", sep=""), sep= "\t")

lncRNAs_bed <- lncRNAs_bed[as.character(lncRNAs_bed$V4) %in% lnc_names,]

write.table(lncRNAs_bed, file = paste(PATH_results,"/HS_genome/Novel_lncRNAs_IPS.bed", sep=""), sep="\t", quote=FALSE, row.name=FALSE, col.names=FALSE)

write.table(lncRNAs, file = paste(PATH_results,"/myRegs_cov/Novel_lncRNAs_IPS.gtf", sep=""), sep="\t", quote=FALSE, row.name=FALSE, col.names=FALSE)

rownames(toc_coff_nc) <- gsub("mRNA", "LncRNA", rownames(toc_coff_nc))

LncRNAs_toc <- toc_coff_nc[rownames(toc_coff_nc) %in% new_ids,]

write.csv(LncRNAs_toc, file = paste(PATH_results,"/LncRNAs_toc_IPS.csv", sep=""))

save(LncRNAs_toc, file = paste(PATH_results,"/LncRNAs_toc_IPS.RData", sep=""))


