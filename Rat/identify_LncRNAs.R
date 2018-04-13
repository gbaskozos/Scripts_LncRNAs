#Load packages
library("data.table")
library(Rsamtools)
library(GenomicFeatures)
library("rtracklayer")
library(biovizBase)
library(GenomicAlignments)
library(BiocParallel)
library(GenomicRanges)
library(biomaRt)
library(parallel)
library("org.Rn.eg.db")
library("doParallel")
library(iterators)
source("/home/george/Desktop/LncRNA_v2/Rat/scripts/dropDetect.R")
source("/home/george/Desktop/LncRNA_v2/Rat/scripts/BAM_to_IOE.R")
####################
# Convert Z-scores #
####################
convert.z.score<-function(z, one.sided=NULL) {
    if(is.null(one.sided)) {
        pval = pnorm(-abs(z));
        pval = 2 * pval
    } else if(one.sided=="-") {
        pval = pnorm(z);
    } else {
        pval = pnorm(-z);                                                                                 
    }
    return(pval);
}   

################################

#############################
# Yield size and nr of CPUs #
#############################



n_cores <- detectCores() - 1

multicoreParam <- MulticoreParam(workers = n_cores, progressbar = TRUE)
options(mc.cores=n_cores)
register(multicoreParam)


filter_SJ_length <- 1
filter_regs_no_SJ <- 1
filter_intronic <- 1

lag <- 31
threshold <- 5
length <- 20
influence <- 0

lanes <- 1

#PATH of input files, bam files etc
PATH <- "/media/george/data3tb/Sequencing_Data/Rat_DRG_snt/star_out"

#PATH to export results
PATH_results <- "/home/george/Desktop/LncRNA_v2/Rat"

########################################
##Create file lists for each condition##
########################################
src_files <- list.files(path=PATH, pattern="*.sortedByCoord.out.bam$", full.names=TRUE)
bai_files <- list.files(path=PATH, pattern="*.sortedByCoord.out.bam.bai$", full.names=TRUE)

file_names <- list.files(path=PATH, pattern="*.sortedByCoord.out.bam$", full.names=FALSE)

sham_names <- c(grep("_229", file_names, value=TRUE), grep("_235", file_names, value=TRUE), grep("_226", file_names, value=TRUE), grep("_232", file_names, value=TRUE))

D21_SNT_names <- c(grep("_231", file_names, value=TRUE), grep("_237", file_names, value=TRUE), grep("_228", file_names, value=TRUE), grep("_234", file_names, value=TRUE)) 

sham <- c(grep("_229", src_files, value=TRUE), grep("_235", src_files, value=TRUE), grep("_226", src_files, value=TRUE), grep("_232", src_files, value=TRUE))

D21_SNT <- c(grep("_231", src_files, value=TRUE), grep("_237", src_files, value=TRUE), grep("_228", src_files, value=TRUE), grep("_234", src_files, value=TRUE))

design_all <- c(sham, D21_SNT)

#bamlist <- BamFileList(c(sham, sni), index = c(sham, sni), yieldSize = 9000000)

bamlist <- BamFileList(design_all, index = design_all, yieldSize = 1000000)

#seqlevels(bamlist, pruning.mode="coarse") <- seqlevels(bamlist)[grepl("^[0-9, X, Y]*$", seqlevels(bamlist))]

bam_names <- c(sham_names, D21_SNT_names)

#bam_names <- paste(substr(bam_names, 0, nchar(bam_names)-29), c(rep("Sham",4), rep("SNT",4)), sep="_")

bam_names <- substr(bam_names, 0, nchar(bam_names)-29)

################################################
##     DOWNLOAD KNOWN GENE SET ANNOTATIONS    ## 
##Create txdb objects, group features by genes## 
################################################
if (!file.exists(file = paste(PATH_results, "/rn6_genome/iGRanges.RData", sep=""))) {

if (!file.exists(file = paste(PATH_results,"/rn6_genome/ENSEMBL_rn6_exonsByGene.RData", sep=""))) {

txdb <- makeTxDbFromBiomart(biomart="ENSEMBL_MART_ENSEMBL", dataset="rnorvegicus_gene_ensembl", transcript_ids=NULL, circ_seqs=DEFAULT_CIRC_SEQS, filter=NULL, id_prefix="ensembl_", host="www.ensembl.org", port=80, taxonomyId=NA, miRBaseBuild=NA)

Rat <- useMart("ensembl",dataset="rnorvegicus_gene_ensembl")

ens_annot <- getBM(attributes=c('ensembl_transcript_id', 'description', 'transcript_biotype'), filters = 'ensembl_transcript_id', values = as.character(transcripts(txdb)$tx_name), mart = Rat)

nc_ensembl <- ens_annot[as.character(ens_annot$transcript_biotype) == "lincRNA" | as.character(ens_annot$transcript_biotype) == "antisense" | as.character(ens_annot$transcript_biotype) == "antisense_RNA" | as.character(ens_annot$transcript_biotype) == "sense_intronic", ]

write.csv(file = paste(PATH_results,"/rn6_genome/ensembl_annot.csv", sep=""), ens_annot)

write.csv(file = paste(PATH_results,"/rn6_genome/ensembl_nc.csv", sep=""), nc_ensembl)

tx_by_gene <- reduce(exonsBy(txdb,"gene"))

seqlevels(tx_by_gene, pruning.mode="coarse") <- seqlevels(tx_by_gene)[grepl("^[0-9, X, Y]*$", seqlevels(tx_by_gene))]

export.bed(tx_by_gene, paste(PATH_results,"/rn6_genome/ENSEMBL_rn6_genes.bed", sep=""))
export(asGFF(tx_by_gene), paste(PATH_results,"/rn6_genome/ENSEMBL_rn6_genes.gtf", sep=""),  format="gtf")

#########################
# Fetch introns by gene #
#########################
introns <- intronsByTranscript(txdb, use.names=TRUE)
ulst <- unlist(introns)
intronsNoDups <- ulst[!duplicated(ulst)]
txnames <- names(intronsNoDups)
map <- select(txdb, keys = txnames, keytype='TXNAME', columns='GENEID')
idx <- map$GENEID[!is.na(map$GENEID)]
intronsByGene <- split(intronsNoDups[!is.na(map$GENEID)], idx)
names(intronsByGene) <- unique(idx)

intronsByGene <- reduce(intronsByGene)

seqlevels(intronsByGene, pruning.mode="coarse") <- seqlevels(intronsByGene)[grepl("^[0-9, X, Y]*$", seqlevels(intronsByGene))]

save(file = paste(PATH_results,"/rn6_genome/ENSEMBL_rn6_intronsByGene.RData", sep=""), intronsByGene)
export(asGFF(intronsByGene), paste(PATH_results,"/rn6_genome/intronsByGene.gtf", sep="") , format="gtf")
export.bed(intronsByGene, paste(PATH_results,"/rn6_genome/intronsByGene.bed", sep=""))
####
tx_by_gene <- reduce(unlist(exonsBy(txdb,"gene")))

save(txdb, tx_by_gene, file = paste(PATH_results,"/rn6_genome/ENSEMBL_rn6_exonsByGene.RData", sep="") )

########################
# Export known lncRNAs #
########################
txdb_lncs <- makeTxDbFromBiomart(biomart="ENSEMBL_MART_ENSEMBL", dataset="rnorvegicus_gene_ensembl", transcript_ids=as.character(nc_ensembl$ensembl_transcript_id), circ_seqs=DEFAULT_CIRC_SEQS, filter=NULL, id_prefix="ensembl_", host="www.ensembl.org", port=80, taxonomyId=NA, miRBaseBuild=NA)

tx_by_lnc <- reduce(exonsBy(txdb_lncs,"gene"))

seqlevels(tx_by_lnc, pruning.mode="coarse") <- seqlevels(tx_by_lnc)[grepl("^[0-9, X, Y]*$", seqlevels(tx_by_lnc))]

export(asGFF(tx_by_lnc), paste(PATH_results,"/rn6_genome/ENSEMBL_rn6_LncRNAs.gtf", sep=""),  format="gtf")
export.bed(tx_by_lnc, paste(PATH_results,"/rn6_genome/ENSEMBL_rn6_LncRNAs.bed", sep=""))

ens_lncs_names <- names(tx_by_lnc)
write.csv(file = paste(PATH_results,"/rn6_genome/ens_lncs_names.csv", sep=""), ens_lncs_names)

ens_lncs_annot <- getBM(attributes=c('ensembl_gene_id', 'rgd_symbol', 'description', 'gene_biotype'), filters = 'ensembl_gene_id', values = as.character(ens_lncs_names), mart = Rat)

antisense_ens <- ens_lncs_annot[as.character(ens_lncs_annot$gene_biotype) == "antisense" | as.character(ens_lncs_annot$gene_biotype) == "antisense_RNA",]

intergenic_ens <- ens_lncs_annot[as.character(ens_lncs_annot$gene_biotype) == "lincRNA",]

antisense_ens_lncs <- tx_by_lnc[as.character(antisense_ens$ensembl_gene_id)]

intergenic_ens_lncs <- tx_by_lnc[as.character(intergenic_ens$ensembl_gene_id)]

export.bed(antisense_ens_lncs, paste(PATH_results,"/rn6_genome/ENSEMBL_rn6_antisense.bed", sep=""))
export(asGFF(antisense_ens_lncs), paste(PATH_results,"/rn6_genome/ENSEMBL_rn6_antisense.gtf", sep=""),  format="gtf")

export.bed(intergenic_ens_lncs, paste(PATH_results,"/rn6_genome/ENSEMBL_rn6_intergenic.bed", sep=""))
export(asGFF(intergenic_ens_lncs), paste(PATH_results,"/rn6_genome/ENSEMBL_rn6_intergenic.gtf", sep=""),  format="gtf")

write.csv(file = paste(PATH_results,"/rn6_genome/ens_lncs_annotation.csv", sep=""), ens_lncs_annot)

save(tx_by_lnc, nc_ensembl, txdb_lncs, ens_lncs_annot, ens_lncs_annot, ens_lncs_names, file = paste(PATH_results,"/rn6_genome/ENSEMBL_rn6_lncRNAs.RData", sep="")) 
}

#######################################################################

load(paste(PATH_results,"/rn6_genome/ENSEMBL_rn6_exonsByGene.RData", sep=""))
load(paste(PATH_results,"/rn6_genome/ENSEMBL_rn6_intronsByGene.RData", sep=""))

if (!file.exists(file = paste(PATH_results,"/rn6_genome/Refseq_rn6_exonsByGene.RData", sep=""))) {

txdb_RG <- makeTxDbFromUCSC(genome='rn6',tablename='refGene')

tx_RG_by_gene <- reduce(unlist(exonsBy(txdb_RG, "gene")))

seqlevels(tx_RG_by_gene) <- sub('chr','',seqlevels(tx_RG_by_gene))

seqlevels(tx_RG_by_gene, pruning.mode="coarse") <- seqlevels(tx_RG_by_gene)[grepl("^[0-9, X, Y]*$", seqlevels(tx_RG_by_gene))]

export.bed(tx_RG_by_gene,paste(PATH_results,"/rn6_genome/rn6_refseq.bed", sep=""))

save(tx_RG_by_gene, file = paste(PATH_results,"/rn6_genome/Refseq_rn6_exonsByGene.RData", sep=""))



####################################

}

load(paste(PATH_results, "/rn6_genome/Refseq_rn6_exonsByGene.RData", sep=""))

#txdb_xRG <- makeTxDbFromUCSC(genome='rn6',tablename='xenoRefGene')

#xenoRefGenes <- reduce(unlist(exonsBy(txdb_xRG, "gene")))

xenoRefGenes <- import.gff(paste(PATH_results,"/rn6_genome/xenoRefGene.gtf", sep=""))

seqlevels(xenoRefGenes) <- sub('chr','',seqlevels(xenoRefGenes))

seqlevels(xenoRefGenes, pruning.mode="coarse") <- seqlevels(xenoRefGenes)[grepl("^[0-9, X, Y]*$", seqlevels(xenoRefGenes))]

save(xenoRefGenes, file = paste(PATH_results, "/rn6_genome/xenoRefGenes.RData", sep=""))

 mcols(xenoRefGenes) <- NULL

gRanges <- reduce(c(tx_by_gene, tx_RG_by_gene, xenoRefGenes))

seqlevels(gRanges, pruning.mode="coarse") <- seqlevels(gRanges)[grepl("^[0-9, X, Y]*$", seqlevels(gRanges))]

save(gRanges, file = paste(PATH_results, "/rn6_genome/rn6_gRanges.RData", sep=""))

################################################################
##Extend genes by 1000 bp to avoid UTRs or non annotated exons##
################################################################
load(paste(PATH_results, "/rn6_genome/rn6_gRanges.RData", sep=""))
gRangesExt <- gRanges
start(ranges(gRangesExt[! start(ranges(gRangesExt)) < 1000])) = start(ranges(gRangesExt[! start(ranges(gRangesExt)) < 1000])) - 1000
maxlengthV <- seqlengths(gRangesExt)[as.character(seqnames(gRangesExt))] # total for that chromosome
end(ranges(gRangesExt))[! maxlengthV - end(ranges(gRangesExt)) <= 1000] <- end(ranges(gRangesExt))[! maxlengthV - end(ranges(gRangesExt)) <= 1000]  + 1000
save(gRangesExt, file = paste(PATH_results, "/rn6_genome/rn6_gRangesExt.RData", sep=""))

#export.bed(gRangesExt, paste(PATH_results, "/rn6_genome/rn6_gRangesExt.bed", sep=""))


#####################################################
##Find gaps, i.e. non annoatted areas of the genome##
#####################################################
igRangesExt <- gaps(gRangesExt)
igRangesExt <- igRangesExt[!strand(igRangesExt)=="*"]


save(igRangesExt, file = paste(PATH_results, "/rn6_genome/iGRanges.RData", sep=""))

#export.bed(igRangesExt, paste(PATH_results, "/rn6_genome/iGRanges.bed", sep=""))
}

load(paste(PATH_results, "/rn6_genome/iGRanges.RData", sep=""))

#####################################################################################

######################################################################################
## Read bam files into Granges lists and discard reads overlapping with known genes ##
######################################################################################
param <- ScanBamParam(flag=scanBamFlag(isProperPair=TRUE, isDuplicate=NA, isSecondaryAlignment=NA))

bplapply(bamlist, FUN = BAM_to_IOE, BPPARAM = MulticoreParam(), PATH=PATH, PATH_results=PATH_results, igRangesExt=igRangesExt, param = param, len=100, dep=2, suffix=29)

##########################################################
#             Read SJ identified from STAR               #
# SJs are represented as genomic coordinates of introns  #
##########################################################
# Filter out known SJ
SJ_path <- paste(PATH, "/SJ", sep="")
system(paste("for file in ", SJ_path,"/*SJ.out.tab; do awk -F \"\t\" '$6==\"0\" {print }' $file > $file.novel.tab; done", sep=""))

SJ_files <- list.files(path= paste(PATH, "/SJ", sep=""), pattern="*.novel.tab$", full.names=TRUE)

SJ_list <- vector('list', lanes * length(bamlist))
SJ_list <- mclapply(SJ_files, readSTARJunctions, mc.cores=n_cores, mc.cleanup=TRUE, mc.allow.recursive=TRUE)

#Filter SJ with < 2 reads and length < 20 or > 100000 
if (filter_SJ_length) {
SJ_list_filtered <- mclapply(SJ_list, function(x) x[((mcols(x)$um_reads + mcols(x)$mm_reads)  >= 2) & (width(x) > 20 & width(x) < 100000)], mc.cores=n_cores, mc.cleanup=TRUE, mc.allow.recursive=TRUE)
} else {
SJ_list_filtered <- mclapply(SJ_list, function(x) x[((mcols(x)$um_reads + mcols(x)$mm_reads)  >= 2)], mc.cores=n_cores, mc.cleanup=TRUE, mc.allow.recursive=TRUE)
}


SJ_regs <- mclapply(SJ_list_filtered, reduce, mc.cores=n_cores, mc.cleanup=TRUE, mc.allow.recursive=TRUE)

rm(SJ_list_filtered, SJ_list)

SJ_regs <- Reduce(c,SJ_regs)

seqlevels(SJ_regs, pruning.mode="coarse") <- seqlevels(SJ_regs)[grepl("^[0-9, X, Y]*$", seqlevels(SJ_regs))]

SJ_regs <- unique(SJ_regs)
SJ_regs <- sortSeqlevels(SJ_regs)
SJ_regs <- sort(SJ_regs, ignore.strand=TRUE)
names(SJ_regs) <- paste(seqnames(SJ_regs), ":",  start(ranges(SJ_regs)), "-", end(ranges(SJ_regs)),"(", as.character(strand(SJ_regs)), ")", sep="")

save(file = paste(PATH_results,"/myRegs_cov/SJ_regs_novel.RData", sep=""), SJ_regs)


###################################  	
## Load the customised annotation##
################################### 
myRegs <- vector('list', length(bam_names))
coverage_F <- NULL
coverage_B <- NULL

i <- 1

for (name in bam_names) {

load(file = paste(PATH_results,"/myRegs_cov/", name,".RData", sep=""))
cat(name, "\n")
myRegs[i] <- bam

if (i == 1) {
coverage_F <- cov_F 
coverage_B <- cov_B
}
else {
coverage_F <- coverage_F + cov_F
coverage_B <- coverage_B + cov_B
}

rm(bam, cov_F, cov_B)

i <- i+1
}

save(file = paste(PATH_results,"/myRegs_cov/coverage_F.RData", sep=""), coverage_F)
save(file = paste(PATH_results,"/myRegs_cov/coverage_B.RData", sep=""), coverage_B)

####################################################		
##   Reconstruct genes outside annotated regions  ##
####################################################


#Create TOC
#Reduce/collapse book-end and overlapping features
#Keep IOE identified in all samples of the experiment

#rRegs_hc <- reduce(c(Reduce(intersect, myRegs[grep("BALB.c", bam_names)]), Reduce(intersect, myRegs[grep("B10.D2", bam_names)])))
rRegs_hc <- reduce(c(Reduce(intersect, myRegs)))
rRegs <- reduce(Reduce(c,myRegs)) 
rm(myRegs)
save(file = paste(PATH_results,"/myRegs_cov/rRegs.RData", sep=""), rRegs, rRegs_hc)

#Identify continuously expressed regions fully contained in introns
load(paste(PATH_results,"/rn6_genome/ENSEMBL_rn6_intronsByGene.RData", sep=""))
rRegs_intronic <- subsetByOverlaps(rRegs, intronsByGene, type="within", ignore.strand=FALSE)
rm(intronsByGene)

#####################
# Ouside gene models #
#####################
#Check overlaps between collapsed features and SJ
SJ_olap <- subsetByOverlaps(SJ_regs, rRegs, maxgap=0L, minoverlap=0L, type="any", ignore.strand=FALSE)
#export(SJ_olap, paste(PATH_results,"/myRegs_cov/SJ_olap.gtf", sep=""), format="gtf")
save(file = paste(PATH_results,"/myRegs_cov/SJ_olap.RData", sep=""), SJ_olap)

#Check overlaps between SJ and collapsed features
rRegs_sj <- subsetByOverlaps(rRegs, SJ_olap, maxgap=0L, minoverlap=0L, type="any", ignore.strand=FALSE)
names(rRegs_sj) <- paste(seqnames(rRegs_sj), ":",  start(ranges(rRegs_sj)), "-", end(ranges(rRegs_sj)),"(", as.character(strand(rRegs_sj)), ")", sep="")
#export(rRegs_sj, paste(PATH_results,"/myRegs_cov/rRegs_sj.gtf", sep=""), format="gtf")
save(file = paste(PATH_results,"/myRegs_cov/rRegs_sj.RData", sep=""), rRegs_sj)

#Identify continuously expressed regions with no SJs, keep only the sub-regions identified in all bam files
rRegs_nosj <- setdiff(rRegs, rRegs_sj, ignore.strand=FALSE)
#rRegs_nosj <- setdiff(rRegs_nosj, rRegs_intronic, ignore.strand=FALSE)
rRegs_nosj <- intersect(rRegs_nosj, rRegs_hc)
#export(rRegs_nosj, paste(PATH_results,"/myRegs_cov/rRegs_nosj.gtf", sep=""), format="gtf")
save(file = paste(PATH_results,"/myRegs_cov/rRegs_nosj.RData", sep=""), rRegs_nosj)

#rm(rRegs_sj, SJ_olap, rRegs_intronic, rRegs_outside_genes, rRegs_nosj)
##############
# Single IOE #
##############
#load(file = paste(PATH_results,"/myRegs_cov/rRegs_nosj.RData", sep=""))
#load(file = paste(PATH_results,"/myRegs_cov/coverage_F.RData", sep=""))

#Filter out mono-exonic with length <= 200bp
rRegs_nosj <- rRegs_nosj[width(rRegs_nosj) > 200]
####################################
#Trim edges and find coverage drops#
####################################

#Forward strand
rRegs_nosj_F <- rRegs_nosj[strand(rRegs_nosj) == "+"] 
coverage_F_nosj <- coverage_F[rRegs_nosj_F]

norm_coverage_F_nosj <- scale(mean(coverage_F_nosj)/unlist(mclapply(coverage_F_nosj, function(x) length(x), mc.cores=n_cores, mc.cleanup=TRUE, mc.allow.recursive=TRUE)), center=TRUE, scale=TRUE) 

rRegs_nosj_F <- rRegs_nosj_F[p.adjust(convert.z.score(norm_coverage_F_nosj, one.sided=TRUE), method = "BH") < 0.1]
coverage_F_nosj <- coverage_F_nosj[p.adjust(convert.z.score(norm_coverage_F_nosj, one.sided=TRUE), method = "BH") < 0.1] 

start_F = start(rRegs_nosj_F)
seqnames_F = names(coverage_F_nosj)

drops <- mcmapply(FUN=dropDetect, coverage_F_nosj, start_F, seqnames_F, MoreArgs = list(strand = "+", lag = lag, threshold = threshold, length = length, influence = influence, intron_identification = 1), mc.cores=n_cores)

drops <- GRangesList(drops)
rRegs_nosj_F_trimmed <- reduce(psetdiff(rRegs_nosj_F, drops, ignore.strand=FALSE))
rRegs_nosj_F_trimmed <- rRegs_nosj_F_trimmed[sum(width(rRegs_nosj_F_trimmed)) > 200]

#Reverse strand
rRegs_nosj_B <- rRegs_nosj[strand(rRegs_nosj) == "-"] 
coverage_B_nosj <- coverage_B[rRegs_nosj_B]

norm_coverage_B_nosj <- scale(mean(coverage_B_nosj)/unlist(mclapply(coverage_B_nosj, function(x) length(x), mc.cores=n_cores, mc.cleanup=TRUE, mc.allow.recursive=TRUE)), center=TRUE, scale=TRUE) 

rRegs_nosj_B <- rRegs_nosj_B[p.adjust(convert.z.score(norm_coverage_B_nosj, one.sided=TRUE), method = "BH") < 0.1]
coverage_B_nosj <- coverage_B_nosj[p.adjust(convert.z.score(norm_coverage_B_nosj, one.sided=TRUE), method = "BH") < 0.1] 

start_B = start(rRegs_nosj_B)
seqnames_B = names(coverage_B_nosj)

drops <- mcmapply(FUN=dropDetect, coverage_B_nosj, start_B, seqnames_B, MoreArgs = list(strand = "-", lag = lag, threshold = threshold, length = length, influence = influence, intron_identification = 1), mc.cores=n_cores)

drops <- GRangesList(drops)
rRegs_nosj_B_trimmed <- reduce(psetdiff(rRegs_nosj_B, drops, ignore.strand=FALSE))
rRegs_nosj_B_trimmed <- rRegs_nosj_B_trimmed[sum(width(rRegs_nosj_B_trimmed)) > 200]

#Merge
rRegs_nosj_trimmed <- reduce(c(rRegs_nosj_F_trimmed, rRegs_nosj_B_trimmed))

rRegs_nosj_trimmed <- mclapply(rRegs_nosj_trimmed, FUN= function(x) x[width(x)>1], mc.cores=n_cores)

rRegs_nosj_trimmed <- GRangesList(rRegs_nosj_trimmed)

names(rRegs_nosj_trimmed) <- paste(unique(seqnames(rRegs_nosj_trimmed)), ":", min(start(rRegs_nosj_trimmed)), "-", max(end(rRegs_nosj_trimmed)), "(", as.character(unique(strand(rRegs_nosj_trimmed))), ")", sep="")

export(asGFF(rRegs_nosj_trimmed), paste(PATH_results,"/myRegs_cov/rRegs_nosj_trimmed.gtf", sep=""), format="gtf")
save(file = paste(PATH_results,"/myRegs_cov/rRegs_nosj_trimmed.RData", sep=""), rRegs_nosj_trimmed)

rm(rRegs_nosj_F_trimmed,rRegs_nosj_B_trimmed,drops, rRegs_nosj_F, rRegs_nosj_B, coverage_F_nosj, coverage_B_nosj, start_F, start_B, seqnames_F, seqnames_B, rRegs_nosj)

####################
# Overlapped by SJ #
####################

SJ_red_regs <- reduce(SJ_olap)

#Calculate disjoined regions
SJ_disjoin <- disjoin(SJ_olap)
names(SJ_disjoin) <- paste(seqnames(SJ_disjoin), ":",  start(ranges(SJ_disjoin)), "-", end(ranges(SJ_disjoin)),"(", as.character(strand(SJ_disjoin)), ")", sep="")

#Keep track of overlapping features after collapsing them
SJ_split_guide <- reduce(SJ_olap, with.revmap=TRUE, ignore.strand=FALSE)
revmap <- mcols(SJ_split_guide)$revmap  # an IntegerList
SJ_regs_grouped <- relist(SJ_olap[unlist(revmap)], revmap)
#export(asGFF(SJ_regs_grouped), paste(PATH_results,"/myRegs_cov/SJ_regs_grouped.gtf", sep=""), format="gtf")
##########################################################################
#Find overlaps between collapsed island of exression and introns/SJs
olap_mapping <- findOverlaps(query = rRegs_sj, subject = SJ_red_regs, type="any", ignore.strand=FALSE, maxgap=-1L, minoverlap=0L)

#Create graph table of collapsed regions connected by introns/SJ
hitsmap <- as.data.frame(olap_mapping)

hitsmap$names <- names(rRegs_sj)[as.numeric(hitsmap$queryHits)] 

AggHits <- aggregate(hitsmap, list(hitsmap$subjectHits), function(x) paste0(unique(x)))

#test <- grep("\\b6861\\b", AggHits$queryHits)

skeleton <- AggHits$names

#Group collapsed regions according to graph table
rRegs_sj_grouped <- relist(rRegs_sj[unlist(AggHits$names)], AggHits$names)

names(SJ_disjoin) <- paste(seqnames(SJ_disjoin), ":",  start(ranges(SJ_disjoin)), "-", end(ranges(SJ_disjoin)),"(", as.character(strand(SJ_disjoin)), ")", sep="")

names(rRegs_sj_grouped) <- paste(unique(seqnames(rRegs_sj_grouped)), ":", min(start(rRegs_sj_grouped)), "-", max(end(rRegs_sj_grouped)), "(", as.character(unique(strand(rRegs_sj_grouped))), ")", sep="")  

#Filter out groups of collapsed regions with width < 200bp
rRegs_sj_grouped <- rRegs_sj_grouped[sum(width(rRegs_sj_grouped)) > 200] #regs width 200bp

##################################################################
# Remove single IoE that are parts of another non-annotated gene #
##################################################################

single_islands <- rRegs_sj_grouped[as.character(names(rRegs_sj_grouped)) %in% as.character(hitsmap$names)]

rRegs_sj_grouped_clean <- rRegs_sj_grouped[!as.character(names(rRegs_sj_grouped)) %in% as.character(hitsmap$names)]

to_remove <- subsetByOverlaps(single_islands, rRegs_sj_grouped_clean, maxgap=0L, minoverlap=100L, type="any", ignore.strand=FALSE)

rRegs_sj_grouped <- rRegs_sj_grouped[!as.character(names(rRegs_sj_grouped)) %in% as.character(names(to_remove))]

export(asGFF(rRegs_sj_grouped), paste(PATH_results,"/myRegs_cov/rRegs_sj_grouped.gtf", sep=""), format="gtf")

name_mapping <- subset(mcols(asGFF(rRegs_sj_grouped)), type=="mRNA", select=c(ID, Name))

rm(to_remove, rRegs_sj_grouped_clean, single_islands)
##########################################################################################
# Claculate consensus intron at the gene level by voting of individual disjoined introns #
#                  Count overlaps of SJs with their disjoined regions                    #
##########################################################################################

SJ_red_olap <- countOverlaps(query = SJ_disjoin, subject = SJ_olap, type="any", ignore.strand=FALSE)

SJ_olap_mapping <- findOverlaps(query = SJ_disjoin, subject = rRegs_sj_grouped, type="any", ignore.strand=FALSE, maxgap=-1L, minoverlap=0L)


#Create graph table of disjoined introns and collapsed regions
SJ_hitsmap <- as.data.frame(SJ_olap_mapping)

SJ_hitsmap$names <- names(SJ_disjoin)[as.numeric(SJ_hitsmap$queryHits)] 

SJ_AggHits <- aggregate(SJ_hitsmap, list(SJ_hitsmap$subjectHits), function(x) paste0(unique(x)))

skeleton <- SJ_AggHits$names

mcols(SJ_disjoin)$score <- SJ_red_olap

#Group according to the graph table
SJ_disjoint_grouped <- relist(SJ_disjoin[unlist(SJ_AggHits$names)], SJ_AggHits$names)

subset_SJ <- function(x){ 

coff <- round(max(as.matrix(mcols(x))) - 2*sd(as.matrix(mcols(x)))) 

if (is.na(coff)) {
coff <- 1 }
x[as.matrix(mcols(x)) >= coff,]

}

SJ_recon <- mclapply(SJ_disjoint_grouped, subset_SJ, mc.cores=n_cores, mc.cleanup=TRUE, mc.allow.recursive=TRUE) 

#SJ_recon <- lapply(SJ_disjoint_grouped, FUN=subset_SJ) 

SJ_recon <- GRangesList(SJ_recon)

#Check SJ_recon falling into known introns
load(paste(PATH_results,"/rn6_genome/ENSEMBL_rn6_intronsByGene.RData", sep=""))
SJ_recon_intronic <- subsetByOverlaps(SJ_recon, intronsByGene, type="within", ignore.strand=FALSE)
#export(asGFF(SJ_recon_intronic), paste(PATH_results,"/myRegs_cov/SJ_recon_intronic.gtf", sep=""), format="gtf")

save(SJ_recon, SJ_disjoint_grouped, rRegs_sj_grouped, SJ_recon_intronic, file = paste(PATH_results,"/myRegs_cov/SJ_recon.RData", sep=""))

#load(paste(PATH_results,"/myRegs_cov/SJ_recon.RData", sep=""))

#rm(SJ_red_olap, SJ_olap_mapping, SJ_hitsmap, SJ_AggHits, skeleton)

#export(asGFF(SJ_recon), paste(PATH_results,"/myRegs_cov/SJ_recon.gtf", sep=""), format="gtf")
###############################################################
# Subtract GRanges introns from exons for the same gene level #
###############################################################
covered_sj <- subsetByOverlaps(rRegs_sj_grouped, SJ_recon, type="within", ignore.strand=FALSE)
to_remove <- subsetByOverlaps(covered_sj, intronsByGene, type="within", ignore.strand=FALSE)
covered_sj <- covered_sj[!as.character(names(covered_sj)) %in% as.character(names(to_remove))]

#export(asGFF(covered_sj), paste(PATH_results,"/myRegs_cov/covered_sj.gtf", sep=""), format="gtf")

#Remove SJ_recon overlapped by intronic SJs
to_remove <- subsetByOverlaps(rRegs_sj_grouped, SJ_recon_intronic, type="any", ignore.strand=FALSE)
covered_sj <- covered_sj[!as.character(names(covered_sj)) %in% as.character(names(to_remove))]


rRegs_grouped_trimmed <- reduce(setdiff(rRegs_sj_grouped, SJ_recon, ignore.strand=FALSE))

############################################################################
# Create non-redundant annottation suitable for counting at the gene level # 
############################################################################
self_olap <- findOverlaps(rRegs_grouped_trimmed, drop.self=TRUE, drop.redundant=TRUE)

hitsmap <- as.data.frame(self_olap)

rRegs_grouped_trimmed_nr <- GRangesList()

for (i in 1:nrow(hitsmap)){

rRegs_grouped_trimmed[hitsmap$queryHits[i]] <- GRangesList(reduce(unlist(c(rRegs_grouped_trimmed[hitsmap$queryHits[i]], rRegs_grouped_trimmed[hitsmap$subjectHits[i]]))))

}

rRegs_grouped_trimmed <- rRegs_grouped_trimmed[ -as.numeric(hitsmap$subjectHits)]

rRegs_grouped_trimmed <- rRegs_grouped_trimmed[sum(width(rRegs_grouped_trimmed)) > 200]

to_remove_nosj <- subsetByOverlaps(rRegs_nosj_trimmed, rRegs_grouped_trimmed, type="any", ignore.strand=FALSE)
rRegs_nosj_trimmed <- rRegs_nosj_trimmed[!as.character(names(rRegs_nosj_trimmed)) %in% as.character(names(to_remove_sj))] 

nr_of_exons <- unlist(mclapply(rRegs_grouped_trimmed, function(x) length(x), mc.cores=n_cores, mc.cleanup=TRUE, mc.allow.recursive=TRUE))
to_remove_sj <- subsetByOverlaps(rRegs_grouped_trimmed[nr_of_exons==1], intronsByGene, type="any", ignore.strand=FALSE)
rRegs_grouped_trimmed <- rRegs_grouped_trimmed[!as.character(names(rRegs_grouped_trimmed)) %in% as.character(names(to_remove_sj))]

nr_of_exons_csj <- unlist(mclapply(covered_sj, function(x) length(x), mc.cores=n_cores, mc.cleanup=TRUE, mc.allow.recursive=TRUE))
to_remove_csj <- subsetByOverlaps(covered_sj[nr_of_exons_csj==1], intronsByGene, type="any", ignore.strand=FALSE)
covered_sj <- covered_sj[!as.character(names(covered_sj)) %in% as.character(names(to_remove_csj))]

rRegs_complete <- reduce(c(rRegs_grouped_trimmed, covered_sj, rRegs_nosj_trimmed))

nr_of_exons_complete <- unlist(mclapply(rRegs_complete, function(x) length(x), mc.cores=n_cores, mc.cleanup=TRUE, mc.allow.recursive=TRUE))

export(asGFF(rRegs_complete), paste(PATH_results,"/myRegs_cov/rRegs_complete.gtf", sep=""), format="gtf")

############################
############################
############################
# Export GFF and BED files #
############################
names(rRegs_complete) <- paste(unique(seqnames(rRegs_complete)), ":", min(start(rRegs_complete)), "-", max(end(rRegs_complete)), "(", as.character(unique(strand(rRegs_complete))), ")", sep="")  

#Remove outlying transcripts with > 42 regions
rRegs_complete <- rRegs_complete[unlist(mclapply(rRegs_complete, function(x) length(x), mc.cores=n_cores, mc.cleanup=TRUE, mc.allow.recursive=TRUE)) <=40 ] 

export(asGFF(rRegs_complete), paste(PATH_results,"/myRegs_cov/rRegs_complete.gtf", sep=""), format="gtf")

export.bed(rRegs_complete, paste(PATH_results,"/myRegs_cov/rRegs_complete.bed", sep=""))

name_mapping <- subset(mcols(asGFF(rRegs_complete)), type=="mRNA", select=c(ID, Name))

#########################################################################
save(file = paste(PATH_results,"/myRegs_cov/transcript_recon.RData", sep=""), olap_mapping, rRegs_sj, SJ_olap, rRegs, AggHits, rRegs_sj_grouped, SJ_disjoin, rRegs_grouped_trimmed, name_mapping, rRegs_complete)
##########################################################################
load(paste(PATH_results,"/myRegs_cov/transcript_recon.RData", sep=""))
##########################################################################
###############################

#Count overlaps
system("Rat_scripts/parallel_count.sh")

#Calculate coding potential
system("Rat_scripts/coding_potential.sh")

####################
# Load count files #
####################
count_files <- list.files(path=PATH, pattern="*.lncs_grouped.counts.tab$", full.names=TRUE)
sham_counts <- c(grep("_229", count_files, value=TRUE), grep("_235", count_files, value=TRUE), grep("_226", count_files, value=TRUE), grep("_232", count_files, value=TRUE))
D21_SNT_counts <- c(grep("_231", count_files, value=TRUE), grep("_237", count_files, value=TRUE), grep("_228", count_files, value=TRUE), grep("_234", count_files, value=TRUE))
count_files_ord <- c(sham_counts, D21_SNT_counts)
toc <- mclapply(count_files_ord,read.table, sep="\t", row.names=1, mc.cores=n_cores, mc.cleanup=TRUE, mc.allow.recursive=TRUE)
toc <- as.data.frame(toc)
names(toc) <- bam_names 
no_counts <- tail(toc)[-1,]
toc <- toc[1:(dim(toc)[1]-5),]

#####################
save(file = paste(PATH_results,"/regToc_no_coff.RData", sep=""), toc, name_mapping)

load(paste(PATH_results,"/regToc_no_coff.RData", sep=""))

## Aplly cut off to condition
cutoff <- 0
toc_coff <- toc[apply(toc, 1, function(x) sum(x[c(1:4)] > cutoff) == 4 | sum(x[c(5:8)] > cutoff) == 4), ]

rownames(name_mapping) <- as.character(name_mapping$ID)

name_mapping_coff <- name_mapping[as.character(rownames(toc_coff)),]

toc_coff$name <- name_mapping_coff[as.character(rownames(toc_coff)),2]

rRegs_complete_coff <- rRegs_complete[names(rRegs_complete) %in% as.character(toc_coff$name)]

#########################
### Read CPAT results ###
#########################
cpat_1 <- read.table(paste(PATH_results,"/myRegs_cov/cpat_1.txt", sep=""), sep="\t", row.name=1, header=TRUE)
cpat_2 <- read.table(paste(PATH_results,"/myRegs_cov/cpat_2.txt", sep=""), sep="\t", row.name=1, header=TRUE)
cpat_3 <- read.table(paste(PATH_results,"/myRegs_cov/cpat_3.txt", sep=""), sep="\t", row.name=1, header=TRUE)
cpat_4 <- read.table(paste(PATH_results,"/myRegs_cov/cpat_4.txt", sep=""), sep="\t", row.name=1, header=TRUE)
cpat_5 <- read.table(paste(PATH_results,"/myRegs_cov/cpat_5.txt", sep=""), sep="\t", row.name=1, header=TRUE)
cpat_6 <- read.table(paste(PATH_results,"/myRegs_cov/cpat_6.txt", sep=""), sep="\t", row.name=1, header=TRUE)
cpat_7 <- read.table(paste(PATH_results,"/myRegs_cov/cpat_7.txt", sep=""), sep="\t", row.name=1, header=TRUE)
cpat_8 <- read.table(paste(PATH_results,"/myRegs_cov/cpat_8.txt", sep=""), sep="\t", row.name=1, header=TRUE)

cpat_out <- rbind(cpat_1, cpat_2, cpat_3, cpat_4, cpat_5, cpat_6, cpat_7, cpat_8)

nc <- cpat_out[cpat_out$Coding.Label=="no", ]

coding <- cpat_out[cpat_out$Coding.Label=="yes", ]

rRegs_complete_coff_nc <- rRegs_complete_coff[names(rRegs_complete_coff) %in% as.character(nc$Sequence.Name)]

rRegs_complete_coff_coding <- rRegs_complete_coff[names(rRegs_complete_coff) %in% as.character(coding$Sequence.Name)]

export.bed(rRegs_complete_coff_nc, paste(PATH_results,"/rn6_genome/rRegs_complete_coff_nc.bed", sep=""))
export(asGFF(rRegs_complete_coff_nc), paste(PATH_results,"/myRegs_cov/rRegs_complete_coff_nc.gtf", sep=""), format="gtf")

export.bed(rRegs_complete_coff_coding, paste(PATH_results,"/rn6_genome/rRegs_complete_coff_coding.bed", sep=""))
export(asGFF(rRegs_complete_coff_coding), paste(PATH_results,"/myRegs_cov/rRegs_complete_coff_coding.gtf", sep=""), format="gtf")

toc_coff_nc <- toc_coff[as.character(toc_coff$name) %in% as.character(nc$Sequence.Name),]

toc_coff_coding <- toc_coff[as.character(toc_coff$name) %in% as.character(coding$Sequence.Name),]

save(toc_coff, rRegs_complete_coff, name_mapping_coff, toc_coff_nc, toc_coff_coding, rRegs_complete_coff_nc, rRegs_complete_coff_coding, coding, nc, file =  paste(PATH_results,"/myRegs_cov/rRegs_complete_coff.RData", sep=""))

###############################################
system("./genomic_context.sh")

source("adjust_lncRNA_ids.R")

source("DE.R")

source("annotate_expression.R", echo=TRUE)

source("closest_gene_lincs.R", echo=TRUE)








