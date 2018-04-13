#PHANTOM TSS (classified as TRUE TSS), lifted over from mm9 to mm10
PATH_results <- "/home/george/Desktop/LncRNA_v2/Mouse"

lncs <- read.table(paste(PATH_results,"/mm10_genome/Novel_lncRNAs_mouse_DRG.bed", sep=""), sep="\t", stringsAsFactors=FALSE)

load(paste(PATH_results,"/DE/expressed_nc_ensembl.RData", sep=""))

lncs_ens <- read.table(paste(PATH_results,"/mm10_genome/ENSEMBL_mm10_LncRNAs.bed", sep=""), sep="\t", stringsAsFactors=FALSE)

lncs_ens_expressed <- lncs_ens[as.character(lncs_ens$V4) %in% rownames(expressed_nc_ensembl),]


ts <- function(x) { 
if(x[[6]] == "+") {
start <- as.numeric(x[2]) - 100
end <- start + 100 
}
else if(x[[6]] == "-") {
end <- as.numeric(x[3]) + 100
start <- end - 100
} 

coord <- data.frame(start=start, end=end)


return(coord)
}


lncs_ts <- apply(lncs, 1, ts)

lncs$V2 <- unlist(sapply(lncs_ts,function(x) x[[1]]))


lncs$V3 <- unlist(sapply(lncs_ts,function(x) x[[2]]))

write.table(file = paste(PATH_results,"/tss/start_sites_all_lncs.bed", sep=""), lncs, sep="\t", row.name=FALSE, col.name=FALSE, quote=FALSE)


lncs_ens_ts <- apply(lncs_ens, 1, ts)

lncs_ens$V2 <- unlist(sapply(lncs_ens_ts,function(x) x[[1]]))


lncs_ens$V3 <- unlist(sapply(lncs_ens_ts,function(x) x[[2]]))

write.table(file = paste(PATH_results,"/tss/start_sites_ens_lncs.bed", sep=""), lncs_ens, sep="\t", row.name=FALSE, col.name=FALSE, quote=FALSE)

lncs_ens_ts_expressed <- apply(lncs_ens_expressed, 1, ts)

lncs_ens_expressed$V2 <- unlist(sapply(lncs_ens_ts_expressed,function(x) x[[1]]))

lncs_ens_expressed$V3 <- unlist(sapply(lncs_ens_ts_expressed,function(x) x[[2]]))

write.table(file = paste(PATH_results,"/tss/start_sites_ens_lncs_expressed.bed", sep=""), lncs_ens_expressed, sep="\t", row.name=FALSE, col.name=FALSE, quote=FALSE)

system("/home/george/Desktop/Desktop/LncRNA_v2/Mouse/scripts/tss_evidence.sh")

####################################
tss_lncs_evidence <- read.table(file = "/home/george/Desktop/LncRNA_v2/Mouse/tss/lncs_tss_evidence.bed", sep="\t", header=FALSE)

tss_lncs_distance <- read.table(file = "/home/george/Desktop/LncRNA_v2/Mouse/tss/lncs_tss_distance.bed", sep="\t", header=FALSE)

ens_lncs_distance <- read.table(file = "/home/george/Desktop/LncRNA_v2/Mouse/tss/ens_lncs_tss_distance.bed", sep="\t", header=FALSE)
###################################
dsBase.iqr <- tss_lncs_distance
 
# Create a variable/vector/collection of the column names you want to remove outliers on.
vars <- c("V22")

Outliers <- c()
 
for(i in vars){
 
  max <- quantile(dsBase.iqr[,i],0.75, na.rm=TRUE) + (IQR(dsBase.iqr[,i], na.rm=TRUE) * 1.5 )
  min <- quantile(dsBase.iqr[,i],0.25, na.rm=TRUE) - (IQR(dsBase.iqr[,i], na.rm=TRUE) * 1.5 )
  
  idx <- which(dsBase.iqr[,i] < min | dsBase.iqr[,i] > max)
  
  print(paste(i, length(idx), sep=''))
  
  Outliers <- c(Outliers, idx) 
}
 
Outliers <- sort(Outliers)
 
dsBase.iqr <- dsBase.iqr[-Outliers,]



pdf(file = paste(PATH_results,"/tss/tss_density.pdf", sep=""))
layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))

d_lncs_no <- density(dsBase.iqr$V22)
plot(d_lncs_no, main="Kernel distances between Novel LncRNAs (no outliers) and TSS", xlab="distance from TSS", cex.main = 0.8, cex.lab=0.7, cex.axis=0.7)
polygon(d_lncs_no, col="red", border="blue")
legend("topleft", title="Distances > 1.5xIQR were removed", legend=c("Median:", median(tss_lncs_distance$V22), "Mean:", mean(tss_lncs_distance$V22), "Median after removal:", median(dsBase.iqr$V22), "Mean after removal:", mean(dsBase.iqr$V22)), cex=0.5, bg="white")

d_lncs <- density(tss_lncs_distance$V22) 
plot(d_lncs, main="Kernel distances Novel LncRNAs and TSS", xlab="distance from TSS", cex.main = 0.8, cex.lab=0.7, cex.axis=0.7)
polygon(d_lncs, col="red", border="red")

d_ens <- density(ens_lncs_distance$V22)
plot(d_ens, main="Kernel distances ENSEMBL LncRNAs and TSS", xlab="distance from TSS", cex.main = 0.8, cex.lab=0.7, cex.axis=0.7)
polygon(d_ens, col="blue", border="blue")


dev.off()
#################################################################


res_sham_vs_sni_B10.D2 <- read.csv("/home/george/Desktop/LncRNA_v2/Mouse/DE/res_sni_vs_sham_B10.D2.csv", row.name=1)
res_sham_vs_sni_BALB.c <- read.csv("/home/george/Desktop/LncRNA_v2/Mouse/DE/res_sni_vs_sham_BALB.c.csv", row.name=1)
res_sni_diff_B10.D2_vs_BALB.c <- read.csv("/home/george/Desktop/LncRNA_v2/Mouse/DE/res_sni_diff_B10.D2_vs_BALB.c.csv", row.name=1)

tss_evidence <- read.table("/home/george/Desktop/LncRNA_v2/Mouse/tss/lncs_tss_evidence.bed", sep="\t", stringsAsFactors=FALSE)


complete_lncs_tss <- lncs[lncs$V4 %in% tss_evidence$V4,]

write.csv(file = "/home/george/Desktop/LncRNA_v2/Mouse/tss/complete_lncs_tss.csv", complete_lncs_tss)

