###################################################################################################
##Find intergenic regions with suficient continuous expression coverage >= dep over a length >len##
###################################################################################################
findRegs <- function(reads, len, dep) {
require(IRanges)
require(GenomicRanges)
  reads <- grglist(reads, use.mcols=TRUE)
  # Separate strands
  forward <- reads[strand(reads) == "+"]
  reverse <- reads[strand(reads) == "-"]
  #Calculate coverage for each strand
  cov_forward <- coverage(forward)
  cov_reverse <- coverage(reverse)
  #Identify island of expression with read coverage of >= dep
  islands.i_forward <- slice(cov_forward, lower = dep) 
  islands.i_reverse <- slice(cov_reverse, lower = dep)
  transcribedRegions.i_forward <- islands.i_forward[width(islands.i_forward) > len]
  transcribedRegions.i_reverse <- islands.i_reverse[width(islands.i_reverse) > len]
  t_forward <- reduce(as(transcribedRegions.i_forward, "GRanges"))
  strand(t_forward) = "+"
  t_reverse <- reduce(as(transcribedRegions.i_reverse, "GRanges"))
  strand(t_reverse) = "-"
  return(list(regions = reduce(c(t_forward, t_reverse)), cov_forward = cov_forward, cov_reverse = cov_reverse))
}  

######################################################################################
## Read bam files into Granges lists and discard reads overlapping with known genes ##
######################################################################################

BAM_to_IOE <- function(bamfile, PATH, PATH_results, igRangesExt, param, len=100, dep=2, suffix=11) {

require(Rsamtools)
require(GenomicFeatures)
require(GenomicAlignments)
require(BiocParallel)
require(GenomicRanges)
require(BiocParallel)

name <- substr(bamfile$path, nchar(PATH) + 2, nchar(bamfile$path)-suffix)

sink(stdout(), type = "message") # sink messages to stdout
message(paste("Processing",name,"\n", sep=" "))
sink(NULL, type="message")


write(paste("Log file for ", name, sep=""), file = paste(PATH_results, "/myRegs_cov/",name,"_log.txt", sep=""), append = FALSE)

open(bamfile)
bam <- NULL
cov_F <- NULL
cov_B <- NULL
chunk_ind <- 1
#Read read-pairs from bam file
while(length(chunk <- readGAlignmentPairs(file = bamfile, index = bamfile, use.names=TRUE, strandMode=2, param=param))) {
#Subset only reads that fall within non-annotated regins
inp <- subsetByOverlaps(chunk, igRangesExt, type="within", ignore.strand=FALSE) 
write(paste("Chunk nr:", chunk_ind,". Number of reads (pairs) overlapping with intergenic regions: ", length(inp), paste = ""), file = paste(PATH_results, "/myRegs_cov/",name,"_log.txt", sep=""), append = TRUE)
res_list <- findRegs(inp, len, dep)
rm(inp)
write(paste("Chunk nr:", chunk_ind,". Number of reads overlapping with islands of expression coverage >= ",dep,  " & width > ", len, " : ", length(res_list$regions), paste = ""), file = paste(PATH_results, "/myRegs_cov/",name,"_log.txt", sep=""), append = TRUE)
if (chunk_ind == 1) {
bam <- res_list$regions
cov_F <- res_list$cov_forward
cov_B <- res_list$cov_reverse
rm(res_list)
rm(chunk) }
else {
bam <- c(bam, res_list$regions)
cov_F <- cov_F + res_list$cov_forward
cov_B <- cov_F + res_list$cov_forward
rm(res_list)
rm(chunk)
}
log_out <- paste("Integrating chunk nr ",chunk_ind, " to bam file ", name,". Nr of reads overlapping with intergenic regions: ", length(bam), "\n", sep="")
write(log_out, file = paste(PATH_results, "/myRegs_cov/",name,"_log.txt", sep=""), append = TRUE)
chunk_ind <- chunk_ind + 1
}
close(bamfile)
save(bam, cov_F, cov_B, file = paste(PATH_results, "/myRegs_cov/", name, ".RData", sep=""))
rm(bam)
rm(cov_F)
rm(cov_B)
rm(chunk)
gc()
flush.console() 

}


