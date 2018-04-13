# Scripts_LncRNAs
Workflow to identify LncRNAs in Mouse DRG, Rat DRG and Human IPSc derived neurons

This is a repository with all the scripts used for the work "Comporehensive analysis of Long non-coding RNA expression in dorsal root ganglion reveals cell type specificity and dysregulation following nerve injury".

There are mostly R scripts, which rely on bioconductor packages.

The script "identify_LncRNAs.R" is the main script of the pipeline.
The script "genomic_context.sh" uses bedtools to classify putative LncRNAs into antisense, intronic and intergenic.
The script "adjust_lncRNA_ids.R" is a helper script to map LncRNA IDs and names in all output files (.bed format, .gtf format, table of counts).
The script "DE_and_WGCNA.R" uses DESeq2 and WGCNA to calculate differential expression and create a weighted co-expression network.
The script "annotate_expression.R" annotates the expression of antisense (and intronic) LncRNAs in the context of the gene on the opposite strand.
The script "closest_gene_lincs.R" annotates the expression of intergenic LncRNAs in the context of the closest neighbouring gene.
The script "calc_start_sites.R" analyses LncRNAs in the context of Phantom annotated transcription start sites.
The script "synteny.sh" identifies LncRNAs overlapping syntenic blocks in Human, Rat and Mouse.
The script "homologs.R" identifies LncRNAs antisense of orthologous genes in Human, Rat and Mouse.

Helper functions:
"BAM_to_IOE.R" reads BAM files and identified continuously expressed regions.
bamfile: list of bamfiles, Rsamtools::BamFileList
PATH: absolute path of BAM files
PATH_results: output path
igRangesExt: Intergenic regions GRanges
param: scan BAM parameters, Rsamtools::ScanBamParam
len: minimum length of intergenic regions
dep: minimum coverage depth of intergenic regions
suffix: length of file suffix

"dropDetect.R" identifies sudden drops in coverage that could indicate introns or trasncription ends.
Arguments:
coverage: RLE vcoverage vector
start: region start bases
seqnames: region chromosome
strand: region strand
lag: the lag of the smoothing window
threshold: drop/peak detection threshold in steps of variance (z-score)
length: minimum intron length
influence: influence of previous drop/peak to current signal
intron_identification: perform intron identification TRUE/FALSE

"novel_SJ.sh" selectcs only de novo identified splicing junctions from STAR aligner
