#!/bin/bash


mkdir $HOME/Desktop/LncRNA_v2/IPS/antisense
mkdir $HOME/Desktop/LncRNA_v2/IPS/non_antisense

bedtools intersect -wo -S -header -split -a $HOME/Desktop/LncRNA_v2/IPS/HS_genome/rRegs_complete_coff_nc.bed -b $HOME/Desktop/LncRNA_v2/IPS/HS_genome/ENSEMBL_HS_genes.bed > $HOME/Desktop/LncRNA_v2/IPS/antisense/antisense_lncs.bed

bedtools intersect -wo -S -header -split -a $HOME/Desktop/LncRNA_v2/IPS/HS_genome/rRegs_complete_coff_nc.bed -b $HOME/Desktop/LncRNA_v2/IPS/HS_genome/HS_refseq.bed > $HOME/Desktop/LncRNA_v2/IPS/antisense/antisense_lncs_refseq.bed

sort -k1,1 -k2,2n $HOME/Desktop/LncRNA_v2/IPS/antisense/antisense_lncs.bed > $HOME/Desktop/LncRNA_v2/IPS/antisense/antisense_lncs_sorted.bed
sort -k1,1 -k2,2n $HOME/Desktop/LncRNA_v2/IPS/antisense/antisense_lncs_refseq.bed > $HOME/Desktop/LncRNA_v2/IPS/antisense/antisense_lncs_refseq_sorted.bed

bedtools subtract -s -A -a $HOME/Desktop/LncRNA_v2/IPS/HS_genome/rRegs_complete_coff_nc.bed -b $HOME/Desktop/LncRNA_v2/IPS/antisense/antisense_lncs_sorted.bed > $HOME/Desktop/LncRNA_v2/IPS/non_antisense/lncs_non_antisense.bed

sort -k1,1 -k2,2n $HOME/Desktop/LncRNA_v2/IPS/non_antisense/lncs_non_antisense.bed > $HOME/Desktop/LncRNA_v2/IPS/non_antisense/lncs_non_antisense.sorted.bed

bedtools intersect -wo -s -header -split -f 1.0 -a $HOME/Desktop/LncRNA_v2/IPS/non_antisense/lncs_non_antisense.sorted.bed -b $HOME/Desktop/LncRNA_v2/IPS/HS_genome/intronsByGene.bed > $HOME/Desktop/LncRNA_v2/IPS/antisense/intronic_all.bed

bedtools intersect -u -s -header -split -F 0.3 -a $HOME/Desktop/LncRNA_v2/IPS/antisense/intronic_all.bed -b $HOME/Desktop/LncRNA_v2/IPS/HS_genome/intronsByGene.bed > $HOME/Desktop/LncRNA_v2/IPS/antisense/intronic_retained.bed

bedtools subtract -s -A -a $HOME/Desktop/LncRNA_v2/IPS/non_antisense/lncs_non_antisense.sorted.bed -b $HOME/Desktop/LncRNA_v2/IPS/antisense/intronic_all.bed > $HOME/Desktop/LncRNA_v2/IPS/non_antisense/intergenic.bed

bedtools subtract -s -A -a $HOME/Desktop/LncRNA_v2/IPS/antisense/intronic_all.bed -b $HOME/Desktop/LncRNA_v2/IPS/antisense/intronic_retained.bed > $HOME/Desktop/LncRNA_v2/IPS/antisense/intronic.bed

sort -k1,1 -k2,2n $HOME/Desktop/LncRNA_v2/IPS/non_antisense/intergenic.bed > $HOME/Desktop/LncRNA_v2/IPS/non_antisense/intergenic.sorted.bed

sort -k1,1 -k2,2n $HOME/Desktop/LncRNA_v2/IPS/antisense/intronic.bed > $HOME/Desktop/LncRNA_v2/IPS/non_antisense/intronic.sorted.bed

sort -k1,1 -k2,2n $HOME/Desktop/LncRNA_v2/IPS/HS_genome/ENSEMBL_HS_genes.bed > $HOME/Desktop/LncRNA_v2/IPS/HS_genome/ENSEMBL_HS_genes.sorted.bed

bedtools closest -s -io -D a -t all -a $HOME/Desktop/LncRNA_v2/IPS/non_antisense/intergenic.sorted.bed -b $HOME/Desktop/LncRNA_v2/IPS/HS_genome/ENSEMBL_HS_genes.sorted.bed >  $HOME/Desktop/LncRNA_v2/IPS/non_antisense/lncs_intergenic_closest.bed

bedtools closest -s -io -D a -t all -id -a $HOME/Desktop/LncRNA_v2/IPS/non_antisense/intergenic.sorted.bed -b $HOME/Desktop/LncRNA_v2/IPS/HS_genome/ENSEMBL_HS_genes.sorted.bed >  $HOME/Desktop/LncRNA_v2/IPS/non_antisense/lncs_intergenic_closest_upstream.bed

bedtools closest -s -io -D a -t all -iu -a $HOME/Desktop/LncRNA_v2/IPS/non_antisense/intergenic.sorted.bed -b $HOME/Desktop/LncRNA_v2/IPS/HS_genome/ENSEMBL_HS_genes.sorted.bed >  $HOME/Desktop/LncRNA_v2/IPS/non_antisense/lncs_intergenic_closest_downstream.bed


####ENSEMBL lncRNAs
sort -k1,1 -k2,2n $HOME/Desktop/LncRNA_v2/IPS/HS_genome/ENSEMBL_HS_antisense.bed > $HOME/Desktop/LncRNA_v2/IPS/antisense/ens_antisense.sorted.bed

sort -k1,1 -k2,2n $HOME/Desktop/LncRNA_v2/IPS/HS_genome/ENSEMBL_HS_intergenic.bed > $HOME/Desktop/LncRNA_v2/IPS/non_antisense/ens_intergenic.sorted.bed

bedtools intersect -wo -S -header -split -a $HOME/Desktop/LncRNA_v2/IPS/antisense/ens_antisense.sorted.bed -b $HOME/Desktop/LncRNA_v2/IPS/HS_genome/ENSEMBL_HS_genes.sorted.bed > $HOME/Desktop/LncRNA_v2/IPS/antisense/lncs_antisense_ENSEMBL.bed

bedtools closest -s -io -D a -t all -a $HOME/Desktop/LncRNA_v2/IPS/non_antisense/ens_intergenic.sorted.bed -b $HOME/Desktop/LncRNA_v2/IPS/HS_genome/ENSEMBL_HS_genes.sorted.bed >  $HOME/Desktop/LncRNA_v2/IPS/non_antisense/lncs_intergenic_closest_ENSEMBL.bed



