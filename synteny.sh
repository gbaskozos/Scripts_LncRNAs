#!/bin/sh


#synetny blocks downloaded from http://bioinfo.konkuk.ac.kr/synteny_portal/
#Resolution: 150000 bp
#Reference species: Human (hg38)
#Target species: Mouse (mm10), Rat (rn5)

PATH_files="$HOME/Desktop/LncRNA_v2"

awk 'BEGIN {FS = "[.: -]"; OFS = "\t" } {if ($1 == "hg38") print $2,$3,$4,"block" (NR-1)/5,0,$5;}' $PATH_files/synteny_blocks_1.txt | sed 's/chr//' > $PATH_files/hg38.synteny_blocks

awk 'BEGIN {FS = "[.: -]"; OFS = "\t" } {if ($1 == "mm10") print $2,$3,$4,"block" (NR-2)/5,0,$5;}' $PATH_files/synteny_blocks_1.txt | sed 's/chr//' > $PATH_files/mm10.synteny_blocks

awk 'BEGIN {FS = "[.: -]"; OFS = "\t" } {if ($1 == "rn5") print $2,$3,$4,"block" (NR-3)/5,0,$5;}' $PATH_files/synteny_blocks_1.txt > $PATH_files/rn5.synteny_blocks


awk 'BEGIN {FS = "\t"; OFS = "\t" } {if ($6 == "") print $1,$2,$3,$4,$5,"-"; else print;}' $PATH_files/hg38.synteny_blocks > $PATH_files/hg38.synteny_blocks.bed

awk 'BEGIN {FS = "\t"; OFS = "\t" } {if ($6 == "") print $1,$2,$3,$4,$5,"-"; else print;}' $PATH_files/mm10.synteny_blocks > $PATH_files/mm10.synteny_blocks.bed

awk 'BEGIN {FS = "\t"; OFS = "\t" } {if ($6 == "") print $1,$2,$3,$4,$5,"-"; else print;}' $PATH_files/rn5.synteny_blocks > $PATH_files/rn5.synteny_blocks.bed

###RN5 regions were lifted over to RN6 using the liftOverTool from UCSC

cat $PATH_files/hglft_genome_2ac2_dfe7d0.bed | sed 's/chr//' > $PATH_files/rn6.synteny_blocks.bed

bedtools intersect -wo -S -header -a $HOME/Desktop/LncRNA_v2/Rat/rn6_genome/Novel_lncRNAs_rat_DRG.bed -b $HOME/Desktop/LncRNA_v2/rn6.synteny_blocks.bed > $HOME/Desktop/LncRNA_v2/Rat_syntenic_lncRNAs.bed

bedtools intersect -wo -S -header -a $HOME/Desktop/LncRNA_v2/Mouse/mm10_genome/Novel_lncRNAs_mouse_DRG.bed -b $HOME/Desktop/LncRNA_v2/mm10.synteny_blocks.bed > $HOME/Desktop/LncRNA_v2/Mouse_syntenic_lncRNAs.bed

bedtools intersect -wo -S -header -a $HOME/Desktop/LncRNA_v2/IPS/HS_genome/Novel_lncRNAs_IPS.bed -b $HOME/Desktop/LncRNA_v2/hg38.synteny_blocks.bed > $HOME/Desktop/LncRNA_v2/Hs_syntenic_lncRNAs.bed
