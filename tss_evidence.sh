#!/bin/bash


#mkdir /home/george/Desktop/LncRNA_v2/Mouse/tss

sort -k1,1 -k2,2n /home/george/Desktop/LncRNA_v2/Mouse/mm10_genome/Novel_lncRNAs_mouse_DRG.bed > /home/george/Desktop/LncRNA_v2/Mouse/mm10_genome/Novel_lncRNAs_mouse_DRG.sorted.bed

sort -k1,1 -k2,2n /home/george/Desktop/LncRNA_v2/Mouse/tss/tss_mm9_to_m10_lift.bed > /home/george/Desktop/LncRNA_v2/Mouse/tss/tss_mm9_to_m10_lift.sorted.bed

sort -k1,1 -k2,2n /home/george/Desktop/LncRNA_v2/Mouse/tss/start_sites_coding.bed > /home/george/Desktop/LncRNA_v2/Mouse/tss/start_sites_coding.sorted.bed

sort -k1,1 -k2,2n /home/george/Desktop/LncRNA_v2/Mouse/tss/start_sites_ens_lncs.bed > /home/george/Desktop/LncRNA_v2/Mouse/tss/start_sites_ens_lncs.sorted.bed

bedtools intersect -s -header -u -a /home/george/Desktop/LncRNA_v2/Mouse/tss/start_sites_all_lncs.bed -b /home/george/Desktop/LncRNA_v2/Mouse/tss/tss_mm9_to_m10_lift.bed > /home/george/Desktop/LncRNA_v2/Mouse/tss/lncs_tss_evidence.bed

bedtools closest -t first -s -header -D a -a /home/george/Desktop/LncRNA_v2/Mouse/mm10_genome/Novel_lncRNAs_mouse_DRG.sorted.bed -b /home/george/Desktop/LncRNA_v2/Mouse/tss/tss_mm9_to_m10_lift.sorted.bed > /home/george/Desktop/LncRNA_v2/Mouse/tss/lncs_tss_distance.bed

bedtools intersect -s -header -u -a /home/george/Desktop/LncRNA_v2/Mouse/tss/start_sites_ens_lncs.bed -b /home/george/Desktop/LncRNA_v2/Mouse/tss/tss_mm9_to_m10_lift.bed > /home/george/Desktop/LncRNA_v2/Mouse/tss/ens_lncs_tss_evidence.bed

bedtools closest -s -header -D a -a /home/george/Desktop/LncRNA_v2/Mouse/tss/start_sites_ens_lncs.sorted.bed -b /home/george/Desktop/LncRNA_v2/Mouse/tss/tss_mm9_to_m10_lift.sorted.bed > /home/george/Desktop/LncRNA_v2/Mouse/tss/ens_lncs_tss_distance.bed

