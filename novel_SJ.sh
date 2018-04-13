#!/bin/sh

PATH_SJ="/media/george/data3tb/Sequencing_Data/Mouse_DRG_sni_jmogil/merged_bams/SJ"

for file in $PATH_SJ/*SJ.out.tab
do
awk -F "\t" '$6=="1" {print }' $file > $file.novel.tab
done



