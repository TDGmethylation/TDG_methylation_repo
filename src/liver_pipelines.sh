#!/bin/bash

#Step1: Concatenate CpH, CHH, CHG files to CN w/ R1 & R2.

declare -i i=0

filepath1='/Volumes/Genome_data/WGBS_liver_ENCSR108ESU/ENCODE_raw/ENCFF099BSV_CHG.bed'
filepath2='/Volumes/Genome_data/WGBS_liver_ENCSR108ESU/ENCODE_raw/ENCFF415EFB_CHH.bed'
filepath3='/Volumes/Genome_data/WGBS_liver_ENCSR108ESU/ENCODE_raw/ENCFF577VGR_CpG.bed'

# filetest1='/Volumes/Genome_data/WGBS_H1/test/test1.bed'
# filetest2='/Volumes/Genome_data/WGBS_H1/test/test2.bed'

# Pipeline 1 start here, we concatenate all files and filter 
# the rows that have either coverage = 0 or metthylation level = 0.
# After that, we do weighted mean of the coverage between R1 & R2.
# The output format is the same as on_hardac folder $1=chr, $2=start, $3=end
# $4=coverage, $5=methylation %, $6 = strand. $7 is pending, waiting for align.

cat $filepath1 $filepath2 $filepath3 |
awk -F '\t' '$10>0 && $11>0 {print $1, $2, $3, $4=$10, $5=$11, $6}' OFS='\t' |
awk -F '\t' '{print $1, $2, $3, $4, $4*$5/100, $6}' OFS='\t' |
awk -F '\t' '{$5=sprintf("%.0f", $5)}1' OFS='\t' |
sort -n |
bedtools groupby -g 1,2,3 -c 4,5,6 -o sum,sum,distinct -i "stdin" |
awk 'BEGIN{OFS="\t"} {print $1, $2, $3, $4, $5=$5/$4, $6}' |
awk -F '\t' 'length($1)<=5 && $1!="chrM" {print $1,$2-=6, $3+=7, $4, $5, $6}' OFS='\t' | 
 bedtools getfasta -fi "/Volumes/Genome_data/Genome_assembly/UCSC/hg38.fa" -bed "stdin" -s -bedOut \
  > "/Volumes/Genome_data/WGBS_liver_ENCSR108ESU/Aligned/Aligned_liver_preupload.bed"
#I will stop here today. Here, the files are weighted-mean-ed. Definitely check it first before 
# moving to the next step: change chr start-end for aligning.


