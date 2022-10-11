#!/bin/bash
set -exuo pipefail 
ANN_DIR=../../../ann/

bedtools coverage -counts -a ${ANN_DIR}/transcript_promoter_merged.bed -b E21-A0_S19__trimmed_R1.sorted.bam > Rad_rep1.txt
bedtools coverage -counts -a ${ANN_DIR}/transcript_promoter_merged.bed -b E121-A0-ab104232_S10__trimmed_R1.sorted.bam > Rad_rep2.txt
bedtools coverage -counts -a ${ANN_DIR}/transcript_promoter_merged.bed -b 1003622_sorted.bam > Swi6_rep1.txt

bedtools coverage -a ${ANN_DIR}/transcript_promoter_merged.bed -b 1003622_sorted.bam  > Swi6_rep1_pe.txt
bedtools coverage -a ${ANN_DIR}/transcript_promoter_merged.bed -b E21-A0_S19__trimmed_R1.sorted.bam > Rad_rep1_pe.txt
bedtools coverage -a ${ANN_DIR}/transcript_promoter_merged.bed -b E121-A0-ab104232_S10__trimmed_R1.sorted.bam > Rad_rep2_pe.txt

samtools view -b -F 4 1003622_sorted.bam > Swi6_uniq.bam
samtools view -b -F 4 E21-A0_S19__trimmed_R1.sorted.bam > Rad_rep1_uniq.bam
samtools view -b -F 4 E121-A0-ab104232_S10__trimmed_R1.sorted.bam > Rad_rep2_uniq.bam
bedtools coverage -a ${ANN_DIR}/transcript_promoter_merged.bed -b Swi6_uniq.bam  > Swi6_rep1_pe_uniq.txt
bedtools coverage -a ${ANN_DIR}/transcript_promoter_merged.bed -b Rad_rep1_uniq.bam > Rad_rep1_pe_uniq.txt
bedtools coverage -a ${ANN_DIR}/transcript_promoter_merged.bed -b Rad_rep2_uniq.bam > Rad_rep2_pe_uniq.txt

Rscript lorenz_curve.R
#samtool -c 1003622_sorted.bam  839002
#E121: 4801789 -> 4196755 + 250020
#E21: 4697588 -> 4202122 + 242632
#22: 839002 -> 473120 + 191002

