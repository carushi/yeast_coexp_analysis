#!/bin/bash

ORF_DIR=../../../ann/
LANG=en_EN sort -k 1 -t $'\t' de_result_sc_genotype.tsv > de_result_sc_genotype_sorted.tsv
LANG=en_EN sort -k 1 -t $'\t' de_result_mc_genotype.tsv > de_result_mc_genotype_sorted.tsv
join -a 1 -1 1 -2 1 de_result_mc_genotype_sorted.tsv ${ORF_DIR}/SGD_orf_id.tab -t' ' > de_result_mc_genotype_sorted_name.tsv
join -a 1 -1 1 -2 1 de_result_sc_genotype_sorted.tsv ${ORF_DIR}/SGD_orf_id.tab -t' ' > de_result_sc_genotype_sorted_name.tsv

LANG=en_EN sort -k 1 -t $'\t' de_result_sc_time.tsv > de_result_sc_time_sorted.tsv
LANG=en_EN sort -k 1 -t $'\t' de_result_sc_each_genot.tsv > de_result_sc_each_genot_sorted.tsv
join -a 1 -1 1 -2 1 de_result_sc_time_sorted.tsv ${ORF_DIR}/SGD_orf_id.tab -t' ' > de_result_sc_time_sorted_name.tsv
join -a 1 -1 1 -2 1 de_result_sc_each_genot_sorted.tsv ${ORF_DIR}/SGD_orf_id.tab -t' ' > de_result_sc_each_genot_sorted_name.tsv
