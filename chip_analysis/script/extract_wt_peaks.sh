#/bin/bash
# extract wt peaks

ANN_DIR="../../../ann/"
if true;
then
    grep -e peak_pe_0_2 -e peak_pe_1_2 -e peak_pe_2_2 -e peak_pe_3_1 -e peak_pe_4_1 -e peak_pe_5_1 all_peaks_pe_narrow_sorted.bed > all_peaks_pe_narrow_wt_sorted.bed
    bedtools merge -i all_peaks_pe_narrow_wt_sorted.bed > all_peaks_pe_narrow_wt_merged.bed
    grep -e peak_pe_0_2 -e peak_pe_3_1  all_peaks_pe_narrow_sorted.bed > all_peaks_pe_narrow_wt1_sorted.bed
    bedtools merge -i all_peaks_pe_narrow_wt1_sorted.bed > all_peaks_pe_narrow_wt1_merged.bed
    grep -e peak_pe_1_2 -e peak_pe_4_1  all_peaks_pe_narrow_sorted.bed > all_peaks_pe_narrow_wt2_sorted.bed
    bedtools merge -i all_peaks_pe_narrow_wt2_sorted.bed > all_peaks_pe_narrow_wt2_merged.bed
    grep -e peak_pe_2_2 -e peak_pe_5_1  all_peaks_pe_narrow_sorted.bed > all_peaks_pe_narrow_wt3_sorted.bed
    bedtools merge -i all_peaks_pe_narrow_wt3_sorted.bed > all_peaks_pe_narrow_wt3_merged.bed
fi

if true;
then
    grep -e peak_pe_3_2 -e peak_pe_4_2 -e peak_pe_5_2 -e peak_pe_6_1 -e peak_pe_7_1 -e peak_pe_8_1 all_peaks_pe_narrow_sorted.bed > all_peaks_pe_narrow_rad_sorted.bed
    bedtools merge -i all_peaks_pe_narrow_rad_sorted.bed > all_peaks_pe_narrow_rad_merged.bed
    grep -e peak_pe_3_2 -e peak_pe_6_1  all_peaks_pe_narrow_sorted.bed > all_peaks_pe_narrow_rad1_sorted.bed
    bedtools merge -i all_peaks_pe_narrow_rad1_sorted.bed > all_peaks_pe_narrow_rad1_merged.bed
    grep -e peak_pe_4_2 -e peak_pe_7_1  all_peaks_pe_narrow_sorted.bed > all_peaks_pe_narrow_rad2_sorted.bed
    bedtools merge -i all_peaks_pe_narrow_rad2_sorted.bed > all_peaks_pe_narrow_rad2_merged.bed
    grep -e peak_pe_5_2 -e peak_pe_8_1  all_peaks_pe_narrow_sorted.bed > all_peaks_pe_narrow_rad3_sorted.bed
    bedtools merge -i all_peaks_pe_narrow_rad3_sorted.bed > all_peaks_pe_narrow_rad3_merged.bed
fi

if true;
then
    grep -e peak_pe_6_2 -e peak_pe_7_2 -e peak_pe_8_2 -e peak_pe_9_1 -e peak_pe_10_1 -e peak_pe_11_1 all_peaks_pe_narrow_sorted.bed > all_peaks_pe_narrow_mrc_sorted.bed
    bedtools merge -i all_peaks_pe_narrow_mrc_sorted.bed > all_peaks_pe_narrow_mrc_merged.bed
    grep -e peak_pe_6_2 -e peak_pe_9_1  all_peaks_pe_narrow_sorted.bed > all_peaks_pe_narrow_mrc1_sorted.bed
    bedtools merge -i all_peaks_pe_narrow_mrc_sorted.bed > all_peaks_pe_narrow_mrc1_merged.bed
    grep -e peak_pe_7_2 -e peak_pe_10_1  all_peaks_pe_narrow_sorted.bed > all_peaks_pe_narrow_mrc2_sorted.bed
    bedtools merge -i all_peaks_pe_narrow_mrc_sorted.bed > all_peaks_pe_narrow_mrc2_merged.bed
    grep -e peak_pe_8_2 -e peak_pe_11_1  all_peaks_pe_narrow_sorted.bed > all_peaks_pe_narrow_mrc3_sorted.bed
    bedtools merge -i all_peaks_pe_narrow_mrc_sorted.bed > all_peaks_pe_narrow_mrc3_merged.bed
fi

for sp in wt rad mrc
do
for tail in '' 1 2 3
do
bedtools intersect -u -wa -a all_peaks_pe_narrow_${sp}${tail}_merged.bed -b ${ANN_DIR}/hyper_chipable_rom.bed > all_peaks_pe_narrow_${sp}${tail}_merged_hcov.bed
bedtools intersect -v -wa -a all_peaks_pe_narrow_${sp}${tail}_merged.bed -b ${ANN_DIR}/histone_mark/hyper_chipable_rom.bed > all_peaks_pe_narrow_${sp}${tail}_merged_hcno.bed
wc -l  all_peaks_pe_narrow_${sp}${tail}_merged_hcov.bed
bedtools intersect -u  -wa -a  all_peaks_pe_narrow_${sp}${tail}_merged_hcov.bed -b ${ANN_DIR}/Saccharomyces_cerevisiae.R64-1-1.98_chr_sorted_transcript.gtf | wc -l
bedtools intersect -u  -wa -a  all_peaks_pe_narrow_${sp}${tail}_merged_hcno.bed -b ${ANN_DIR}/Saccharomyces_cerevisiae.R64-1-1.98_chr_sorted_transcript.gtf | wc -l
bedtools intersect -u -f 0.5 -wa -a  all_peaks_pe_narrow_${sp}${tail}_merged_hcov.bed -b ${ANN_DIR}/ann/Saccharomyces_cerevisiae.R64-1-1.98_chr_sorted_transcript.gtf | wc -l
bedtools intersect -u -f 0.5 -wa -a  all_peaks_pe_narrow_${sp}${tail}_merged_hcno.bed -b ${ANN_DIR}/ann/Saccharomyces_cerevisiae.R64-1-1.98_chr_sorted_transcript.gtf | wc -l
done
done

# bedtools intersect -wa -a all_peaks_pe_narrow_wt1_merged.bed -b ../../../histone_mark/hyper_chipable_rom.bed | wc -l
# bedtools intersect -wa -a all_peaks_pe_narrow_wt2_merged.bed -b ../../../histone_mark/hyper_chipable_rom.bed | wc -l
# bedtools intersect -wa -a all_peaks_pe_narrow_wt3_merged.bed -b ../../../histone_mark/hyper_chipable_rom.bed | wc -l


# 238 hyper chipable region
# 214/2504 all
# 182/1867 t1
# 184/1833 t2 
# 202/2394 t3
