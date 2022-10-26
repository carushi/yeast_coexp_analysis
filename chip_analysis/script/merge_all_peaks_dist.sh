#!/bin/bash
set -euxo pipefail 

for end in pe se;
do
for peak in broad narrow ;
do
truncate -s 0 all_peaks_${end}_${peak}.bed
done
done


for i in $( seq 0 8 );
do
OFFSET=$(( 3+${i} ))
for end in pe se;
do
for peak in broad narrow ;
do
echo $OFFSET
sort -k1,1 -k2,2n peak_${end}_${OFFSET}_1_peaks.${peak}Peak >> all_peaks_${end}_${peak}.bed
sort -k1,1 -k2,2n peak_${end}_${i}_2_peaks.${peak}Peak >> all_peaks_${end}_${peak}.bed
done
done
done

for end in pe se ;
do
for peak in broad narrow ;
do
	sort -k1,1 -k2,2n all_peaks_${end}_${peak}.bed > all_peaks_${end}_${peak}_sorted.bed
	bedtools merge -i all_peaks_${end}_${peak}_sorted.bed -d 0  > merged_peaks_${end}_${peak}_0.bed
for i in $( seq 10 10 100 );
do
	bedtools merge -i all_peaks_${end}_${peak}_sorted.bed -d $i > merged_peaks_${end}_${peak}_${i}.bed
done
for i in 200 300 400 500 750 1000 ;
do
	bedtools merge -i all_peaks_${end}_${peak}_sorted.bed -d $i > merged_peaks_${end}_${peak}_${i}.bed
done

done
done
