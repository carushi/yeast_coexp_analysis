#!/bin/bash
set -exuo pipefail 

DATA_DIR=../data/
GENOME=../../genome/chromFaMasked/
GENOME_DIR=${GENOME}
ANN_DIR=../../ann/

DIRS="300600 301136 301417 301760 302393 303277 305186 305218 305704 305826 306560 306595 306600"

ulimit -n 2048

if [ -d "${GENOME}" ] 
then
    echo ${GENOME}
else
    mkdir ${GENOME}
fi

STAR --genomeSAindexNbases 10 --runThreadN 4 --runMode genomeGenerate --genomeDir ${GENOME} \
--genomeFastaFiles ${GENOME_DIR}/*.masked --sjdbGTFfile ${ANN_DIR}/Saccharomyces_cerevisiae.R64-1-1.98_chr.gtf \
--sjdbOverhang 50

cd ${DATA_DIR}
for dir in $DIRS 
do
if [ -d "LID${dir}" ] 
then
    cd LID${dir}
    # ls  *.gz | xargs -I{} fastqc {}
    find ./ -name "*.fastq.gz" -exec gunzip {} \;
    cd ../
fi
done

for dir in $DIRS
do
if [ -d "LID${dir}" ] 
then
    cd LID${dir}
    FILES=(*_R1_001.fastq)
    for file in "${FILES[@]}"
    do
    echo base
    echo $file
    file=$(basename "$file")
    bfname="${file%R1_001.fastq}"
    echo $bfname
    echo $file
    if test -f "${bfname}Aligned.sortedByCoord.out.bam"; then
        continue
    fi
    cutadapt --cores=5 --minimum-length 5 --nextseq-trim=20 -a AGATCGGAAGAG -A AGATCGGAAGAG -o ${bfname}_trimmed_R1.fastq -p ${bfname}_trimmed_R2.fastq ${bfname}R1_001.fastq ${bfname}R2_001.fastq 
    ls ../${GENOME}
    STAR --outSAMtype BAM SortedByCoordinate \
    --quantMode GeneCounts \
    --limitBAMsortRAM 1145611175 \
    --twopassMode Basic \
    --twopass1readsN -1 \
    --runThreadN 5  --genomeDir ../${GENOME} --readFilesIn ${bfname}_trimmed_R1.fastq ${bfname}_trimmed_R2.fastq --outFileNamePrefix ${bfname} 
    # bzip2 ${bfname}_trimmed_R*.fastq
    done
fi
done
