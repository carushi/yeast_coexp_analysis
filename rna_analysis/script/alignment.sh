#!/bin/bash
set -exuo pipefail 

DATA_DIR=../data/
GENOME=../../genome/chromFaMasked/
GENOME_DIR=${GENOME}

ANN_DIR=../../ann/


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



for dir in LID304275 LID304917
do
cd ${dir}
find R*R1_*.fastq.gz | while read file
do
bfname="${file%R1_001.fastq.gz}"
echo $bfname
echo $file
if test -f "${bfname}_trimmed_R1.fastq.bz2"; then
	continue
fi
cutadapt --cores=6 --minimum-length 5 --nextseq-trim=20 -a AGATCGGAAGAG -A AGATCGGAAGAG -o ${bfname}_trimmed_R1.fastq -p ${bfname}_trimmed_R2.fastq ${bfname}R1_001.fastq.gz ${bfname}R2_001.fastq.gz 
bunzip2 ${bfname}_trimmed_R*fastq.bz2
STAR --outSAMtype BAM SortedByCoordinate \
--quantMode GeneCounts \
--twopassMode Basic \
--twopass1readsN -1 \
--outTmpDir /home/rkawaguc/ipython/ChIP/startmp --runThreadN 5 --genomeDir data/SGD/star --readFilesIn ${bfname}_trimmed_R1.fastq ${bfname}_trimmed_R2.fastq --outFileNamePrefix ${bfname}
bzip2 ${bfname}_trimmed_R*.fastq
done
cd ../
done

#find $TOP_DIR/ -name "*.bam" | xargs -I{} echo "samtools index {}" | bash
#find $TOP_DIR/ -name "*.sam" -exec rm {} \;

