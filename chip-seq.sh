#!/bin/bash

basedir=$1

if [ ! ${basedir} ]; then
	basedir="."
fi

if [ ! -d ${basedir} ]; then
	echo "This is not a directory!!"
fi

echo "     Using basedir: ${basedir}/"

current_dir=`pwd`

#########################lists##########################################################################

list0=(
"C-3E"
"INP-C-3E"
"INP-SIHJ-2E"
"INP-SIHJ-3E"
"SIHJ-2E"
"SIHJ-3E"
)

list1=(
"L001"
"L002"
)

num=( 0 
1 
2
3
4
5
)

########################################################################################################################################

echo "################################    trimming      ################################"

echo "################################      and         ################################"

echo "########################    fastq quality control    #########################"


mkdir -p ${basedir}/fastqc/raw
mkdir -p ${basedir}/fastqc/post-trim
mkdir -p ${basedir}/trimmed_data

for smp in ${list0[@]}
do
	for l in ${list1[@]}; do
		fastqc -t 12 -o fastqc/raw/ raw_data/${smp}_*_${l}_R1_*.fastq.gz
		fastqc -t 12 -o fastqc/raw/ raw_data/${smp}_*_${l}_R2_*.fastq.gz
		java -jar ~/programs/trimmomatic.jar PE -threads 12 -phred33 raw_data/${smp}_*_${l}_R1_*.fastq.gz raw_data/${smp}_*_${l}_R2_*.fastq.gz trimmed_data/${smp}_${l}_1.fastq trimmed_data/${smp}_${l}_1_unpaired.fastq.gz trimmed_data/${smp}_${l}_2.fastq trimmed_data/${smp}_${l}_2_unpaired.fastq.gz ILLUMINACLIP:${current_dir}/../../ref_data/adapters/truseq_all.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:3:20 MINLEN:36
		fastqc -t 12 -o fastqc/post-trim/ trimmed_data/${smp}_${l}_1.fastq
		fastqc -t 12 -o fastqc/post-trim/ trimmed_data/${smp}_${l}_2.fastq
	done
done

echo "################################     Aligning      ################################"

echo "################################     mapping        ################################"

echo "##############################    and assembling      ################################"

for sp in ${list0[@]}
do
	mkdir -p ${basedir}/bowtie2/${sp}
	bowtie2 -p 12 -q -t -x ${current_dir}/../../ref_data/GRCh38/bt-index/bt -1 trimmed_data/${sp}_L001_1.fastq,trimmed_data/${sp}_L002_1.fastq -2 trimmed_data/${sp}_L001_2.fastq,trimmed_data/${sp}_L002_2.fastq -S bowtie2/${sp}/${sp}.sam
	cd bowtie2/${sp}
	samtools view -bS -T ${current_dir}/../../ref_data/GRCh38/GCF_000001405.39_GRCh38.p13_genomic.fna -o ${sp}.bam ${sp}.sam
	samtools sort -@ 10 ${sp}.bam -o ${sp}.sorted.bam
done

for nm in ${num[@]}; do
	mkdir -p ${basedir}/bedtools/${list2[${nm}]}
	mkdir -p ${basedir}/bedtools/${list3[${nm}]}
	bedtools bamtobed -i bowtie2/sorted-chip/${list2[${nm}]}.sorted.bam > bedtools/${list2[${nm}]}/${list2[${nm}]}.bed
	bedtools bamtobed -i bowtie2/sorted-control/${list3[${nm}]}.sorted.bam > bedtools/${list3[${nm}]}/${list3[${nm}]}.bed
done

mkdir -p ${basedir}/diffReps
cd diffReps

diffReps.pl --treatment ../bedtools/SIHJ-2E_L001/SIHJ-2E_L001.bed ../bedtools/SIHJ-2E_L002/SIHJ-2E_L002.bed ../bedtools/SIHJ-3E_L001/SIHJ-3E_L001.bed ../bedtools/SIHJ-3E_L002/SIHJ-3E_L002.bed --control ../bedtools/C-3E_L001/C-3E_L001.bed ../bedtools/C-3E_L002/C-3E_L002.bed --btr ../bedtools/INP-SIHJ-2E_L001/INP-SIHJ-2E_L001.bed ../bedtools/INP-SIHJ-2E_L002/INP-SIHJ-2E_L002.bed ../bedtools/INP-SIHJ-3E_L001/INP-SIHJ-3E_L001.bed ../bedtools/INP-SIHJ-3E_L002/INP-SIHJ-3E_L002.bed --bco ../bedtools/INP-C-3E_L001/INP-C-3E_L001.bed ../bedtools/INP-C-3E_L002/INP-C-3E_L002.bed --nproc 12 --mode p --meth nb --report diff.w.o.C2.nb.txt --chrlen ../../../ref_data/GRCh38/GRCh38.len.txt

cd ${current_dir}

exit
