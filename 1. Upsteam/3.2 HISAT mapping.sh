#!/bin/bash
#$ -pe smp 32
#$ -cwd
#$ -V
#$ -l "h_vmem=8g,h_rt=96:00:00"
#$ -o HISAT_output.my
#$ -e HISAT_error.my
#$ -q large.q

for f1 in *_R1.fastq
do
f2=${f1%%_R1.fastq}"_R2.fastq"
f3=${f1%%_R1.fastq}".sam"
f4=${f1%%_R1.fastq}".bam"
echo $f1
echo $f2
echo $f3
echo $f4
hisat2 -p 32 -x /share/data2/public/huihong.li/Star/HISAT2/genome -1 $f1 -2 $f2 -S /share/data2/public/huihong.li/test1/EP/clean/HISAT/$f3
samtools view -@ 32 -S HISAT/$f3 -1b -o HISAT/$f4
rm -rf /share/data2/public/huihong.li/test1/EP/clean/HISAT/$f3
done