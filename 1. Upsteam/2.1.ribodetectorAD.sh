#!/bin/bash
#$ -pe smp 30
#$ -cwd
#$ -V
#$ -l "h_vmem=10g,h_rt=96:00:00"
#$ -o RNA_output.my
#$ -e RNA_error.my
#$ -q large.q


for f1 in *_R1.fastq
do
f2=${f1%%_1.fastq}"_R2.fastq"
echo $f1
echo $f2
ribodetector_cpu -t 32 \
  -l 101 \
  -i $f1 $f2 \
  -e rrna \
  --chunk_size 1000 \
  -o /share/data2/public/huihong.li/test1/AD/clean/$f1 /share/data2/public/huihong.li/test1/AD/clean/$f2
done