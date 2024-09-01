#!/bin/bash
#$ -pe smp 32
#$ -cwd
#$ -V
#$ -l "h_vmem=8g,h_rt=96:00:00"
#$ -o RNA_output2.my
#$ -e RNA_error2.my
#$ -q large.q

for f1 in *_R1.fastq
do
f2=${f1%%_R1.fastq}"_R2.fastq"
if [ ! -f "clean/$f1" ]; then
echo $f1
echo $f1
ribodetector_cpu -t 64 \
  -l 58 \
  -i $f1 $f2 \
  -e rrna \
  --chunk_size 2000 \
  -o /share/data2/public/huihong.li/test1/EP/clean/$f1 /share/data2/public/huihong.li/test1/EP/clean/$f2
fi
done