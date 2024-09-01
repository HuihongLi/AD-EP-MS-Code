#!/bin/bash
#$ -pe smp 16
#$ -cwd
#$ -V
#$ -l "h_vmem=4g,h_rt=96:00:00"
#$ -o output.my
#$ -e trimerror.my
#$ -q long.q

for f1 in *_1.fastq.gz
do
f2=${f1%%_1.fastq.gz}"_2.fastq.gz"
f3=${f1%%_1.fastq.gz}"_R1.fastq"
f4=${f1%%_1.fastq.gz}"_R2.fastq"
echo $f1
echo $f2
echo $f3
echo $f4
fastp -i $f1 -I $f2 -o $f3 -O $f4 --thread 16 --qualified_quality_phred 30 --n_base_limit 5 --detect_adapter_for_pe --cut_tail --trim_front1 15 --dedup --trim_poly_x --trim_poly_g
done