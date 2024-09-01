#!/bin/bash
#$ -pe smp 16
#$ -cwd
#$ -V
#$ -l "h_vmem=4g,h_rt=96:00:00"
#$ -o fastqc_output.my
#$ -e fastqc_error.my
#$ -q large.q

mkdir qc
fastqc -o qc -t 16  *.fastq.gz