#!/bin/bash
#$ -pe smp 4
#$ -cwd
#$ -V
#$ -l "h_vmem=4g,h=compute-0-3,h_rt=96:00:00"
#$ -o qc_output.my
#$ -e qc_error.my
#$ -q large.q

multiqc -f qc
