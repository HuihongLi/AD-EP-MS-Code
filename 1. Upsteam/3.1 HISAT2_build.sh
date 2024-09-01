#!/bin/bash
#$ -pe smp 8
#$ -cwd
#$ -V
#$ -l "h_vmem=8g,h_rt=96:00:00"
#$ -o hisat2_output.my
#$ -e hisat2_error.my
#$ -q large.q

hisat2-build Mus_musculus.GRCm39.dna_sm.primary_assembly.fa genome