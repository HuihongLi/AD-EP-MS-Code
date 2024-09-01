#!/bin/bash
#$ -pe smp 32
#$ -cwd
#$ -V
#$ -l "h_vmem=8g,h_rt=96:00:00"
#$ -o count2_output.my
#$ -e count2_error.my
#$ -q large.q


featureCounts -T 32 \
-p \
-t  exon -g gene_id \
-a /share/data2/public/huihong.li/Star/HISAT2/Mus_musculus.GRCm39.111.gtf \
-o counts2unpair.txt \
*.bam \
