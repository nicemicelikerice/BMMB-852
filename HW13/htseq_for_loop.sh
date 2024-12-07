#!/bin/bash

conda activate rnaseq_workshop

cd ~/jk/work/hw13/02-alignments

GTF="~/jk/work/hw13/Genome/GCF_000001635.27_GRCm39_genomic.gtf"

for sample in $(cat samples.txt); do

htseq-count -f bam -r pos -s reverse -n 4 \
${sample}/Aligned.sortedByCoord.out.bam ${GTF} \
 > ../03-readcounts/${sample}_batch.txt

done