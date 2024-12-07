#!/bin/bash

conda activate rnaseq_workshop

cd ~/jk/work/hw13/01-trimmedreads

# "sample.txt" has sample names
# "sample" directories should be premade ("mkdir" or "cat samples.txt | xargs mkdir")

for sample in $(cat samples.txt); do

STAR --runThreadN 4 \
--genomeDir ~/jk/work/hw13/Genome/NCBI_STAR_index \
--sjdbGTFfile ~/jk/work/hw13/Genome/GCF_000001635.27_GRCm39_genomic.gtf \
--readFilesIn ${sample}_trimmed.fq.gz \
--readFilesCommand zcat \
--outSAMtype BAM SortedByCoordinate \
--outFileNamePrefix ../02-alignments/${sample}/

done