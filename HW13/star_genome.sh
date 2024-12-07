#!/bin/bash

conda activate rnaseq_workshop

cd ~/jk/work/hw13/Genome

mkdir NCBI_STAR_index

STAR --runMode genomeGenerate \
--genomeFastaFiles GCF_000001635.27_GRCm39_genomic.fna \
--sjdbGTFfile GCF_000001635.27_GRCm39_genomic.gtf \
--runThreadN 4 \
--genomeSAindexNbases 13 \
--genomeDir NCBI_STAR_index
