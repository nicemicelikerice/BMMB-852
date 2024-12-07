#!/bin/bash

conda activate rnaseq_workshop

cd ~/jk/work/hw13/00-rawreads

for sample in $(cat samples.txt); do
	trim_galore -o ~/jk/work/hw13/01-trimmedreads -j 4 --length 70 --gzip ${sample}.fastq.gz
done