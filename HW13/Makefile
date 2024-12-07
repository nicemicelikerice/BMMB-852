# This makefile automates the following:
# indexes a genome, trims + aligns + analyzes rna reads, and creates a table of the read counts

ACC=GCF_000001635.27
FULLNAME=GCF_000001635.27_GRCm39_genomic
GENOME = ncbi_dataset/data/${ACC}/${ACC}*genomic.fna
GFF = ncbi_dataset/data/${ACC}/genomic.gff
GTF = ncbi_dataset/data/${ACC}/genomic.gtf
SRR=SRR17145315
N=100

downloadgenome:
	mkdir Genome
	cd Genome
	datasets download genome accession ${ACC} --include genome,gff3,gtf
	@echo "unzipping genome"
	unzip -n ncbi_dataset.zip
	mv ncbi_dataset/data/${GFF}/* .
	rm -r ncbi_dataset

downloadsrr:
	@echo "downloading srr files"
	fastq-dump -X ${N} -F --outdir 00-rawreads --split-files ${SRR}

index:
	@echo "indexing genome"
	sh star_genome.sh

trim: 
	@echo "trimmming reads"
	sh trim_for_loop.sh

align: 
	@echo "aligning reads"
	sh star_for_loop.sh

htseq:
	@echo "analyzing reads"
	sh htseq_for_loop.sh

counttable:
	@echo "creating read count table"
	micromamba run -n stats Rscript deseq.R
.PHONY: usage

