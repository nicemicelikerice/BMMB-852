# This makefile automates the following:
# RNAseq starting from scratch

ACC = GCF_000001635.27
GENOME = ${ACC}*genomic.fna
GFF = genomic.gff
GTF = genomic.gtf
SRR = SRR17145315
N = 100000
Samples = samples.txt

usage:
	@echo "downloadgenome: downloads a genome"
	@echo "downloadsrr:    downloads a srr file. requires a shell script for multiple downloads"
	@echo "index:          indexes the genome file"
	@echo "trim:           trims raw reads"
	@echo "align:          aligns trimmed reads against genome"
	@echo "htseq:          analyzes reads"
	@echo "deseq2:         runs deseq2 for DEG analysis"

downloadgenome:
	@echo "downloading genome"
	mkdir -p Genome
	cd Genome && \
	datasets download genome accession ${ACC} --include genome,gff3,gtf && \
	unzip -n ncbi_dataset.zip && \
	mv ncbi_dataset/data/${ACC}/* . && \
	rm -r ncbi_dataset

downloadsrr:
	@echo "downloading srr files"
	mkdir -p 00-rawreads
	fastq-dump -X ${N} -F --outdir 00-rawreads --split-files ${SRR}

index:
	@echo "indexing genome"
	cd Genome && \
	mkdir -p NCBI_STAR_index && \
	STAR --runMode genomeGenerate \
	--genomeFastaFiles ${GENOME} \
	--sjdbGTFfile ${GTF} \
	--runThreadN 4 \
	--genomeSAindexNbases 13 \
	--genomeDir NCBI_STAR_index

trim: 
	@echo "trimmming reads"
	mkdir -p 01-trimmedreads
	cd 00-rawreads && \
	for sample in $(cat ../${Samples}); do
		trim_galore -o ../01-trimmedreads -j 4 --length 70 --gzip ${sample}.fastq.gz
	done

align: 
	@echo "aligning reads"
	mkdir -p 02-alignments
	cd 01-trimmedreads && \
	for sample in $(cat ../${Samples}); do	
		STAR --runThreadN 4 \
		--genomeDir ../Genome/NCBI_STAR_index \
		--sjdbGTFfile ../Genome${GTF} \
		--readFilesIn ${sample}_trimmed.fq.gz \
		--readFilesCommand zcat \
		--outSAMtype BAM SortedByCoordinate \
		--outFileNamePrefix ../02-alignments/${sample}/	
	done

htseq:
	@echo "analyzing reads"
	mkdir -p 03-readcounts
	cd 02-alignments && \
	for sample in $(cat ../${Samples}); do
		htseq-count -f bam -r pos -s reverse -n 4 \
		${sample}/Aligned.sortedByCoord.out.bam ../Genome/${GTF} \
		 > ../03-readcounts/${sample}_batch.txt
	done

deseq2:
	@echo "creating count matrix and running deseq2 on R. generates PCA plot, MA plot, heatmap"
	mkdir results
	mkdir plots
	micromamba run -n stats Rscript deseq.R

.PHONY: usage downloadgenome downloadsrr index trim align htseq deseq2