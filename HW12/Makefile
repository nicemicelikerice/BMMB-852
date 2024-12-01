# This makefile automates the following:
# dowloads a genome, simulates reads, downloads SRA files, trims and aligns them, and calls variants

# Accession number
ACC = GCF_000002035.6

# Genome file
GENOME = ncbi_dataset/data/${ACC}/${ACC}*genomic.fna

# GFF file
GFF = ncbi_dataset/data/${ACC}/genomic.gff

# GTF file
GTF = ncbi_dataset/data/${ACC}/genomic.gtf

# SRR number
SRR = SRR31379322

# Number of reads
N = 100000

# Output read names
R1 = reads/${SRR}_1.fastq

# Trimmed read names
T1 = reads/${SRR}_1.trimmed.fastq

# Adapter sequence
ADAPTER = AGATCGGAAGAGCACACGTCTGAACTCCAGTCA

# Read directory
RDIR = reads

# Report directory
PDIR = reports

# Reference genome
REF = ncbi_dataset/data/${ACC}/${ACC}*genomic.fna

# SAM file
SAM = ${SRR}_align.sam

# BAM file
BAM = ${SRR}_align.bam

# VCF file
VCF = variants/${SRR}_variants.vcf

# BCF flags for piling
PILE_FLAGS = -d 100 --annotate 'INFO/AD,FORMAT/DP,FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/SP'

# BCF flags for calling
CALL_FLAGS = --ploidy 2 --annotate 'FORMAT/GQ'

#Print the help
usage:
	@echo "make genome    # Download the genome"
	@echo "make download  # Download reads from SRA"
	@echo "make trim      # Trim reads"
	@echo "make index     # Index the reference genome"
	@echo "make align     # Align the reads to the reference genome"
	@echo "make variant   # Call variants on the bam file"
	@echo "make all       # Runs everything except downloading the genome and simulating reads"

#Download the genome and get smallest chromosome
genome:
	mkdir -p refs
	@echo "downloading genome"
	datasets download genome accession ${ACC} --include genome,gff3,gtf
	@echo "unzipping genome"
	unzip -n ncbi_dataset.zip
	cp ${GFF} refs/${ACC}.gff

#Download reads from SRA
download:
	mkdir -p ${RDIR} ${PDIR}
	@echo "downloading ${SRR} file"
	fastq-dump -X ${N} -F --outdir reads --split-files ${SRR}
	@echo "running fastqc for ${SRR}"
	fastqc -q ${R1} -o ${PDIR}

#QC Trimming
trim:
	@echo "trimming ${SRR}"
	fastp --adapter_sequence=${ADAPTER} --cut_right --cut_mean_quality=20 -i ${R1} -o ${T1}
	@echo "running fastqc for trimmed ${SRR}"
	fastqc -q -o ${PDIR} ${T1}

#Index the reference genome
index:
	@echo "indexing genome"
	bwa index ${REF}

#Align the reads to the reference genome
align:
	mkdir -p bam
	@echo "aligning ${SRR} to reference"
	bwa mem ${REF} ${T1} > bam/${SAM}
	cat bam/${SAM} | samtools sort > bam/${BAM}
	samtools index bam/${BAM}

#Call variants
variant:
	mkdir -p variants
	@echo "calling variants for ${SRR}"
	bcftools mpileup ${PILE_FLAGS} -O u -f ${GENOME} bam/${BAM} | \
		bcftools call ${CALL_FLAGS} -mv -O u | \
		bcftools norm -f ${GENOME} -d all -O u | \
		bcftools sort -O z > ${VCF}

all: download trim align variant

.PHONY: usage