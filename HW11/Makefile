# This makefile dowloads a genome, simulates reads, downloads SRA files, and trims them

# Accession number
ACC = GCF_000001735.3

# Genome file
GENOME = ncbi_dataset/data/${ACC}/${ACC}*genomic.fna

# GFF file
GFF = ncbi_dataset/data/${ACC}/genomic.gff

# GTF file
GTF = ncbi_dataset/data/${ACC}/genomic.gtf

# SRR number
SRR = SRR31238423

# Number of reads
N = 1000000

# Output read names
R1 = reads/${SRR}_1.fastq
R2 = reads/${SRR}_2.fastq

# Trimmed read names
T1 = reads/${SRR}_1.trimmed.fastq
T2 = reads/${SRR}_2.trimmed.fastq

# Adapter sequence
ADAPTER = AGATCGGAAGAGCACACGTCTGAACTCCAGTCA

# Read directory
RDIR = reads

# Report directory
PDIR = reports

# Reference genome
REF = ncbi_dataset/data/${ACC}/${ACC}*genomic.fna

# SAM file
SAM = align.sam

# BAM file
BAM = align.bam

# BCF flags for piling
PILE_FLAGS = -d 100 --annotate 'INFO/AD,FORMAT/DP,FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/SP'

# BCF flags for calling
CALL_FLAGS = --ploidy 2 --annotate 'FORMAT/GQ'

# VCF file
VCF = variants/variants.vcf

#Print the help
usage:
	@echo "make genome    # Download the genome"
	@echo "make simulate  # Simulate reads for the genome"
	@echo "make download  # Download reads from SRA"
	@echo "make trim      # Trim reads"
	@echo "make index     # Index the reference genome"
	@echo "make align     # Align the reads to the reference genome"
	@echo "make variant   # Call variants on the bam file"
	@echo "make gffsort   # Sorts, compresses, and indexes the GFF file"
	@echo "make vep       # Predicts variant effects"

#Download the genome and get smallest chromosome
genome:
	mkdir -p refs
	datasets download genome accession ${ACC} --include genome,gff3,gtf
	unzip -n ncbi_dataset.zip
	cp ${GFF} refs/${ACC}.gff

#Simulate sequencing
simulate:
	wgsim -e 0 -r 0 -R 0 -1 100 -2 100 -N ${N} ${REF} reads/sim_read1.fq reads/sim_read2.fq

#Download reads from SRA
download:
	mkdir -p ${RDIR} ${PDIR}
	fastq-dump -X ${N} -F --outdir reads --split-files ${SRR}
	echo $(seqkit stats ${SRR}_1.fastq ${SRR}_2.fastq)
	fastqc -q ${R1} ${R2} -o ${PDIR}

#QC Trimming
trim:
	fastp --adapter_sequence=${ADAPTER} --cut_tail \      -i ${R1} -I ${R2} -o ${T1} -O ${T2}
	fastqc -q -o ${PDIR} ${T1} ${T2}

#Index the reference genome
index:
	bwa index ${REF}

#Align the reads to the reference genome
align:
	mkdir -p bam
	bwa mem ${REF} ${T1} ${T2} > bam/${SAM}
	cat bam/${SAM} | samtools sort > bam/${BAM}
	samtools index bam/${BAM}

#Call variants
variant:
	mkdir -p variants
	bcftools mpileup ${PILE_FLAGS} -O u -f ${GENOME} bam/${BAM} | \
		bcftools call ${CALL_FLAGS} -mv -O u | \
		bcftools norm -f ${GENOME} -d all -O u | \
		bcftools sort -O z > ${VCF}

# VEP needs a sorted and compressed GFF file
gffsort:
	cat ${GFF} | sort -k1,1 -k4,4n -k5,5n -t$$'\t' | bgzip -c > ${GFF}.gz
	tabix -p gff ${GFF}.gz

# VEP is installed in the environment called vep
vep: ${GFF}.gz
	mkdir -p results
	micromamba run -n vep \
		~/src/ensembl-vep/vep \
		-i ${VCF} \
		-o results/vep.txt \
		--gff ${GFF}.gz \
		--fasta ${REF} \
		--species tetraodon_nigroviridis
		--force_overwrite 

	# Show the resulting files
	ls -lh results/*

all: genome download trim index align variant

.PHONY: usage