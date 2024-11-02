# This makefile dowloads a genome, simulates reads, downloads SRA files, and trims them

# Accession number
ACC = GCF_018350195.1

# Genome file
GENOME = ncbi_dataset/data/${ACC}/${ACC}*genomic.fna

# Chromosome name
CHR = NC_056673.1

# SRR number
SRR = SRR25597804

# Number of reads
N = 3000000

#Set output read names
R1=reads/${SRR}_1.fastq
R2=reads/${SRR}_2.fastq

#Set trimmed read names
T1=reads/${SRR}_1.trimmed.fastq
T2=reads/${SRR}_2.trimmed.fastq

#Set adapter sequence
ADAPTER=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA

#Set read directory
RDIR=reads

#Set report directory
PDIR=reports

#Set reference genome
REF=refs/${CHR}.fa

#Set SAM file name
SAM=align.sam

#Set BAM file name
BAM=align.bam

#Print the help
usage:
	@echo "make genome    # Download the genome"
	@echo "make simulate  # Simulate reads for the genome"
	@echo "make download  # Download reads from SRA"
	@echo "make trim      # Trim reads"
	@echo "make index     # Index the reference genome"
	@echo "make align     # Align the reads to the reference genome"
	@echo "make variant   # Call variants on the bam file"

#Download the genome and get smallest chromosome
genome:
	mkdir -p refs
	datasets download genome accession ${ACC}
	unzip -n ncbi_dataset.zip
	cat ${GENOME} | sed -n '/^>${CHR}/,/^>/p' | sed '$d' > ${REF}

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
	bcftools mpileup -f ${GENOME} bam/${BAM} | bcftools call -mv -Oz -o variants/variants.vcf.gz
	bcftools filter -s LOWQUAL -e 'QUAL<20 || DP<10' variants/variants.vcf.gz > variants/filtered_variants.vcf


.PHONY: usage