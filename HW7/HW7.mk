# This makefile dowloads a genome, simulates reads, downloads SRA files, and trims them

# Accession number
ACC = GCF_018350195.1

# Genome file
GENOME = ncbi_dataset/data/${ACC}/${ACC}*genomic.fna

# Chromosome name
CHR = NC_056673.1

# SRR number
SRR = SRR639755

# Number of reads
N = 3000000

#Set output read names
r1=reads/${SRR}_1.fastq
r2=reads/${SRR}_2.fastq

#Set trimmed read names
t1=reads/${SRR}_1.trimmed.fastq
t2=reads/${SRR}_2.trimmed.fastq

#Set adapter sequence
ADAPTER=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA

#Set read directory
rdir=reads

#Set report directory
pdir=reports

#Print the help
usage:
	@echo "make genome    # Download the genome"
	@echo "make simulate  # Simulate reads for the genome"
	@echo "make download  # Download reads from SRA"
	@echo "make trim      # Trim reads"
	@echo "make fastqc    # Run fastqc on the reads"

#Download the genome and get smallest chromosome
genome:
	datasets download genome accession ${ACC}
	unzip -n ncbi_dataset.zip
	cat ${GENOME} | sed -n '/^>${CHR}/,/^>/p' | sed '$d' > ${CHR}.fa

#Simulate sequencing
simulate:
	wgsim -e 0 -r 0 -R 0 -1 100 -2 100 -N ${N} ${CHR}.fa read1.fq read2.fq
	echo $(seqkit stats read1.fq read2.fq)

#Download reads from SRA
download:
	mkdir -p ${rdir} ${pdir}
	fastq-dump -X ${N} -F --outdir reads --split-files ${SRR}
	echo $(seqkit stats ${SRR}_1.fastq ${SRR}_2.fastq)
	fastqc -q ${r1} ${r2} -o ${pdir}

#QC Trimming
trim:
	fastp --adapter_sequence=${ADAPTER} --cut_tail \      -i ${r1} -I ${r2} -o ${t1} -O ${t2}
	fastqc -q -o ${pdir} ${t1} ${t2}

.PHONY: usage