# In which we make a makefile to combine HW5 and HW6

Before we make the makefile (lol), I had to decide which variables to make.

The basic variables are for the common tiger, panthera tigris.

```bash
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
```

The following commands were used.

Usage explains which commands the makefile contains.

```bash
usage:
	@echo "make genome    # Download the genome"
	@echo "make simulate  # Simulate reads for the genome"
	@echo "make download  # Download reads from SRA"
	@echo "make trim      # Trim reads"
	@echo "make fastqc    # Run fastqc on the reads"
```

The genome command downloads the genome and gets the smallest chromosome
```bash
genome:
	datasets download genome accession ${ACC}
	unzip -n ncbi_dataset.zip
	cat ${GENOME} | sed -n '/^>${CHR}/,/^>/p' | sed '$d' > ${CHR}.fa
```

The simulate command simulates sequencing for the previous command's output
```bash
simulate:
	wgsim -e 0 -r 0 -R 0 -1 100 -2 100 -N ${N} ${CHR}.fa read1.fq read2.fq
	echo $(seqkit stats read1.fq read2.fq)
```

The download command downloads reads from the SRA database and creates initial fastqc files
```bash
download:
	mkdir -p ${rdir} ${pdir}
	fastq-dump -X ${N} -F --outdir reads --split-files ${SRR}
	echo $(seqkit stats ${SRR}_1.fastq ${SRR}_2.fastq)
	fastqc -q ${r1} ${r2} -o ${pdir}
```

The trim command takes the previous command's output and trims from the right to create hopefully better fastqc files
```bash
trim:
	fastp --adapter_sequence=${ADAPTER} --cut_right \      -i ${r1} -I ${r2} -o ${t1} -O ${t2}
	fastqc -q -o ${pdir} ${t1} ${t2}
```

And I add a final line that makes sure the following commands are always executed
```bash
.PHONY: usage
```