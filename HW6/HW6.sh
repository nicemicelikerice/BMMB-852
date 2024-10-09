#This script downloads a sra file and performs quality control
#Make sure to run this script in the right environment and folder

#Set trace
set -uex

##--------Change This Section to Update the File--------##

#Set SRR
srr="SRR639755"

#Set number of reads
N=10000

#Set output read names
r1=reads/${srr}_1.fastq
r2=reads/${srr}_2.fastq

#Set trimmed read names
t1=reads/${srr}_1.trimmed.fastq
t2=reads/${srr}_2.trimmed.fastq

#Set adapter sequence
ADAPTER=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA

#Set read directory
rdir=reads

#Set report directory
pdir=reports

##--------Nothing Changed Below This Line--------##

#Make directories
mkdir -p ${rdir} ${pdir}

#Download data from SRA
fastq-dump -X ${N} -F --outdir reads --split-files ${srr}

#Run FastQC
fastqc -q ${r1} ${r2} -o ${pdir}

#QC Trimming
fastp --adapter_sequence=${ADAPTER} --cut_tail \
      -i ${r1} -I ${r2} -o ${t1} -O ${t2}

#Run FastQC on trimmed reads
fastqc -q -o ${pdir} ${t1} ${t2}

#Assuming that the environment menv has been created
#Create multiqc report
micromamba run -n menv multiqc -o ${pdir} ${pdir}