# Simulating FASTQ files

**Part 1: The Basics**
Download, unzip, and link file
```bash
datasets download genome accession GCF_018350195.1
unzip ncbi_dataset.zip
ln -sf ncbi_dataset/data/GCF_018350195.1/GCF_018350195.1_P.tigris_Pti1_mat1.1_genomic.fna tiger.fa
```

File size: 2.3 Gb
-L is necessary to dereference link
```bash
ls -lhL tiger.fa
```

Genome size: 2.4 billion
```bash
cat tiger.fa | grep -v '^>' | wc -m
```

Number of chromosomes: 19
*this wont work for other organisms!
```bash
cat tiger.fa | grep 'Pti1 chromosome' | grep -v 'unlocalized' | wc -l
```

ID and length 
This command grabs all the chromosomes, lists the number of basepairs, and filters out unlocalized or mitochondrial sections
```bash
cat tiger.fa | awk '$0 ~ ">" {if (NR > 1) {print c;} c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' | grep -v 'unlocalized' | grep 'Pti1 chromosome' 
```

**Part 2: Electric Boogaloo**

Extracting the smallest chromosome to save my sanity and pc
*insert picture of Ash Ketchum*
"I choose you, NC_056673.1 (E1)!"
```bash
cat tiger.fa | sed -n '/^>NC_056673.1/,/^>/p' | sed '$d' > chromosome_1.fa
```

check it has the same number of characters as it should...
```bash
seqkit stats chromosome_1.fa
```

E1 has about 61 million bp length so a 10x coverage with 100 bp read length means I want my number of reads to be 3 million
*Praying that my pc doesnt melt*

using wgsim for 10x coverage and checking the FASTQ files
```bash
wgsim -e 0 -r 0 -R 0 -1 100 -2 100 -N 3000000 chromosome_1.fa read1.fq read2.fq
seqkit stats read1.fq read2.fq
```
The results:
```bash
file      format  type   num_seqs      sum_len  min_len  avg_len  max_len
read1.fq  FASTQ   DNA   3,000,000  300,000,000      100      100      100
read2.fq  FASTQ   DNA   3,000,000  300,000,000      100      100      100
```

Number of reads: 3 million
Average read length: 100 bp
FASTQ file size: 731 Mb
```bash
ls -lH read1.fq read2.fq
```
Compressed FASTQ file size: 138 Mb
```bash
gzip read1.fq read2.fq
ls -lH read1.fq.gz read2.fq.gz
```

Yes, I could get the same coverage by manipulating the read length to be longer and the read number to be less but to my knowledge, technical limitations make long reads hard to do.

**Part 3: The End**

Assumptions for the following estimations
    Yeast genome size: 12 million bp
    Drosophila genome size: 180 million bp
    Human genome size: 3 billion bp

As the tiger genome is 2.4 billion bp and has a FASTA file size of 2.3Gb, I can estimate the following

Size of the FASTA file that holds the genome
    Yeast: 11.7 Mb
    Drosophila: 176.6 Mb
    Human: 2.9 Gb

If the read length parameters are equal to my data,

Number of FASTQ reads needed for 30x
    Yeast: 1,800,000
    Drosophila: 27,000,000
    Human: 450,000,000

My chromosome 1 is 60 Mb, read1.fq and read2.fq were 10x coverage and 731 Mb
So, if they were 30x they should be 3 times bigger at 2.1 Gb
That means that the read files would be around 36x the original data
And as my compression rate was around 81%,

Size of the FASTQ files before and after compression
    Yeast: 421.2 Mb => 80 Mb
    Drosophila: 6.2 Gb => 1.2 Gb
    Human: 104.4 Gb => 19.8 Gb