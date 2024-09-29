#This script downloads a genome file, unzips it, gives some basic information about the file, and simulates some sequencing
#Make sure to run this script in the right environment and folder

##--------Change This Section to Update the Genome File--------##

#download, unzip, and link genome file
datasets download genome accession GCF_018350195.1
unzip ncbi_dataset.zip
ln -sf ncbi_dataset/data/GCF_018350195.1/GCF_018350195.1_P.tigris_Pti1_mat1.1_genomic.fna tiger.fa

##--------Nothing Changed Below This Line--------##

#get file size
#-L is necessary to dereference link
echo $(ls -lhL tiger.fa)

#get genome size
echo "Genome size: $(cat tiger.fa | grep -v '^>' | wc -m)"

#get number of chromosomes
echo "Number of chromosomes: $(cat tiger.fa | grep 'Pti1 chromosome' | grep -v 'unlocalized' | wc -l)"
echo "This wont work for other organisms!"

#get chromosome ids and lengths
echo "Chromosome IDs and lengths: $(cat tiger.fa | awk '$0 ~ ">" {if (NR > 1) {print c;} c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' | grep -v 'unlocalized' | grep 'Pti1 chromosome')"
echo "This wont work for other organisms!"

#get smallest chromosome
cat tiger.fa | sed -n '/^>NC_056673.1/,/^>/p' | sed '$d' > chromosome_1.fa

#simulate sequencing at 10x coverage
wgsim -e 0 -r 0 -R 0 -1 100 -2 100 -N 3000000 chromosome_1.fa read1.fq read2.fq
echo $(seqkit stats read1.fq read2.fq)

#compress the fastq files and compare sizes
echo "Uncompressed FASTQ file size: $(ls -lh read1.fq read2.fq)"
gzip read1.fq read2.fq
echo "Compressed FASTQ file size: $(ls -lh read1.fq.gz read2.fq.gz)"