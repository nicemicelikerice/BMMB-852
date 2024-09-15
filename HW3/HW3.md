#   Downloading genome and files
(bioinfo)
jk@JiwooHPbar ~/work/hw3
$ datasets download genome accession GCF_018350195.1 --include gff3,cds,protein,rna,genome
New version of client (16.28.0) available at https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/LATEST/linux-amd64/datasets.
Collecting 1 genome record [================================================] 100% 1/1
Downloading: ncbi_dataset.zip    848MB valid zip structure -- files not checked
Validating package [================================================] 100% 9/9
(bioinfo)
jk@JiwooHPbar ~/work/hw3
$ ls
README.md  md5sum.txt  ncbi_dataset  ncbi_dataset.zip
#   Unzipping
(bioinfo)
jk@JiwooHPbar ~/work/hw3
$ unzip ncbi_dataset.zip
Archive:  ncbi_dataset.zip
replace README.md? [y]es, [n]o, [A]ll, [N]one, [r]ename: y
  inflating: README.md
replace ncbi_dataset/data/assembly_data_report.jsonl? [y]es, [n]o, [A]ll, [N]one, [r]ename: A
  inflating: ncbi_dataset/data/assembly_data_report.jsonl
  inflating: ncbi_dataset/data/GCF_018350195.1/GCF_018350195.1_P.tigris_Pti1_mat1.1_genomic.fna
  inflating: ncbi_dataset/data/GCF_018350195.1/genomic.gff
  inflating: ncbi_dataset/data/GCF_018350195.1/cds_from_genomic.fna
  inflating: ncbi_dataset/data/GCF_018350195.1/protein.faa
  inflating: ncbi_dataset/data/GCF_018350195.1/rna.fna
  inflating: ncbi_dataset/data/dataset_catalog.json
  inflating: md5sum.txt
#   Attempting to open genomic.gff in IGV forced to sort, creating genomic.sorted.gff
#   Extracting genes
(bioinfo)
jk@JiwooHPbar ~/work/hw3
$ cat ncbi_dataset/data/GCF_018350195.1/genomic.sorted.gff | awk ' $3=="gene" { print $0 }' > ncbi_dataset/data/GCF_018350195.1/gene.gff
#   Extracting CDS
(bioinfo)
jk@JiwooHPbar ~/work/hw3
$ cat ncbi_dataset/data/GCF_018350195.1/genomic.sorted.gff | awk ' $3=="CDS" { print $0 }' > ncbi_datase
t/data/GCF_018350195.1/cds.gff