https://github.com/nicemicelikerice/BMMB-852/blob/main

Tell us a bit about the organism.
A: Mus musculus is more commonly known as the house mouse and is perhaps the most widely used animal model in research. It is a small rodent that shares many genes with humans, and breeds extremely fast.

How many features does the file contain?
$ conda activate bioinfo
$ mkdir hw2
$ cd hw2
$ wget https://ftp.ensembl.org/pub/current_gff3/mus_musculus/Mus_musculus.GRCm39.112.gff3.gz
$ gunzip Mus_musculus.GRCm39.112.gff3.gz
$ cat Mus_musculus.GRCm39.112.gff3 | grep -v '#' > clean.gff3
$ cat clean.gff3 | cut -f 3 | wc -l
A: 1850382

How many sequence regions (chromosomes) does the file contain? 
$ cat clean.gff3 | cut -f 1 | sort | uniq
1
10
11
12
13
14
15
16
17
18
19
2
3
4
5
6
7
8
9
GL456210.1
GL456211.1
GL456212.1
GL456219.1
GL456221.1
GL456233.2
GL456239.1
GL456354.1
GL456359.1
GL456360.1
GL456366.1
GL456367.1
GL456368.1
GL456370.1
GL456372.1
GL456378.1
GL456379.1
GL456381.1
GL456382.1
GL456383.1
GL456385.1
GL456387.1
GL456389.1
GL456390.1
GL456392.1
GL456394.1
GL456396.1
JH584295.1
JH584296.1
JH584297.1
JH584298.1
JH584299.1
JH584300.1
JH584301.1
JH584302.1
JH584303.1
JH584304.1
MT
MU069434.1
MU069435.1
X
Y
A: 20 (19 + 1 sex chromosome)

How many genes are listed for this organism?
$ cat clean.gff3 | cut -f 3 | grep -w 'gene' | wc -l
A: 25425

What are the top-ten most annotated feature types (column 3) across the genome?
$ cat clean.gff3 | cut -f 3 | sort-uniq-count-rank | head
864295  exon
522293  CDS
95372   five_prime_UTR
86987   three_prime_UTR
74994   biological_region
65883   mRNA
45359   lnc_RNA
25425   gene
18004   ncRNA_gene
14744   transcript
A: exons, CDS, 5'-UTRs, 3'-UTRs, biological regions, mRNAs, lncRNAs, genes, ncRNA genes, and transcripts

Having analyzed this GFF file, does it seem like a complete and well-annotated organism?
Yes, it has a large amount of annotated features. Which it should, as it is the common lab mouse.