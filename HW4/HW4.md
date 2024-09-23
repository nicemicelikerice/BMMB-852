# This markdown file shows how I created a script based off of HW2

**Create a .sh file**

Make sure to move to the folder you want to before all of this and set up the right environment!

```bash
touch HW4.sh
```
Now edit the contents with a code editor like VSC

**Make a variable that will store the link to the file you want and the exact data of the organism**

```bash
gff3="https://ftp.ensembl.org/pub/current_gff3/mus_musculus/Mus_musculus.GRCm39.112.gff3.gz"
org="Mus_musculus.GRCm39.112"
```

Cut off a section of your script so that you don't have to change anything below that line

**Download that file and unzip it**

```bash
wget $gff3
gunzip *.gz
```

We should also create a clean file to work out of

```bash
cat $org.gff3 | grep -v '#' > clean.gff3
```

This should give a file named clean.gff

**Make your script say the things you want it to**

1. The script should get the number of features
2. The script should get the number of sequence regions
3. The script should get the number of genes listed for this organism
4. The script should list the top ten most annotated feature types for this organism
5. The script should be snarky and question the decision-making capability of the user

```bash
echo "Number of features in the file: $(cat clean.gff3 | cut -f 3 | wc -l)"
echo "Number of sequence regions in the file: $(cat clean.gff3 | cut -f 1 | sort | uniq | wc -l)"
echo "Number of genes in the file: $(cat clean.gff3 | cut -f 3 | grep -w 'gene' | wc -l)"
echo -e "Top ten most annotated feature types in the file:\n$(cat clean.gff3 | cut -f 3 | sort-uniq-count-rank | head)"
echo "Does it seem like a complete and well-annotated organism? Well, I can't do everything for you! Use your brain!"
```

The echo -e is because I wanted to change lines just for the aesthetics

Unfortunately I had to add in an extra line because the mus musculus file contains a bunch of fragmented chromosomes that get counted by the above code

```bash
echo "For the mus musculus genome, the number of sequence regions in the file is different due to the chromosome names (this counts X and Y as separate chromosomes): $(cat clean.gff3 | cut -f 1 | sort | uniq | grep -E '^[0-9]+$|^[xXyY]$' | wc -l)"
```

But now your script is done! Save it and run it in the folder of choice

**Practicing on other people's data**

This is on Marti Albuja's data

When I ran my data on their data and changed the variables to 

```bash
gff3="https://ftp.ensembl.org/pub/current_gff3/prolemur_simus/Prolemur_simus.Prosim_1.0.112.gff3.gz"
org="Prolemur_simus.Prosim_1.0.112"
```

I was able to replicate her results!

```bash
Number of features in the file: 1104650
Number of sequence regions in the file: 128596
For the mus musculus genome, the number of sequence regions in the file is different due to the chromosome names (this counts X and Y as separate chromosomes): 0
Number of genes in the file: 20354
Top ten most annotated feature types in the file:
400695  exon
393577  CDS
128596  region
111882  biological_region
37979   mRNA
20354   gene
5290    ncRNA_gene
1902    transcript
1467    snRNA
924     rRNA
Does it seem like a complete and well-annotated organism? Well, I can't do everything for you! Use your brain!
```