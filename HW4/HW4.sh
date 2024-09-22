#This script automates HW2
#It downloads a gff3 file, unzips it, cleans it, and then provides some basic statistics about the file
#It also makes fun of the user for not being able to make decisions for themselves
#Make sure to run this script in the right environment and folder

##--------Change This Section to Update the GFF3 File--------##

#link to the gff3 file and data name
gff3="https://ftp.ensembl.org/pub/current_gff3/mus_musculus/Mus_musculus.GRCm39.112.gff3.gz"
org="Mus_musculus.GRCm39.112"

##--------Nothing Changed Below This Line--------##

#download gff3 file
wget $gff3

#unzip the file
gunzip *.gz

#clean up gff3 file
cat $org.gff3 | grep -v '#' > clean.gff3

#get number of features from clean file
echo "Number of features in the file: $(cat clean.gff3 | cut -f 3 | wc -l)"

#get number of sequence regions from clean file (may not work for other names)
echo "Number of sequence regions in the file: $(cat clean.gff3 | cut -f 1 | sort | uniq | wc -l)"

#backup for previous answer
echo "For the mus musculus genome, the number of sequence regions in the file is different due to the chromosome names (this counts X and Y as separate chromosomes): $(cat clean.gff3 | cut -f 1 | sort | uniq | grep -E '^[0-9]+$|^[xXyY]$' | wc -l)"

#number of genes listed for this organism
echo "Number of genes in the file: $(cat clean.gff3 | cut -f 3 | grep -w 'gene' | wc -l)"

#top ten most annotated feature types for this organism
echo -e "Top ten most annotated feature types in the file:\n$(cat clean.gff3 | cut -f 3 | sort-uniq-count-rank | head)"

#question decision-making capability of user
echo "Does it seem like a complete and well-annotated organism? Well, I can't do everything for you! Use your brain!"