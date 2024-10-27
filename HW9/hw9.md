# In which we use the .bam file from HW8 to look at statistics

First, I re-created a bam file that used a good quality SRA dataset and adjusted the makefile to align with both read files.

Then I generated some basic statistics so that I can cross-check against our queries for the questions.
```bash
samtools flagstat bam/p_align.bam
5941425 + 0 in total (QC-passed reads + QC-failed reads)
5795336 + 0 primary
0 + 0 secondary
146089 + 0 supplementary
0 + 0 duplicates
0 + 0 primary duplicates
2139359 + 0 mapped (36.01% : N/A)
1993270 + 0 primary mapped (34.39% : N/A)
5795336 + 0 paired in sequencing
2897668 + 0 read1
2897668 + 0 read2
1559448 + 0 properly paired (26.91% : N/A)
1798548 + 0 with itself and mate mapped
194722 + 0 singletons (3.36% : N/A)
0 + 0 with mate mapped to a different chr
0 + 0 with mate mapped to a different chr (mapQ>=5)
```

1. How many reads did not align with the reference genome?
```bash
samtools view -c -f 0x4 bam/p_align.bam
3802066
```
Unaligned: 5941425-2139359 = 3802066 reads (63.99%)

2-1. How many primary, secondary, and supplementary alignments are in the BAM file?
```bash
samtools view -c -F 0x900 bam/p_align.bam
5795336
```
Primary alignments: 5795336

2-2. 
```bash
samtools view -c -f 0x100 bam/p_align.bam
0
```
Secondary alignments: 0

2-3. 
```bash
samtools view -c -f 0x800 bam/p_align.bam
146089
```
Supplementary alignments: 146089

3. How many properly paired alignments on the reverse strand are formed by reads contained in the first pair (read1) file?
```bash
samtools view -c -f 66 bam/p_align.bam
797665
```
Properly paired alignments on the reverse strand formed by reads contained in the first pair (read1) file: 797665 

4. Making a new BAM file that contains only the properly paired primary alignments with a mapping quality of over 10.
```bash
samtools view -b -h -f 0x2 -F 0x900 -q 10 -o bam/q10_ppp_p_align.bam bam/p_align.bam
```

5. Compare the flagstats for your original and your filtered BAM file.
```bash
samtools flagstat bam/q10_ppp_p_align.bam 
1015967 + 0 in total (QC-passed reads + QC-failed reads)
1015967 + 0 primary
0 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
0 + 0 primary duplicates
1015967 + 0 mapped (100.00% : N/A)
1015967 + 0 primary mapped (100.00% : N/A)
1015967 + 0 paired in sequencing
508106 + 0 read1
507861 + 0 read2
1015967 + 0 properly paired (100.00% : N/A)
1015967 + 0 with itself and mate mapped
0 + 0 singletons (0.00% : N/A)
0 + 0 with mate mapped to a different chr
0 + 0 with mate mapped to a different chr (mapQ>=5)
```
If we compare this with the unfiltered flagstat results at the top, 1015967 reads made the quality cut. All of these reads are mapped (100%), compared to the 36.01% of the unfiltered results. And all of the filtered reads were able to be properly paired as well.