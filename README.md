# SplitASM

SplitASM is a tool that corrects contigs using linked reads.
The connecting errors made by contig assemblers are detected by analyzing linked reads molecules profiles for each contig and splitting contigs where outlier values are identified.

## Dependancies
- gcc-7.2.0 or higher

## Installation

Clone this repository and run make in the folder

    git clone https://forgemia.inra.fr/seqoccin/axis1-assembly-methodo.git
    cd splitASM
    make

## Creating input files

1st step: Align linked reads to contigs to create a bam file (example with LongRanger on genotoul clster)


    module load bioinfo/bwa-0.7.17
    module load bioinfo/longranger-2.2.2

    longranger mkref contigFileName.fasta
    longranger align --id=outputFolderName --fastqs=fastqPATH --sample=sampleNAME --reference=refdata-contigFileName --localcores=20 --localmem=100

The bam file would be in the outputFolderName/out folder with the name possorted_bam.bam

2nd step: Create molecule tsv file

Using the modified LongRanger version from this repository:
    module load bioinfo/bwa-0.7.17
    module load bioinfo/samtools-1.10

    samtools sort -t BX -@ 20 -m 2G -o bcsorted_bam.bam possorted_bam.bam
    longranger-2.2.2_reportMolecules/longranger reportMolecules --id=outputFolderName \
     --fastqs=fastqPATH \
     --sample=sampleNAME \
     --reference=refdata-contigFileName \
     --vcmode=disable \
     --bam_bc=bcsorted_bam.bam \
     --bam_pos=possorted_bam.bam  \
     --localcores=20 --localmem=100

The molecule file would be in the outputFolderName/out folder with the name
fragments_tsv.tsv

3rd step: Create contig size files

    module load bioinfo/samtools-1.10

    samtools faidx contigs.fasta
    cut -f1,2 contigs.fasta.fai > contigs_size.tsv

4th step: Sort both molecules and contig size files with

    sort -k1,1 -s contigs_size.tsv > contigs_size_sort_by_name.tsv
    sort -k1,1 -s fragments_tsv.tsv > molecules_sort_by_name.tsv


## Usage
SplitASM requires two files containing the contigs size and the molecules information, both sorted by contig name. See **Creating input file** for how to obtain them.


    module load compiler/gcc-7.2.0

    splitASM -w 10000 -t 20  -c contigs_sort_by_name.tsv -o contigs_splited.bed molecules_sort_by_name.tsv

The output is a bed file contigs_splited.bed containing the new contig structure.

To obtain the fasta file you can use bedtools getfasta:

    module load bioinfo/bedtools2-2.29.0

    bedtools getfasta -fi contigs.fasta -bed contigs_splited.bed > contigs_splited.fasta


## SplitASM parameters


    Usage: ./splitASM <option(s)> MOLECULE FILE
    Options:
        -h,--help		Show this help message
        -t,--threads INT	Number of threads (default 1)
        -w, --window INT	Window size for outliers detection (default 10kb)
        -c, --contigs FILE	Contig size file name
        -o, --output FILE	Output bed file name
        -s, --sampleSize INT	Sample size for outlier detection, if zero or not specified then the sample is the input size
