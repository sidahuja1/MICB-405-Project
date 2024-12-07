#!/bin/bash

# Make sure you put in the stuff
if [ "$#" -ne 5 ]; then
    echo "Usage: $0 <reference_dir> <gtf_file> <read1.fastq> <read2.fastq> <output_prefix>"
    exit 1
fi

# Input arguments in order, so call this script and the paths for the following and a prefix for your output files
REFERENCE=$1      # Reference genome directory
GTF=$2            # GTF annotation file
READ1=$3          # FASTQ file for read 1
READ2=$4          # FASTQ file for read 2
OUT_PREFIX=$5     # Prefix for output files

# Check if the reference directory contains STAR index files
if [ ! -f "${REFERENCE}/SA" ]; then # Looks for SA file that is made by STAR for a reference index
    echo "No STAR index found in ${REFERENCE}. Creating index now..."
    
    # Create STAR index if no index from above (No SA file found in reference directory)
    STAR --runThreadN 8 --runMode genomeGenerate --genomeDir $REFERENCE \
         --genomeFastaFiles $REFERENCE/*.fna \  # if your reference is .fa or .fna, you have to change this
         --sjdbGTFfile $GTF --sjdbOverhang 100
else
    echo "STAR index already exists in ${REFERENCE}. Skipping index creation."
fi

# Run STAR alignment
if [[ $READ1 == *.gz ]] && [[ $READ2 == *.gz ]]; then
     echo "Fastq read files are zipped, will use zcat"
     UNZIP_COMMAND="--readFilesCommand zcat"
else
     echo "Fastq read files are unzipped, will continue without zcat"
     UNZIP_COMMAND=""
fi

echo "Running STAR with reference: $REFERENCE, GTF: $GTF, Read1: $READ1, Read2: $READ2, Output prefix: $OUT_PREFIX" # Making sure we put in the right things
STAR --runThreadN 8 --genomeDir $REFERENCE --sjdbGTFfile $GTF \
     --readFilesIn $READ1 $READ2 \
     --outFileNamePrefix $OUT_PREFIX \
     --outSAMtype BAM SortedByCoordinate \
     --quantMode GeneCounts \
     $UNZIP_COMMAND

# Rename the BAM file to include the output prefix
mv "${OUT_PREFIX}Aligned.sortedByCoord.out.bam" "${OUT_PREFIX}.bam"

# Index BAM file with Samtools
samtools index "${OUT_PREFIX}.bam"

echo "Alignment and indexing finished very swag. Output files are ${OUT_PREFIX}.bam and ${OUT_PREFIX}.bam.bai"

