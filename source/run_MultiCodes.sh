#!/usr/bin/env bash

# Project name
PROJ=$1

# Input fastq.gz directory
IDIR="data/projects/${PROJ}"

# Intermediate output directory ("codes" and "counts" files from MultiCodes)
ODIR="intermediate/projects/${PROJ}"

# Path to FEBA MultiCodes.pl script
MC="source/feba/bin/MultiCodes.pl"

# List all input fastq.gz files
ls ${IDIR} | grep ".fastq.gz$" | while read Infile;
  do
    # Extract ID of fastq.gz file
    ID=`echo ${Infile} | cut -f 1 -d \.`
    # Set up custom MultiCodes command to run in parallel
    echo "zcat ${IDIR}/${ID}.fastq.gz | \
      ${MC} -out ${ODIR}/${ID} -index ${ID} -minQuality 0 \
      -preseq CAGCGTACG -postseq AGAGACCTC -nPreExpected 11:17 \
      >& ${ODIR}/${ID}.log"
  done | \
parallel --no-notice --jobs 16
