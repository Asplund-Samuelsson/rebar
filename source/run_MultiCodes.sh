#!/usr/bin/env bash
#
# Script to trim, map and summarize sequencing reads obtained from a barcoded transposon library.
# Authors: Johannes Asplund-Samuelsson, Michael Jahn

# optional input parameters
input_dir=${input_dir:-"./"}
output_dir=${output_dir:-"./"}
pattern=${pattern:-".fastq.gz"}
mc_path=${mc_path:-"./source/feba/bin/MultiCodes.pl"}

# assign optional parameters that were passed with "--"
while [ $# -gt 0 ]; do
  if [[ $1 == *"--"* ]]; then
    param="${1/--/}"
    declare $param="$2"
  fi
  shift
done


# if in and output folders are not present, throw error
for dir in $input_dir $output_dir
do
	if [ -d ${dir} ]; then
	  echo "Input/Output directory: ${dir} exists"
  else
    echo "ERROR: Input/Output directory: ${dir} does not exist, stopping."
      exit 9
  fi
done


# List all input fastq.gz files
ls ${input_dir} | grep ${pattern} | while read fastq;
  do
    # extract name of fastq.gz file
    filename=`echo ${fastq} | cut -f 1 -d \.`
    # Set up custom MultiCodes command to run in parallel
    echo "zcat ${input_dir}/${filename}.fastq.gz | \
      ${mc_path} -out ${output_dir}/${filename} -index ${filename} -minQuality 0 \
      -preseq CAGCGTACG -postseq AGAGACCTC -nPreExpected 11:17 \
      >& ${output_dir}/${filename}.log"
  done | \
parallel --no-notice --jobs 16
