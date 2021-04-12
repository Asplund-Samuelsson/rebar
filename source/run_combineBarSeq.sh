#!/usr/bin/env bash
#
# Script to trim, map and summarize sequencing reads obtained from a barcoded transposon library.
# Authors: Johannes Asplund-Samuelsson, Michael Jahn

# optional input parameters
input_dir=${input_dir:-"./"}
output_dir=${output_dir:-"./"}
poolfile=${poolfile:-"./ref/poolfile.tsv"}
cb_path=${cb_path:-"./source/feba/bin/combineBarSeq.pl"}

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

# Run combineBarSeq.pl script on "codes" files
${cb_path} ${output_dir}/result ${poolfile} ${input_dir}/*.codes
