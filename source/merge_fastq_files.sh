#!/usr/bin/env bash
#
# Script to merge gzipped sequencing result files ("fastq.gz") for different 
# lanes into one master file
# based on an example on https://medium.com/ngs-sh/merging-illumina-lanes-with-a-bash-script-112e0e0e0224
# by Andrea Telatin 2017, Bash Training for Bioinformatics, Quadram Institute

# Project name
PROJ=$1

# Input fastq.gz directory (Illuminas default dir structure)
input_dir="data/projects/${PROJ}/Data/Intensities/BaseCalls/"

# outout folder for merged files
output_dir="data/projects/${PROJ}"

# default file extension
ext='fastq.gz'

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


# Loop files in {input_dir} with extension {ext}
for sample_file in ${input_dir}/*_*.${ext};
do
  sample_name=$(basename "$sample_file"   | sed -e "s#_S[1-9]*_L00[1-4]_R._001.fastq.gz##")
  sample_index=$(basename "$sample_file"  | grep -o "_S[1-9]*_L00[1-4]_R._001.fastq.gz" | grep -o "_S[0-9]*")
  sample_strand=$(basename "$sample_file" | grep -o "_S[1-9]*_L00[1-4]_R._001.fastq.gz" | grep -o "_R[1-4]")
  
  echo " > Adding $sample_file to ${sample_name}${sample_index}${sample_strand}.${ext}";
  cat $sample_file >> ${output_dir}/${sample_name}${sample_index}${sample_strand}.${ext};
done
