#!/usr/bin/env bash
#
# Script to merge gzipped sequencing result files ("fastq.gz") for different 
# lanes into one master file
# based on an example on https://medium.com/ngs-sh/merging-illumina-lanes-with-a-bash-script-112e0e0e0224
# by Andrea Telatin 2017, Bash Training for Bioinformatics, Quadram Institute

# optional input parameters
input_dir=${input_dir:-"./"}
output_dir=${output_dir:-"./"}
file_ext=${file_ext:-"fastq.gz"}

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
      exit
  fi
done

# make list of target files
filenames=(`ls ${input_dir} | grep ${file_ext}`)
#echo ${filenames[@]}
echo "Input files matching pattern: ${#filenames[@]}"
if [ ${#filenames[@]} == 0 ]; then
  echo "Found no files to process. Exiting."
  exit
fi

# loop through target files and merge to new files
for sample_file in ${filenames[@]};
do
  sample_name=$(basename "$sample_file"   | sed -e "s#_S[0-9]*_L00[1-4]_R[1-2]_001.${file_ext}##")
  sample_index=$(basename "$sample_file"  | grep -o  "_S[0-9]*_L00[1-4]_R[1-2]_001.${file_ext}" | grep -o "_S[0-9]*")
  sample_strand=$(basename "$sample_file" | grep -o  "_S[0-9]*_L00[1-4]_R[1-2]_001.${file_ext}" | grep -o "_R[1-2]")
  echo " > Adding $sample_file to ${sample_name}${sample_index}${sample_strand}.${file_ext}";
  cat ${input_dir}/$sample_file >> ${output_dir}/${sample_name}${sample_index}${sample_strand}.${file_ext};
done