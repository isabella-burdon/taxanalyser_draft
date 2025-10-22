#!/bin/bash

# Define variables and directory paths
file_name=$1
analysis_output_path=$2

# Chopper params inputs
chopper_lengths_dir="$analysis_output_path/$run_name/data/QCMetrics/chopper_lengths"
lengths_file="$chopper_lengths_dir"/"$file_name"_nanoQClengths.txt

# Construct the input fastq.gz path
rebarcoded_dir="$analysis_output_path/data/rebarcoded"
input_fastq="$rebarcoded_dir/$file_name.fastq.gz"

# Construct the output fastq.gz path
final_fastqs_dir="$analysis_output_path/data/final_fastqs"
output_fastq="$final_fastqs_dir/$file_name.fastq.gz"

# Read head_length and tail_length from the lengths file
head_length=$(head -n 1 "$lengths_file")
echo "head_length: $head_length"
tail_length=$(sed -n '2p' "$lengths_file")
echo "tail_length: $tail_length"

# Run chopper with head_length and tail_length
gunzip -c "$input_fastq" | chopper --headcrop "$head_length" --tailcrop "$tail_length" | gzip > "$output_fastq"