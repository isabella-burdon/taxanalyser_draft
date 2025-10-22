#!/bin/bash

# Variables
data_root_path="$1"
run_name="$2"
file_basename="$3"  # Passed as the first argument
fastq_file="$data_root_path/$run_name/data/host_depleted/$file_basename.fastq.gz"
output_dir="$data_root_path/$run_name/output/assembly/$file_basename"

mkdir -p "$output_dir"

# -- (1) assemble MAGs with flye (default minimum overlap set to 1000)
echo ""
echo "Assembling contigs/MAGs for $file_basename..."
echo ""
mamba run -n flye flye --nano-hq $fastq_file --out-dir $output_dir/flye_ouput --meta 