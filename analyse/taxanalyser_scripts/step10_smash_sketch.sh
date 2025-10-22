#!/bin/bash
set -euo pipefail

# -------------------- Parse Input Arguments --------------------
file_basename="$1"
analysis_output_path="$2"

# -------------------- Define Paths --------------------
fastq_file="$analysis_output_path/data/final_fastqs/${file_basename}.fastq.gz"
temp_dir="$analysis_output_path/output/tax_profiling/smash_tmp"
sig_file="$temp_dir/${file_basename}.sig.gz"

# -------------------- Validate Inputs --------------------
if [ ! -f "$fastq_file" ]; then
    echo "Error: FASTQ file not found: $fastq_file" >&2
    exit 1
fi

if [ -f "$sig_file" ]; then
    echo "Skipping: $file_basename (signature already exists at $sig_file)"
    exit 0
fi

# -------------------- Run Sourmash Gather --------------------
mamba run -n smash sourmash sketch dna \
    -p k=31,abund "$fastq_file" \
    -o "$sig_file" \
    --name "$file_basename"