#!/bin/bash

# Positional args (from caller)
data_root_path="$1"
file_basename="$2"

# Derive sample name from filename
sample="${file_basename%.gz}"
sample="${sample%.fastq}"
sample="${sample%.fq}"

bakta_db="/Users/issyburdon/taxanalyser_app_MB/reference_dbs/bakta_db"

# Primary Flye contig path
fasta_file="$data_root_path/output/assembly/flye/$sample/assembly.fasta"

# Check contig file; fallback to hybracter incomplete if needed
if [ -f "$fasta_file" ]; then
    echo "Found fasta file: $fasta_file"
else
    echo "No fasta file found at: $fasta_file"
    genome_completeness='incomplete'
    fasta_file="$data_root_path/output/hybracter/$sample/FINAL_OUTPUT/incomplete/${sample}_final.fasta"
    if [ -f "$fasta_file" ]; then
        echo "Found fallback fasta file: $fasta_file"
    else
        echo "No fallback fasta found at: $fasta_file"
        exit 1
    fi
fi

output_dir="$data_root_path/output/bakta/$sample"
mkdir -p "$data_root_path/output/bakta" "$output_dir"

# Run Bakta
bakta --db "$bakta_db" \
      --output "$output_dir" \
      --force --verbose --keep-contig-headers \
      "$fasta_file"
