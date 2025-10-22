#!/bin/bash
set -euo pipefail

# -------------------- Parse Input Arguments --------------------
file_basename="$1"
analysis_output_path="$2"
db_name="$3"
threshold_bp="$4"
lineage_ident_level="$5"


# -------------------- Define Paths --------------------
genome_db="reference_dbs/$db_name/${db_name}_files"
taxonomy_db="reference_dbs/$db_name/$db_name.taxonomy.sqldb"
fastq_file="$analysis_output_path/data/final_fastqs/${file_basename}.fastq.gz"
output_dir="$analysis_output_path/output/tax_profiling/smash_output"
temp_dir="$analysis_output_path/output/tax_profiling/smash_tmp"
sig_file="$temp_dir/${file_basename}.sig.gz"

sample_name="${file_basename}_${threshold_bp}bp_${db_name}"
gather_csv="$temp_dir/${sample_name}.gather.k31.csv"
matches_file="$temp_dir/matches_${sample_name}"
profile_csv="$output_dir/${sample_name}.profile.csv"

# -------------------- Create Directories --------------------
mkdir -p "$output_dir" "$temp_dir"

# -------------------- Validate Inputs --------------------
if [ ! -f "$sig_file" ]; then
    echo "Error: Signature file not found: $sig_file" >&2
    exit 1
fi

if [ -f "$profile_csv" ]; then
    echo "Profile already exists: $profile_csv — skipping."
    exit 0
fi

# -------------------- Run Sourmash Gather --------------------
mamba run -n smash sourmash gather "$sig_file" "$genome_db" \
    --threshold-bp "$threshold_bp" \
    --save-matches "$matches_file" \
    --output "$gather_csv"

# -------------------- Run Sourmash Taxonomy Profiling --------------------
mamba run -n smash sourmash tax metagenome -g "$gather_csv" \
    -o "$profile_csv" \
    -t "$taxonomy_db" -F krona -r "$lineage_ident_level"