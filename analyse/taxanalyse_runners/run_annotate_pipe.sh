# chmod +x analyse/run_annotate_pipe.sh
# bash +x analyse/run_annotate_pipe.sh

################# PARAMS #################
data_root_path="/Volumes/SSD06/initial_mbiome_project/new_analysis"
input_data_directory="$data_root_path/data/final_fastqs"
output_dir="$data_root_path/output/assembly/flye"

# Use --nano-hq for HAC/SUP ONT reads; switch to --nano-raw for raw
flye_mode="--nano-hq"

threads_per_job=4
n_jobs=8

mkdir -p "$output_dir"

echo "Scanning FASTQ files in: $input_data_directory"

# Build a list of FASTQ basenames (no readarray; macOS-safe)
file_basename_list=()
while IFS= read -r -d '' f; do
  file_basename_list+=( "$(basename "$f")" )
done < <(find "$input_data_directory" -mindepth 1 -maxdepth 1 -type f \
  \( -name "*.fastq" -o -name "*.fastq.gz" -o -name "*.fq" -o -name "*.fq.gz" \) -print0)

# Show what we'll process
printf '  - %s\n' "${file_basename_list[@]}"

export output_dir threads_per_job flye_mode input_data_directory

################# RUN FLYE #################
run_flye_meta() {
  file_basename="$1"
  fastq_file="$input_data_directory/$file_basename"

  # derive per-file output directory (avoids clobbering)
  sample="${file_basename%.gz}"
  sample="${sample%.fastq}"
  sample="${sample%.fq}"
  out="$output_dir/$sample"
  mkdir -p "$out"

  echo "Running flye-assembly of $file_basename → $out"
  mamba run -n flye flye $flye_mode "$fastq_file" \
    --out-dir "$out" \
    --threads "$threads_per_job" \
    --meta
}
export -f run_flye_meta
parallel --color --tag --will-cite -j "$n_jobs" run_flye_meta ::: "${file_basename_list[@]}"

################ PROCESS MAGS ################
database_name='sinus_alldomain_database'
run_process_mags() {
    file_basename="$1"
    echo "Processing MAGs for $file_basename against $database_name"
    mamba run -n smash python analyse/scripts_parallel/process_metagenome_assemblies.py --data_root_path "$data_root_path" --db_name "$database_name" --file_basename "$file_basename"
}
export -f run_process_mags
parallel --color --tag run_process_mags ::: "${file_basename_list[@]}"

################ RUN BAKTA ################
run_bakta() {
    file_basename="$1"
    # pass positional args to the script
    mamba run -n bakta analyse/scripts_parallel/run_bakta.sh "$data_root_path" "$file_basename"
}
export -f run_bakta
export data_root_path
parallel --color --tag --will-cite -j 3 run_bakta ::: "${file_basename_list[@]}"