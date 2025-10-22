#!/bin/bash

# TO RUN (ENSURE IN ROOT DIRECTORY): 
# chmod +x taxanalyser_scripts/run_taxanalyser.sh
#  
# bash ./taxanalyser_scripts/run_smash_anystage.sh --smash_stage host_depleted_nofilt --barcode_to_sample_sheet "/Volumes/SSD06/initial_mbiome_project/new_analysis/sample_id_sheet.xlsx" --input_data_path "/Volumes/SSD06/initial_mbiome_project/raw_fastqs" --analysis_output_path "/Volumes/SSD06/initial_mbiome_project/new_analysis" --ont_kit_name "SQK-NBD114-24" --threshold_windows 5000 1000 --reference_databases "sinus_alldomain_database" --lineage_ident_level "strain"
# options to start from: , nanoqc, chopper, pipe_check, smash, flye_meta, process_mags
# smash stages: concatenated, hqfiltered, host_depleted, host_depleted_nofilt, final_fastqs
################ Parse command-line arguments ################
while [[ $# -gt 0 ]]; do
    case $1 in
        --smash_stage)
            smash_stage="$2"
            shift 2
            ;;
        --barcode_to_sample_sheet)
            barcode_to_sample_sheet="$2"
            shift 2
            ;;
        --input_data_path)
            input_data_path="$2"
            shift 2
            ;;
        --analysis_output_path)
            analysis_output_path="$2"
            shift 2
            ;;
        --ont_kit_name)
            ont_kit_name="$2"
            shift 2
            ;;
        --threshold_windows)
            shift # Skip the flag itself
            threshold_windows=() # Initialize an array
            while [[ $# -gt 0 && $1 != --* ]]; do
                threshold_windows+=("$1") # Add each threshold value to the array
                shift
            done
            ;;
        --reference_databases)
            shift # Skip the flag itself
            reference_databases=() # Initialize an array
            while [[ $# -gt 0 && $1 != --* ]]; do
                reference_databases+=("$1") # Add each database name to the array
                shift
            done
            ;;
        --lineage_ident_level)
            lineage_ident_level="$2"
            shift 2
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done
export smash_stage barcode_to_sample_sheet input_data_path analysis_output_path ont_kit_name threshold_windows reference_databases lineage_ident_level

# Use the variables as needed
echo " "
echo "Running TAXANALYSER with the following parameters:"
echo " - Input Data Path: $input_data_path"
echo " - Data Output Path: $analysis_output_path"
echo " - ONT Kit Name: $ont_kit_name"
echo " - Lineage Identification Level: $lineage_ident_level"
echo " - Reference Databases: ${reference_databases[@]}"
echo " - Threshold Windows: ${threshold_windows[@]}"
echo " "

################# RUN SOURMASH #################
final_fastq_dir="$analysis_output_path/data/$smash_stage"
export final_fastq_dir

file_basename_list=()
while IFS= read -r line; do
    [[ "$(basename "$line")" == ._* ]] && continue
    file_basename_list+=("$line")
done < <(find "$final_fastq_dir" -mindepth 1 -maxdepth 1 -type f -name "*.fastq.gz" -exec basename {} \; | sed 's/\.fastq\.gz$//')

smash_sketch() {
    file_basename="$1"
    echo "Running sourmash sketch for: $file_basename"
    bash ./taxanalyser_scripts/step09_smash_sketch_anystage.sh $file_basename $analysis_output_path $smash_stage
    echo "***"
}
export -f smash_sketch
chmod +x ./taxanalyser_scripts/step09_smash_sketch_anystage.sh
parallel --will-cite --color --tag --eta -j 24 smash_sketch ::: ${file_basename_list[@]}

# -- run host depleted files against various databases with the human reference genome
smash_gather() {
    file_basename="$1"
    database_name="$2"
    threshold_window="$3"
    # don't run if file exists
    output_dir="$analysis_output_path/output/tax_profiling/${smash_stage}/smash_output"
    sample_name="${file_basename}_${threshold_bp}bp_${db_name}"
    profile_csv="$output_dir/${sample_name}.profile.csv"

    if [ ! -f "$profile_csv" ]; then
        echo "Running sourmash on $file_basename from stage $smash_stage with a threshold of $threshold_window bp using database: $database_name"
        bash ./taxanalyser_scripts/step10_smash_taxgather_anystage.sh $file_basename $analysis_output_path $database_name $threshold_window $lineage_ident_level $smash_stage
        echo "***"
    else
        echo "Smash analysis of $sample_name completed, skipping..."
    fi
}
export -f smash_gather
chmod +x ./taxanalyser_scripts/step10_smash_taxgather_anystage.sh
parallel --will-cite --color --tag --eta -j 12 smash_gather \
    ::: "${file_basename_list[@]}" \
    ::: "${reference_databases[@]}" \
    ::: "${threshold_windows[@]}"