#!/bin/bash

# TO RUN (ENSURE IN ROOT DIRECTORY): 
# chmod +x analyse/run_taxanalyser.sh
#  
# bash ./analyse/run_taxanalyser.sh --start_from "concatenate_B96" --barcode_to_sample_sheet "analyse/barcode_to_sample_sheets/ISMS_RUN13.xlsx" --input_data_path "/Volumes/SSD07/ISMS_rawdata/L13_ISMS_B96" --analysis_output_path "/Volumes/SSD07/ISMS_analysis/L13_ISMS_B96" --ont_kit_name "SQK-NBD114-96" --threshold_windows 50000 15000 1000 --reference_databases "sinus_alldomain_database" --lineage_ident_level "strain"
# bash ./analyse/run_taxanalyser.sh --start_from "qual_length_splits" --barcode_to_sample_sheet "analyse/barcode_to_sample_sheets/ISMS_RUN15.xlsx" --input_data_path "/Volumes/SSD07/ISMS_rawdata/L15_ISMS_B96" --analysis_output_path "/Volumes/SSD07/ISMS_analysis/L15_ISMS_B96" --ont_kit_name "SQK-NBD114-96" --threshold_windows 50000 15000 1000 --reference_databases "sinus_alldomain_database" --lineage_ident_level "strain"

# ./analyse/clear_metapipe.sh L5_ISMS_B96 "/Volumes/SSD05/ISMS_rawdata"

# options to start from: concatenate_B96, concatenate_B24, hqfilter, host_deplete, double_end_barcode, nanoqc, chopper, pipe_check, smash, flye_meta, process_mags

################ Parse command-line arguments ################
while [[ $# -gt 0 ]]; do
    case $1 in
        --start_from)
            start_from="$2"
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
export start_from barcode_to_sample_sheet input_data_path analysis_output_path ont_kit_name threshold_windows reference_databases lineage_ident_level

# Set default values if not provided
start_from="${start_from:-concatenate}" # Default to concatenation if not specified

# Use the variables as needed
echo " "
echo "Running TAXANALYSER with the following parameters:"
echo " - Start From: $start_from"
echo " - Input Data Path: $input_data_path"
echo " - Data Output Path: $analysis_output_path"
echo " - ONT Kit Name: $ont_kit_name"
echo " - Lineage Identification Level: $lineage_ident_level"
echo " - Reference Databases: ${reference_databases[@]}"
echo " - Threshold Windows: ${threshold_windows[@]}"
echo " "

################ (STEP 01) SET UP ################

#### make new folders ####
chmod +x analyse/taxanalyser_scripts/step01_mkdirs_permissions.sh
bash ./analyse/taxanalyser_scripts/step01_mkdirs_permissions.sh $analysis_output_path
cat_output_dir="$analysis_output_path/data/concatenated"
export cat_output_dir 

#### save output_log.txt ####
exec > >(tee -a "$analysis_output_path/output/pipeline_stats/output_log.txt") 2>&1
set -xe
trap 'echo "Error occurred at line $LINENO"; exit 1;' ERR

# tell the script where to start based on start_from variable
if [[ "$start_from" == "concatenate" ]]; then
    echo "Starting from concatenation step... for SQKNBD96 kit"
    steps_to_do=("concatenate" "hqfilter" "host_deplete" "double_end_barcode" "nanoqc" "chopper" "pipe_check" "qual_length_splits" "smash" "flye_meta" "process_mags" "annotate_function")
    echo "Steps to do: ${steps_to_do[@]}"
elif [[ "$start_from" == "concatenate_B96" ]]; then
    echo "Starting from concatenation step... for SQKNBD96 kit"
    steps_to_do=("concatenate_B96" "hqfilter" "host_deplete" "double_end_barcode" "nanoqc" "chopper" "pipe_check" "qual_length_splits" "smash" "flye_meta" "process_mags" "annotate_function")
    echo "Steps to do: ${steps_to_do[@]}"
elif [[ "$start_from" == "concatenate_B24" ]]; then
    echo "Starting from concatenation step... for SQKNBD24 kit"
    steps_to_do=("concatenate_B24" "hqfilter" "host_deplete" "double_end_barcode" "nanoqc" "chopper" "pipe_check" "qual_length_splits" "smash" "flye_meta" "process_mags" "annotate_function")
    echo "Steps to do: ${steps_to_do[@]}"
elif [[ "$start_from" == "hqfilter" ]]; then
    echo "Starting from hqfilter step."
    steps_to_do=("hqfilter" "host_deplete" "double_end_barcode" "nanoqc" "chopper" "pipe_check" "qual_length_splits" "smash" "flye_meta" "process_mags" "annotate_function")
    echo "Steps to do: ${steps_to_do[@]}"
elif [[ "$start_from" == "host_deplete" ]]; then
    echo "Starting from host_deplete step."
    steps_to_do=("host_deplete" "double_end_barcode" "nanoqc" "chopper" "pipe_check" "qual_length_splits" "smash" "flye_meta" "process_mags" "annotate_function")
    echo "Steps to do: ${steps_to_do[@]}"
elif [[ "$start_from" == "double_end_barcode" ]]; then
    echo "Starting from double_end_barcode step."
    steps_to_do=("double_end_barcode" "nanoqc" "chopper" "pipe_check" "qual_length_splits" "smash" "flye_meta" "process_mags" "annotate_function")
    echo "Steps to do: ${steps_to_do[@]}"
elif [[ "$start_from" == "nanoqc" ]]; then
    echo "Starting from nanoqc step."
    steps_to_do=("nanoqc" "chopper" "pipe_check" "qual_length_splits" "smash" "flye_meta" "process_mags" "annotate_function")
    echo "Steps to do: ${steps_to_do[@]}"
elif [[ "$start_from" == "chopper" ]]; then
    echo "Starting from chopper step."
    steps_to_do=("chopper" "pipe_check" "qual_length_splits" "smash" "flye_meta" "process_mags" "annotate_function")
    echo "Steps to do: ${steps_to_do[@]}"
elif [[ "$start_from" == "pipe_check" ]]; then
    echo "Starting from pipe_check step."
    steps_to_do=("pipe_check" "qual_length_splits" "smash" "flye_meta" "process_mags" "annotate_function")
    echo "Steps to do: ${steps_to_do[@]}"
elif [[ "$start_from" == "qual_length_splits" ]]; then
    echo "Starting from qual_length_splits step."
    steps_to_do=("qual_length_splits" "smash" "flye_meta" "process_mags" "annotate_function")
    echo "Steps to do: ${steps_to_do[@]}"
elif [[ "$start_from" == "smash" ]]; then
    echo "Starting from smash step."
    steps_to_do=("smash" "flye_meta" "process_mags" "annotate_function")
    echo "Steps to do: ${steps_to_do[@]}"
elif [[ "$start_from" == "flye_meta" ]]; then
    echo "Starting from flye_meta step."
    steps_to_do=("flye_meta" "process_mags" "annotate_function")
    echo "Steps to do: ${steps_to_do[@]}"
elif [[ "$start_from" == "process_mags" ]]; then
    echo "Starting from process_mags step."
    steps_to_do=("process_mags" "annotate_function")
    echo "Steps to do: ${steps_to_do[@]}"
elif [[ "$start_from" == "annotate_function" ]]; then
    echo "Starting from process_mags step."
    steps_to_do=("annotate_function")
    echo "Steps to do: ${steps_to_do[@]}"
else
    echo "Invalid start_from option. Exiting."
    exit 1
fi

################ (STEP 02) CONCATENATION ################

### 96 Barcode concat (4 barcodes per sample)
if [[ " ${steps_to_do[@]} " =~ " concatenate_B96 " ]]; then
    echo "Starting concatenation step."
    bash ./analyse/taxanalyser_scripts/step02_concatenate_B96.sh $input_data_path $cat_output_dir
else
    echo "No concatenate B96 step."
fi

### 24 Barcode concat (1 barcodes per sample)
if [[ " ${steps_to_do[@]} " =~ " concatenate_B24 " ]]; then
    echo "Starting concatenation step."
    bash ./analyse/taxanalyser_scripts/step02_concatenate_B24.sh $input_data_path $cat_output_dir
else
    echo "No concatenate B24 step."
fi

################ (STEP 03) HQFILTER ################
if [[ " ${steps_to_do[@]} " =~ " hqfilter " ]]; then
    echo "Starting hqfilter step."

    concatenated_dir="$analysis_output_path/data/concatenated"
    export concatenated_dir

    small_file_basename_list=()
    while IFS= read -r line; do
        [[ "$line" == ._* ]] && continue
        small_file_basename_list+=("$line")
    done < <(find "$concatenated_dir" -mindepth 1 -maxdepth 1 -type f -name "*.fastq.gz" -size -5G -exec basename {} \; | sed 's/\.fastq\.gz$//')

    big_file_basename_list=()
    while IFS= read -r line; do
        [[ "$line" == ._* ]] && continue
        big_file_basename_list+=("$line")
    done < <(find "$concatenated_dir" -mindepth 1 -maxdepth 1 -type f -name "*.fastq.gz" -size +5G -exec basename {} \; | sed 's/\.fastq\.gz$//')

    #### -- (5) filter out low quality reads  the chopped reads (default is phred 20)
    # Define and export filter function
    run_hq_filter_stringent() {
        file_basename="$1"
        mamba run -n SeqIO python analyse/taxanalyser_scripts/step03_hq_filter_stringent.py --data_root_path "$analysis_output_path" --file_name "$file_basename"
    }
    export -f run_hq_filter_stringent
    # Run in parallel
    parallel --will-cite --color --tag --eta -j 24 run_hq_filter_stringent ::: "${small_file_basename_list[@]}"
    parallel --will-cite --color --tag --eta -j 6 run_hq_filter_stringent ::: "${big_file_basename_list[@]}"

else
    echo "Skipping hqfilter step."
fi

################ (STEP 04) HOST DEPLETE ################
if [[ " ${steps_to_do[@]} " =~ " host_deplete " ]]; then
    echo "Starting host deplete step."

    concatenated_dir="$analysis_output_path/data/concatenated"
    export concatenated_dir

    small_file_basename_list=()
    while IFS= read -r line; do
        [[ "$line" == ._* ]] && continue
        small_file_basename_list+=("$line")
    done < <(find "$concatenated_dir" -mindepth 1 -maxdepth 1 -type f -name "*.fastq.gz" -size -5G -exec basename {} \; | sed 's/\.fastq\.gz$//')

    big_file_basename_list=()
    while IFS= read -r line; do
        [[ "$line" == ._* ]] && continue
        big_file_basename_list+=("$line")
    done < <(find "$concatenated_dir" -mindepth 1 -maxdepth 1 -type f -name "*.fastq.gz" -size +5G -exec basename {} \; | sed 's/\.fastq\.gz$//')
    #### -- (6) deplete reads that map to chm13
    # Define and export host depletion function
    host_deplete() {
        file_basename="$1"
        echo "Human (host) reads filtered out of $file_basename"
        mamba run -n minimap2 python analyse/taxanalyser_scripts/step04_map_deplete.py --analysis_output_path "$analysis_output_path" --file_basename "$file_basename"
    }
    export -f host_deplete
    # Run in parallel
    parallel --will-cite --color --tag --eta -j 8 host_deplete ::: "${small_file_basename_list[@]}" 
    parallel --will-cite --color --tag --eta -j 4 host_deplete ::: "${big_file_basename_list[@]}" 

else
    echo "Skipping host deplete step."
fi


################ (STEP 05) DOUBLE ENDED BARCODING ################
if [[ " ${steps_to_do[@]} " =~ " double_end_barcode " ]]; then
    echo "Starting double end barcode step."
    mamba run -n dorado python analyse/taxanalyser_scripts/step05_rebarcode.py --input_data_path "$input_data_path" --analysis_output_path "$analysis_output_path" --ont_kit_name "$ont_kit_name" --barcode_to_sample_sheet "$barcode_to_sample_sheet"
    gzip $analysis_output_path/data/rebarcoded/*.fastq
else
    echo "Skipping double end barcode step."
fi

################ (STEP 06) NANOQC ################
if [[ " ${steps_to_do[@]} " =~ " nanoqc " ]]; then
    echo "Starting nanoQC step."

    # Define and export nanoQC inputs/functions
    rebarcoded_dir="$analysis_output_path/data/rebarcoded"
    export rebarcoded_dir

    file_basename_list=()
    while IFS= read -r line; do
        [[ "$(basename "$line")" == ._* ]] && continue
        file_basename_list+=("$line")
    done < <(find "$rebarcoded_dir" -mindepth 1 -maxdepth 1 -type f -name "*.fastq.gz" -exec basename {} \; | sed 's/\.fastq\.gz$//')
    
    qcmetric_output_dir="$analysis_output_path/data/QCMetrics"
    export qcmetric_output_dir 

    run_nanoqc() {
        file_basename="$1"
        input_fastq="$rebarcoded_dir/$file_basename.fastq.gz"
        output_dir="$qcmetric_output_dir/$file_basename"
        mamba run -n nanoqc nanoQC -o "$output_dir" "$input_fastq"
        echo "Completed quality check of $file_basename"
    }
    export -f run_nanoqc

    # Run in parallel
    parallel --will-cite --colour --tag --eta -j 24 run_nanoqc ::: "${file_basename_list[@]}"

else
    echo "Skipping nanoQC step."
fi

################ (STEP 07) CHOPPER ################
if [[ " ${steps_to_do[@]} " =~ " chopper " ]]; then
    echo "Starting chopper step."

    # Define and export nanoQC inputs/functions
    rebarcoded_dir="$analysis_output_path/data/rebarcoded"
    export rebarcoded_dir

    file_basename_list=()
    while IFS= read -r line; do
        [[ "$(basename "$line")" == ._* ]] && continue
        file_basename_list+=("$line")
    done < <(find "$rebarcoded_dir" -mindepth 1 -maxdepth 1 -type f -name "*.fastq.gz" -exec basename {} \; | sed 's/\.fastq\.gz$//')

    #### -- (3) get parameters for chopper (fast no parallellisation required)
    mamba run -n chopper python analyse/taxanalyser_scripts/step06_chopper_params.py --analysis_output_path "$analysis_output_path"
    
    #### -- (4) trim low qual bases from end of reads 
    # Define and export trim function
    run_chopper() {
        file_basename="$1"
        mamba run -n chopper bash ./analyse/taxanalyser_scripts/step07_run_chopper.sh "$file_basename" "$analysis_output_path"
        echo "Completed trimming of $file_basename"
    }
    export -f run_chopper
    # Run in parallel
    parallel --will-cite --color --tag --eta -j 24 run_chopper ::: "${file_basename_list[@]}"

else
    echo "Skipping chopper step."
fi

################ (STEP 08) ANALYSE PIPE EFFECTS ################
if [[ " ${steps_to_do[@]} " =~ " pipe_check " ]]; then
    echo "Starting pipe_check."

    #### -- (7) analyse read counts and quals at key stages in the pipeline
    # Define and export pipe check function
    pipe_check() {
        file_basename="$1"
        check_directory="$2"
        check_directory_path="$analysis_output_path/data/$check_directory"
        echo "Analysed read counts and quals from $file_basename"
        mamba run -n SeqIO python analyse/taxanalyser_scripts/step08_pipe_check.py --check_directory_path "$check_directory_path" --file_basename "$file_basename"
        echo "***"
    }
    export -f pipe_check

    file_basename_list=()
    while IFS= read -r line; do
        [[ "$line" == ._* ]] && continue
        file_basename_list+=("$line")
    done < <(find "$analysis_output_path/data/concatenated" -mindepth 1 -maxdepth 1 -type f -name "*.fastq.gz" -exec basename {} \; | sed 's/\.fastq\.gz$//')
    parallel --will-cite --color --tag --eta -j 24 pipe_check {} concatenated ::: "${file_basename_list[@]}"

    file_basename_list=()
    while IFS= read -r line; do
        [[ "$line" == ._* ]] && continue
        file_basename_list+=("$line")
    done < <(find "$analysis_output_path/data/hqfiltered" -mindepth 1 -maxdepth 1 -type f -name "*.fastq.gz" -exec basename {} \; | sed 's/\.fastq\.gz$//')
    parallel --will-cite --color --tag --eta -j 24 pipe_check {} hqfiltered ::: "${file_basename_list[@]}"


    file_basename_list=()
    while IFS= read -r line; do
        [[ "$line" == ._* ]] && continue
        file_basename_list+=("$line")
    done < <(find "$analysis_output_path/data/host_depleted" -mindepth 1 -maxdepth 1 -type f -name "*_bacterial.fastq.gz" -exec basename {} \; | sed 's/\.fastq\.gz$//')
    parallel --will-cite --color --tag --eta -j 24 pipe_check {} host_depleted ::: "${file_basename_list[@]}"

    file_basename_list=()
    while IFS= read -r line; do
        [[ "$line" == ._* ]] && continue
        file_basename_list+=("$line")
    done < <(find "$analysis_output_path/data/final_fastqs" -mindepth 1 -maxdepth 1 -type f -name "*.fastq.gz" -exec basename {} \; | sed 's/\.fastq\.gz$//')
    parallel --will-cite --color --tag --eta -j 24 pipe_check {} final_fastqs ::: "${file_basename_list[@]}"
    
else
    echo "Completed pipe_check step."
fi

################ (STEP 09) SPLIT FILE BY QUALS AND LENGTHS ################
if [[ " ${steps_to_do[@]} " =~ " qual_length_splits " ]]; then
    echo "Starting qual_length_splits step."

    # run rebarcoding
    mamba run -n SeqIO python analyse/taxanalyser_scripts/step09_qual_length_splits.py --analysis_output_path "$analysis_output_path"

else
    echo "Skipping double end barcode step."
fi

################ (STEP 10-A) RUN SMASH ON FINAL FASTQs ################
if [[ " ${steps_to_do[@]} " =~ " smash " ]]; then
    echo "Starting sourmash mash screen step."

    smash_fastq_dir="$analysis_output_path/data/final_fastqs"
    export smash_fastq_dir

    file_basename_list=()
    while IFS= read -r line; do
        [[ "$(basename "$line")" == ._* ]] && continue
        file_basename_list+=("$line")
    done < <(find "$smash_fastq_dir" -mindepth 1 -maxdepth 1 -type f -name "*.fastq.gz" -exec basename {} \; | sed 's/\.fastq\.gz$//')

    smash_sketch() {
        file_basename="$1"
        echo "Running sourmash sketch for: $file_basename"
        bash ./analyse/taxanalyser_scripts/step10_smash_sketch.sh $file_basename $analysis_output_path
        echo "***"
    }
    export -f smash_sketch
    parallel --will-cite --color --tag --eta -j 24 smash_sketch ::: ${file_basename_list[@]}
    
    # -- run host depleted files against various databases with the human reference genome
    smash_gather() {
        file_basename="$1"
        database_name="$2"
        threshold_window="$3"
        echo "Running sourmash on $file_basename with a threshold of $threshold_window bp using databse: $database_name"
        bash ./analyse/taxanalyser_scripts/step10_smash_taxgather.sh $file_basename $analysis_output_path $database_name $threshold_window $lineage_ident_level
        echo "***"
    }
    export -f smash_gather
    parallel --will-cite --color --tag --eta -j 12 smash_gather \
        ::: "${file_basename_list[@]}" \
        ::: "${reference_databases[@]}" \
        ::: "${threshold_windows[@]}"
else
    echo "Skipping smash step."
fi

################ (STEP 10-B) RUN SMASH ON LENGTH QUAL SPLIT FASTQs ################
if [[ " ${steps_to_do[@]} " =~ " smash " ]]; then
    echo "Starting sourmash mash screen step."

    smash_fastq_dir="$analysis_output_path/data/length_qual_splits"
    export smash_fastq_dir

    file_basename_list=()
    while IFS= read -r line; do
        [[ "$(basename "$line")" == ._* ]] && continue
        file_basename_list+=("$line")
    done < <(find "$smash_fastq_dir" -mindepth 1 -maxdepth 1 -type f -name "*.fastq.gz" -exec basename {} \; | sed 's/\.fastq\.gz$//')

    smash_sketch() {
        file_basename="$1"
        echo "Running sourmash sketch for: $file_basename"
        bash ./analyse/taxanalyser_scripts/step10_smash_sketch.sh $file_basename $analysis_output_path
        echo "***"
    }
    export -f smash_sketch
    parallel --will-cite --color --tag --eta -j 24 smash_sketch ::: ${file_basename_list[@]}
    
    # -- run host depleted files against various databases with the human reference genome
    smash_gather() {
        file_basename="$1"
        database_name="$2"
        threshold_window="$3"
        echo "Running sourmash on $file_basename with a threshold of $threshold_window bp using databse: $database_name"
        bash ./analyse/taxanalyser_scripts/step10_smash_taxgather.sh $file_basename $analysis_output_path $database_name $threshold_window $lineage_ident_level
        echo "***"
    }
    export -f smash_gather
    parallel --will-cite --color --tag --eta -j 12 smash_gather \
        ::: "${file_basename_list[@]}" \
        ::: "${reference_databases[@]}" \
        ::: "${threshold_windows[@]}"
else
    echo "Skipping smash step."
fi

# ################# RUN FLYE #################
# if [[ " ${steps_to_do[@]} " =~ " flye_meta " ]]; then
#     echo "Starting flye-assembly step."

#     hostdep_dir="$data_root_path/$run_name/data/host_depleted"
#     export hostdep_dir

#     file_basename_list=()
#     while IFS= read -r line; do
#         [[ "$(basename "$line")" == ._* ]] && continue
#         file_basename_list+=("$line")
#     done < <(find "$hostdep_dir" -mindepth 1 -maxdepth 1 -type f -name "*.fastq.gz" -exec basename {} \; | sed 's/\.fastq\.gz$//')

#     run_flye_meta() {
#         file_basename="$1"
#         mamba run -n flye bash ./analyse/scripts_parallel/run_flye_meta.sh "$data_root_path" "$run_name" "$file_basename"
#         echo "Running flye-assembly of $file_basename"
#     }
#     export -f run_flye_meta
#     # Run in parallel
#     parallel --will-cite --color --tag --eta -j 6 run_flye_meta ::: "${file_basename_list[@]}"

#     break

# else
#     echo "Skipping flye-assembly step."
# fi

# ################ PROCESS MAGS ################
# if [[ " ${steps_to_do[@]} " =~ " process_mags " ]]; then
#     echo "Starting MAG processing step."

#     hostdep_dir="$data_root_path/$run_name/data/host_depleted"
#     export hostdep_dir

#     file_basename_list=()
#     while IFS= read -r line; do
#         [[ "$(basename "$line")" == ._* ]] && continue
#         file_basename_list+=("$line")
#     done < <(find "$hostdep_dir" -mindepth 1 -maxdepth 1 -type f -name "*.fastq.gz" -exec basename {} \; | sed 's/\.fastq\.gz$//')

#     run_process_mags() {
#         database_name="$1"
#         file_basename="$2"
#         echo "Processing MAGs for $file_basename against $database_name"
#         mamba run -n smash python analyse/scripts_parallel/process_metagenome_assemblies.py --run_name "$run_name" --data_root_path "$data_root_path" --db_name "$database_name" --file_basename "$file_basename"
#     }
#     export -f run_process_mags
#     # Run in parallel
#     parallel --will-cite --color --tag --eta run_process_mags ::: ${reference_databases[@]} ::: "${file_basename_list[@]}"

#     break

# else
#     echo "Skipping MAG processing step."
# fi

# # tidy pipe
# rm -rf $data_root_path/$run_name/data/concatenated 
# rm -rf $data_root_path/$run_name/data/hqfiltered 
# rm -rf $data_root_path/$run_name/data/host_depleted 
# rm -rf $data_root_path/$run_name/data/QCMetrics 