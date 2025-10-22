#!/bin/bash
#!/bin/bash

input_data_path=$1
cat_output_dir=$2

# ---------- (1) Define functions ----------

fix_corrupted_files_concat() {
    folder_name="$1"
    folder_path="${fastq_pass_dir}/${folder_name}"
    mkdir -p "${cat_output_dir}/${sub_runname}"
    output_file="${cat_output_dir}/${sub_runname}/${folder_name}.fastq.gz"

    tmp_concat=$(mktemp)

    # Try decompressing all files safely
    for fq in "${folder_path}"/*.fastq.gz; do
        if gzip -t "$fq" 2>/dev/null; then
            gunzip -c "$fq" >> "$tmp_concat"
        else
            echo "Skipping corrupted file: $fq"
        fi
    done

    # Recompress the concatenated, uncompressed fastqs
    gzip -c "$tmp_concat" > "$output_file"
    rm -f "$tmp_concat"
    echo "Concatenated and recompressed ${folder_name} to ${output_file}"
}
export -f fix_corrupted_files_concat

four_barcode_concat() {
    b1=$(printf "%02d" "$1")
    b2=$(printf "%02d" "$2")
    b3=$(printf "%02d" "$3")
    b4=$(printf "%02d" "$4")
    output="${cat_output_dir}/${sub_runname}/barcode${b1}-${b4}.fastq.gz"

    gunzip -c \
        "$cat_output_dir/$sub_runname/barcode${b1}.fastq.gz" \
        "$cat_output_dir/$sub_runname/barcode${b2}.fastq.gz" \
        "$cat_output_dir/$sub_runname/barcode${b3}.fastq.gz" \
        "$cat_output_dir/$sub_runname/barcode${b4}.fastq.gz" | gzip > "$output"

    rm -f \
        "$cat_output_dir/$sub_runname/barcode${b1}.fastq.gz" \
        "$cat_output_dir/$sub_runname/barcode${b2}.fastq.gz" \
        "$cat_output_dir/$sub_runname/barcode${b3}.fastq.gz" \
        "$cat_output_dir/$sub_runname/barcode${b4}.fastq.gz"

    echo "Concatenation of $output complete"
}
export -f four_barcode_concat

final_concat() {
    file_basename="$1"
    all_paths=$(find $concat_folders -type f -name "${file_basename}.fastq.gz")
    gunzip -c $all_paths | gzip > "$cat_output_dir/${file_basename}.fastq.gz"
    gzip -t "$cat_output_dir/${file_basename}.fastq.gz" && echo "Final concat of ${file_basename}.fastq.gz complete"
}
export -f final_concat

# ---------- (2) Process each subrun directory ----------

ont_directory="$input_data_path"
ont_directory_list=$(find "$ont_directory" -mindepth 1 -maxdepth 1 -type d)

echo "ONT directories LISTED: $ont_directory_list"

for full_ont_dir_path in $ont_directory_list; do
    fastq_pass_dir=$(find "$full_ont_dir_path" -mindepth 2 -maxdepth 2 -type d -name "fastq_pass" | head -n 1)
    echo "Processing directory: $fastq_pass_dir"
    sub_runname=$(basename "$full_ont_dir_path")
    export fastq_pass_dir sub_runname

    file_basename_list=($(find "$fastq_pass_dir" -mindepth 1 -maxdepth 1 -type d -exec basename {} \; | grep -v '^unclassified$' | sort))
    echo "Found data directories: ${file_basename_list[*]}"

    mkdir -p "$cat_output_dir/$sub_runname"
    parallel --will-cite --color --tag --eta -j 96 fix_corrupted_files_concat ::: "${file_basename_list[@]}"

    if (( ${#file_basename_list[@]} >= 4 )); then
        parallel --will-cite --color --tag --eta -j 24 four_barcode_concat \
        ::: $(seq 1 4 $((${#file_basename_list[@]} - 3))) \
        :::+ $(seq 2 4 $((${#file_basename_list[@]} - 2))) \
        :::+ $(seq 3 4 $((${#file_basename_list[@]} - 1))) \
        :::+ $(seq 4 4 $((${#file_basename_list[@]})))
    else
        echo "Skipping barcode concat: only ${#file_basename_list[@]} barcodes found"
    fi
done

# ---------- (3) Final concat across all runs ----------

concat_folders=$(find "$cat_output_dir" -mindepth 1 -maxdepth 1 -type d)
first_concat_folder=$(echo "$concat_folders" | head -n 1)
file_basename_list=($(find "$first_concat_folder" -type f -name "*.fastq.gz" -exec basename {} .fastq.gz \; | grep -v '^\._'))

export concat_folders file_basename_list
parallel --will-cite --color --tag --eta -j 24 final_concat ::: "${file_basename_list[@]}"

# ---------- (4) Clean up ----------
for concat_folder in $concat_folders; do
    rm -rf "$concat_folder"
done