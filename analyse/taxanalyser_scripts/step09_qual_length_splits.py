# Import dependencies
import os
import gzip
import numpy as np
from Bio import SeqIO
import pandas as pd
import argparse

parser = argparse.ArgumentParser(description='Quality check for the run.')
parser.add_argument('--analysis_output_path', type=str, help='Name of the file')
args = parser.parse_args()
analysis_output_path = args.analysis_output_path

input_dir_path = f'{analysis_output_path}/data/final_fastqs'
output_dir_path = f'{analysis_output_path}/data/length_qual_splits'

# Thresholds
param_LR = 5000
param_MR = 2000
param_SR = 500
param_HQ = 30
param_MQ = 20

def process_file(input_dir_path, output_dir_path, file_base_name):
    input_file = f'{input_dir_path}/{file_base_name}.fastq.gz'

    output_suffixes = {
        "lrhq": 'LR-HQ', "lrmq": 'LR-MQ',
        "mrhq": 'MR-HQ', "mrmq": 'MR-MQ',
        "srhq": 'SR-HQ', "srmq": 'SR-MQ',
    }

    output_handles = {}
    written_counts = {}  # Track number of reads written per file

    # Open all output files
    for label, suffix in output_suffixes.items():
        out_path = f"{output_dir_path}/{file_base_name}_{suffix}.fastq.gz"
        output_handles[label] = gzip.open(out_path, "wt")
        written_counts[label] = 0

    try:
        with gzip.open(input_file, "rt") as input_handle:
            for record in SeqIO.parse(input_handle, "fastq"):
                mean_quality = np.mean(record.letter_annotations['phred_quality'])
                read_length = len(record.seq)

                # Assign tiers
                if read_length >= param_LR:
                    tier = "lr"
                elif read_length >= param_MR:
                    tier = "mr"
                elif read_length >= param_SR:
                    tier = "sr"
                else:
                    tier = "vsr"

                if mean_quality >= param_HQ:
                    qual = "hq"
                elif mean_quality >= param_MQ:
                    qual = "mq"
                else:
                    qual = "lq"

                label = f"{tier}{qual}"
                handle = output_handles.get(label)
                if handle is not None:
                    SeqIO.write(record, handle, "fastq")
                    written_counts[label] += 1
    finally:
        # Close all handles
        for handle in output_handles.values():
            handle.close()
        # Delete any empty files (no reads written)
        for label, count in written_counts.items():
            if count == 0:
                suffix = output_suffixes[label]
                file_path = f"{output_dir_path}/{file_base_name}_{suffix}.fastq.gz"
                if os.path.exists(file_path):
                    os.remove(file_path)
                    print(f"Removed empty file: {file_path}")

file_basename_list = sorted([x.split('.')[0] for x in os.listdir(input_dir_path) if x.endswith('.fastq.gz') and x.startswith('classified-single-end')])

for file_base_name in file_basename_list:
    print(f'Processing file: {file_base_name}')
    process_file(input_dir_path, output_dir_path, file_base_name)