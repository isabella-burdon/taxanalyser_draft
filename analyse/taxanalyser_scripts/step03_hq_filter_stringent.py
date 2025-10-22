# Import dependencies
import os
import gzip
import numpy as np
import argparse
import time
import signal

# Parse arguments
parser = argparse.ArgumentParser(description='Quality check for the run.')
parser.add_argument('--data_root_path', type=str, required=True, help='Root path to data directory')
parser.add_argument('--file_name', type=str, required=True, help='Name of the file')
args = parser.parse_args()
data_root_path = args.data_root_path
file_base_name = args.file_name

# Define paths
input_dir_path = f'{data_root_path}/data/concatenated'
output_path = f'{data_root_path}/data/hqfiltered'
faulty_dir = f'{data_root_path}/output/pipeline_stats/faulty_fastqs'
txt_path = f'{faulty_dir}/hqfilter_badfastq_summary.txt'

# Create required directories
os.makedirs(output_path, exist_ok=True)
os.makedirs(faulty_dir, exist_ok=True)

# Filtering parameters
param_LR = 2000
param_SR = 500
param_HQ = 30
param_MQ = 20
record_timeout = 10  # seconds before a record is declared "stuck"

# Timeout signal handler
class TimeoutException(Exception):
    pass

def timeout_handler(signum, frame):
    raise TimeoutException("Record processing timed out.")

# Set the signal handler
signal.signal(signal.SIGALRM, timeout_handler)

# Function to process reads and filter
def process_file(file_base_name):
    input_file = f'{input_dir_path}/{file_base_name}.fastq.gz'
    output_file_lr = f"{output_path}/{file_base_name}_LR.fastq.gz"
    output_file_sr = f"{output_path}/{file_base_name}_SR.fastq.gz"
    malformed_fastq_path = f"{faulty_dir}/{file_base_name}_MALFORMED.fastq"

    lr_handle = gzip.open(output_file_lr, "wt")
    sr_handle = gzip.open(output_file_sr, "wt")
    malformed_handle = open(malformed_fastq_path, "w")

    bad_read_count = 0
    bad_read_indices = []

    with gzip.open(input_file, "rt") as handle:
        record_index = 0

        while True:
            try:
                # Start the timeout timer for a single record
                signal.alarm(record_timeout)

                lines = [next(handle) for _ in range(4)]
                signal.alarm(0)  # Cancel alarm on success

                if not lines[0].startswith("@") or not lines[2].startswith("+"):
                    raise ValueError("Invalid FASTQ structure")

                seq = lines[1].strip()
                qual = lines[3].strip()

                if len(seq) != len(qual):
                    raise ValueError("Sequence and quality length mismatch")

                phred_scores = [ord(ch) - 33 for ch in qual]
                mean_quality = sum(phred_scores) / len(phred_scores)
                read_length = len(seq)

                # Tier classification
                if read_length >= param_LR:
                    tier = "lr"
                elif read_length >= param_SR:
                    tier = "sr"
                else:
                    tier = "vsr"

                if mean_quality >= param_HQ:
                    qual_tier = "hq"
                elif mean_quality >= param_MQ:
                    qual_tier = "mq"
                else:
                    qual_tier = "lq"

                label = f"{tier}{qual_tier}"

                if label in ["lrhq", "lrmq"]:
                    lr_handle.writelines(lines)
                # elif label == "srhq":
                elif label in ["srhq", "srmq"]:
                    sr_handle.writelines(lines)

            except StopIteration:
                signal.alarm(0)
                break

            except Exception as e:
                signal.alarm(0)
                bad_read_count += 1
                bad_read_indices.append(record_index)
                malformed_handle.writelines(lines if 'lines' in locals() else [])
                print(f"[Warning] Skipping malformed record #{record_index} in {file_base_name}: {e}")

            # Progress update
            if record_index > 0 and record_index % 100_000 == 0:
                print(f"[{file_base_name}] Processed {record_index} reads...")

            record_index += 1

    lr_handle.close()
    sr_handle.close()
    malformed_handle.close()

    txt_lines = [
        f"File name: {file_base_name}",
        f"  Total number of bad reads: {bad_read_count}",
    ]
    if bad_read_count > 0:
        txt_lines.append(f"  Bad read indices: {', '.join(map(str, bad_read_indices))}")
    txt_lines.append("\n")

    # print("\n".join(txt_lines))

    with open(txt_path, 'a') as f:
        f.write("\n".join(txt_lines) + "\n")

    return bad_read_indices

# Function to extract bad records based on known indices
def extract_bad_fastq_records(file_base_name, bad_indices):
    start_time = time.time()

    input_file = f'{input_dir_path}/{file_base_name}.fastq.gz'
    output_file = f'{faulty_dir}/{file_base_name}_bad.fastq'

    bad_indices_set = set(bad_indices)  # faster lookup

    count_extracted = 0
    record_index = 0

    try:
        with gzip.open(input_file, "rt") as handle, open(output_file, "w") as output_handle:
            while True:
                try:
                    line1 = next(handle)
                    line2 = next(handle)
                    line3 = next(handle)
                    line4 = next(handle)
                except StopIteration:
                    break

                if record_index in bad_indices_set:
                    output_handle.writelines([line1, line2, line3, line4])
                    count_extracted += 1

                record_index += 1

    except (KeyboardInterrupt, TimeoutError) as e:
        print(f"[Error] Extraction interrupted: {e}")
    except Exception as e:
        print(f"[Error] Unexpected error during extraction: {e}")

    duration = time.time() - start_time
    print(f"Extracted {count_extracted} bad records to: {output_file} in {duration:.2f} seconds.")

# Run the process with timing
start_time = time.time()

bad_read_indices = process_file(file_base_name)

if bad_read_indices:
    extract_bad_fastq_records(file_base_name, bad_read_indices)
else:
    print(f"{file_base_name} processed.")

end_time = time.time()
elapsed = end_time - start_time
print(f"Processing of {file_base_name} completed in {elapsed:.2f} seconds.")