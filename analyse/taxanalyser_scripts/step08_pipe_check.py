# Dependancies
import os
import pandas as pd
import numpy as np
import gzip
import argparse
from Bio import SeqIO

# parse arguments
parser = argparse.ArgumentParser(description='Read counts and quals at pipeline stages.')
parser.add_argument('--check_directory_path', type=str)
parser.add_argument('--file_basename', type=str)
args = parser.parse_args()
check_directory_path = args.check_directory_path
file_basename = args.file_basename
analysis_root_path = '/'.join(check_directory_path.split('/')[:-2])

prefix_dict = {
    'concatenated': 'nopreprocessing',
    'hqfiltered': 'qualityfiltered',
    'host_depleted': 'hostdepleted',
    'final_fastqs': 'finalfastq'
}

dir_to_check = check_directory_path.split('/')[-1]
prefix = prefix_dict[dir_to_check]

sample_path = f'{check_directory_path}/{file_basename}.fastq.gz'
output_path = f'{analysis_root_path}/output/pipeline_stats/read_counts_quals'
output_file_basename = f'{prefix}_{file_basename}'

seq_lens = []
avg_quals = []

if sample_path.endswith('.gz'):
    with gzip.open(sample_path, "rt") as input_file:
        for record in SeqIO.parse(input_file, "fastq"):
            seq_lens.append(len(record.seq))
            avg_quals.append(np.mean(record.letter_annotations['phred_quality']))

print(f'Median seq len: {np.median(seq_lens)}')
print(f'Median avg quality score: {np.median(avg_quals)}')

# save the data named after its sample 'sample[idx]'
data = pd.DataFrame({'seq_len': seq_lens, 'avg_qual': avg_quals})
data.to_csv(f'{output_path}/{output_file_basename}.csv', index=False)