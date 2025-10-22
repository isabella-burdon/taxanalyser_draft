import os
import subprocess as sp
import argparse
import sys
import pysam

# parse args
parser = argparse.ArgumentParser(description='Quality check for the run.')
parser.add_argument('--analysis_output_path', type=str, help='Name of the run')
parser.add_argument('--file_basename', type=str, help='Name of the file')
args = parser.parse_args()
analysis_output_path = args.analysis_output_path
file_basename = args.file_basename

# Define functions
def extract_read_ids(bam_path):
    mapped_ids = set()
    unmapped_ids = set()
    with pysam.AlignmentFile(bam_path, "rb") as bamfile:
        for read in bamfile:
            if read.is_unmapped:
                unmapped_ids.add(read.query_name)
            else:
                mapped_ids.add(read.query_name)
    return mapped_ids, unmapped_ids
        
# set data paths
input_data_directory = f'{analysis_output_path}/data/hqfiltered'
output_data_directory = f'{analysis_output_path}/data/host_depleted'
chm13_fasta_filepath = "reference_dbs/deplete_genomes/chm13.fasta"

input_lr_fastq_file = f"{input_data_directory}/{file_basename}_LR.fastq.gz"
input_sr_fastq_file = f"{input_data_directory}/{file_basename}_SR.fastq.gz"
output_file_prefix = f"{output_data_directory}/{file_basename}"

#### MAP READS AGAINST CHM13 ####

# Run minimap2 to generate BAM files with unmapped reads
lr_bam = f"{output_file_prefix}_LR_chm13.bam"
sr_bam = f"{output_file_prefix}_SR_chm13.bam"

# Commands
lr_command = f"""
minimap2 -t 16 -ax map-hifi {chm13_fasta_filepath} {input_lr_fastq_file} |
samtools view -b -o {lr_bam}
"""

sr_command = f"""
minimap2 -t 16 -ax sr {chm13_fasta_filepath} {input_sr_fastq_file} |
samtools view -b -o {sr_bam}
"""

# Execute alignment + BAM
try:
    sp.run(lr_command, shell=True, check=True, stderr=sp.PIPE, text=True)
    sp.run(sr_command, shell=True, check=True, stderr=sp.PIPE, text=True)
except sp.CalledProcessError:
    print("An error occurred during mapping. Please check minimap2 and samtools.")
    sys.exit(1)

#### EXTRACT FASTQ IDS - (UNMAPPED = BACTERIAL) (MAPPED = HUMAN) ####
lr_chm13_ids, lr_unmapped_ids = extract_read_ids(lr_bam)
sr_chm13_ids, sr_unmapped_ids = extract_read_ids(sr_bam)
read_id_lists = [lr_unmapped_ids, sr_unmapped_ids, lr_chm13_ids, sr_chm13_ids]

# txt files for read IDs
lr_bacterial_txt = f"{output_file_prefix}_LR_unmapped_ids.txt"
sr_bacterial_txt = f"{output_file_prefix}_SR_unmapped_ids.txt"
lr_chm13_txt = f"{output_file_prefix}_LR_chm13_ids.txt"
sr_chm13_txt = f"{output_file_prefix}_SR_chm13_ids.txt"
read_txt_files = [lr_bacterial_txt, sr_bacterial_txt, lr_chm13_txt, sr_chm13_txt]

for read_ids, txt_file in zip(read_id_lists, read_txt_files):
    with open(txt_file, 'w') as f:
        for read_id in read_ids:
            f.write(f"{read_id}\n")

#### EXTRACT FASTQ READS - MAPPED (HUMAN) ####

# Set output filenames
output_lr_bacterial_fastq = f"{output_file_prefix}_LR_unmapped.fastq.gz"
output_sr_bacterial_fastq = f"{output_file_prefix}_SR_unmapped.fastq.gz"
output_lr_chm13_fastq = f"{output_file_prefix}_LR_chm13.fastq.gz"
output_sr_chm13_fastq = f"{output_file_prefix}_SR_chm13.fastq.gz"

# Output fastqs with original headers
print(f"Extracting bacterial long reads:")
sp.run(f"seqkit grep -f {lr_bacterial_txt} {input_lr_fastq_file} -o {output_lr_bacterial_fastq}", shell=True, check=True)
print(f"Extracting bacterial short reads:")
sp.run(f"seqkit grep -f {sr_bacterial_txt} {input_sr_fastq_file} -o {output_sr_bacterial_fastq}", shell=True, check=True)
print(f"Extracting human long reads:")
sp.run(f"seqkit grep -f {lr_chm13_txt} {input_lr_fastq_file} -o {output_lr_chm13_fastq}", shell=True, check=True)
print(f"Extracting human short reads:")
sp.run(f"seqkit grep -f {sr_chm13_txt} {input_sr_fastq_file} -o {output_sr_chm13_fastq}", shell=True, check=True)

### CLEAN UP ###
final_bacterial_fastq = f"{output_file_prefix}_bacterial.fastq.gz"
sp.run(f"gzcat {output_lr_bacterial_fastq} {output_sr_bacterial_fastq} | gzip > {final_bacterial_fastq}", shell=True, check=True)

final_chm13_fastq = f"{output_file_prefix}_chm13.fastq.gz"
sp.run(f"gzcat {output_lr_chm13_fastq} {output_sr_chm13_fastq} | gzip > {final_chm13_fastq}", shell=True, check=True)

# final_hqfiltered_fastq = f"{input_data_directory}/{file_basename}.fastq.gz"
# sp.run(f"gzcat {input_lr_fastq_file} {input_sr_fastq_file} | gzip > {final_hqfiltered_fastq}", shell=True, check=True)

final_hqfiltered_fastq = f"{input_data_directory}/{file_basename}.fastq.gz"
sp.run(f"cat {input_lr_fastq_file} {input_sr_fastq_file} > {final_hqfiltered_fastq}",
       shell=True, check=True)


# Remove intermediate files
os.remove(lr_bam)
os.remove(sr_bam)
os.remove(lr_bacterial_txt)
os.remove(sr_bacterial_txt)
os.remove(lr_chm13_txt)
os.remove(sr_chm13_txt)
os.remove(output_lr_bacterial_fastq)
os.remove(output_sr_bacterial_fastq)
os.remove(output_lr_chm13_fastq)
os.remove(output_sr_chm13_fastq)
os.remove(input_lr_fastq_file)
os.remove(input_sr_fastq_file)