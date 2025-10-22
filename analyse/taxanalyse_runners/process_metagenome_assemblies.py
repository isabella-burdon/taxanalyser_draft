import os
import subprocess as sp
import argparse
import sys
import pandas as pd

### PARSE ARGUMENTS ###
parser = argparse.ArgumentParser(description='Quality check for the run.')
parser.add_argument('--data_root_path', type=str, help='Path to the data directory')
parser.add_argument('--db_name', type=str, help='Database name')
parser.add_argument('--file_basename', type=str, help='Name of the file')
args = parser.parse_args()
data_root_path = args.data_root_path
db_name = args.db_name
folder = args.file_basename

### INITIALISE MANIFEST DF ###
manifest_add_df = pd.DataFrame(columns=['internal_location', 'md5', 'md5short', 'ksize', 'moltype', 'num',
                                   'scaled', 'n_hashes', 'with_abundance', 'name', 'filename'])

manifest_human_df = pd.DataFrame(columns=['internal_location', 'md5', 'md5short', 'ksize', 'moltype', 'num',
                                   'scaled', 'n_hashes', 'with_abundance', 'name', 'filename'])

### GET DATA ###
assembly_directory = f'{data_root_path}/output/assembly'
pipestats_directory = f"{data_root_path}/output/pipeline_stats"
genome_db=f"reference_dbs/{db_name}/{db_name}_files"
taxonomy_db=f"reference_dbs/{db_name}/{db_name}.taxonomy.sqldb"

txt_filepath = f"{assembly_directory}/flye/{folder}/assembly_info.txt"
assembly_info_df = pd.read_csv(txt_filepath, sep="\t")
all_contigs = assembly_info_df['#seq_name'].tolist()

long_circular_contigs = []
short_circular_contigs = []
verylong_contigs = []
long_contigs = []
short_contigs = []

for contig in assembly_info_df['#seq_name']:
    subset = assembly_info_df[assembly_info_df['#seq_name'] == contig]
    if int(subset['length'].values[0]) >= 100000 and str(subset['circ.'].values[0]) == 'Y':
        long_circular_contigs.append(contig)
    if int(subset['length'].values[0]) < 100000 and str(subset['circ.'].values[0]) == 'Y':
        short_circular_contigs.append(contig)
    if int(subset['length'].values[0]) >= 100000 and str(subset['circ.'].values[0]) == 'N':
        verylong_contigs.append(contig)
    if 10000 <= int(subset['length'].values[0]) < 100000 and str(subset['circ.'].values[0]) == 'N':
        long_contigs.append(contig)
    if int(subset['length'].values[0]) < 10000 and str(subset['circ.'].values[0]) == 'N':
        short_contigs.append(contig)

### CREATE SOURMASH SIGNATURES FOR ASSEMBLED CONTIGS ###
assembly_fasta_file = f"{assembly_directory}/flye/{folder}/assembly.fasta"
sourmash_output_dir = f'{assembly_directory}/{folder}/smash'
sig_dir = f'{assembly_directory}/{folder}/smash/sigs'

command = f'''mkdir -p {sourmash_output_dir}
mkdir -p {sig_dir}
sourmash sketch dna {assembly_fasta_file} -p scaled=1000 --force --singleton --outdir {sig_dir}
sourmash sig split {sig_dir}/assembly.fasta.sig --outdir {sig_dir} 
rm {sig_dir}/assembly.fasta.sig'''
sp.run(command, shell=True)

sig_filenames = os.listdir(sig_dir)
sig_filenames = [x for x in sig_filenames if x.endswith('.sig')]
md5_contig_dict = {}

for sig_filename in sig_filenames:

    full_file_path = f'{sig_dir}/{sig_filename}'
    command = f'sourmash sig describe "{full_file_path}" --csv "{full_file_path}.csv"'
    sp.run(command, shell=True)

    sig_info_df = pd.read_csv(f"{full_file_path}.csv")
    this_md5 = str(sig_info_df['md5'][0])
    this_contig_name = str(sig_info_df['name'][0])
    n_hashes = int(sig_info_df['n_hashes'][0])
    md5_contig_dict[this_md5] = this_contig_name, n_hashes

    command = f'mv "{full_file_path}" "{sig_dir}/{this_md5}.sig"; ' \
              f'mv "{full_file_path}.csv" "{sig_dir}/{this_md5}.sig.csv"'
    sp.run(command, shell=True)

### COMPARE CONTIGS ###
command = f'''mkdir -p {sourmash_output_dir}/sig_comparisons
sourmash compare {sig_dir}/*.sig -o {sourmash_output_dir}/sig_comparisons/compare_all.mat
sourmash plot {sourmash_output_dir}/sig_comparisons/compare_all.mat --output-dir {sourmash_output_dir}/sig_comparisons'''
sp.run(command, shell=True)

### ASSIGN TAXONOMY TO CONTIGS ###
sig_filenames = os.listdir(sig_dir)
sig_filenames = [x for x in sig_filenames if x.endswith('.sig')]

command = f'mkdir -p {sourmash_output_dir}/smash_gather_output'
sp.run(command, shell=True)

for key in md5_contig_dict.keys():
    n_hashes = md5_contig_dict[key][1]
    if n_hashes == 0:
        command = f'rm {sig_dir}/{key}.sig'
        sp.run(command, shell=True)

sig_filenames = os.listdir(sig_dir)
sig_filenames = [x for x in sig_filenames if x.endswith('.sig')]

for sig_filename in sig_filenames:

    this_md5 = sig_filename.split('.')[0]
    this_contig_name = md5_contig_dict[this_md5][0]

    if this_contig_name in long_circular_contigs or this_contig_name in verylong_contigs:
        threshold_bp = 50000
    elif this_contig_name in long_contigs:
        threshold_bp = 10000
    else:
        threshold_bp = 1000

    output_file_basepath = f'{sourmash_output_dir}/smash_gather_output/{this_contig_name}_{this_md5[:8]}_{threshold_bp}bp'

    command = f'sourmash gather {sig_dir}/{sig_filename} {genome_db} --threshold-bp {threshold_bp} --save-matches {output_file_basepath}_matches'
    sp.run(command, shell=True)

    if os.path.exists(f'{output_file_basepath}_matches'):
        command = f'''sourmash gather {sig_dir}/{sig_filename} {genome_db} --threshold-bp {threshold_bp} -o {output_file_basepath}.gather.gtdb.csv 
            sourmash tax genome -g {output_file_basepath}.gather.gtdb.csv -t {taxonomy_db} -o {output_file_basepath} -r species'''
        sp.run(command, shell=True)

        species = pd.read_csv(f'{output_file_basepath}.classifications.csv')['lineage'].values[0].split(';s__')[1]
        ident= f'BHI_{species.replace(" ", "_")}_{output_file_basepath.split("/")[-3]}_{output_file_basepath.split("/")[-1]}'

        command = f'sourmash signature rename {sig_dir}/{sig_filename} {ident} -o {sig_dir}/{sig_filename}'
        sp.run(command, shell=True)

        sig_info_df = pd.read_csv(f"{sig_dir}/{sig_filename}.csv")
        manifest_row = [
            f"signatures/{sig_info_df['md5'][0]}.sig", 
            sig_info_df['md5'][0], 
            sig_info_df['md5'][0][:8], 
            sig_info_df['ksize'][0], 
            sig_info_df['moltype'][0], 
            sig_info_df['num'][0], 
            sig_info_df['scaled'][0], 
            sig_info_df['n_hashes'][0], 
            sig_info_df['with_abundance'][0], 
            ident, 
            f"{sig_dir}/{sig_filename}"
        ]
        print(manifest_row)
        if species == 'Homo sapiens':
            manifest_human_df = pd.concat([manifest_human_df, pd.DataFrame([manifest_row], columns=manifest_human_df.columns)], ignore_index=True)
        else:
            manifest_add_df = pd.concat([manifest_add_df, pd.DataFrame([manifest_row], columns=manifest_add_df.columns)], ignore_index=True)

manifest_add_df = manifest_add_df.sort_values(by="n_hashes", ascending=False)
manifest_add_df.to_csv(f'{sourmash_output_dir}/manifest_bacterialcontigs.csv', index=False)
manifest_human_df = manifest_human_df.sort_values(by="n_hashes", ascending=False)
manifest_human_df.to_csv(f'{sourmash_output_dir}/manifest_humancontigs.csv', index=False)