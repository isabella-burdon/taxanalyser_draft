import os
import argparse
import pandas as pd
import re

parser = argparse.ArgumentParser(description='Quality check for the run.')
parser.add_argument('--input_data_path', type=str, help='Path to the input data directory')
parser.add_argument('--analysis_output_path', type=str, help='Name of the file')
parser.add_argument('--ont_kit_name', type=str, help='Name of the ONT kit used for sequencing')
parser.add_argument('--barcode_to_sample_sheet', type=str, help='Path to the barcode to sample mapping sheet (optional)', default=None)
args = parser.parse_args()
input_data_path = args.input_data_path
analysis_output_path = args.analysis_output_path
ont_kit_name = args.ont_kit_name
barcode_to_sample_sheet = args.barcode_to_sample_sheet

# Read in the barcode to sample mapping sheet
barcode_to_sample_sheet_df = pd.read_excel(barcode_to_sample_sheet)

# get all run identifiers from the raw fastq files
raw_folders = [x for x in os.listdir(input_data_path) if not x.startswith('.')]

all_files_in_raw_folders = []

for folder in raw_folders:
    print(folder)
    base_dir = f"{input_data_path}/{folder}"
    subdirs = [d for d in os.listdir(base_dir) if os.path.isdir(os.path.join(base_dir, d))]
    if len(subdirs) == 1:
        folder_path = os.path.join(base_dir, subdirs[0], "fastq_pass")
    else:
        raise ValueError(f"Expected 1 subdir in {base_dir}, found {len(subdirs)}: {subdirs}")
    print(folder_path)
    for barcode_dir in os.listdir(folder_path):
        barcode_path = os.path.join(folder_path, barcode_dir)
        if os.path.isdir(barcode_path):
            file_list = [re.sub(r'(_\d+)?\.fastq(?:\.gz)?$', '', f) for f in os.listdir(barcode_path) if f.endswith(('.fastq', '.fastq.gz')) and not f.startswith('.')]
            all_files_in_raw_folders.extend(file_list)
all_files_in_raw_folders = sorted(set(all_files_in_raw_folders))

print(all_files_in_raw_folders)

df_rows = []
for this_file in all_files_in_raw_folders:
    row = this_file.split('_')
    df_rows = df_rows + [row]

header_names = ['flow_cell_number', 'read_status', 'ONT_barcode', 'runid_prefix', 'runid_suffix']
df_run_info = pd.DataFrame(df_rows, columns=header_names)
df_run_info['runid'] = df_run_info['runid_prefix'] + '_' + df_run_info['runid_suffix']

df_run_sample_info = pd.merge(
    df_run_info[['flow_cell_number', 'read_status', 'ONT_barcode', 'runid']],
    barcode_to_sample_sheet_df[['flow_cell_number', 'ONT_barcode', 'runid', 'sample_id']],
    left_on=['flow_cell_number', 'runid', 'ONT_barcode'],
    right_on=['flow_cell_number', 'runid', 'ONT_barcode'],
    how='inner'
)
df_run_sample_info['runid_prefix'] = [x.split('_')[0] for x in df_run_sample_info['runid']]

runid_to_flowcell_dict = (
    df_run_sample_info
    .groupby('runid_prefix')['flow_cell_number'].unique()
    .apply(list)
    .to_dict()
)

# helper to normalize '04' -> 'barcode04' to match your Excel
def _normalize_barcode(bc: str) -> str:
    bc = str(bc)
    return bc if bc.startswith('barcode') else f'barcode{bc}'

def rebarcode_files(analysis_output_path, ont_kit_name):

    host_depleted_dir=f"{analysis_output_path}/data/host_depleted"
    file_basename_list = sorted([f.split('/')[-1].split('_bacterial')[0] for f in os.listdir(host_depleted_dir) if f.endswith('_bacterial.fastq.gz')])
    output_dir=f"{analysis_output_path}/data/rebarcoded"

    # initialise shell script for dorado commands
    with open(f"dorado_commands_tmp.sh", 'w') as f:
        f.write("#!/bin/bash\n\n")

    for file_basename in file_basename_list:

        # DOUBLE END BARCODE
        working_dir=f"{output_dir}/both_ends/{file_basename}"
        os.makedirs(working_dir, exist_ok=True)
        command = f'dorado demux "{host_depleted_dir}/{file_basename}_bacterial.fastq.gz" -v -o "{working_dir}" --kit-name "{ont_kit_name}" --barcode-both-ends --emit-fastq --emit-summary'
        with open(f"dorado_commands_tmp.sh", 'a') as f:
            f.write(f"{command}\n")

        # SINGLE END BARCODE
        working_dir=f"{output_dir}/single_end/{file_basename}"
        os.makedirs(working_dir, exist_ok=True)
        command = f'dorado demux "{host_depleted_dir}/{file_basename}_bacterial.fastq.gz" -v -o "{working_dir}" --kit-name "{ont_kit_name}" --emit-fastq --emit-summary'
        with open(f"dorado_commands_tmp.sh", 'a') as f:
            f.write(f"{command}\n")

def process_rebarcoded_files(analysis_output_path, which_end):

    rebarcoded_dir=f"{analysis_output_path}/data/rebarcoded"

    which_path = f'{rebarcoded_dir}/{which_end}'
    if not os.path.isdir(which_path):
        print(f"{which_path} does not exist. Skipping {which_end}.")
        return

    file_basename_list = sorted([f for f in os.listdir(which_path) if not f.startswith('.')])

    for file_basename in file_basename_list:
        
        working_dir=f"{rebarcoded_dir}/{which_end}/{file_basename}"
        if not os.path.isdir(working_dir):  # FIX: guard stray files
            continue

        # FIX: accept both .fastq and .fastq.gz
        input_file_names = sorted([
            f for f in os.listdir(working_dir)
            if (f.endswith('.fastq') or f.endswith('.fastq.gz'))
        ])
        
        # extract and save info from the info txt file
        info_txt_path = f"{working_dir}/barcoding_summary.txt"
        info_txt_contents = []
        if os.path.exists(info_txt_path):
            with open(info_txt_path, 'r') as f:
                info_txt_contents = f.readlines()
        else:
            print(f"{info_txt_path} does not exist. Skipping info parsing.")
            # don't continue; still move/concat files

        loop_df = []
        for line in info_txt_contents[1:]:
            barcode_identity = line.split('.fastq')[0].split('_')[-1]
            read_identity = line.split('\t')[1]
            row = [read_identity, barcode_identity]
            loop_df.append(row)
        if loop_df:
            loop_df = pd.DataFrame(loop_df, columns=['read_identity', 'barcode_identity'])
            stats_dir = f"{analysis_output_path}/output/pipeline_stats/barcoding"
            os.makedirs(stats_dir, exist_ok=True)
            if which_end == 'both_ends':
                output_stats_dir = f"{stats_dir}/double_barcoding"
            else:
                output_stats_dir = f"{stats_dir}/single_barcoding"
            os.makedirs(output_stats_dir, exist_ok=True)
            loop_df.to_csv(f"{output_stats_dir}/{file_basename}_barcode_counts.csv", index=False)

        # rename and move unclassified files
        unclassified_file_names = [f for f in input_file_names if 'unclassified' in f]
        if len(unclassified_file_names) == 0:
            print(f"No unclassified files found for {file_basename}. Skipping concatenation.")
        else:
            run_id = unclassified_file_names[0].split('-')[0]
            flowcell_id = runid_to_flowcell_dict.get(run_id, [])[0]
            print(flowcell_id)
            raw_bc = file_basename.split('barcode')[-1].split('-')[-1].split('.')[0]
            print(raw_bc)
            barcode = _normalize_barcode(raw_bc)
            print(barcode)

            matches = df_run_sample_info.loc[
                (df_run_sample_info['flow_cell_number'] == flowcell_id) & 
                (df_run_sample_info['ONT_barcode'] == barcode), 'sample_id'
            ].values
            print(matches)

            if len(matches) > 0 and pd.notna(matches[0]):
                sample_id = str(matches[0])
            else:
                sample_id = f"unknown_{flowcell_id}-{barcode}"
                print(f"No sample_id match for flow_cell_number={flowcell_id}, barcode={barcode}. Using '{sample_id}'.")

            # FIX: preserve extension (.fastq or .fastq.gz) and concat safely
            out_ext = '.fastq.gz' if any(n.endswith('.gz') for n in unclassified_file_names) else '.fastq'
            out_path = f"{rebarcoded_dir}/unclassified-{which_end.replace('_', '-')}_{sample_id.replace('_', '-')}{out_ext}"
            if out_ext.endswith('.gz'):
                cmd = f"cat {' '.join([f'{working_dir}/{x}' for x in unclassified_file_names])} > {out_path}"
            else:
                cmd = f"cat {' '.join([f'{working_dir}/{x}' for x in unclassified_file_names])} > {out_path}"
            print(cmd)
            os.system(cmd)

        # rename and move classified files
        classified_file_names = [f for f in input_file_names if 'unclassified' not in f]
        if len(classified_file_names) == 0:
            print(f"No classified files found for {file_basename}. Skipping concatenation.")
        else:
            print(f"Processing classified files for {file_basename}...")
            for idx, input_file_name in enumerate(classified_file_names):
                run_id = input_file_name.split('_')[0][:8]
                flowcell_id_list = runid_to_flowcell_dict.get(run_id, [])
                flowcell_id = flowcell_id_list[0] if flowcell_id_list else None

                # last underscore token before extension -> 'barcode04'
                stem = input_file_name.rsplit('.', 1)[0]
                barcode = stem.split('_')[-1]
                barcode = _normalize_barcode(barcode)

                matches = df_run_sample_info.loc[
                    (df_run_sample_info['flow_cell_number'] == flowcell_id) & 
                    (df_run_sample_info['ONT_barcode'] == barcode), 'sample_id'
                ].values

                if len(matches) > 0 and pd.notna(matches[0]):
                    sample_id = str(matches[0])
                else:
                    sample_id = f"unknown_{flowcell_id}-{barcode}"
                    print(f"No sample_id match for flow_cell_number={flowcell_id}, barcode={barcode}. Using '{sample_id}'.")

                # FIX: preserve original extension
                ext = '.fastq.gz' if input_file_name.endswith('.gz') else '.fastq'
                output_file_name = f"classified_{which_end.replace('_', '-')}_{sample_id.replace('_', '-')}_{flowcell_id}-{barcode}_n{idx}.fastq"
                input_file_path = f"{working_dir}/{input_file_name}"
                output_file_path = f"{rebarcoded_dir}/{which_end}/{output_file_name}"
                command = f"mv {input_file_path} {output_file_path}"
                print(command)
                os.system(command)

    # make dictionary of existing files to match by sample ID
    classified_existing_files = [
        f for f in os.listdir(f"{rebarcoded_dir}/{which_end}/")
        if not f.startswith('.') and 'unclassified' not in f and (f.endswith('.fastq') or f.endswith('.fastq.gz'))]

    classified_existing_files_dict = {}
    for existing_file in classified_existing_files:
        sample_id = existing_file.split('_')[2]
        classified_existing_files_dict.setdefault(sample_id, []).append(existing_file)

    # concatenate classified files based on sample ID
    for sample_id, items_files in classified_existing_files_dict.items():
        out_ext = '.fastq.gz' if any(x.endswith('.gz') for x in items_files) else '.fastq'
        output_file_name = f"classified-{which_end.replace('_', '-')}_{sample_id}{out_ext}"
        cmd = f"cat {' '.join([f'{rebarcoded_dir}/{which_end}/{x}' for x in items_files])} > {rebarcoded_dir}/{output_file_name}"
        print(cmd)
        os.system(cmd)

    command = f"rm -rf {rebarcoded_dir}/{which_end}"
    print(command)
    os.system(command)

rebarcode_files(analysis_output_path, ont_kit_name)
shell_script_path = 'dorado_commands_tmp.sh'
os.system(f"chmod +x {shell_script_path}")
os.system(f"./{shell_script_path}")
os.system(f"rm {shell_script_path}")

process_rebarcoded_files(analysis_output_path, 'both_ends')
process_rebarcoded_files(analysis_output_path, 'single_end')

# remove all unclassified single-end files
single_unclassified_files = [f for f in os.listdir(f"{analysis_output_path}/data/rebarcoded") if f.startswith('unclassified-single-end') and (f.endswith('.fastq') or f.endswith('.fastq.gz')) and not f.startswith('.')]
for f in single_unclassified_files:
    os.remove(f"{analysis_output_path}/data/rebarcoded/{f}")