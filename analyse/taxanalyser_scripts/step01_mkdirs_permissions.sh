analysis_output_path=$1

# make directories
mkdir -p $analysis_output_path
mkdir -p $analysis_output_path/data
mkdir -p $analysis_output_path/data/concatenated
mkdir -p $analysis_output_path/data/hqfiltered
mkdir -p $analysis_output_path/data/host_depleted
mkdir -p $analysis_output_path/data/rebarcoded
mkdir -p $analysis_output_path/data/rebarcoded/single_end
mkdir -p $analysis_output_path/data/rebarcoded/both_ends
mkdir -p $analysis_output_path/data/QCMetrics
mkdir -p $analysis_output_path/data/QCMetrics/chopper_lengths
mkdir -p $analysis_output_path/data/final_fastqs
mkdir -p $analysis_output_path/data/length_qual_splits
mkdir -p $analysis_output_path/output
mkdir -p $analysis_output_path/output/assembly
mkdir -p $analysis_output_path/output/tax_profiling
mkdir -p $analysis_output_path/output/tax_profiling/smash_output
mkdir -p $analysis_output_path/output/tax_profiling/smash_tmp
mkdir -p $analysis_output_path/output/pipeline_stats
mkdir -p $analysis_output_path/output/pipeline_stats/read_counts_quals

# confirm permissions
chmod +x analyse/taxanalyser_scripts/step01_mkdirs_permissions.sh
chmod +x analyse/taxanalyser_scripts/step02_concatenate_B24.sh
chmod +x analyse/taxanalyser_scripts/step02_concatenate_B96.sh
chmod +x analyse/taxanalyser_scripts/step07_run_chopper.sh
chmod +x analyse/taxanalyser_scripts/step09_smash_sketch.sh
chmod +x analyse/taxanalyser_scripts/step10_smash_taxgather.sh