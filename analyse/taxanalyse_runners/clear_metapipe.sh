# chmod +x analyse/clear_metapipe.sh
# ./analyse/clear_metapipe.sh <run_name>
# ./analyse/clear_metapipe.sh mock9strain data_root_path

run_name=$1
data_root_path=$2 

# delete folders
rm -rf $data_root_path/$run_name/data
rm -rf $data_root_path/$run_name/output