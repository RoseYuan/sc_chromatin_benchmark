#!/bin/bash

# Function to extract all rows from a bed that contains cell ID in a barcode file
function extract_rows_by_barcodes() {
  # Get the arguments
  local bed_file="$1"
  local cell_ids_file="$2"
  # Use grep to extract the rows that contain cell IDs from the cell IDs file
  grep -Fwf "$cell_ids_file" "$bed_file" 

}

input_dir=/home/siluo/projects/sc_chromatin_benchmark/raw_data
output_dir=/home/siluo/projects/sc_chromatin_benchmark/cleaned_data
bed_file_gz="${input_dir}"/Chen_2019/fragments.sort.bed.gz
bed_file="${input_dir}"/Chen_2019/fragments.sort.bed
cell_ids_file="Cells_Chen_2019.txt"


# gzip -dk $bed_file_gz

# extract_rows_by_barcodes $bed_file $cell_ids_file > "${output_dir}"/Chen_2019_filtered_fragments.bed

# sort according to coordinate
sort -k1,1 -k2,2n -k3,3n "${output_dir}"/Chen_2019_filtered_fragments.bed > "${output_dir}"/Chen_2019_filtered_fragments_sorted.bed
# bgzip
bgzip "${output_dir}"/Chen_2019_filtered_fragments_sorted.bed
# tabix index
tabix -p bed "${output_dir}"/Chen_2019_filtered_fragments_sorted.bed.gz