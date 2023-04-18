#!/bin/bash

# Function to extract all rows from a bed that contains cell ID in a barcode file
function extract_rows_by_barcodes() {
  # Get the arguments
  local bed_file="$1"
  local cell_ids_file="$2"
  # Use grep to extract the rows that contain cell IDs from the cell IDs file
  grep -Fwf "$cell_ids_file" "$bed_file" 

}

fragment_file=/home/siluo/public/SiyuanLuo/projects/benchmark/raw_data/PBMC_multiomics/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz
fragment_file2=/home/siluo/public/SiyuanLuo/projects/benchmark/raw_data/PBMC_multiomics/pbmc_granulocyte_sorted_10k_atac_fragments.tsv
cell_ids_file=/home/siluo/public/SiyuanLuo/projects/benchmark/scripts/data_cleaning/PBMC_multiomics/PBMC_multiomics_Cells.txt
output_dir=/home/siluo/public/SiyuanLuo/projects/benchmark/cleaned_data/PBMC_multiomics

# mkdir $output_dir
# gzip -dk $fragment_file
extract_rows_by_barcodes $fragment_file2 $cell_ids_file > "${output_dir}"/pbmc_granulocyte_sorted_10k_atac_fragments.filtered.tsv

# sort according to coordinate
sort -k1,1 -k2,2n -k3,3n "${output_dir}"/pbmc_granulocyte_sorted_10k_atac_fragments.filtered.tsv > "${output_dir}"/pbmc_granulocyte_sorted_10k_atac_fragments.filtered.sorted.tsv
# bgzip
bgzip "${output_dir}"/pbmc_granulocyte_sorted_10k_atac_fragments.filtered.sorted.tsv
# tabix index
tabix -p bed "${output_dir}"/pbmc_granulocyte_sorted_10k_atac_fragments.filtered.sorted.tsv.gz