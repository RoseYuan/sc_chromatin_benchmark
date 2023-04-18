#!/bin/bash

# gzip -dk GSM4156597_skin.late.anagen.atac.fragments.bed.gz
# sort -k1,1 -k2,2n -k3,3n GSM4156597_skin.late.anagen.atac.fragments.bed > "${output_dir}"/GSM4156597_skin.late.anagen.atac.fragments.sorted.bed
# bgzip "${output_dir}"/GSM4156597_skin.late.anagen.atac.fragments.sorted.bed
# tabix -p bed "${output_dir}"/GSM4156597_skin.late.anagen.atac.fragments.sorted.bed.gz



# Function to extract all rows from a bed that contains cell ID in a barcode file
function extract_rows_by_barcodes() {
  # Get the arguments
  local bed_file="$1"
  local cell_ids_file="$2"
  # Use grep to extract the rows that contain cell IDs from the cell IDs file
  grep -Fwf "$cell_ids_file" "$bed_file" 

}

output_dir=/home/siluo/projects/sc_chromatin_benchmark/cleaned_data/Ma_2019
bed_file_gz="${output_dir}"/GSM4156597_skin.late.anagen.atac.fragments.sorted.bed.gz
bed_file="${output_dir}"/GSM4156597_skin.late.anagen.atac.fragments.sorted.bed
cell_ids_file="/home/siluo/public/SiyuanLuo/projects/benchmark/scripts/data_cleaning/Ma_2019/Cells_Ma_2019.txt"

gzip -dk $bed_file_gz
echo $bed_file
echo $cell_ids_file
extract_rows_by_barcodes $bed_file $cell_ids_file > "${output_dir}"/GSM4156597_skin.late.anagen.atac.fragments.sorted.filtered.bed

# sort according to coordinate
sort -k1,1 -k2,2n -k3,3n "${output_dir}"/GSM4156597_skin.late.anagen.atac.fragments.sorted.filtered.bed > "${output_dir}"/GSM4156597_skin.late.anagen.atac.fragments.sorted.filtered.sorted.bed
# add count column (as 10X fragment file looks like)
awk 'BEGIN{FS=OFS="\t"} {print $0, 1}' "${output_dir}"/GSM4156597_skin.late.anagen.atac.fragments.sorted.filtered.sorted.bed > "${output_dir}/GSM4156597_skin.late.anagen.atac.fragments.sorted.filtered.sorted.1.bed"
# bgzip
bgzip "${output_dir}/GSM4156597_skin.late.anagen.atac.fragments.sorted.filtered.sorted.1.bed"
# tabix index
tabix -p bed "${output_dir}/GSM4156597_skin.late.anagen.atac.fragments.sorted.filtered.sorted.1.bed.gz"

cp "${output_dir}"/GSM4156597_skin.late.anagen.atac.fragments.sorted.filtered.sorted.1.bed.gz* /home/siluo/public/SiyuanLuo/projects/benchmark/cleaned_data/Ma_2019