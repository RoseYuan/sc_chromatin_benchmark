#!/bin/bash

# Fragment file is already filtered, this will disable snaptools to generate quality control metrics.
# decompress the gz file

fragfile=$1
genome=$2
genome_sizefile=$3
binsize=$4
data_dir=$5

#file_prefix=$(basename $fragfile .tsv.gz)
#file_prefix=$(basename $fragfile .bed.gz)
file_prefix=$(basename $fragfile .gz)

suffix=".gz"
unziped_fragfile=${fragfile%"$suffix"}

sorted_file=${data_dir}/${file_prefix}_sorted.bed
compressed_sorted_file=${data_dir}/${file_prefix}_sorted.bed.gz
output_snapfile=${data_dir}/${file_prefix}.snap

# decompress
if [ ! -f "$unziped_fragfile" ]; then 
    gunzip -k $fragfile || gunzip $fragfile; 
fi


if [ ! -f "$$sorted_file" ]; then 
    # sort the tsv file using the 4th column (barcode column)
    sort --parallel=4 -k4,4 $unziped_fragfile > $sorted_file
    # replace ":" in barcode
    sed -i 's/:/+/g' $sorted_file  # TODO
fi


# compress the bed file 
if [ ! -f "$sorted_file.gz" ]; then 
    bgzip -k $sorted_file || gzip $sorted_file; 
    tabix -p bed "$sorted_file.gz"
fi

if [ ! -f "$output_snapfile" ]; then 
    # run snap files using the bed file
    snaptools snap-pre  \
      --input-file=$compressed_sorted_file  \
      --output-snap=$output_snapfile  \
      --genome-name=$genome  \
      --genome-size=$genome_sizefile  \
      --min-mapq=30  \
      --min-flen=0  \
      --max-flen=1000  \
      --keep-chrm=TRUE  \
      --keep-single=TRUE  \
      --keep-secondary=False  \
      --overwrite=True  \
      --verbose=False
    
    # create the cell-by-bin matrix. Can specify multiple bin size. The cell-by-bin matrices will be added to the snap file without creating another file. Same with snap-add-pmat and snap-add-gmat.
    snaptools snap-add-bmat  \
    --snap-file=$output_snapfile  \
    --bin-size-list $binsize  \
    --verbose=False
fi


  
