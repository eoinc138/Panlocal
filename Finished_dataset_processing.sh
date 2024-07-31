#!/bin/bash

#One EGGNOG groups file required for the header, and any subsequent number of files for further merging.
if [ "$#" -lt 2 ]; then
  echo "Usage: $0 <header_file> <tsv_file1> [<tsv_file2> ... <tsv_fileN>]"
  exit 1
fi

#Keep the header file separate from the subsequent merging
header_origin="$1"
shift

#Merge the files, skipping the header for all but the first file
awk 'FNR==1 && NR!=1 {next} {print}' "$@" > master_dataset.tsv

#Sorting the master dataset by the biosample column to ensure genes are in correct order.
sort -k1 master_dataset.tsv > sorted_dataset.tsv

#An additional column comprising the Biosample ID in GTDB format is added.
awk -F'\t' '{
  split($1, a, ".");
  new_col = a[2];
  for (i=3; i<=length(a); i++) {
    new_col = new_col "." a[i];
  }
  $0 = $0 "\t" new_col;
  print
}' sorted_dataset.tsv > dataset_noheader.tsv

#The header is produced from the Euryarchaeota eggNOG file, it is ran through a similar process 
# and has an additional two columns added for the Contig ID and GTDB biosample.
head -n 1 "$header_origin" > dataset_header.tsv
cut -f1-6,9-13,15,18-20 dataset_header.tsv > dataset_headerv1.tsv
awk -F'\t' -v OFS='\t' '{print $0, "CONTIG_ID", "GTDB_BIOSAMPLE"}' dataset_headerv1.tsv > dataset_headerv2.tsv

#The header is then added to complete dataset.
cat dataset_headerv2.tsv dataset_noheader.tsv > Finished_dataset.tsv

#Remove interim files to clean the directory.
rm sorted_dataset.tsv
rm master_dataset.tsv
rm dataset_noheader.tsv
rm dataset_header.tsv
rm dataset_headerv1.tsv
rm dataset_headerv2.tsv
