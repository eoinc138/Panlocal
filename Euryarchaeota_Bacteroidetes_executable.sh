#!/bin/bash

#Checking for the number of arguments provided.
if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <Biosample_IDs_filtered> <eggNOG_groups> <gene_annotations> <output_prefix>"
    exit 1
fi

#The arguments/files are assigned to variables for later calling.
Biosample_IDs_filtered=$1
eggNOG_groups=$2
gene_annotations=$3
output_prefix=$4

#Filtering the eggNOG_groups using the Biosample IDs.
awk 'NR==FNR {Biosample_IDs_filtered[$1]; next} $1 in Biosample_IDs_filtered' "$Biosample_IDs_filtered" "$eggNOG_groups" > "${output_prefix}.eggNOG_groups_filtered.tsv"

#Removing unneccessary columns from the dataset eggNOG groups.
cut -f1-6,9-13,15,18-20 "${output_prefix}.eggNOG_groups_filtered.tsv" > "${output_prefix}.dataset.tsv"

#Removing unneccessary columns from the gene annotation files.
cut -f1-2 "$gene_annotations" > "${output_prefix}.gene_annotations.filtered.tsv"

#Joining the gene annotations file to the dataset file to give our final file.
awk -F'\t' 'NR==FNR{a[$1]=$2; next} {print $0, ($2 in a ? a[$2] : "NA")}' OFS='\t' "${output_prefix}.gene_annotations.filtered.tsv" "${output_prefix}.dataset.tsv" > "${output_prefix}.tsv"

#Removing temporary files from the directory
rm "${output_prefix}.gene_annotations.filtered.tsv"
rm "${output_prefix}.eggNOG_groups_filtered.tsv"
rm "${output_prefix}.dataset.tsv"
