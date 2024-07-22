#!/bin/bash

#Ensuring the arguments required for the script are present.
if [ "$#" -ne 5 ]; then
    echo "Usage: $0 <Biosample_IDs_filtered> <eggNOG_groups> <gene_annotations> <GeneIDs> <output_prefix>"
    exit 1
fi

#Required arguments are assigned to variables.
Biosample_IDs_filtered=$1
eggNOG_groups=$2
gene_annotations=$3
GeneIDs=$4
output_prefix=$5

#Filtering the EGGNOG file based on the Biosample.
awk 'NR==FNR {Biosample_IDs_filtered[$1]; next} $1 in Biosample_IDs_filtered' "$Biosample_IDs_filtered" "$eggNOG_groups" > "${output_prefix}.eggNOG_groups_filtered.tsv"

#Unneccessary columns are removed from the EGGNOG file.
cut -f1-6,9-13,15,18-20 "${output_prefix}.eggNOG_groups_filtered.tsv" > "${output_prefix}.dataset.tsv"

#Similarly unneccessary columns are removed from the gene_annotations file.
cut -f1-2 "$gene_annotations" > "${output_prefix}.gene_annotations.filtered.tsv"

#To reduce size the gene annotations file is trimmed using the Gene IDs file to reduce size.
awk 'NR==FNR {GeneIDS_formatted[$1]; next} $1 in GeneIDS_formatted' "$GeneIDs" "${output_prefix}.dataset.tsv" > "${output_prefix}.dataset_parsed.tsv"

#The gene annotation file is then joined to the EGGNOG annotation file for its' contig ID.
awk -F'\t' 'NR==FNR {a[$1]=$2; next} {print $0, ($2 in a ? a[$2] : "NA")}' OFS='\t' "${output_prefix}.gene_annotations.filtered.tsv" "${output_prefix}.dataset_parsed.tsv" > "${output_prefix}.tsv"

#Removing temporary files from the Directory.
rm "${output_prefix}.gene_annotations.filtered.tsv"
rm "${output_prefix}.eggNOG_groups_filtered.tsv"
rm "${output_prefix}.dataset_parsed.tsv"
