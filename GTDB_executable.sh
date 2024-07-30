#!/bin/bash

gzip -d bac120_metadata.tsv.gz
gzip -d ar53_metadata.tsv.gz

cut -f 20,53 bac120_metadata.tsv > bac_GTDB.tsv
cut -f 20,53 ar53_metadata.tsv > arc_GTDB.tsv

cat bac_GTDB.tsv arc_GTDB.tsv > GTDB.tsv

grep -v "gtdb_taxonomy" GTDB.tsv > GTDB_noheader.tsv

awk -F'\t' '{
 n = split($1, arr, ";")
 printf $2
 for (i = 1; i <= n; i++) {
     printf "\t" arr[i]
 }
 for (i = n+1; i <= NF; i++) {
     printf "\t" $(i+1)
  }
  print ""
}' GTDB_noheader.tsv > GTDB_split_taxonomy_data.tsv

awk 'BEGIN{FS=OFS="\t"} {for(i=2; i<=8; i++) $i=substr($i,4)} 1' GTDB_split_taxonomy_data.tsv > GTDB_taxonomy_noheader.tsv

echo -e "ncbi_biosample\tDomain\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies" | cat - GTDB_taxonomy_noheader.tsv > GTDB_taxonomy.tsv

rm bac_GTDB.tsv
rm arc_GTDB.tsv
rm GTDB.tsv
rm GTDB_noheader.tsv
rm GTDB_split_taxonomy_data.tsv
rm GTDB_taxonomy_noheader.tsv
