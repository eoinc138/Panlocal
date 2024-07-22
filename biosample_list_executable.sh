#!/bin/bash

#Read the representative proteins file and only extract the FASTA headers.
bzcat progenomes3.proteins.representatives.fasta.bz2 | pv | bzgrep '^>' > FASTA_HEADERS.txt>

#Remove everything after the gene ID in the FASTA headers.
sed 's/.GCA*//' FASTA_HEADERS.txt 

#Sort and remove all duplicates from the Biosample text file.
sort FASTA_HEADERS.txt | uniq > Biosample_IDs_filtered.txt

#Remove intermediate files from the directory.
rm FASTA_HEADERS.txt

