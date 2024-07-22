#!/bin/bash

#Read the representative gene file and only extract the FASTA headers.
bzcat progenomes3.genes.representatives.fasta.bz2 | pv | bzgrep '^>' > FASTA_HEADERS.txt

#Remove everything after the < character in the FASTA headers.
sed 's/<*//' FASTA_HEADERS.txt 

#Remove the first character and everything the the < onwards from each entry
cut -c2- FASTA_HEADERS.TXT | sed 's/<.*//' > GeneIDS.txt

#Remove intermediate files from the directory
rm FASTA_HEADERS.txt

