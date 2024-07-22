#!/bin/bash

#Unzip the file into your directory
gzip -d 1_annotations.tsv.gz

#Remove unneccessary columns from file
cut -f2,4 1_annotations.tsv > EGGNOG_freetext_description.tsv

#Remove unneeded files from directory
rm 1_annotations.tsv
