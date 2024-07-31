#!/bin/bash

#Unzip the file into your directory
gzip -d NOG.annotations.tsv.gz

#Remove unneccessary columns from file
cut -f2,6 NOG_annotations.tsv > EGGNOG_IDS.tsv

#Remove unneeded files from directory
rm NOG_annotations.tsv
