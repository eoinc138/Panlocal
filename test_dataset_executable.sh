#!/bin/bash

#Script exits in the event of failure.
set -e

#First column is removed from the master dataset.
cut -f 1 Final_dataset.tsv > Sample.txt

#Duplicate entries are removed from the first column text file
uniq Sample.txt > Sample1.txt

#These 3 commands count the total number of lines, calculate 5% of them and then shuffle the text file to select 5% at random.
total_lines=$(wc -l < Sample1.txt)
num_lines_to_select=$(awk "BEGIN {print int($total_lines * 0.05)}")
shuf -n $num_lines_to_select Sample1.txt > output.txt

#The dataset is then parsed using the random 5% to create a smaller test dataset
awk 'NR==FNR {output[$1]; next} $1 in output' output.txt Final_dataset.tsv > Final_test_dataset.tsv

#Interim files are deleted to clear directory
rm Sample.txt
rm Sample1.txt
rm output.txt
