#!/bin/bash

#Download the Cog and arCog data
wget https://ftp.ncbi.nlm.nih.gov/pub/COG/COG2020/data/cog-20.def.tab
wget https://ftp.ncbi.nlm.nih.gov/pub/wolf/COGs/arCOG/ar14.arCOGdef.tab

#Combine the two files
cat ar14.arCOGdef.tab cog-20.def.tab > combinedfile.tab

#Rename file for easier use later
mv combinedfile.tab combinedfile.tsv

#Cut unneccessary columns
cut -f1,3,4 combinedfile.tsv > COGnarCOG_descriptors.tsv

#Remove unneeded files
rm cog-20.def.tab ar14.arCOGdef.tab

