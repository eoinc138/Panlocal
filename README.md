Execution of the Panlocal script requires initial assembly of the master dataset for either the individual phylums or the full concatenated dataset.
## Creating representative genome annotation database
1) Obtain the reference biosamples list by downloading the reference proteins file from Progenomes, available here:
```
wget "https://progenomes.embl.de/data/repGenomes/progenomes3.proteins.representatives.fasta.bz2"
```
3) In your working directory execute biosample_list_executable.sh to obtain the text file of reference biosamples.
4) Obtain the reference gene list by downloading the reference gene file from Progenomes, available here:
```
wget "https://progenomes.embl.de/data/repGenomes/progenomes3.genes.representatives.fasta.bz2"
```
6) In your working directory execute
```
sh gene_list_executable.sh
```
to obtain the text file of reference biosample genes.
10) Ensure both the reference biosamples and reference gene list text files are present in your working directory for the subsequent steps.

11) Depending on the dataset you wish to assemble download the following EGGNOG and Gene annotation lists below. Ensure a space of at least 3TB is available in your directory if you intend to use all 5 phylums.
  The gene annotations are available at the following links
  Proteobacteria: https://progenomes.embl.de/dumpAnnotation.cgi?p=Proteobacteria&t=ga&a=phylum
  Firmicutes: https://progenomes.embl.de/dumpAnnotation.cgi?p=Firmicutes&t=ga&a=phylum
  Bacteroidetes: https://progenomes.embl.de/dumpAnnotation.cgi?p=Bacteroidetes&t=ga&a=phylum
  Actinobacteria: https://progenomes.embl.de/dumpAnnotation.cgi?p=Actinobacteria&t=ga&a=phylum
  Euryarchaeota: https://progenomes.embl.de/dumpAnnotation.cgi?p=Euryarchaeota&t=ga&a=phylum

  The EGGNOG annotation files are available from the following links
  Proteobacteria: https://progenomes.embl.de/dumpAnnotation.cgi?p=Proteobacteria&t=ae&a=phylum
  Firmicutes: https://progenomes.embl.de/dumpAnnotation.cgi?p=Firmicutes&t=ae&a=phylum
  Bacteroidetes: https://progenomes.embl.de/dumpAnnotation.cgi?p=Bacteroidetes&t=ae&a=phylum
  Actinobacteria: https://progenomes.embl.de/dumpAnnotation.cgi?p=Actinobacteria&t=ae&a=phylum
  Euryarchaeota: https://progenomes.embl.de/dumpAnnotation.cgi?p=Euryarchaeota&t=ae&a=phylum

6) In the working directory with the annotation files for Euryarchaeota and Bacteroidetes run the Euryarchaeota_Bacteroidetes_executable.sh as the following input
```
sh Euryarchaeota_Bacteroidetes_executable.sh Biosample_IDs_filtered.txt Euryarchaeota.eggNOG_groups.tsv phylum-Euryarchaeota.gene_annotations.tsv Euryarchaeota.tsv
```
Alter the input as appropriate for Bacteroidetes
8) In the working directory with the annotation file for Firmicutes, Actinobacteria and Proteobacteria run the Proteobacteria_Actinobacteria_Firmicutes_executable.sh as the following input
```
./Proteobacteria_Actinobacteria_Firmicutes_executable.sh Biosample_IDs_filtered.txt Proteobacteria.eggNOG_groups.tsv phylum-Proteobacteria.gene_annotations.tsv GeneIDS.txt Proteobacteria.tsv
```
10) Alter the input as appropriate for Firmicutes and Actinobacteria.
11) To concatenate the phylums execute the Finished_dataset_processing.sh script in the working directory with any of the 5 phyla present or only one if you intend to examine one phylum.
```
sh Finished_dataset_processing.sh Euryarchaeota.eggNOG_groups.tsv Euryarchaeota.tsv
```
12) To prepare the GTDB file obtain the most recent Archaea and Bacterial taxonomy data from GTDB, available here:
```
wget "https://data.gtdb.ecogenomic.org/releases/latest/bac120_metadata.tsv.gz"
wget "https://data.gtdb.ecogenomic.org/releases/latest/ar53_metadata.tsv.gz"
```
15) In the same directory as these files execute
```
GTDB_executable.sh
```
17) To prepare the additional EGGNOG file obtain the annotations file for EGGNOG V4.5 here 
```
wget "http://eggnog45.embl.de/download/eggnog_4.5/data/NOG/NOG.annotations.tsv.gz"
```
18) In the same directory as this file execute
```
sh EGGNOG_descriptions_executable.sh
```
19) Additional python modules are required to run Panlocal. These are as follows; numpy, matplotlib, collections, csv, os, re, pandas, base64, io, anytree and pdfkit. It is recommended to install pdfkit using !pip as it is known to run into errors when installed in a conda environment. 
Lastly before running the Panlocal script as an executable ensure the directory you are using contains the Master dataset tsv, EGGNOG_freetext_description.tsv and the GTDB_taxonomy.tsv files.
## Usage
```
python Panlocal.py
```
