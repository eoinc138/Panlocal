Execution of the Panlocal script requires initial assembly of the master dataset for either the individual phylums or the full concatenated dataset.
1) Obtain the reference biosamples list by downloading the reference proteins file from Progenomes, available here: https://progenomes.embl.de/data/repGenomes/progenomes3.proteins.representatives.fasta.bz2.
2) In your working directory execute biosample_list_executable.sh to obtain the text file of reference biosamples.
3) Obtain the reference gene list by downloading the reference gene file from Progenomes, available here: https://progenomes.embl.de/data/repGenomes/progenomes3.genes.representatives.fasta.bz2.
4) In your working directory execute gene_list_executable.sh to obtain the text file of reference biosample genes.
5) Ensure both the reference biosamples and reference gene list text files are present in your working directory for the subsequent steps.

6) Depending on the dataset you wish to assemble download the following EGGNOG and Gene annotation lists below. Ensure a space of at least 3TB is available in your directory if you intend to use all 5 phylums.
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
   ./Euryarchaeota_Bacteroidetes_executable.sh Biosample_IDs_filtered.txt Euryarchaeota.eggNOG_groups.tsv phylum-Euryarchaeota.gene_annotations.tsv Euryarchaeota.tsv
   Alter the input as appropriate for Bacteroidetes
7) In the working directory with the annotation file for Firmicutes, Actinobacteria and Proteobacteria run the Proteobacteria_Actinobacteria_Firmicutes_executable.sh as the following input
8) ./Proteobacteria_Actinobacteria_Firmicutes_executable.sh Biosample_IDs_filtered.txt Proteobacteria.eggNOG_groups.tsv phylum-Proteobacteria.gene_annotations.tsv GeneIDS.txt Proteobacteria.tsv
9) Alter the input as appropriate for Firmicutes and Actinobacteria.
