# Step by step guide on how to run the *region_size_and_abundance* Python script
This Python script takes a bowtie aligned SAM file (single-end reads) and determines the length of mapped reads for each feature:
1. 5'-UTR
2. start codon
3. open-reading frame (ORF)
4. stop codon
5. 3'-UTR

The *window* parameter can be adjusted (standard setting = 4). This setting determines the regions around the start and stop codon in which a mapped read needs to overlap: 
![alt text](https://github.com/kyrakerkhofs/figures/blob/main/region_size_and_abundance.png)

The script includes a number of normalization options described below.
If a list of genes is provided, this script can be run on a subset of genes.


# How to run script guide:
We provide an example on how to create a virtual python environment using anaconda [here](https://github.com/kyrakerkhofs/ribofootprinter_python_guide/tree/main). Other options are also available. This package has dependencies on:
1. biopython
2. pandas
3. matplotlib (optional for generating plot)

Create and navigate to the ribofootprinter folder. 
This folder contains:
1.	The reduced transcriptome file MANEv1.4_longnames.fasta, download [here](https://github.com/kyrakerkhofs/MANE_v1.4_Preparation)
2.	Transcriptome-aligned SAM files
3.	OPTIONAL: A genelist containing a selection of genes of interest. For example a selection of uORF identified transcripts. 

## Download python script to ribofootprinter folder
You can either download the Python script manually from this page or copy the entire folder into the directory as follows:
```unix
git clone https://github.com/kyrakerkhofs/region_size_and_abundance/
```
Navigate to this folder in the terminal
```unix
cd ./region_size_and_abundance
```

## Run script in terminal
The script is run in the terminal by calling python and providing the python script and settings as described below. An overview of the different settings is in Table format below.

This code will analyze all transcripts. 
```unix
python region_size_and_abundance.py ../MANEv1.4_longnames.fasta ../80S_subset.SAM 80S_region_size 25 34 4 none
```
This code will analyze a subset of the dataset for predicted uORF containing transcripts. 
```unix
python region_size_and_abundance.py ../MANEv1.4_longnames.fasta ../80S_subset.SAM 80S_region_size_subset 25 34 4 ./subset_list.xlsx
```

![alt text](https://github.com/kyrakerkhofs/figures/blob/main/settings_region_size.png)


