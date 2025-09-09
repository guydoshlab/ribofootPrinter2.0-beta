# Step by step guide on how to run the *3D_metagene* Python script
This script takes a bowtie aligned SAM file (single-end reads) and generates start or stop 3D metagene plots from ribosome profiling datasets.
![alt text](https://github.com/kyrakerkhofs/figures/blob/main/3D_metagene.png)

The script includes a number of normalization options described below.
If a list of genes is provided, this script can be run on a subset of genes. An example xlsx file is provided in this Github page.


# How to run script guide:
We provide an example on how to create a virtual python environment using anaconda [here](https://github.com/kyrakerkhofs/ribofootprinter_python_guide/tree/main). Other options are also available. This package has dependencies on:
1. biopython
2. pandas
3. openpyxl
4. matplotlib (optional for generating plot)

Create and navigate to the ribofootprinter folder. 
This folder contains:
1.	The reduced transcriptome file MANEv1.4_longnames.fasta
2.	Transcriptome-aligned SAM files
3.	OPTIONAL: A genelist containing a selection of genes of interest. For example a selection of uORF identified transcripts. 

## Download python script to ribofootprinter folder
You can either download the Python script manually from this page or copy the entire folder into the directory as follows:
```unix
git clone https://github.com/kyrakerkhofs/3D_metagene/
```
Navigate to this folder in the terminal
```unix
cd ./3D_metagene
```


## Run script in terminal
The script is run in the terminal by calling python and providing the python script and settings as described below. An overview of the different settings is in Table format below.

This code will analyze the start codon metagene. 
```unix
python 3D_metagene.py ../MANEv1.4_longnames.fasta ../80S_subset.SAM ./80S_3D_start.csv none 25 34 100 100 1
```
This code will analyze the stop codon metagene. 
```unix
python 3D_metagene.py ../MANEv1.4_longnames.fasta ../80S_subset.SAM ./80S_3D_stop.csv none 25 34 100 100 2
```


![alt text](https://github.com/kyrakerkhofs/figures/blob/main/settings_3D_metagene.png)


## OPTIONAL: Generate plots
The Python script outputs a csv file which can be used to generate plots using the package matplotlib. Only a single argument is needed for this script.
```unix
python 3D_metagene_plot.py ./80S_3D_start.csv
```
