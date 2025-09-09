# Step by step guide on how to run the *builddense* Python script to generate ROCC files
This Python script converts SAM files into ROCC files which are required for the following packages:
1. [writegene2](https://github.com/kyrakerkhofs/writegene2
)
2. [metagene](https://github.com/kyrakerkhofs/metagene
)
3. [genelist](https://github.com/kyrakerkhofs/genelist
)
4. [smorflist](https://github.com/kyrakerkhofs/smorflist
)
5. [posavg](https://github.com/kyrakerkhofs/posavg
)
6. [posstats](https://github.com/kyrakerkhofs/posstats
)

![alt text](https://github.com/kyrakerkhofs/figures/blob/main/ribofootprinter.png)

Preparation of the MANE transcriptome for alignment and ribofootprinter can be found [here](https://github.com/kyrakerkhofs/MANE_v1.4_Preparation).

A guide on how to view your aligned reads in IGV can be found [here](https://github.com/kyrakerkhofs/MANE_v1.4_IGV).


# How to run script guide:
We provide an example on how to create a virtual python environment using anaconda [here](https://github.com/kyrakerkhofs/ribofootprinter_python_guide/tree/main). Other options are also available.

Navigate to the ribofootPrinter2.0 folder. 
This folder contains:
1. SAM file containing bowtie aligned reads (single-end or paired-end)
2. Reduced transcriptome file MANEv1.4_longnames.fasta, download [here](https://github.com/kyrakerkhofs/MANE_v1.4_Preparation)

## Download python script to ribofootprinter folder
You can either download the Python script manually from this page or copy the entire folder into the directory as follows:
```unix
git clone https://github.com/kyrakerkhofs/builddense/
```
Navigate to this folder in the terminal
```unix
cd ./builddense
```

## Run script in terminal
The script is run in the terminal by calling python and providing the python script and settings as described below. An overview of the different settings is in Table format below.

This code will generate 5'-end mapped ROCC files from single-end reads
```unix
python builddense.py "./MANEv1.4_longnames.fasta" "./80S_subset.SAM" "80S_subset" -1 25 34 1 > 80S_builddense_metadata.txt
```
```unix
python builddense.py "./MANEv1.4_longnames.fasta" "./40S_subset.SAM" "40S_subset" -1 20 80 1 > 40S_builddense_metadata.txt
```
 

![alt text](https://github.com/kyrakerkhofs/figures/blob/main/settings_builddense.png)

If paired-end mRNA-seq data is used, we suggest using normalization setting 0 for coverage.

