# Step by step guide on how to run the *genelist* Python script
This Python script takes a ROCC file and counts reads that map to ORF and UTR regions of gene models. It also calculates open reading frames on individual transcripts. ROCC files are generated from SAM files [here](https://github.com/kyrakerkhofs/mammalian_builddense_edits). 
![alt text](https://github.com/kyrakerkhofs/figures/blob/main/genelist.png)

The 5'-end aligned data is shifted to accomodate the P-site of the ribosome for example. We used a shift of 12 for our riboseq dataset in this study. 
![alt text](https://github.com/kyrakerkhofs/figures/blob/main/shift.png)


# How to run script guide:
We provide an example on how to create a virtual python environment using anaconda [here](https://github.com/kyrakerkhofs/ribofootprinter_python_guide/tree/main). Other options are also available.

Create and navigate to the ribofootprinter folder. 
This folder contains:
1.	The reduced transcriptome file MANEv1.4_longnames.fasta
2.	ROCC files (generated from SAM files)

## Download python script to ribofootprinter folder
You can either download the Python script manually from this page or copy the entire folder into the directory as follows:
```unix
git clone https://github.com/kyrakerkhofs/genelist/
```
Navigate to this folder in the terminal
```unix
cd ./genelist
```


## Run script in terminal
The script is run in the terminal by calling python and providing the python script and settings as described below. An overview of the different settings is in Table format below.

```unix
python genelist.py "../80S_subset.rocc" 12 1 "80S_genelist" > 80S_genelist_metadata.txt
```
```unix
python genelist.py "../40S_subset.rocc" 12 1 "40S_genelist" > 40S_genelist_metadata.txt
```



![alt text](https://github.com/kyrakerkhofs/figures/blob/main/settings_genelist.png)

Note that genes with 5'-UTRs shorter than the shift value are excluded from the analysis.
