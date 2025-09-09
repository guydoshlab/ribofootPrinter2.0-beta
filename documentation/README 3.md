
# Step by step guide on how to run the *metagene* Python script
This Python script takes a ROCC file and calculates the average around start or stop codons. ROCC files are generated from SAM files [here](https://github.com/kyrakerkhofs/mammalian_builddense_edits). 
![alt text](https://github.com/kyrakerkhofs/figures/blob/main/metagene.png)

A metagene plot averages profiling data around start or stop codons:
![alt text](https://github.com/kyrakerkhofs/figures/blob/main/metagene_script.png)


# How to run script guide:
We provide an example on how to create a virtual python environment using anaconda [here](https://github.com/kyrakerkhofs/ribofootprinter_python_guide/tree/main). Other options are also available.

Create and navigate to the ribofootprinter folder. 
This folder contains:
1.	The reduced transcriptome file MANEv1.4_longnames.fasta
2.	ROCC files (generated from SAM files)
3.	OPTIONAL: A genelist containing a selection of genes of interest. For example a selection of uORF identified transcripts. 

## Download python script to ribofootprinter folder
You can either download the Python script manually from this page or copy the entire folder into the directory as follows:
```unix
git clone https://github.com/kyrakerkhofs/metagene/
```
Navigate to this folder in the terminal
```unix
cd ./metagene
```


## Run script in terminal
The script is run in the terminal by calling python and providing the python script and settings as described below. An overview of the different settings is in Table format below.


This code will generate a start codon metagene for all transcripts. 
```unix
python metagene.py "../80S_subset.rocc" 1 1 5 50 300 "none" "./80S_metagene_start_equalweight" > 80S_metagene_start_metadata.txt
```

This code will generate a stop codon metagene for all transcripts. 
```unix
python metagene.py "../80S_subset.rocc" 2 1 5 300 50 "none" "./80S_metagene_stop_equalweight" > 80S_metagene_stop_metadata.txt
```

This code will generate a start codon metagene for a subset of the dataset containing predicted uORF containing transcripts. 
```unix
python metagene.py "../80S_subset.rocc" 1 1 5 50 300 "../subset_list.xlsx" "./80S_metagene_start_equalweight_subset" > 80S_metagene_start_subset_metadata.txt
```

![alt text](https://github.com/kyrakerkhofs/figures/blob/main/settings_metagene.png)

Note that if the UTRs are shorter than the 5' or 3' range settings, these transcripts will be excluded from the analysis. Information of number of transcripts included in the analysis can be found as an output in the terminal.
