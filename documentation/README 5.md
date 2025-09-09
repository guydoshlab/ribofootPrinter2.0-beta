
# Step by step guide on how to run the *posstats* Python script
This Python script takes a ROCC file and computes pause scores for every individual occurrence of a sequence feature of interest. ROCC files are generated from SAM files [here](https://github.com/kyrakerkhofs/mammalian_builddense_edits). 
![alt text](https://github.com/kyrakerkhofs/figures/blob/main/posstats.png)

# How to run script guide:
We provide an example on how to create a virtual python environment using anaconda [here](https://github.com/kyrakerkhofs/ribofootprinter_python_guide/tree/main). Other options are also available.

Create and navigate to the ribofootprinter folder. 
This folder contains:
1.	The reduced transcriptome file MANEv1.4_longnames.fasta
2.	ROCC files (generated from SAM files)

## Download python script to ribofootprinter folder
You can either download the Python script manually from this page or copy the entire folder into the directory as follows:
```unix
git clone https://github.com/kyrakerkhofs/posstats/
```
Navigate to this folder in the terminal
```unix
cd ./posstats
```

## Run script in terminal
The script is run in the terminal by calling python and providing the python script and settings as described below. An overview of the different settings is in Table format below.

```unix
python posstats.py "../80S_subset.rocc" "PPG" 1 0 10 1 12 1 "80S_posstats_PPG" > 80S_posstats_PPG_metadata.txt
```
```unix
python posstats.py "../80S_subset.rocc" "AAA" 1 0 10 1 12 1 "80S_posstats_AAA" > 80S_posstats_AAA_metadata.txt
```

![alt text](https://github.com/kyrakerkhofs/figures/blob/main/settings_posstats.png)
