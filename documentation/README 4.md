# Step by step guide on how to run the *writegene2* Python script
This Python script takes a ROCC file and determines the abundance of 5'- or 3'-end mapped reads at a certain location within the transcript. ROCC files are generated from SAM files [here](https://github.com/kyrakerkhofs/mammalian_builddense_edits). 
![alt text](https://github.com/kyrakerkhofs/figures/blob/main/writegene2.png)

The script tolerates genenames (*e.g.* ENSG00000111640.15) or aliases (*e.g.* GAPDH) as inputs.

# How to run script guide:
We provide an example on how to create a virtual python environment using anaconda [here](https://github.com/kyrakerkhofs/ribofootprinter_python_guide/tree/main). Other options are also available.

Create and navigate to the ribofootprinter folder. 
This folder contains:
1.	The reduced transcriptome file MANEv1.4_longnames.fasta
2.	ROCC files (generated from SAM files)


## Download python script to ribofootprinter folder
You can either download the Python script manually from this page or copy the entire folder into the directory as follows:
```unix
git clone https://github.com/kyrakerkhofs/writegene2/
```
Navigate to this folder in the terminal
```unix
cd ./writegene2
```

## Run script in terminal
The script is run in the terminal by calling python and providing the python script and settings as described below. An overview of the different settings is in Table format below.

This script will output data for ACTB.
```unix
python writegene2.py "../80S_subset.rocc" "ACTB" "80S_writegene2_ACTB" > 80S_writegene2_ACTB_metadata.txt
```
This script will output data for uORF-containing EIF4G2.
```unix
python writegene2.py "../80S_subset.rocc" "EIF4G2" "80S_writegene2_EIF4G2" > 80S_writegene2_EIF4G2_metadata.txt
```


![alt text](https://github.com/kyrakerkhofs/figures/blob/main/settings_writegene2.png)

