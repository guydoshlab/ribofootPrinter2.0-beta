# Optional guide on how to run *ribofootPrinter's* Python scripts locally
This repository provides an optional step-by-step guide on how to run the Python toolbox ribofootPrinter locally. This includes setting up your python environment and installation of dependencies. Multiple other methods are also available.
A summary of the different ribofootPrinter packages are listed below and hyperlinks to their corresponding Github page with how to run these specific packages are provided.

The following packages use SAM files as input files (ROCC file generation not required):
1. [region_size_and_abundance](https://github.com/kyrakerkhofs/region_size_and_abundance
)
2. [metagene_3D](https://github.com/kyrakerkhofs/3D_metagene
)

All the other packages require a ROCC file as input file. The ROCC file is generated from a SAM file using the [builddense](https://github.com/kyrakerkhofs/builddense
) script.
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

![alt text](https://github.com/kyrakerkhofs/ribofootprinter_guide/blob/main/ribofootprinter.png)

Other steps summarized in the figure above are found elsewhere: 

Preparation of the MANE transcriptome for alignment and ribofootprinter can be found [here](https://github.com/kyrakerkhofs/MANE_v1.4_Preparation).

A guide on how to view your aligned reads in IGV can be found [here](https://github.com/kyrakerkhofs/MANE_v1.4_IGV).


# Setting up the environment:
This can be done locally using a virtual environment, such as Anaconda, but not required. 

Reference for Anaconda: https://github.com/conda/conda?tab=readme-ov-file.

Details for Anaconda follow.

## Generate conda environment (only once)
```unix
conda create --name python-env
```
## Activate the environment (each time)
```unix
conda activate python-env
```
## Install Python (only once)
This installs version 3.9. 
```unix
conda install python=3.9
```
# Install BioPython (only once)
```unix
pip install biopython
```
# Install openpyxl (only once)
```unix
pip install openpyxl
```

# OPTIONAL: Install pandas, numpy and matplotlib (only once)
Two ribofootPrinter packages (region_size_and_abundance.py and 3D_metagene) also have dependencies on pandas and numpy. If figures need to be generated, matplotlib should also be installed.

```unix
pip install pandas
pip install matplotlib
```

## Create folder where *ribofootPrinter's* Github packages containing python scripts will be downloaded. 
This folder also contains the following files:
1. SAM files which will serve as input files.
2. The fully annotated transcriptome longnames FASTA file, download [here](https://github.com/kyrakerkhofs/MANE_v1.4_Preparation).

```unix
mkdir -p ribofootprinter
```

Navigate to the ribofootprinter folder. 
```unix
cd ./ribofootprinter
```


Each package in ribofootPrinter contains an individual Github page with a how to run guide (see links above). The code is run by passing relevant parameters to the Python script using the command line. We recommend capturing the metadata output including error messages and the read-in parameters in a separate txt file. 
