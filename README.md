### Epidemiologically-informed estimates of mutation rates in cancer
#### Adith S. Arun, David Liarakos, Gaurav Mendiratta, Edward C. Stites

There are two main parts to generating these estimates: 1) generating mutation by ROSETTA code count matrices, and 2) perform mutation rate reweighting with SEER

#### How to Use

##### Generating mutation by ROSETTA code matrices

cBioPortal exome sequencing data is found in `data/cbioportal_raw_EXOME139` which is downloadable using `data/cbioportal_raw_EXOME139/Raw_Studies_Downloader.py` and the associated clinical data with ROSETTA code mapped for each patient is found for each study in the folder `data/clinical`. 

For all python scripts, make sure that the following packages are installed: `os, multiprocess, pandas, numpy, argparse`. 

To generate a count matrix (gene by ROSETTA code) for all mutations, navigate to this directory and run `python3 generate-entire-matrix.py`. The count matrix is output at `results/all-muts/all-counts-matrix.txt`. Note that the folder `data/reference_data_files` contains helper files necessary for counting mutations. 

Now, to generate a count matrix (mutation by ROSETTA code) for mutations for a specific amino acid `Xyz`, run `python3 generate-aa-specific-matrix.py --aa XYZ` where `Xyz` is the 3 letter amino acid code with first letter capitalized. This creates a count matrix at `results/counts-Xyz.txt`. 
