### Epidemiologically-informed estimates of mutation rates in cancer
#### Adith S. Arun, David Liarakos, Gaurav Mendiratta, Edward C. Stites

There are two main parts to generating these estimates: 1) generating mutation by ROSETTA code count matrices, and 2) perform mutation rate reweighting with SEER

#### How to Use

##### Generating mutation by ROSETTA code matrices

The `data/` folder is not tracked on GitHub because it is too large (~20GB). It is available as a .zip file on Zenodo (here)[]. Download it from the link and unzip it in the `aa-mut-counter` folder. 

cBioPortal exome sequencing data is found in `data/cbioportal_raw_EXOME139` which is downloadable using `data/cbioportal_raw_EXOME139/Raw_Studies_Downloader.py` and the associated clinical data with ROSETTA code mapped for each patient is found for each study in the folder `data/clinical`. 

For all python scripts, make sure that the following packages are installed: `os, multiprocess, pandas, numpy, argparse`. 

To generate a count matrix (gene by ROSETTA code) for all mutations, navigate to this directory and run `python3 generate-entire-matrix.py`. The count matrix is output at `results/all-muts/all-counts-matrix.txt`. Note that the folder `data/reference_data_files` contains helper files necessary for counting mutations. 

Now, to generate a count matrix (mutation by ROSETTA code) for mutations for a specific amino acid `Xyz`, run `python3 generate-aa-specific-matrix.py --aa XYZ` where `Xyz` is the 3 letter amino acid code with first letter capitalized. This creates a count matrix at `results/counts-Xyz.txt`. 

Now, suppose you want to generate counts for a list of specific mutations. Define this list as a txt file in `code/tmp/{filename}.txt` with each line being a mutation in the format `p.{Xyz}{#}{Abc}` (e.g., p.Val600Glu) separated by a new line with no other punctuation. Then, go to `code/` and run `python3 generate-specific-muts-matrix.py --mutations_file tmp/{filename}.txt` which will generate a count matrix at `results/counts-{filename}.txt`. 


##### Performing Epidemio-Genomic reweighting

Go to `code/` and run `Rscript reweighting.R` and it will create a folder `../results/Cys/` with the results, and master tables. Figures and other materials can be run by `Rscript fig1.R`, `Rscript fig2.R`, `Rscript fig3.R` and `Rscript create-tables.R`. 

Please reach out at adith.3.arun@gmail.com with questions or comments. 
