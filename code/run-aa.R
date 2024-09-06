#stitch analyses, figures, and tables from scratch

library(tidyverse)

lookuptable <- "amino-acid-lookup.csv" %>% read_csv()

lookuptable <- lookuptable %>% filter(`amino acid` == "Cys")

for (x in 1:nrow(lookuptable)){

    amino_acid <- lookuptable$`amino acid`[x]

    print(amino_acid)

    #system(paste0("python3 generate-aa-specific-matrix.py --aa ", amino_acid))

    system(paste0("Rscript stitch-figures.R -a ", amino_acid))

    system(paste0("Rscript subset-weighting.R ", amino_acid, " subset-file.csv"))

    system(paste0("Rscript subset-weighting.R ", amino_acid, " lung-subset-file.csv"))
    system(paste0("Rscript subset-weighting.R ", amino_acid, " ac-subset-file.csv"))
    system(paste0("Rscript subset-weighting.R ", amino_acid, " mm-subset-file.csv"))
    system(paste0("Rscript subset-weighting.R ", amino_acid, " scc-subset-file.csv"))
    system(paste0("Rscript subset-weighting.R ", amino_acid, " tcc-subset-file.csv"))

}

