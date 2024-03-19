#stitch analyses, figures, and tables from scratch

library(tidyverse)

lookuptable <- "amino-acid-lookup.csv" %>% read_csv()

for (x in 1:nrow(lookuptable)){

  generated_status <- lookuptable$generated[x]
  amino_acid <- lookuptable$`amino acid`[x]
  amino_acid_long <- lookuptable$long[x]

  if (is.na(generated_status)){

    system(paste0("python3 generate-aa-specific-matrix.py --aa ", amino_acid))
    system(paste0("Rscript stitch-figures.R -a ", amino_acid, " -l ", amino_acid_long))

  }

}


