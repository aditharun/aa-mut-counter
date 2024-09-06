library(tidyverse)

aa <-  "amino-acid-lookup.csv" %>% read_csv()

for (x in 1:nrow(aa)){

	system(paste0("Rscript stitch-figures.R -a ", aa$`amino acid`[x], " -l ", aa$long[x]))

}