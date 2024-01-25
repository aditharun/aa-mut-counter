#compute analysis values used in paper 

library(tidyverse)
library(readxl)

set.seed(123)

outdir <- "../analyses/Cys"
if (!dir.exists(outdir)){
	dir.create(outdir, recursive = TRUE)
}

