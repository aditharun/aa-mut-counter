library(tidyverse)
library(readxl)

set.seed(123)

outdir <- "../analyses"
if (!dir.exists(outdir)){
	dir.create(outdir)
}

