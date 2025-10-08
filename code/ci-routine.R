library(tidyverse)

set.seed(123)

args = commandArgs(trailingOnly=TRUE)
amino_acid <- args[1]

tcga.muts <- readRDS(file.path("../analyses", amino_acid, "tcga-muts.rds"))


quanum <- tcga.muts %>% select(cancer, rosetta, incidence_frac, tcount) %>% mutate(quanum = incidence_frac / tcount) %>% select(cancer, quanum) %>% distinct()

tcga.muts <- tcga.muts %>% select(id, count, cancer)
tcga.muts <- tcga.muts %>% pivot_wider(names_from = cancer, values_from = count)

write_csv(tcga.muts, file.path("../analyses", amino_acid, "tcga-muts.csv"))
write_csv(quanum, file.path("../analyses", amino_acid, "quanum.csv"))

system(paste0("python3 generate-ci.py ", amino_acid))

tcgamutsfile <- file.path("../analyses", amino_acid, "tcga-muts.csv")
outputfile <- file.path("../analyses", amino_acid, "main-ci.csv")

mutratefile <- file.path("../analyses", amino_acid, "mut-rates.rds")

cifile <- file.path("../analyses", amino_acid, "main-ci.csv")


ci <- paste0(cifile) %>% read_csv(., col_names = FALSE)

mut <- paste0(mutratefile) %>% readRDS()

colnames(ci) <- c("hugo", "lb", "ub")

ci <- ci %>% mutate(pct_lb = lb*100, pct_ub = ub*100) %>% select(-c(ub, lb))

mut <- mut %>% select(-c(pct_lb, pct_ub))

df <- mut %>% left_join(ci, by=c("hugo"="hugo"))

saveRDS(object = df, file = mutratefile)

outfile <- file.path("../analyses", amino_acid, "mutation-rates.csv")

write_csv(x = df, file = outfile)
