library(tidyverse)
library(readxl)

set.seed(123)

args = commandArgs(trailingOnly=TRUE)
amino_acid <- args[1]

outdir <- paste0("../results/")

if (!dir.exists(outdir)){
	dir.create(outdir, recursive = TRUE)
}

if (!dir.exists(file.path("../analyses", amino_acid))){
	dir.create(file.path("../analyses", amino_acid), recursive = TRUE)
}

df <- paste0("../results/counts-", amino_acid, ".txt") %>% read_tsv()

df <- df %>% mutate(aa_change = gsub("(.*)\\.", "", Hugo_Symbol)) %>% relocate(aa_change)

indels <- df %>% filter(grepl("ins|del", aa_change)) 

df <- df %>% filter(!grepl("ins|del", aa_change)) 

df <- df %>% select(-aa_change)

dffile <- paste0(outdir, "counts-", amino_acid, ".txt")

indelfile <- paste0("../analyses/", amino_acid, "/indel-", amino_acid, ".csv")

df %>% write_tsv(., file = dffile)

indels %>% write_csv(., file = indelfile)