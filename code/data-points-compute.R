#compute analysis values used in paper 

library(tidyverse)
library(readxl)

set.seed(123)

args = commandArgs(trailingOnly=TRUE)
amino_acid <- args[1]

outdir <- paste0("../analyses/", amino_acid)

if (!dir.exists(outdir)){
	dir.create(outdir, recursive = TRUE)
}

mut <- file.path(outdir, "mut-rates.rds") %>% readRDS()

seer.df <- "../data/seer-abundances.xlsx" %>% read_excel() %>% magrittr::set_colnames(c("rosetta", "cancer", "incidence")) %>% mutate(incidence_frac = incidence / 100)

all_cases <- "../results/all-muts/all-counts-matrix.txt" %>% read_table() %>% filter(Hugo_Symbol == "Total")

total_sequenced_cases <- all_cases$All
case_counts <- all_cases %>% pivot_longer(-Hugo_Symbol) %>% filter(name != "All")

case.df <- case_counts %>% mutate(seq_frac = value / total_sequenced_cases) %>% left_join(., seer.df %>% select(rosetta, cancer, incidence), by=c("name"="rosetta")) %>% select(name, cancer, incidence, seq_frac) %>% mutate(seq_pct = seq_frac * 100) %>% select(-seq_frac) %>% mutate(diff = abs(incidence - seq_pct)) %>% arrange(desc(diff)) 
