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

#percentage of driver genes outside of 95% CI and percentage of non-drivers outside of CI?

genetiers <- file.path("../data", "Lists_CT1and2.txt") %>% read_tsv()

tierdf <- mut %>% mutate(driver = ifelse(gene %in% paste0(genetiers$Consensus_Tier1, "p"), TRUE, FALSE))
tierdf <- tierdf %>% mutate(inside = ((pct_tcga <= pct_ub) & (pct_tcga >= pct_lb)))

round((1 - (tierdf %>% filter(driver) %>% pull(inside) %>% mean()))*100, 2) #rate of driver genes outside CI 

round((1 - (tierdf %>% filter(!driver) %>% pull(inside) %>% mean()))*100, 2) #rate of nondriver genes outside CI



# mutation level
 
n.outside.tier1 <- (tierdf %>% filter(driver) %>% nrow()) - (tierdf %>% filter(driver) %>% pull(inside) %>% sum())

n.inside.tier1 <-  (tierdf %>% filter(driver) %>% pull(inside) %>% sum())

n.outside.complement <- (tierdf %>% filter(!driver) %>% nrow()) - (tierdf %>% filter(!driver) %>% pull(inside) %>% sum())
n.inside.complement <-  (tierdf %>% filter(!driver) %>% pull(inside) %>% sum())

df <- data.frame(outside = c(n.outside.complement, n.outside.tier1), inside = c(n.inside.complement, n.inside.tier1))

fisher.test(df)


# gene level
gene <- file.path(outdir, "gene-rates.rds") %>% readRDS()

tierdf <- gene %>% mutate(driver = ifelse(gene %in% paste0(genetiers$Consensus_Tier1, "p"), TRUE, FALSE))
tierdf <- tierdf %>% mutate(inside = ((pct_tcga <= pct_ub) & (pct_tcga >= pct_lb)))

round((1 - (tierdf %>% filter(driver) %>% pull(inside) %>% mean()))*100, 2) #rate of driver genes outside CI 

round((1 - (tierdf %>% filter(!driver) %>% pull(inside) %>% mean()))*100, 2) #rate of nondriver genes outside CI


n.outside.tier1 <- (tierdf %>% filter(driver) %>% nrow()) - (tierdf %>% filter(driver) %>% pull(inside) %>% sum())

n.inside.tier1 <-  (tierdf %>% filter(driver) %>% pull(inside) %>% sum())

n.outside.complement <- (tierdf %>% filter(!driver) %>% nrow()) - (tierdf %>% filter(!driver) %>% pull(inside) %>% sum())
n.inside.complement <-  (tierdf %>% filter(!driver) %>% pull(inside) %>% sum())

df <- data.frame(outside = c(n.outside.complement, n.outside.tier1), inside = c(n.inside.complement, n.inside.tier1))

fisher.test(df)



















##

seer.df <- "../data/seer-abundances.xlsx" %>% read_excel() %>% magrittr::set_colnames(c("rosetta", "cancer", "incidence")) %>% mutate(incidence_frac = incidence / 100)

all_cases <- "../results/all-muts/all-counts-matrix.txt" %>% read_table() %>% filter(Hugo_Symbol == "Total")

total_sequenced_cases <- all_cases$All
case_counts <- all_cases %>% pivot_longer(-Hugo_Symbol) %>% filter(name != "All")

case.df <- case_counts %>% mutate(seq_frac = value / total_sequenced_cases) %>% left_join(., seer.df %>% select(rosetta, cancer, incidence), by=c("name"="rosetta")) %>% select(name, cancer, incidence, seq_frac) %>% mutate(seq_pct = seq_frac * 100) %>% select(-seq_frac) %>% mutate(diff = abs(incidence - seq_pct)) %>% arrange(desc(diff)) 
