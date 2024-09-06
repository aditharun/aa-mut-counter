library(tidyverse)
library(readxl)

args = commandArgs(trailingOnly=TRUE)
amino_acid <- args[1]

set.seed(123)

outdir <- paste0("../analyses/", amino_acid)

if (!dir.exists(outdir)){
	dir.create(outdir, recursive = TRUE)
}


mut_cancer_table <- file.path(outdir, "mut-cancer-contributions.rds") %>% readRDS()

top3_cancers_per_mutation <- mut_cancer_table %>% select(label_text, cancer, frac, gene, informed_rate)  %>% mutate(indiv_imp = paste0(cancer, " (", round(frac*100, 2), "%)")) %>% select(-c(cancer, frac)) %>% group_by(informed_rate, label_text) %>% mutate(row_id = row_number()) %>% pivot_wider(names_from = row_id, values_from = indiv_imp, names_prefix = "V") %>% mutate(cancers = paste0(V1, ", ", V2, ", ", V3)) %>% select(-starts_with("V")) %>% arrange(desc(informed_rate)) %>% relocate(label_text) %>% ungroup() %>% select(-gene) %>% magrittr::set_colnames(c("Point Mutation", "Mutation Proportion (%)", "The most common cancers to harbor a mutation in this gene (ROSETTA classification, % of all instances of this mutation)"))

writexl::write_xlsx(top3_cancers_per_mutation %>% slice(1:25), file.path(outdir, "Table1.xlsx"))


mut <- file.path(outdir, "mut-rates.rds") %>% readRDS()
gene <- file.path(outdir, "gene-rates.rds") %>% readRDS()

muts.df <- mut %>% mutate(diff=(pct_us - pct_tcga)) %>% arrange(desc(pct_us)) %>% mutate(gene = str_sub(gene, 1, -2)) %>% mutate(aa_change = gsub("(.*)\\.", "", hugo)) %>% mutate(mutation = paste0(gene, ".", aa_change)) %>% select(-hugo) %>% relocate(mutation, gene, aa_change, pct_tcga, pct_us, pct_lb, pct_ub, diff) 

genes.df <- gene %>% mutate(diff=(pct_us - pct_tcga)) %>% arrange(desc(abs(diff))) %>% mutate(gene = str_sub(gene, 1, -2)) %>% relocate(gene, pct_tcga, pct_us, pct_lb, pct_ub) 


muts.df <- muts.df %>% filter(pct_us > 0 | pct_tcga > 0)

muts.df <- muts.df %>% arrange(desc(pct_us))

muts.df <- muts.df %>% magrittr::set_colnames(c("Mutation", "Gene", "Amino Acid Change", "Naive Pan-Cancer (NPC) Rate (%)", "Epidemio-Genomic (EG) Rate (%)", "Lower Bound EG (%)", "Upper Bound EG (%)", "EG - NPC"))


genes.df <- genes.df %>% filter(pct_us > 0 | pct_tcga > 0)

genes.df <- genes.df %>% arrange(desc(pct_us))

genes.df <- genes.df %>% magrittr::set_colnames(c("Gene", "Naive Pan-Cancer (NPC) Rate (%)", "Epidemio-Genomic (EG) Rate (%)", "Lower Bound EG (%)", "Upper Bound EG (%)", "EG - NPC"))

writexl::write_xlsx(muts.df, file.path(outdir, "Supplemental Table 1 (no rounding).xlsx"))


muts.df <- muts.df %>% mutate_if(is.numeric, ~signif(., 4))
genes.df <- genes.df %>% mutate_if(is.numeric, ~signif(., 4))


writexl::write_xlsx(muts.df, file.path(outdir, "Supplemental Table 1.xlsx"))

writexl::write_xlsx(genes.df, file.path(outdir, "supplemental-gene-rates.xlsx"))
