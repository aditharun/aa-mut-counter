library(tidyverse)
library(readxl)

set.seed(123)

args = commandArgs(trailingOnly=TRUE)
amino_acid <- args[1]

outdir <- paste0("../analyses/", amino_acid)

if (!dir.exists(outdir)){
	dir.create(outdir, recursive = TRUE)
}

datadir <- "../data/seer-abundances.xlsx"

#load Seer data
#same as SuppTable1_ROSETTA_Abundance_added12.xlsx in Stites 2021
seer.df <- datadir %>% read_excel() %>% magrittr::set_colnames(c("rosetta", "cancer", "incidence")) %>% mutate(incidence_frac = incidence / 100)


#by mutation
aa_tcga <- paste0("../results/counts-", amino_acid, ".txt") %>% read_table()

#by gene (not important for getting total cases)
all_cases <- "../results/all-muts/all-counts-matrix.txt" %>% read_table() %>% filter(Hugo_Symbol == "Total")

#get rid of <10 cases
less_than_10_cases <- all_cases[-c(1:2)] %>% t() %>% as.data.frame() %>% filter(V1 < 10) %>% rownames(.)

seer.df <- seer.df %>% filter(!rosetta %in% less_than_10_cases)
all_cases <- all_cases[which(!colnames(all_cases) %in% less_than_10_cases)]
aa_tcga <- aa_tcga[which(!colnames(aa_tcga) %in% less_than_10_cases)]


all_cases %>% select(-2) %>% pivot_longer(-Hugo_Symbol) %>% select(-Hugo_Symbol) %>% left_join(seer.df %>% select(rosetta, cancer), by=c("name"="rosetta")) %>% write_csv(., "../count-per-code.csv")

total_cases <- aa_tcga %>% filter(Hugo_Symbol == "Total")

all_cases_subset <- all_cases %>% select(intersect(colnames(all_cases), colnames(aa_tcga))) %>% relocate(colnames(aa_tcga))

#we use this sum as our total because we threw away < 10 cases
npc_total <- (all_cases[,-c(1:2)] %>% t() %>% as.numeric() %>% sum())

all_cases_subset <- all_cases_subset %>% mutate(All = npc_total)

aa_tcga <- rbind(aa_tcga %>% filter(Hugo_Symbol != "Total"), all_cases_subset)


aa_tcga_by_mut <- aa_tcga %>% select(All)

aa_tcga_by_rosetta <- aa_tcga %>% filter(Hugo_Symbol == "Total")


aa_tcga.df <- aa_tcga %>% filter(Hugo_Symbol != "Total") %>% select(-All)


### WRITE OUT FILE ####
#What fraction of cases for a given ROSETTA code have at least 1 mutation of interest?
total_cases$Hugo_Symbol <- paste0("cases with >= 1 acquired ", amino_acid)
all_cases_subset$Hugo_Symbol <- "total sequenced cases"
colnames(all_cases_subset)[1] <- "status"
colnames(total_cases)[1] <- "status"
rbind(total_cases, all_cases_subset) %>% as_tibble() %>% writexl::write_xlsx(., file.path(outdir, "cases-by-rosetta.xlsx"))


#renormalize weight vector to 1
total <- aa_tcga_by_rosetta %>% pivot_longer(-c(Hugo_Symbol, All)) %>% left_join(seer.df, by=c("name"="rosetta")) 


total <- total %>% mutate(incidence = incidence / sum(incidence), incidence_frac = incidence_frac / sum(incidence_frac))


total <- total %>% magrittr::set_colnames(c("hugo", "all", "rosetta", "count", "cancer", "incidence", "incidence_frac"))

saveRDS(total, file.path(outdir, "total-inc-frac.rds"))

tcga.muts <- aa_tcga.df %>% pivot_longer(-Hugo_Symbol) %>% left_join(total %>% select(-c(count, hugo, all)), by=c("name"="rosetta"))

tcga.muts <- tcga.muts %>% mutate(gene = gsub("\\.(.*)", "", Hugo_Symbol)) %>% magrittr::set_colnames(c("hugo", "rosetta", "count", "cancer", "incidence", "incidence_frac", "gene"))

#tcga.gene <- tcga.muts %>% group_by(gene, cancer) %>% summarize(net_count=sum(count), incidence_frac = unique(incidence_frac), gene = unique(gene), rosetta = unique(rosetta)) %>% select(rosetta, incidence_frac, gene, net_count)

#

tcga.muts <- tcga.muts %>% left_join(total %>% select(rosetta, count) %>% magrittr::set_colnames(c("rosetta", "tcount")), by=c("rosetta"="rosetta"))

tcga.muts <- tcga.muts %>% mutate(frac_count = count / tcount)

#tcga.gene <- tcga.gene %>% left_join(total %>% select(rosetta, count) %>% magrittr::set_colnames(c("rosetta", "tcount")), by=c("rosetta"="rosetta"))

#tcga.gene <- tcga.gene %>% mutate(frac_count = net_count / tcount)

#tcga.gene <- tcga.gene %>% mutate(weight.gene = incidence_frac * frac_count)
tcga.muts <- tcga.muts %>% mutate(weight.mut=incidence_frac * frac_count) 

tcga.muts <- tcga.muts %>% arrange(desc(count))
#tcga.gene <- tcga.gene %>% arrange(desc(net_count))

tcga.muts <- tcga.muts %>% mutate(idx = 1:n())
#tcga.gene <- tcga.gene %>% ungroup() %>% mutate(idx = 1:n())

nonzeromut <- tcga.muts %>% filter(count > 0) %>% pull(idx) %>% max()
#nonzerogene <- tcga.gene %>% filter(net_count > 0) %>% pull(idx) %>% max()

colnames(tcga.muts)[which(colnames(tcga.muts)=="hugo")] <- "id"

#tcga.gene <- tcga.gene %>% magrittr::set_colnames(c("rosetta", "incidence_frac", "id", "count", "tcount", "frac_count", "weight.gene", "idx"))

saveRDS(tcga.muts, file.path(outdir, "tcga-muts.rds"))



mut <- tcga.muts %>% left_join(total %>% select(rosetta, all), by=c("rosetta"="rosetta")) %>% group_by(id) %>% summarize(gene=unique(gene), pct_us = sum(weight.mut), pct_tcga = sum(count) / unique(all)) %>% mutate(across(where(is.double), ~ . * 100))

# %>% left_join(confint.muts, by=c("id"="id")) %>% mutate(across(where(is.double), ~ . * 100))

#gene <- tcga.gene %>% left_join(total %>% select(rosetta, all), by=c("rosetta"="rosetta")) %>% group_by(id) %>% summarize(pct_us = sum(weight.gene), pct_tcga = sum(count) / unique(all)) %>% left_join(confint.gene, by=c("id"="id")) 

colnames(mut)[which(colnames(mut) == "id")] <- "hugo"
#colnames(gene)[which(colnames(gene) == "id")] <- "gene"

mut$lb <- NA
mut$ub <- NA

colnames(mut)[which(colnames(mut)=="lb")] <- "pct_lb"
colnames(mut)[which(colnames(mut)=="ub")] <- "pct_ub"

#colnames(gene)[which(colnames(gene)=="lb")] <- "pct_lb"
#colnames(gene)[which(colnames(gene)=="ub")] <- "pct_ub"

saveRDS(mut, file.path(outdir, "mut-rates.rds"))
#saveRDS(gene, file.path(outdir, "gene-rates.rds"))

colnames(tcga.muts)[which(colnames(tcga.muts)=="id")] <- "hugo"
#which cancers contribute to the mutation rates for each mutation? which cancers have the most mutations of a given type (e.g., kras g12c)

cancer_contr <- tcga.muts %>% select(hugo, rosetta, cancer, incidence, count) %>% distinct() %>% left_join(., total %>% select(cancer, count, rosetta) %>% magrittr::set_colnames(c("cancer", "count_cancer", "rosetta")), by=c("cancer"="cancer", "rosetta"="rosetta"))

#get number of mutations / cases w at least 1 mutation for each cancer code
#the total column is the same as aa_tcga_by_mut variable defined above which is correct (for each mutation, how many cases have at least 1 mutation of that type)
ranked_pct_contr <- cancer_contr %>% group_by(hugo) %>% mutate(total = sum(count), frac = count / total) %>% arrange(desc(frac))

#rank mutations by naive pan-cancer mutation rates
ranked_muts <- mut %>% mutate(diff=abs(pct_tcga - pct_us)) %>% arrange(desc(pct_us)) %>% mutate(us_rank = 1:n()) 

mut_cancer_df <- left_join(ranked_pct_contr, ranked_muts %>% select(pct_us, hugo, gene), by=c("hugo"="hugo")) %>% filter(!is.na(pct_us)) 

### WRITE OUT FILE ###
#for each mutation, in which cancer do we find mutations of this sort?
mut_cancer_table <- mut_cancer_df %>% ungroup() %>% select(hugo, rosetta, cancer, count, count_cancer, total, pct_us, gene, frac) %>% magrittr::set_colnames(c("mutation", "rosetta", "cancer", "count_mutation_for_given_rosetta", "total_cancer_sequenced_cases", "total_mutation_for_mutation_across_rosetta", "informed_rate", "gene", "frac")) %>% mutate(cancer = str_replace_all(cancer, "_", " ")) %>% mutate(gene_fix = str_sub(gene, 1, -2), mut_fix = paste0(" (", gsub("(.*)\\.", "", mutation) %>% str_sub(., 1, -4) %>% str_replace(., "([[:alpha:]])(\\d)", "\\1 \\2"), ")")) %>% mutate(label_text = paste0(gene_fix, mut_fix) ) %>% select(-c(mutation, gene_fix, mut_fix))

mut_cancer_table <- mut_cancer_table %>% mutate(cancer = ifelse(cancer == "Lung Small Cell Carcionoma", "Lung Small Cell Carcinoma", cancer))

saveRDS(mut_cancer_table, file.path(outdir, "mut-cancer-contributions.rds"))




#

