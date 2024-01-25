library(tidyverse)
library(readxl)

set.seed(123)

outdir <- "../analyses/Cys"
if (!dir.exists(outdir)){
	dir.create(outdir, recursive = TRUE)
}

#functions
#perform confidence interval simulations

conf_int <- function(lambda, n=2000){

	values <- rpois(n=n, lambda=lambda)

	bounds <- quantile(values, c(0.025, 0.975)) %>% unname()

	return(list(lb = bounds[1], ub = bounds[2]))

}

seer.df <- "../data/seer-abundances.xlsx" %>% read_excel() %>% magrittr::set_colnames(c("rosetta", "cancer", "incidence")) %>% mutate(incidence_frac = incidence / 100)

all_cases <- "../results/all-muts/all-counts-matrix.txt" %>% read_table() %>% filter(Hugo_Symbol == "Total")


indiv_muts <- "../results/counts-v600-g12c-mutations.txt" %>% read_table()

#cases with at least 1 mutation counted amongst the set of mutations supplied initially (not used)
total_cases <- indiv_muts %>% filter(Hugo_Symbol == "Total")
#


indiv_muts <- indiv_muts %>% filter(Hugo_Symbol %in% c("BRAFp.Val600Glu", "BRAFp.Val600Lys", "KRASp.Gly12Cys")) 


all_cases_subset <- all_cases %>% select(intersect(colnames(all_cases), colnames(indiv_muts))) %>% relocate(colnames(indiv_muts))

indiv_muts <- rbind(indiv_muts, all_cases_subset)

indiv_muts_by_rosetta <- indiv_muts %>% filter(Hugo_Symbol == "Total")

indiv_muts.df <- indiv_muts %>% filter(Hugo_Symbol != "Total") %>% select(-All)

total <- indiv_muts_by_rosetta %>% pivot_longer(-c(Hugo_Symbol, All)) %>% left_join(seer.df, by=c("name"="rosetta")) 

total <- total %>% mutate(incidence = incidence / sum(incidence), incidence_frac = incidence_frac / sum(incidence_frac))

total <- total %>% magrittr::set_colnames(c("hugo", "all", "rosetta", "count", "cancer", "incidence", "incidence_frac"))


tcga.muts <- indiv_muts.df %>% pivot_longer(-Hugo_Symbol) %>% left_join(total %>% select(-c(count, hugo, all)), by=c("name"="rosetta"))

tcga.muts <- tcga.muts %>% mutate(gene = gsub("\\.(.*)", "", Hugo_Symbol)) %>% magrittr::set_colnames(c("hugo", "rosetta", "count", "cancer", "incidence", "incidence_frac", "gene"))

tcga.gene <- tcga.muts %>% group_by(gene, cancer) %>% summarize(net_count=sum(count), incidence_frac = unique(incidence_frac), gene = unique(gene), rosetta = unique(rosetta)) %>% select(rosetta, incidence_frac, gene, net_count)


tcga.gene <- tcga.gene %>% mutate(weight.gene = incidence_frac * net_count)
total <- total %>% mutate(weight.tot = incidence_frac * count)
tcga.muts <- tcga.muts %>% mutate(weight.mut=incidence_frac * count) 

tcga.muts <- tcga.muts %>% mutate(idx = 1:n())
tcga.gene <- tcga.gene %>% ungroup() %>% mutate(idx = 1:n())


tcga.muts.confint <- tcga.muts %>% filter(count > 0)  %>% mutate(result = purrr::map(count, conf_int), lb = purrr::map_dbl(result, "lb"), ub = purrr::map_dbl(result, "ub")) %>% select(-result)
tcga.muts <- tcga.muts %>% left_join(tcga.muts.confint %>% select(idx, lb, ub), by=c("idx"="idx"))


tcga.gene.confint <- tcga.gene %>% filter(net_count > 0)  %>% mutate(result = purrr::map(net_count, conf_int), lb = purrr::map_dbl(result, "lb"), ub = purrr::map_dbl(result, "ub")) %>% select(-result)

tcga.gene <- tcga.gene %>% left_join(tcga.gene.confint %>% select(idx, lb, ub), by=c("idx"="idx"))

tcga.gene <- tcga.gene %>% mutate(weight.lb = lb * incidence_frac, weight.ub = ub * incidence_frac)

tcga.muts <- tcga.muts %>% mutate(weight.lb = lb * incidence_frac, weight.ub = ub * incidence_frac)

tcga.gene <- tcga.gene %>% mutate(weight.lb = ifelse(is.na(weight.lb), 0, weight.lb), weight.ub = ifelse(is.na(weight.ub), 0, weight.ub))

tcga.muts <- tcga.muts %>% mutate(weight.lb = ifelse(is.na(weight.lb), 0, weight.lb), weight.ub = ifelse(is.na(weight.ub), 0, weight.ub))


mut <- tcga.muts %>% left_join(total %>% select(rosetta, weight.tot, all), by=c("rosetta"="rosetta")) %>% group_by(hugo) %>% summarize(gene=unique(gene), pct_us = sum(weight.mut) / sum(weight.tot), pct_tcga = sum(count) / unique(all) , pct_lb = sum(weight.lb) / sum(weight.tot), pct_ub = sum(weight.ub) / sum(weight.tot)) %>% mutate(across(where(is.double), ~ . * 100))


gene <- tcga.gene %>% left_join(total %>% select(rosetta, weight.tot, all), by=c("rosetta"="rosetta")) %>% group_by(gene) %>% summarize(pct_us = sum(weight.gene) / sum(weight.tot), pct_tcga = sum(net_count) / unique(all), pct_lb = sum(weight.lb) / sum(weight.tot), pct_ub = sum(weight.ub) / sum(weight.tot) ) %>% mutate(across(where(is.double), ~ . * 100))


df <- mut %>% mutate(gene = str_sub(gene, 1, -2)) %>% mutate(mutation = gsub("(.*)\\.", "", hugo)) %>% relocate(gene, mutation) %>% select(-hugo) %>% mutate(label_text = paste0(gene, "\n(", mutation, ")")) %>% select(-c(gene, mutation))

df <- df %>% select(label_text, pct_us, pct_tcga) %>% arrange(desc(pct_us)) %>% magrittr::set_colnames(c("mutation", "EG", "NPC")) %>% pivot_longer(-mutation) 

mut.levels <- df$mutation %>% unique()
df$mutation <- factor(df$mutation, levels = mut.levels)

fig1b <- df %>% ggplot(aes(x=mutation, y=value, fill=name, group=name)) + geom_bar(stat = "identity", position="dodge") +  scale_fill_manual(values = c("EG"="maroon", "NPC"="navy"), labels = c("EG"="EG (epidemiologic-genomic) Rate", "NPC"="NPC (naive pan-cancer) Rate"), guide = guide_legend(nrow = 2)) + theme_minimal() + theme(panel.background = element_blank(), panel.border = element_rect(color = "black", fill = "transparent"), panel.grid = element_blank()) + theme(axis.text = element_text(size =12), axis.title = element_text(size = 16)) + xlab("") + ylab("Mutation Rate (%)") + theme(axis.ticks =  element_line(color = "black")) + theme(legend.title = element_blank()) + theme(legend.text = element_text(size = 14), legend.position = "bottom")  + scale_y_continuous(breaks = scales::pretty_breaks(n=5))



#Fig A
total_sequenced_cases <- all_cases$All
case_counts <- all_cases %>% pivot_longer(-Hugo_Symbol) %>% filter(name != "All")

case.df <- case_counts %>% mutate(seq_frac = value / total_sequenced_cases) %>% left_join(., seer.df %>% select(rosetta, cancer, incidence), by=c("name"="rosetta")) %>% select(name, cancer, incidence, seq_frac) %>% mutate(seq_pct = seq_frac * 100) %>% select(-seq_frac) %>% mutate(diff = abs(incidence - seq_pct)) %>% arrange(desc(diff)) %>% slice(1:20) %>% arrange(desc(incidence))

cancer_names <- case.df %>% pull(cancer)
cancer_names[12] <- "DLBCL"
cancer_names[13] <- "CLL/SLL"
cancer_names[17] <- "Renal Cell Carcinoma (Chromophobe)"
cancer_names[18] <- "Neuroblastoma (NOS)"
cancer_names[14] <- "Glioblastoma (NOS)"
cancer_names[15] <- "AML (NOS)"
cancer_names[19] <- "Ewing Sarcoma"

case.df <- case.df %>% mutate(cancer = cancer_names)
case.df$cancer <- factor(case.df$cancer, levels = case.df$cancer)

fig1a <- case.df %>% select(cancer, incidence, seq_pct) %>% magrittr::set_colnames(c("cancer", "SEER Incidence", "Sequenced Cases")) %>% pivot_longer(-cancer) %>% mutate(value = as.numeric(value)) %>% ggplot(aes(x=cancer, y=value, fill=name, group=name)) + geom_bar(stat = "identity", position="dodge") + scale_fill_manual(values = c("#E69F00", "#56B4E9"), guide = guide_legend(nrow = 1)) + theme_minimal() + theme(panel.background = element_blank(), panel.border = element_rect(color = "black", fill = "transparent"), panel.grid = element_blank()) + theme(axis.text = element_text(size =13), axis.title = element_text(size = 16)) + xlab("") + ylab("Proportion (%)") + theme(axis.ticks =  element_line(color = "black")) + theme(legend.title = element_blank()) + theme(legend.text = element_text(size = 14), legend.position = "bottom") + scale_y_continuous(breaks = scales::pretty_breaks(n=8)) + theme(axis.text.x = element_text(size = 12, angle = 60, hjust = 1))




cowplot::plot_grid(fig1a, fig1b, nrow = 1, rel_widths = c(2.4, 1), labels = c("A", "B"), label_size = 25) %>% ggsave(., filename = file.path(outdir, "fig1.pdf"), height = 8, width = 13, device = cairo_pdf, units = "in")






#