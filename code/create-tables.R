library(tidyverse)
library(readxl)

set.seed(123)

outdir <- "../analyses/Cys"
if (!dir.exists(outdir)){
	dir.create(outdir, recursive = TRUE)
}


mut_cancer_table <- file.path(outdir, "mut-cancer-contributions.rds")

top3_cancers_per_mutation <- mut_cancer_table %>% select(label_text, cancer, frac, gene, informed_rate)  %>% mutate(indiv_imp = paste0(cancer, " (", round(frac*100, 2), "%)")) %>% select(-c(cancer, frac)) %>% group_by(informed_rate, label_text) %>% mutate(row_id = row_number()) %>% pivot_wider(names_from = row_id, values_from = indiv_imp, names_prefix = "V") %>% mutate(cancers = paste0(V1, ", ", V2, ", ", V3)) %>% select(-starts_with("V")) %>% arrange(desc(informed_rate)) %>% relocate(label_text) %>% ungroup() %>% select(-gene) %>% magrittr::set_colnames(c("Point Mutation", "Mutation Proportion (%)", "The most common cancers to harbor a mutation in this gene (ROSETTA classification, % of all instances of this mutation)"))

writexl::write_xlsx(top3_cancers_per_mutation %>% slice(1:25), file.path(outdir, "table1.xlsx"))

