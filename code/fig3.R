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

points_mut <- mut %>% mutate(diff=abs(pct_tcga - pct_us)) %>% arrange(desc(pct_us)) %>% mutate(us_rank = 1:n()) %>% mutate(gene_fix = str_sub(gene, 1, -2), mut_fix = paste0(" (", gsub("(.*)\\.", "", hugo) %>% str_sub(., 1, -4) %>% str_replace(., "([[:alpha:]])(\\d)", "\\1 \\2"), ")")) %>% mutate(label_text = paste0(gene_fix, mut_fix) ) %>% mutate(label = ifelse(us_rank <= 7, label_text, NA)) %>% select(-label_text) %>% filter(!is.na(label)) %>% mutate(position = ifelse(us_rank == 5, "right", "left")) %>% mutate(nx = ifelse(us_rank == 5, -0.02, 0.02)) %>% mutate(ny = ifelse(us_rank == 7, -0.01, 0.02)) %>% mutate(label = ifelse(us_rank == 5, paste0(gene_fix, "\n", mut_fix), label) ) 

mut_plot <- mut %>% ggplot(aes(x=pct_tcga, y=pct_us)) + theme_minimal() + geom_errorbar(aes(x=pct_tcga, ymin = pct_lb, ymax = pct_ub), size = 0.02, width = 0, color = "grey70") + geom_point(alpha=0.8, size=1.75, color="black")  + geom_abline(color="grey70", size=1, alpha=0.8, linetype="dashed", slope = 1, intercept = 0) + theme(panel.grid = element_blank(), panel.border = element_rect(color="black", fill="transparent")) + scale_x_continuous(limits=c(0,2), breaks=seq(0, 2, 0.5)) + scale_y_continuous(limits=c(0,2), breaks=seq(0, 2, 0.5)) + ylab("Estimed Epidemio-Genomic (EG)\nMutation Rate (%)") + xlab("Estimated Naive Pan-Cancer (NPC) Mutation Rate (%)") + theme(axis.title=element_text(size=16), axis.text=element_text(size=14), axis.ticks=element_line(color="black")) + geom_text(data=points_mut, mapping=aes(x=pct_tcga+nx, y=pct_us+ny, label=label, hjust = position), lineheight = 0.65)

no_ci_mut_plot <- mut %>% ggplot(aes(x=pct_tcga, y=pct_us)) + theme_minimal()  + geom_point(alpha=0.8, size=1.75, color="black")  + geom_abline(color="grey70", size=1, alpha=0.8, linetype="dashed", slope = 1, intercept = 0) + theme(panel.grid = element_blank(), panel.border = element_rect(color="black", fill="transparent")) + scale_x_continuous(limits=c(0,2), breaks=seq(0, 2, 0.5)) + scale_y_continuous(limits=c(0,2), breaks=seq(0, 2, 0.5)) + ylab("Estimed Epidemio-Genomic (EG)\nMutation Rate (%)") + xlab("Estimated Naive Pan-Cancer (NPC) Mutation Rate (%)") + theme(axis.title=element_text(size=16), axis.text=element_text(size=14), axis.ticks=element_line(color="black")) + geom_text(data=points_mut, mapping=aes(x=pct_tcga+nx, y=pct_us+ny, label=label, hjust = position), lineheight = 0.65)


mut_plot %>% ggsave(., filename = file.path(outdir, "fig3.pdf"), height = 8, width = 8, device = cairo_pdf, units = "in")

no_ci_mut_plot %>% ggsave(., filename = file.path(outdir, "fig3_no_ci.pdf"), height = 8, width = 8, device = cairo_pdf, units = "in")

#0.772 correlation w/ 95% CI: 0.769 to 0.7745, p-value < 2.2e-16
#spearman: 0.42, p-val < 2.2e-16