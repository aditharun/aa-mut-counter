library(tidyverse)
library(readxl)

set.seed(123)

args = commandArgs(trailingOnly=TRUE)
amino_acid <- args[1]

outdir <- paste0("../analyses/", amino_acid)

if (!dir.exists(outdir)){
	dir.create(outdir, recursive = TRUE)
}

gene <- file.path(outdir, "gene-rates.rds") %>% readRDS()


df <- gene %>% mutate(gene = str_sub(gene, 1, -2)) %>% filter(gene %in% c("BRAF", "FGFR3", "TP53", "IDH1", "GNAS", "FBXW7", "CTNNB1", "DNMT3A")) %>% magrittr::set_colnames(c("gene", "EG", "NPC", "LB", "UB"))

fig2a <- df %>% select(gene, EG, NPC) %>% pivot_longer(-gene) %>% mutate(value = as.numeric(value), gene = factor(gene, levels = df %>% mutate(EG = as.numeric(EG)) %>% arrange(desc(EG)) %>% pull(gene))) %>% ggplot(aes(x=gene, y=value, fill=name, group=name)) + geom_bar(stat = "identity", position="dodge") + scale_fill_manual(values = c("EG"="maroon", "NPC"="navy"), labels = c("EG"="EG (epidemiologic-genomic) Rate", "NPC"="NPC (naive pan-cancer) Rate"), guide = guide_legend(nrow = 2)) + theme_minimal() + theme(panel.background = element_blank(), panel.border = element_rect(color = "black", fill = "transparent"), panel.grid = element_blank()) + theme(axis.text = element_text(size =13), axis.title = element_text(size = 16)) + xlab("Gene") + ylab("Mutation Rate (%)") + theme(axis.ticks =  element_line(color = "black")) + theme(legend.title = element_blank()) + theme(legend.text = element_text(size = 14), legend.position = "bottom") + scale_y_continuous(breaks = scales::pretty_breaks(n = 5))


fig2b <- df %>% type.convert(as.is = TRUE) %>% ggplot(aes(x=NPC, y=EG)) + geom_point(size = 2.5)+ theme_minimal() + theme(panel.background = element_blank(), panel.border = element_rect(color = "black", fill = "transparent"), panel.grid = element_blank()) + theme(axis.text = element_text(size =13), axis.title = element_text(size = 16)) + xlab("NPC Mutation Rate (%)") + ylab("EG Mutation Rate (%)") + theme(axis.ticks =  element_line(color = "black")) + theme(legend.title = element_blank()) + theme(legend.text = element_text(size = 14), legend.position = "top") + geom_abline(slope = 1, intercept = 0, color = "grey70", linetype = "dashed", linewidth = 0.8) + scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) + scale_y_continuous(breaks = scales::pretty_breaks(n = 5))

cowplot::plot_grid(fig2a, fig2b, labels = c("A", "B"), label_size = 25, nrow = 1, rel_widths = c(1.5, 1)) %>% ggsave(., filename = file.path(outdir, "genes-fig.pdf"), device = cairo_pdf, units = "in", height = 5, width = 11)

df <- df %>% type.convert(as.is = TRUE)
cor.test(df$NPC, df$EG)
#p-value = 2.308e-7 and cor = 0.995 and 95% CI = 0.97 - 0.99 for Cys
