library(tidyverse)
library(readxl)

set.seed(123)

outdir <- "../analyses/Cys"
if (!dir.exists(outdir)){
	dir.create(outdir, recursive = TRUE)
}

mut <- file.path(outdir, "mut-rates.rds") %>% readRDS()

mut.bar <- mut %>% arrange(desc(pct_tcga)) %>% slice(1:50) 
mut.bar <- mut.bar %>% mutate(label_text = paste0(str_sub(gene, 1, -2), " (", gsub("(.*)\\.", "", hugo) %>% str_sub(., 1, -4) %>% str_replace(., "([[:alpha:]])(\\d)", "\\1 \\2"), ")"))

mut.bar$label_text <- factor(mut.bar$label_text, levels = mut.bar$label_text[order(mut.bar$pct_us, decreasing = TRUE)])

x2c_mut_bar <- mut.bar %>% ggplot(aes(x=label_text, y=pct_us)) + geom_bar(position = position_dodge(width = 5), stat="identity",color="grey50", fill="grey50", alpha = 0.75) + theme_minimal() + theme(panel.grid = element_blank(), panel.border = element_rect(color="black", fill="transparent")) + geom_hline(yintercept = c(0.5, 1, 1.5), linetype="dashed", color="grey60") + scale_y_continuous(name = "Estimated Mutation Proportion\nfor U.S. Population (%)", limits=c(0,1.8), breaks=seq(0,2,0.5), sec.axis = sec_axis(~ . * 0.01 * 1958310, name = "Estimated Number of New Cases In 2023\nWith Acquired Cysteine Mutations In U.S.", breaks = seq(0, 30000, 30000/6))) + theme(axis.ticks = element_line(color="black")) + xlab("Top 50 Point Mutations")  + theme(axis.title = element_text(size=18), axis.text.y = element_text(size = 16), axis.text.x = element_text(size = 12, angle = 70, hjust = 1)) 


ggsave(plot = x2c_mut_bar, filename = file.path(outdir, "fig3.pdf"), device = cairo_pdf, units = "in", height = 8, width = 15)
