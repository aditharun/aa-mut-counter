library(tidyverse)
library(readxl)

set.seed(123)

args = commandArgs(trailingOnly=TRUE)
amino_acid <- args[1]
amino_acid_long <- args[2]

outdir <- paste0("../analyses/", amino_acid)

if (!dir.exists(outdir)){
	dir.create(outdir, recursive = TRUE)
}

mut <- file.path(outdir, "mut-rates.rds") %>% readRDS()

mut.bar <- mut %>% arrange(desc(pct_tcga)) %>% slice(1:50) 
mut.bar <- mut.bar %>% mutate(label_text = paste0(str_sub(gene, 1, -2), " (", gsub("(.*)\\.", "", hugo) %>% str_sub(., 1, -4) %>% str_replace(., "([[:alpha:]])(\\d)", "\\1 \\2"), ")"))

mut.bar$label_text <- factor(mut.bar$label_text, levels = mut.bar$label_text[order(mut.bar$pct_us, decreasing = TRUE)])

mutbar_second_y_axis_label <- paste0("Estimated Number of New Cases In 2023\nWith Acquired ", amino_acid_long ," Mutations In U.S.")

#mut.bar <- mut.bar %>% arrange(desc(pct_us)) %>% mutate(idx = 1:n())


x2c_mut_bar <- mut.bar %>% ggplot(aes(x=label_text, y=pct_us)) + geom_bar(position = position_dodge(width = 5), stat="identity",color="grey50", fill="grey50", alpha = 0.75) + geom_errorbar(aes(x= label_text, ymin = pct_lb, ymax = pct_ub), color = "grey55", size = 0.2, width = 0.4) + theme_minimal() + theme(panel.grid = element_blank(), panel.border = element_rect(color="black", fill="transparent")) + geom_hline(yintercept = c(0.5, 1, 1.5), linetype="dashed", color="grey60") + scale_y_continuous(name = "Estimated Mutation Proportion\nfor U.S. Population (%)", limits=c(0,1.8), breaks=seq(0,2,0.5), sec.axis = sec_axis(~ . * 0.01 * 1958310, name = mutbar_second_y_axis_label, breaks = seq(0, 30000, 30000/6))) + theme(axis.ticks = element_line(color="black")) + xlab("Top 50 Point Mutations")  + theme(axis.title = element_text(size=18), axis.text.y = element_text(size = 16), axis.text.x = element_text(size = 12, angle = 70, hjust = 1)) 


non_ci_x2c_mut_bar <- mut.bar %>% ggplot(aes(x=label_text, y=pct_us)) + geom_bar(position = position_dodge(width = 5), stat="identity",color="grey50", fill="grey50", alpha = 0.75) + theme_minimal() + theme(panel.grid = element_blank(), panel.border = element_rect(color="black", fill="transparent")) + geom_hline(yintercept = c(0.5, 1, 1.5), linetype="dashed", color="grey60") + scale_y_continuous(name = "Estimated Mutation Proportion\nfor U.S. Population (%)", limits=c(0,1.8), breaks=seq(0,2,0.5), sec.axis = sec_axis(~ . * 0.01 * 1958310, name = mutbar_second_y_axis_label, breaks = seq(0, 30000, 30000/6))) + theme(axis.ticks = element_line(color="black")) + xlab("Top 50 Point Mutations")  + theme(axis.title = element_text(size=18), axis.text.y = element_text(size = 16), axis.text.x = element_text(size = 12, angle = 70, hjust = 1)) 


ggsave(plot = x2c_mut_bar, filename = file.path(outdir, "fig2.pdf"), device = cairo_pdf, units = "in", height = 8, width = 15)

ggsave(plot = non_ci_x2c_mut_bar, filename = file.path(outdir, "fig2_no_ci.pdf"), device = cairo_pdf, units = "in", height = 8, width = 15)
