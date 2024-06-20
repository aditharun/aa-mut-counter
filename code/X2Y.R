library(tidyverse)
library(readxl)
library(Biostrings)

set.seed(123)

outdir <- "../x2y"

if (!dir.exists(outdir)){
	dir.create(outdir)
}

cancercontr <- list.files("../analyses", pattern = "*cancer-contr", full.names = TRUE, recursive = TRUE)

mutrates <- list.files("../analyses", full.names = TRUE, recursive = TRUE, pattern = "*mut-rates.rds")
generates <- list.files("../analyses", full.names = TRUE, recursive = TRUE, pattern = "*gene-rates.rds")

dfmut <- do.call(rbind, lapply(mutrates, function(x) readRDS(x) %>% mutate(aa = gsub(".*p.", "", hugo))))

dfmut <- dfmut %>% mutate(diff = pct_us - pct_tcga) %>% arrange(desc(pct_us))

t25mut <- dfmut %>% arrange(desc(pct_us)) %>% dplyr::slice(1:25)

t25mut <- t25mut %>% mutate(gene_label = str_sub(gene, 1, -2), init_aa = str_sub(gsub("(.*)\\.", "", hugo), 1, 3), number = str_extract_all(gsub("(.*)\\.", "", hugo), "\\d+") %>% unlist(), final_aa = str_sub(hugo, -3, -1))

t25mutlong <- t25mut

aa_converter <- AMINO_ACID_CODE %>% as.data.frame() %>% magrittr::set_colnames("three_letter") %>% mutate(one_letter = rownames(.)) %>% as_tibble()

t25mut <- t25mut %>% left_join(aa_converter, by=c("init_aa" = "three_letter")) %>% left_join(aa_converter, by = c("final_aa"="three_letter")) %>% mutate(label_text = paste0(gene_label, " (", one_letter.x, number, one_letter.y, ")")) %>% select(-c(gene_label, init_aa, number, final_aa, one_letter.x, one_letter.y))

mutbar_second_y_axis_label <- paste0("Estimated Number of New\nCases in the US, 2023")

label_order <- t25mut %>% arrange(desc(pct_us)) %>% pull(label_text)

fig1a <- t25mut %>% mutate(label_text = factor(label_text, levels = label_order)) %>% ggplot(aes(x=pct_us, y=label_text)) + geom_bar(position = position_dodge(width = 5), stat="identity",color="grey50", fill="grey50", alpha = 0.75) + theme_minimal() + theme(panel.grid = element_blank(), panel.border = element_rect(color="black", fill="transparent")) + geom_vline(xintercept = c(1), linetype="dashed", color="grey60") + scale_x_continuous(name = "Estimated Mutation Proportion\nfor U.S. Population (%)", limits=c(0,5), breaks=seq(0,5,0.5), sec.axis = sec_axis(~ . * 0.01 * 1958310, name = mutbar_second_y_axis_label, breaks = seq(0, 90000, 90000/9))) + theme(axis.ticks = element_line(color="black")) + ylab("Top 25 Point Mutations")  + theme(axis.title = element_text(size=18), axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 14, hjust = 1)) 

t1p <- dfmut %>% arrange(desc(pct_us)) %>% mutate(idx = 1:n()) %>% filter(idx < (nrow(dfmut)*0.01))

points_mut <- t1p %>% arrange(desc(pct_us)) %>% mutate(us_rank = 1:n()) %>% filter(abs(diff) > 0.5 | us_rank < 12) %>% mutate(gene_label = str_sub(gene, 1, -2), init_aa = str_sub(gsub("(.*)\\.", "", hugo), 1, 3), number = str_extract_all(gsub("(.*)\\.", "", hugo), "\\d+") %>% unlist(), final_aa = str_sub(hugo, -3, -1)) %>% left_join(aa_converter, by=c("init_aa" = "three_letter")) %>% left_join(aa_converter, by = c("final_aa"="three_letter")) %>% mutate(label_text = paste0(gene_label, " (", one_letter.x, number, one_letter.y, ")")) %>% select(-c(gene_label, init_aa, number, final_aa, one_letter.x, one_letter.y)) %>% mutate(nx = -0.04) %>% mutate(ny = 0.02) %>% mutate(position = "right") %>% mutate(ny = ifelse(us_rank == 10, -0.02, ny)) %>% mutate(position = ifelse(us_rank %in% c(11, 39), "left", position)) %>% mutate(nx = ifelse(us_rank %in% c(11, 39), 0.03, nx))

fig1b <- t1p %>% ggplot(aes(x=pct_tcga, y=pct_us)) + geom_point(alpha=0.8, size=1.75, color="black") + theme_minimal() + geom_abline(color="grey70", size=1, alpha=0.8, linetype="dashed", slope = 1, intercept = 0) + theme(panel.grid = element_blank(), panel.border = element_rect(color="black", fill="transparent")) + scale_x_continuous(limits=c(0,4.5), breaks=seq(0, 4.5, 0.5)) + scale_y_continuous(limits=c(0,4.5), breaks=seq(0, 4.5, 0.5)) + ylab("Estimed Epidemio-Genomic (EG)\nMutation Rate (%)") + xlab("Estimated Naive Pan-Cancer (NPC) Mutation Rate (%)") + theme(axis.title=element_text(size=16), axis.text=element_text(size=14), axis.ticks=element_line(color="black")) + geom_text(data=points_mut, mapping=aes(x=pct_tcga+nx, y=pct_us+ny, label=label_text, hjust = position), lineheight = 0.5)


cowplot::plot_grid(fig1a, fig1b, nrow = 1, ncol = 2, labels = c("A", "B"), label_size = 25) %>% ggsave(plot = ., filename = file.path(outdir, "fig1.pdf"), height = 8, width = 18, device = cairo_pdf, units = "in")

cor.test(t1p$pct_us, t1p$pct_tcga) #Pearson correlation = 0.91 and p-val < 2.2e-16


cowplot::plot_grid(fig1a, fig1b, nrow = 1, labels = c("A", "B"), label_size = 25) %>% ggsave(plot = ., filename = file.path(outdir, "fig1.pdf"), units = "in", device = cairo_pdf, width = 15, height = 8)

#Top 25 mutations with the largest difference between US estimate and TCGA estimate - the data show that we can just reference this in the paper and need not specifically have a figure on this

nuc <- dfmut %>% mutate(gene_label = str_sub(gene, 1, -2), final_aa = str_sub(hugo, -3, -1)) %>% filter(final_aa %in% c("His", "Lys", "Tyr", "Thr", "Ser", "Arg")) 

nuc <- nuc %>% mutate(gene_label = str_sub(gene, 1, -2), init_aa = str_sub(gsub("(.*)\\.", "", hugo), 1, 3), number = (str_extract_all(gsub("(.*)\\.", "", hugo), "\\d+")), final_aa = str_sub(hugo, -3, -1))

nuc$number <- lapply(nuc$number, function(x) x[[1]]) %>% unlist()

nuc <- nuc %>% left_join(aa_converter, by=c("init_aa" = "three_letter")) %>% left_join(aa_converter, by = c("final_aa"="three_letter")) %>% mutate(label_text = paste0(gene_label, " (", one_letter.x, number, one_letter.y, ")")) %>% select(-c(gene_label, init_aa, number, final_aa, one_letter.x, one_letter.y))

nuc_bar <- nuc %>% arrange(desc(pct_us)) %>% dplyr::slice(1:25)

mutbar_second_y_axis_label <- paste0("Estimated Number of New\nCases in the US, 2023")

label_order <- nuc_bar %>% arrange(desc(pct_us)) %>% pull(label_text)

nuc_bar <- nuc_bar %>% mutate(cases = round(pct_us * 0.01 * 1958310))

fig3 <- nuc_bar %>% mutate(label_text = factor(label_text, levels = label_order)) %>% ggplot(aes(x=pct_us, y=label_text)) + geom_bar(position = position_dodge(width = 5), stat="identity",color="grey50", fill="grey50", alpha = 0.75) + theme_minimal() + theme(panel.grid = element_blank(), panel.border = element_rect(color="black", fill="transparent")) + geom_vline(xintercept = c(1), linetype="dashed", color="grey60") + scale_x_continuous(name = "Estimated Mutation Proportion\nfor U.S. Population (%)", limits=c(0,5), breaks=seq(0,5,0.5), sec.axis = sec_axis(~ . * 0.01 * 1958310, name = mutbar_second_y_axis_label, breaks = seq(0, 90000, 90000/9))) + theme(axis.ticks = element_line(color="black")) + ylab("Top 25 Acquired Nucleophilic Residues")  + theme(axis.title = element_text(size=18), axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12, hjust = 1))  + geom_text(aes(x = pct_us + 0.2, y = label_text, label = cases), size = 3.5)

ggsave(plot = fig3, filename = file.path(outdir, "fig3.pdf"), units = "in", device = cairo_pdf, width = 10, height = 8)


#nuc_bar %>% select(label_text, pct_us, pct_tcga) %>% mutate(cases = pct_us * 0.01 * 1958310) %>% arrange(desc(cases)) %>% mutate(pct_us = round(pct_us, 2), pct_tcga = round(pct_tcga, 2), cases = round(cases)) %>% magrittr::set_colnames(c("Mutation", "EG Rate (%)", "NPC Rate (%)", "New cases in 2023")) %>% write_csv(., file = file.path(outdir, "supp-table-2.csv"), col_names = TRUE)

#FIG3: mut-cancer-contribution.rds 
	#can make this a table too of the top 3 cancers for each disease 
	#take top 20 cancers and top 25 mutations and show the distirbution of % mutation found in these cancers

#switch in order so that figure 3 comes after figure 1 and figure 2 comes after figure 3

incidence_df <- "../data/seer-abundances.xlsx" %>% read_excel() %>% filter(`CUSTOM SITES` != "Other Adenocarcinoma") %>% arrange(desc(Incidence)) %>% mutate(c = cumsum(Incidence)) 
incidence_df <- incidence_df %>% filter(Incidence > 1) %>% magrittr::set_colnames(c("rosetta", "name", "incidence"))

parse_cancer_contribution_file <- function(t, t25mutlong, incidence_df){

	filepath <- t
	final_aa_x <- filepath %>% dirname() %>% basename()
	x <- filepath %>% readRDS()

	top3table <- x %>% select(label_text, cancer, frac, gene, informed_rate)  %>% mutate(indiv_imp = paste0(cancer, " (", round(frac*100, 2), "%)")) %>% select(-c(cancer, frac)) %>% group_by(informed_rate, label_text) %>% mutate(row_id = row_number()) %>% pivot_wider(names_from = row_id, values_from = indiv_imp, names_prefix = "V") %>% mutate(cancers = paste0(V1, ", ", V2, ", ", V3)) %>% select(-starts_with("V")) %>% arrange(desc(informed_rate)) %>% relocate(label_text) %>% ungroup() %>% select(-gene) %>% mutate(final_aa = final_aa_x)

	top3table <- top3table %>% mutate(label_text = label_text %>% str_sub(., 1, -2) %>% paste0(., " ", final_aa, ")")) %>% select(-final_aa)

	#%>% magrittr::set_colnames(c("Point Mutation", "Mutation Proportion (%)", "The most common cancers to harbor a mutation in this gene (ROSETTA classification, % of all instances of this mutation)"))

	t25subset <- t25mutlong %>% filter(final_aa == final_aa_x)

	if (nrow(t25subset) > 0){

		t25subset <- t25subset %>% mutate(matchstring = paste0(gene_label, " (", init_aa, " ", number, ")"))

		y <- x %>% select(label_text, rosetta, frac, gene, informed_rate) %>% filter(rosetta %in% incidence_df$rosetta) %>% select(label_text, rosetta, frac) %>% pivot_wider(values_from = frac, names_from = rosetta)

		y <- y %>% right_join(t25subset %>% select(hugo, matchstring), by=c("label_text" = "matchstring")) %>% mutate_if(is.numeric, ~round(.*100, 2)) 

		return(list(t25_df = y, top3 = top3table))

	}

	return(list(top3 = top3table))

}

cancercontr_df <- lapply(cancercontr, function(p) parse_cancer_contribution_file(p, t25mutlong, incidence_df))

t25cancercontr_df <- lapply(cancercontr_df, function(p) p$t25_df) %>% do.call(rbind, .) 

df_25_cc <- t25cancercontr_df %>% select(-label_text) %>% relocate(hugo) %>% left_join(.,t25mut %>% select(hugo, label_text, pct_us)) %>% relocate(label_text, pct_us) %>% select(-hugo) 

cc_factor_levels <- df_25_cc %>% select(label_text, pct_us) %>% arrange(desc(pct_us)) %>% pull(label_text)

df_25_cc <- df_25_cc %>% mutate(label_text = factor(label_text, levels = cc_factor_levels))

df_25_cc <- df_25_cc  %>% select(-pct_us) %>% pivot_longer(-label_text) %>% left_join(., incidence_df %>% select(1,2), by=c("name" = "rosetta"))

cancer_types <- c(
  "Breast Carcinoma", 
  "Prostate Cancer",
  "Colorectal Adenocarcinoma",
  "Lung Adenocarcinoma",
  "Malignant Melanoma",
  "Urothelial Cancer",
  "Lung Squamous Cell Carcinoma",
  "Head Neck Squamous Cell Carcinoma",
  "Renal Cell Carcinoma",
  "Endometrial Cancer",
  "Pancreatic Adenocarcinoma",
  "Thyroid Carcinoma",
  "Lung Small Cell Carcionoma",
  "Hepatocellular Carcinoma",
  "Gastric Adenocarcinoma",
  "Diffuse_large_B_cell_lymphoma",
  "Ovarian Cancer",
  "Plasma_cell_myeloma",
  "Chronic_lymphocytic_leukemiasmall_lymphocytic_lymphoma"
)

cancer_types_short <- c(
  "Breast\nCarcinoma", 
  "Prostate\nCancer",
  "Colorectal\nAdenocarcinoma",
  "Lung\nAdenocarcinoma",
  "Malignant\nMelanoma",
  "Urothelial\nCancer",
  "Lung Squamous\nCell Carcinoma",
  "Head Neck Squamous\nCell Carcinoma",
  "Renal Cell\nCarcinoma",
  "Endometrial\nCancer",
  "Pancreatic\nAdenocarcinoma",
  "Thyroid\nCarcinoma",
  "Lung Small\nCell Carcionoma",
  "Hepatocellular\nCarcinoma",
  "Gastric\nAdenocarcinoma",
  "DLBCL",
  "Ovarian\nCancer",
  "Plasma cell\nmyeloma",
  "CLL/SLL"
)



short_levels <- data.frame(cancer_types = cancer_types, cancer_types_short = cancer_types_short) %>% mutate(cancer_types = factor(cancer_types, levels = (incidence_df %>% select(2,3) %>% arrange(desc(incidence)) %>% pull(name)) ) ) %>% as_tibble() %>% pull(cancer_types_short) %>% rev()

df_25_cc <- df_25_cc %>% left_join(., data.frame(cancer_types = cancer_types, cancer_types_short = cancer_types_short) %>% mutate(cancer_types_short  = factor(cancer_types_short, levels = short_levels)), by=c("name.y" = "cancer_types")) 

#11 hieght by 8 width
fig2 <- ggplot(df_25_cc, aes(y= cancer_types_short, x= label_text, fill = value)) + geom_tile() +  scale_fill_gradient(low = "lightgrey", high = "red", na.value = "lightgrey", breaks = scales::pretty_breaks(n = 6))  + labs(x = "", y="", fill = "% of samples") + theme_minimal() + theme(plot.title = element_text(size = 9, hjust = 0.5), axis.text.y = element_text(size = 9), axis.text.x = element_text(size = 10, angle = 45, hjust = 1),legend.title = element_text(size = 10), legend.text = element_text(size = 8), legend.key.size = unit(0.25, "cm")) + theme(panel.background = element_blank(), panel.grid = element_blank()) + theme(axis.ticks = element_line(color = "grey70"))  + xlab("Top 25 Point Mutations") + ylab("Cancers with >1% Incidence") + theme(axis.title = element_text(size = 14))


ggsave(plot = fig2, filename = file.path(outdir, "fig2.pdf"), units = "in", width = 11, height = 8, device = cairo_pdf)


#df_25_cc %>% group_by(label_text) %>% summarize(v = sum(value)) %>% pull(v)
#percent of mutations covered for each of the top25 by the cancers w/ incidence > 1 (n = 19)


df_25_cc_table <- df_25_cc %>% select(label_text, name.y, value) %>% pivot_wider(names_from = name.y, values_from = value) 
colnames(df_25_cc_table)[which(colnames(df_25_cc_table) == "label_text")] <- "Mutation"

df_25_cc_table %>% write_csv(., file = file.path(outdir, "supp-table-2.csv"), col_name = TRUE)

top3 <- lapply(cancercontr_df, function(p) p$top3) %>% do.call(rbind, .) 
dfmut2 <- dfmut %>% mutate(label_text = paste0(str_sub(gene, 1, -2), " (", str_sub(gsub("(.*)\\.", "", hugo), 1, 3), " ", str_extract_all(gsub("(.*)\\.", "", aa), "\\d+"), " ", str_sub(hugo, -3, -1), ")")) %>% mutate(acquired_nucleophile = ifelse(str_sub(hugo, -3, -1) %in% c("His", "Lys", "Tyr", "Thr", "Ser", "Arg"), TRUE, FALSE))

dfmut2 <- dfmut2 %>% left_join(top3, by=c("label_text"="label_text")) %>% mutate(gene = str_sub(gene, 1, -2)) %>% select(label_text, gene, pct_us, pct_lb, pct_ub, pct_tcga, acquired_nucleophile, cancers) %>% magrittr::set_colnames(c("Mutation", "Gene", "EG Rate", "EG lower bound (2.5%) rate", "EG upper bound (97.5%) rate", "NPC Rate", "Acquired Nucleophile", "Top 3 Cancers"))

dfmut2 %>% write_csv(., file = file.path(outdir, "supp-table-1.csv"), col_names = TRUE)


#Correlation between tcga and us rate for top x% of US cases
df <- dfmut %>% arrange(desc(pct_us)) %>% mutate(idx = 1:n())


#### supp-table-3.xlsx

pct_keep <- 10
df2 <- df %>% filter(idx < (nrow(df) * pct_keep * (1/100)) )
cor.test(df2$pct_us, df2$pct_tcga)

###




#FIG5: most common individual mutations in the most common/important cancer driver genes (e.g., tp53, kras, pik3ca) 
#Can get Ed's thoughts on this before building it
	#panel A: fig2 cysteine
	#panel B: comparison of us and tcga rates in stacked bar plot
	#panel C/table: which cancers do these mutations most commonly occur in?

	#most important cancer driver genes: KRAS, TP53, LRP1B, PTEN, PIK3CA, KMT2D, NRAS, ARID1A (https://www.nature.com/articles/s41568-020-0290-x)
	#we actually end up covering most of this stuff, so I'm not sure how much utility there is here

#FIG6: most common individual changes that are not in a gene known to be a cancer driver gene - we can select these based on US and TCGA percentages 
#Can get Ed's thoughts on this before building it
	#look at mut-cancer-contribution for further cancer-specific gene specificity
	#panel A: fig 2 cysteine 
	#panel B: 
	#but can we just comment on it in the results? Is there are reason to show this?



#Supplementary table: table 1 - all data + top3 cancer contributions; table 2 - nucleophilic data (?)
#Supplementary figure / table: cancer contribution data for top 25
#supplementary table: correlation between tcga and us in values and rank for top 1%, 10%, 50%, 100% of data, correlation and p-value


#Supplementary Figure (optional at this point): Top 100 mutations by Top 20 cancers, cancer contribution heatmap
#Shiny app for the mutations!! Like the Allen shiny atlas I made for Kathy. 
















#