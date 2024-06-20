#setwd("~/Documents/GitHub/aa-mut-counter/code/")
library(tidyverse)
library(readxl)
library(Biostrings)

set.seed(123)

standard_aa <- AMINO_ACID_CODE[1:20] %>% unname()

conf_int <- function(lambda, n=2000){

  values <- rpois(n=n, lambda=lambda)

  bounds <- quantile(values, c(0.025, 0.975)) %>% unname()

  return(list(lb = bounds[1], ub = bounds[2]))

}

align_columns_generic <- function(all_cases_cols, df){

    master_cols <- all_cases_cols

    s_cols <- df

    cols_to_add <- setdiff(master_cols, colnames(s_cols))

    if (length(cols_to_add) > 0){

      for (x in cols_to_add){
        s_cols[[x]] <- 0
      }
    }

    s_cols %>% relocate(master_cols)

  }

outdir <- "../x2nuc"

if (!dir.exists(outdir)){
	dir.create(outdir)
}

aatheme <- theme_bw() + theme(panel.grid.minor=element_blank(), panel.grid.major = element_blank()) + theme(axis.text = element_text(size=12), axis.title=element_text(size=16), legend.text=element_text(size=10), legend.title=element_text(size=14), plot.title=element_text(size=18, hjust=0.5))

aatheme_small <- theme_bw() + theme(panel.grid.minor=element_blank(), panel.grid.major = element_blank()) + theme(axis.text = element_text(size=9), axis.title=element_text(size=12), legend.text=element_text(size=9), legend.title=element_text(size=10), plot.title=element_text(size=14, hjust=0.5))

cpal  <- c("#374E5599", "#Df8F4499", "#00A1D5FF", "#B24745FF", "#79AF97FF", "#6A6599F")

source("generic_reweight.R")

nucleophile_aa <- c("His", "Lys", "Tyr", "Thr", "Ser", "Arg", "Cys")

aa_vec <- (AMINO_ACID_CODE %>% unlist() %>% unname())[1:20]
#unweighted: % of all nucleophilic mutations that are Cys mutations?
#            % of all mutations that are nucleopphilic mutations?

#the sum I gather here is the sum() of all mutations for a given amino acid
total_muts <- data.frame(aa = aa_vec, n.muts = paste0("../results/counts-", aa_vec, ".txt") %>% lapply(., function(x) x %>% read_table() %>% filter(Hugo_Symbol != "Total") %>% pull(All) %>% sum()) %>% unlist())

total_muts <- total_muts %>% mutate(nucleophile = ifelse(aa %in% nucleophile_aa, TRUE, FALSE))

total_muts <- total_muts %>% mutate(n.nuc = sum(n.muts[nucleophile])) %>% mutate(n.total = sum(n.muts))


fig1_a <- total_muts %>% mutate(frac = round( (n.muts / n.total) * 100, 2)) %>% ggplot(aes(x=reorder(aa, -frac), y=frac, color = nucleophile, fill = nucleophile)) + geom_col(alpha = 0.6) + aatheme + ylab("Fraction of Total\nSequenced Mutations (%)") + scale_fill_manual(name = "Acquired\nNucleophile", values = c("TRUE" = cpal[1], "FALSE" = cpal[2])) + xlab("Acquired Amino Acid") + scale_color_manual(name = "Acquired\nNucleophile", values = c("TRUE" = cpal[1], "FALSE" = cpal[2])) 

fig1_b <- total_muts %>% filter(nucleophile) %>% mutate(frac = round((n.muts / n.nuc) * 100, 2)) %>% ggplot(aes(x=reorder(aa, -frac), y = frac), alpha = 0.75, color = "grey70") + geom_col() + aatheme + ylab("Fraction of Acquired\nNucleophile Mutations (%)") + xlab("Acquired Nucleophile")


#Weighted version is done for all cancer
#Of all mutations, what % is cysteine?

cancercontr <- list.files("../analyses", pattern = "*cancer-contr", full.names = TRUE, recursive = TRUE)

df.list <- lapply(cancercontr, function(x) x %>% readRDS() %>% group_by(rosetta) %>% summarize(count = sum(count_mutation_for_given_rosetta)) %>% mutate(aa = dirname(x) %>% basename()) %>% pivot_wider(names_from = rosetta, values_from = count) )

df.list <- df.list %>% lapply(., function(x) x %>% mutate(All = 0) %>% relocate(aa, All) %>% dplyr::rename(Hugo_Symbol = aa))

all_cases_cols <- "../results/all-muts/all-counts-matrix.txt" %>% read_table() %>% colnames()

total <- lapply(df.list, function(t) align_columns_generic( all_cases_cols, t)) %>% do.call(rbind, .) %>% as_tibble() %>% summarize_if(is.numeric, sum) %>% mutate(Hugo_Symbol = "Total") %>% relocate(Hugo_Symbol)

reweighted_all_muts <- generic_reweight(df.list, total)


allmuts.l <- reweighted_all_muts %>% arrange(desc(pct_us)) %>% mutate(idx = 1:n()) %>% pivot_longer(-c(hugo, pct_lb, pct_ub, idx)) %>% mutate(pct_lb = ifelse(name == "pct_us", pct_lb, NA), pct_ub = ifelse(name == "pct_us", pct_ub, NA)) 

fig1_a <- allmuts.l %>% ggplot(., aes(fill=name, y=value, color = name, x=idx, group=name)) + geom_bar(stat="identity",position=position_dodge(width = 0.8), width = 0.8, alpha = 0.4) + geom_errorbar(aes( x = (idx + 0.2), ymax = pct_ub, ymin = pct_lb, group = name),width=0.15) + aatheme + xlab("Acquired Amino Acid") + ylab("Mutation Rate (%)") + scale_color_manual(name = "", values = c("pct_tcga" = cpal[3], "pct_us" = cpal[4]), labels = c("pct_tcga" = "NPC rate", "pct_us" = "EG rate")) +  scale_fill_manual(name = "", values = c("pct_tcga" = cpal[3], "pct_us" = cpal[4]), labels = c("pct_tcga" = "NPC rate", "pct_us" = "EG rate")) + scale_x_continuous(breaks = allmuts.l$idx, labels = allmuts.l$hugo) + scale_y_continuous(breaks = scales::pretty_breaks(n = 8))


allmuts.l <- allmuts.l %>% filter(hugo %in% nucleophile_aa) %>% arrange(hugo) %>% mutate(idx = ((1:n())/2) %>% ceiling() )
fig1_anuc <- allmuts.l %>% ggplot(., aes(fill=name, y=value, color = name, x=idx, group=name)) + geom_bar(stat="identity",position=position_dodge(width = 0.8), width = 0.8, alpha = 0.4) + geom_errorbar(aes( x = (idx + 0.2), ymax = pct_ub, ymin = pct_lb, group = name),width=0.15) + aatheme_small + xlab("Acquired Nucleophile") + ylab("Mutation Rate (%)") + scale_color_manual(name = "", values = c("pct_tcga" = cpal[3], "pct_us" = cpal[4]), labels = c("pct_tcga" = "NPC rate", "pct_us" = "EG rate")) +  scale_fill_manual(name = "", values = c("pct_tcga" = cpal[3], "pct_us" = cpal[4]), labels = c("pct_tcga" = "NPC rate", "pct_us" = "EG rate")) + scale_x_continuous(breaks = allmuts.l$idx, labels = allmuts.l$hugo) + scale_y_continuous(breaks = scales::pretty_breaks(n = 8), limits = c(0, 10)) + theme(legend.position = "none")



#Weighted version for nucleophile
#Of all nucleophile mutations, what % is cysteine?

cancercontr <- list.files("../analyses", pattern = "*cancer-contr", full.names = TRUE, recursive = TRUE)

#cancercontr <- cancercontr[which(dirname(cancercontr) %>% basename() %in% nucleophile_aa)]


#same as the df.list loaded above
#df.list <- lapply(cancercontr, function(x) x %>% readRDS() %>% group_by(rosetta) %>% summarize(count = sum(count_mutation_for_given_rosetta)) %>% mutate(aa = dirname(x) %>% basename()) %>% pivot_wider(names_from = rosetta, values_from = count) )

#df.list <- df.list %>% lapply(., function(x) x %>% mutate(All = 0) %>% relocate(aa, All) %>% dplyr::rename(Hugo_Symbol = aa))

all_cases_cols <- "../results/all-muts/all-counts-matrix.txt" %>% read_table() %>% colnames()

df.list2 <- df.list[which(dirname(cancercontr) %>% basename() %in% nucleophile_aa)]

total <- lapply(df.list2, function(t) align_columns_generic( all_cases_cols, t)) %>% do.call(rbind, .) %>% as_tibble() %>% summarize_if(is.numeric, sum) %>% mutate(Hugo_Symbol = "Total") %>% relocate(Hugo_Symbol)

reweighted_nuc_muts <- generic_reweight(df.list2, total)

nucmuts.l <- reweighted_nuc_muts %>% arrange(desc(pct_us)) %>% mutate(idx = 1:n()) %>% pivot_longer(-c(hugo, pct_lb, pct_ub, idx)) %>% mutate(pct_lb = ifelse(name == "pct_us", pct_lb, NA), pct_ub = ifelse(name == "pct_us", pct_ub, NA)) 

nucmuts.l <- nucmuts.l %>% arrange(hugo) %>% mutate(idx = ((1:n())/2) %>% ceiling() )

fig1_b <- nucmuts.l %>% ggplot(., aes(fill=name, y=value, color = name, x=idx, group=name)) + geom_bar(stat="identity",position=position_dodge(width = 0.8), width = 0.8, alpha = 0.4) + geom_errorbar(aes( x = (idx + 0.2), ymax = pct_ub, ymin = pct_lb, group = name),width=0.15) + aatheme_small + xlab("Acquired Nucleophile") + ylab("Nucleophile Mutation Rate (%)") + scale_color_manual(name = "", values = c("pct_tcga" = cpal[3], "pct_us" = cpal[4]), labels = c("pct_tcga" = "NPC rate", "pct_us" = "EG rate")) +  scale_fill_manual(name = "", values = c("pct_tcga" = cpal[3], "pct_us" = cpal[4]), labels = c("pct_tcga" = "NPC rate", "pct_us" = "EG rate")) + scale_x_continuous(breaks = nucmuts.l$idx, labels = nucmuts.l$hugo) + scale_y_continuous(breaks = scales::pretty_breaks(n = 8), limits = c(0, 20))



#916,075 nucleophile mutations and 2,226,189 total mutations, makes up 41.1% of mutations. It is 7 aa and 7/20 is 35%, so a mild positive bias relative to uniform distribution

#If I am trying to answer, "% of all patients w/ at least 1 nucleophilic mutation that have at least 1 Cys mutation?" then I should use the total row


#unweighted/weighted: % of patients with a nucleophile mutation that have a cys mutation?
#                     % of patients that have a cys mutation
#                     I need to use the total vector for each amino acid and look at it reweighted and unreweighted


total_df <- lapply(seq_along(aa_vec), function(x) paste0("../results/counts-", aa_vec[x], ".txt") %>% read_table() %>% filter(Hugo_Symbol == "Total") %>% mutate(Hugo_Symbol = ifelse(Hugo_Symbol == "Total", aa_vec[x], NA))) 

seer.df <- "../data/seer-abundances.xlsx" %>% read_excel() %>% magrittr::set_colnames(c("rosetta", "cancer", "incidence")) %>% mutate(incidence_frac = incidence / 100)
all_cases <- "../results/all-muts/all-counts-matrix.txt" %>% read_table() %>% filter(Hugo_Symbol == "Total")

#need to align column names to the all_cases 
#get unweighted values 
#then reweight according to reweighting.R and get reweighted values

align_columns <- function(all_cases, df){

  master_cols <- all_cases %>% colnames()

  s_cols <- df

  cols_to_add <- setdiff(master_cols, colnames(s_cols))

  if (length(cols_to_add) > 0){

    for (x in cols_to_add){
      s_cols[[x]] <- 0
    }
  }

  s_cols %>% relocate(master_cols)

}



total_df2 <- lapply(total_df, function(x) align_columns(all_cases, x))

total_df <- total_df2 %>% do.call(rbind, .) %>% as_tibble()

#weighting
total <- all_cases %>% pivot_longer(-c(Hugo_Symbol, All)) %>% left_join(seer.df, by=c("name"="rosetta")) 

total <- total %>% mutate(incidence = incidence / sum(incidence), incidence_frac = incidence_frac / sum(incidence_frac))

total <- total %>% magrittr::set_colnames(c("hugo", "all", "rosetta", "count", "cancer", "incidence", "incidence_frac"))

tcga.aa <- total_df %>% select(-All) %>% pivot_longer(-Hugo_Symbol) %>% left_join(total %>% select(-c(count, hugo, all)), by=c("name"="rosetta"))

tcga.aa <- tcga.aa %>% magrittr::set_colnames(c("hugo", "rosetta", "count", "cancer", "incidence", "incidence_frac"))

total <- total %>% mutate(weight.tot = incidence_frac * count)
tcga.aa <- tcga.aa %>% mutate(weight.aa=incidence_frac * count) 
tcga.aa <- tcga.aa %>% mutate(idx = 1:n())
tcga.aa.confint <- tcga.aa %>% filter(count > 0)  %>% mutate(result = purrr::map(count, conf_int), lb = purrr::map_dbl(result, "lb"), ub = purrr::map_dbl(result, "ub")) %>% select(-result)
tcga.aa <- tcga.aa %>% left_join(tcga.aa.confint %>% select(idx, lb, ub), by=c("idx"="idx"))
tcga.aa <- tcga.aa %>% mutate(weight.lb = lb * incidence_frac, weight.ub = ub * incidence_frac)
tcga.aa <- tcga.aa %>% mutate(weight.lb = ifelse(is.na(weight.lb), 0, weight.lb), weight.ub = ifelse(is.na(weight.ub), 0, weight.ub))
aa <- tcga.aa %>% left_join(total %>% select(rosetta, weight.tot, all), by=c("rosetta"="rosetta")) %>% group_by(hugo) %>% summarize(pct_us = sum(weight.aa) / sum(weight.tot), pct_tcga = sum(count) / unique(all) , pct_lb = sum(weight.lb) / sum(weight.tot), pct_ub = sum(weight.ub) / sum(weight.tot)) %>% mutate(across(where(is.double), ~ . * 100)) 
aa <- aa %>% mutate(nucleophile = ifelse(hugo %in% nucleophile_aa, TRUE, FALSE)) 

#what % of patients newly diagnosed with cancer each year have at least 1 acquired XYZ amino acid?

 plt.df <- aa %>% arrange(desc(pct_us)) %>% mutate(idx = 1:n()) %>% pivot_longer(-c(idx, hugo, nucleophile, pct_ub, pct_lb)) %>% mutate(pct_ub = ifelse(name == "pct_tcga", NA, pct_ub), pct_lb = ifelse(name == "pct_tcga", NA, pct_lb)) 

 fig1c <- plt.df %>% ggplot(., aes(fill=name, y=value, color = name, x=idx, group=name)) + geom_bar(stat="identity",position=position_dodge(width = 0.8), width = 0.8, alpha = 0.4) + geom_errorbar(aes( x = (idx + 0.2), ymax = pct_ub, ymin = pct_lb, group = name),width=0.15) + aatheme + xlab("Acquired Amino Acid") + ylab("Patients with at least 1 mutation (%)") + scale_color_manual(name = "", values = c("pct_tcga" = cpal[3], "pct_us" = cpal[4]), labels = c("pct_tcga" = "NPC rate", "pct_us" = "EG rate")) +  scale_fill_manual(name = "", values = c("pct_tcga" = cpal[3], "pct_us" = cpal[4]), labels = c("pct_tcga" = "NPC rate", "pct_us" = "EG rate")) + scale_x_continuous(breaks = plt.df$idx, labels = plt.df$hugo) + scale_y_continuous(breaks = scales::pretty_breaks(n = 8))


 plt.df <- aa %>% filter(nucleophile) %>% arrange(desc(pct_us)) %>% mutate(idx = 1:n()) %>% pivot_longer(-c(idx, hugo, nucleophile, pct_ub, pct_lb)) %>% mutate(pct_ub = ifelse(name == "pct_tcga", NA, pct_ub), pct_lb = ifelse(name == "pct_tcga", NA, pct_lb)) 


plt.df <- plt.df %>% arrange(hugo) %>% mutate(idx = ((1:n())/2) %>% ceiling() )


 fig1d <- plt.df %>% ggplot(., aes(fill=name, y=value, color = name, x=idx, group=name)) + geom_bar(stat="identity",position=position_dodge(width = 0.8), width = 0.8, alpha = 0.4) + geom_errorbar(aes( x = (idx + 0.2), ymax = pct_ub, ymin = pct_lb, group = name),width=0.15) + aatheme_small + xlab("Acquired Nucleophile") + ylab("Patients with at least 1 mutation (%)") + scale_color_manual(name = "", values = c("pct_tcga" = cpal[3], "pct_us" = cpal[4]), labels = c("pct_tcga" = "NPC rate", "pct_us" = "EG rate")) +  scale_fill_manual(name = "", values = c("pct_tcga" = cpal[3], "pct_us" = cpal[4]), labels = c("pct_tcga" = "NPC rate", "pct_us" = "EG rate")) + scale_x_continuous(breaks = plt.df$idx, labels = plt.df$hugo) + scale_y_continuous(breaks = scales::pretty_breaks(n = 8), limits = c(0, 100))


#Figure 1 for X2Nuc paper
cowplot::plot_grid(fig1_anuc, fig1d, nrow = 1, labels = c("A", "B"), label_size = 18, rel_widths = c(1, 1.1)) %>% ggsave(filename = "../x2nuc/fig1.pdf", ., units = "in", device = cairo_pdf, height = 5, width = 11)


#Supplemental Figure 1 for X2Nuc paper
ggsave(filename = "../x2nuc/sfig2.pdf", fig1_b, units = "in", device = cairo_pdf, height = 5, width = 5.5)




#I am counting the number of mutations for the given rosetta code, and I am counting the total number of mutations across all rosetta codes for that mutation

cancercontr <- list.files("../analyses", pattern = "*cancer-contr", full.names = TRUE, recursive = TRUE)

nuc_cancercontr <- cancercontr[which(dirname(cancercontr) %>% basename() %in% nucleophile_aa)]

incidence_df <- "../data/seer-abundances.xlsx" %>% read_excel() %>% filter(`CUSTOM SITES` != "Other Adenocarcinoma") %>% arrange(desc(Incidence)) %>% mutate(c = cumsum(Incidence)) 
incidence_df <- incidence_df %>% filter(Incidence > 1) %>% magrittr::set_colnames(c("rosetta", "name", "incidence"))

parse_cancer_contribution_file <- function(t, incidence_df){

  x <- t %>% readRDS()
  acquired_aa <- dirname(t) %>% basename()

  total <- x$count_mutation_for_given_rosetta %>% sum()

  y <- x %>% filter(rosetta %in% incidence_df$rosetta) %>% group_by(rosetta) %>% summarize(n.mut = sum(count_mutation_for_given_rosetta), rosetta = unique(rosetta)) %>% left_join(., incidence_df %>% magrittr::set_colnames(c("rosetta", "name", "incidence", "cum_incidence"))) 

  y %>% mutate(pct = round((n.mut / total)*100, 2)) %>% mutate(acquired_aa = acquired_aa)

}

nuc_cancercontr_df <- lapply(nuc_cancercontr, function(p) parse_cancer_contribution_file(p, incidence_df))

df_cc <- do.call(rbind, nuc_cancercontr_df) 

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

df_cc2 <- df_cc %>% left_join(., data.frame(cancer_types = cancer_types, cancer_types_short = cancer_types_short) %>% mutate(cancer_types_short  = factor(cancer_types_short, levels = short_levels)), by=c("name" = "cancer_types")) 

#incidence bar plot for >1% cancers

#11 hieght by 8 width
fig1e <- ggplot(df_cc2, aes(x=cancer_types_short, y= acquired_aa, fill = pct)) + geom_tile() +  scale_fill_gradient(low = "lightgrey", high = "red", na.value = "lightgrey", breaks = scales::pretty_breaks(n = 6), limits = c(0, 35))  + labs(x = "", y="", fill = "% of\nsamples") + theme_minimal() + theme(plot.title = element_text(size = 9, hjust = 0.5), axis.text.y = element_text(size = 9), axis.text.x = element_text(size = 10)) + theme(panel.background = element_blank(), panel.grid = element_blank())   + ylab("Acquired Nucleophile") + xlab("Cancers with >1% Incidence") + theme(axis.title = element_text(size = 14)) + coord_flip() + theme(legend.position = "left", legend.text = element_text(size = 8, vjust = 0.5), legend.key.size = unit(0.8, "cm"), legend.title = element_text(size = 10, hjust = 0.5)) 

fig1f <- incidence_df %>% left_join(., data.frame(cancer_types = cancer_types, cancer_types_short = cancer_types_short) %>% mutate(cancer_types_short  = factor(cancer_types_short, levels = short_levels)), by=c("name" = "cancer_types")) %>% select(3,5) %>% ggplot(aes(y=cancer_types_short, x = incidence)) + geom_col(alpha = 0.6, color = "grey40") + aatheme + xlab("Incidence (%)") + ylab("") + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) + theme(panel.border = element_blank()) + theme(axis.line.x = element_line(color = "grey30")) + scale_x_continuous(limits = c(0, 17), expand = c(0,0), breaks = c(0, 5, 10, 15))

fig1heatmap <- cowplot::plot_grid(fig1e, NULL, fig1f, align = "h", axis = "l", nrow = 1, rel_widths = c(2.9, -0.05, 0.8)) 



#Supplemental Table 1 for X2Nuc paper
df_cc2 %>% select(name, acquired_aa, pct) %>% pivot_wider(values_from = pct, names_from = acquired_aa) %>% write_csv(x = ., file = "../x2nuc/stable2.csv")



#I am not making statments about % of patients with at least 1 such transistion. I could, I would look back at the dfmutfiles. But there seems to be limited utility in this. Squeeze not worth the juice. 

# If this is a X2Nuc paper then I ought to be showing primarily for nucleophiles. 

# In XYZ aa mutation, what is the most common aa it switches to (different total row for each aa) DEFINETLY DO THIS and show as figure 1) for nucleophile, 2) for all aa

#can make {heatmap} or {bar chart of the top 10} for the two below. Bar chart likely more useful. 

# In all aa mutation, what is the most common aa it switch to (same total row) 

# In all nuc mutation, what is the most common aa it switches to (same total row, configured once for nuc)

x <-  lapply(seq_along(aa_vec), function(t) paste0("../results/counts-", aa_vec[t], ".txt") %>% read_table() %>% filter(Hugo_Symbol != "Total") %>% mutate(final_aa = str_sub(Hugo_Symbol, -3, -1), init_aa = str_sub(gsub("(.*)\\.", "", Hugo_Symbol), 1, 3)) %>% relocate(final_aa, init_aa) %>% group_by(init_aa, final_aa) %>% summarize_if(is.numeric, sum) %>% ungroup() %>% mutate(Hugo_Symbol = paste0(init_aa, "-", final_aa)) %>% select(-c(init_aa, final_aa)) %>% relocate(Hugo_Symbol)) 

dfl <- x %>% lapply(., function(p) p %>% relocate(Hugo_Symbol, All))

#need to think about how to summarize this to get the total matrix
#each individual matrix seems in dfl seems like it should gets its own total row; 
all_cases_cols <- "../results/all-muts/all-counts-matrix.txt" %>% read_table() %>% colnames()

total_aa_specific <- lapply(dfl, function(p) align_columns_generic(all_cases_cols, p) %>% summarize_if(is.numeric, sum) %>% mutate(Hugo_Symbol = "Total") %>% relocate(Hugo_Symbol))

total_each_aa <- lapply(seq_along(dfl), function(t) align_columns_generic(all_cases_cols, dfl[[t]]) %>% summarize_if(is.numeric, sum) %>% mutate(Hugo_Symbol = "Total") %>% mutate(aa = aa_vec[t])) %>% do.call(rbind, .) %>% as_tibble()

total_nuc  <- total_each_aa %>% filter(aa %in% nucleophile_aa) %>% summarize_if(is.numeric, sum) %>% mutate(Hugo_Symbol = "Total") %>% relocate(Hugo_Symbol)

total_aa <- total_each_aa %>% summarize_if(is.numeric, sum) %>% mutate(Hugo_Symbol = "Total") %>% relocate(Hugo_Symbol)

#reweight for each aa
#max difference is ~2% (75 vs 78%) so not much of a difference between US and NPC rate

trans.aa.specific <- lapply(seq_along(x), function(idx) generic_reweight(x[[idx]], total_aa_specific[[idx]]) %>% mutate(aa = aa_vec[idx])) %>% do.call(rbind, .) %>% as_tibble()


trans.aa.specific <- trans.aa.specific %>% mutate(init_aa = str_sub(hugo, 1, 3)) %>% dplyr::rename(final_aa = aa) %>% filter(init_aa %in% standard_aa)

plt.aa.specific <- ggplot(trans.aa.specific %>% complete(init_aa, final_aa), aes(x=init_aa, y= final_aa, fill = pct_us)) + geom_tile() +  scale_fill_gradient(low = "lightgrey", high = "red", na.value = "lightgrey", breaks = scales::pretty_breaks(n = 8), limits = c(0, 85))  + labs(x = "", y="", fill = "Acquired\nMutation (%)") + theme_minimal() + theme(plot.title = element_text(size = 9, hjust = 0.5), axis.text.y = element_text(size = 9), axis.text.x = element_text(size = 10, angle = 45, hjust = 1)) + theme(panel.background = element_blank(), panel.grid = element_blank()) + theme(axis.ticks = element_line(color = "grey70"))  + ylab("Acquired Amino Acid") + xlab("Mutated Amino Acid") + theme(axis.title = element_text(size = 14)) + theme(legend.position = "right", legend.text = element_text(size = 8, vjust = 0.5), legend.key.size = unit(0.8, "cm"), legend.title = element_text(size = 10, hjust = 0.5)) 

#transistion probability to a specific nucleophile across all possible acquired amino acids for a given initial amino acid

plt.aa.specific.nuc <- ggplot(trans.aa.specific %>% complete(init_aa, final_aa) %>% filter(final_aa %in% nucleophile_aa), aes(x=init_aa, y= final_aa, fill = pct_us)) + geom_tile() +  scale_fill_gradient(low = "lightgrey", high = "red", na.value = "lightgrey", breaks = scales::pretty_breaks(n = 8), limits = c(0, 75))  + labs(x = "", y="", fill = "EG Mutation\nRate (%)") + theme_minimal() + theme(plot.title = element_text(size = 9, hjust = 0.5), axis.text.y = element_text(size = 9), axis.text.x = element_text(size = 10, hjust = 0.5)) + theme(panel.background = element_blank(), panel.grid = element_blank())  + ylab("Acquired Nucleophile") + xlab("Mutated Amino Acid") + theme(axis.title = element_text(size = 14)) + theme(legend.position = "left", legend.text = element_text(size = 8, vjust = 0.5), legend.key.size = unit(0.8, "cm"), legend.title = element_text(size = 10, hjust = 0.5)) 

#Supplemental Figure 2 X2Nuc
sf2 <- trans.aa.specific %>% complete(init_aa, final_aa) %>% filter(final_aa %in% nucleophile_aa)

sf2 <- sf2 %>% mutate(diff = pct_us - pct_tcga)


sf2a <- sf2 %>% ggplot(., aes(x=init_aa, y= final_aa, fill = diff)) + geom_tile() +  scale_fill_gradient(low = "lightgrey", high = "red", na.value = "lightgrey", breaks = scales::pretty_breaks(n = 8), limits = c(-2, 4))  + labs(x = "", y="", fill = "Difference\n(EG% - NPC%)") + theme_minimal() + theme(plot.title = element_text(size = 16, hjust = 0.5), axis.text.y = element_text(size = 9), axis.text.x = element_text(size = 10, hjust = 0.5)) + theme(panel.background = element_blank(), panel.grid = element_blank())  + ylab("Acquired Nucleophile") + xlab("Mutated Amino Acid") + theme(axis.title = element_text(size = 14)) + theme(legend.position = "left", legend.text = element_text(size = 8, vjust = 0.5), legend.key.size = unit(0.8, "cm"), legend.title = element_text(size = 10, hjust = 0.5)) + ggtitle("Transitions to a Specific Nucleophile Across All Possible Acquired Amino Acids\nFor a Given Initial Amino Acid") 


sf2b <- sf2 %>% ggplot(aes(x=final_aa, y = diff, color= final_aa, group = final_aa, fill = final_aa)) + geom_boxplot(alpha = 0.3, outlier.shape = NA) + theme_bw() + theme(legend.position = "none") + geom_jitter(width = 0.15, alpha = 0.5) + xlab("Acquired Nucleophile") + ylab("EG% - NPC%") + ylim(c(-2, 2))  + theme(axis.text = element_text(size=9), axis.title=element_text(size=12))


cowplot::plot_grid(sf2a, cowplot::plot_grid(NULL, sf2b, nrow = 2, rel_heights = c(0.1, 1)), nrow = 1, rel_widths = c(1, 0.4), labels = c("A", "B"), label_size = 18) %>% ggsave(plot = ., filename = "../x2nuc/sfig3.pdf", device = cairo_pdf, units = "in", width = 17, height = 7)

#Figure 2 X2Nuc

fig2a <- plt.aa.specific.nuc + coord_flip() 
fig2b <- fig1heatmap 
fig2 <- cowplot::plot_grid(fig2a, fig2b, nrow = 1, labels = c("A", "B"), rel_widths = c(1, 1.65), label_size = 20)
ggsave(filename = "../x2nuc/fig2.pdf", plot = fig2, units = "in", height = 10, width = 14, device = cairo_pdf)

#Supplemental Table 2 X2Nuc
trans.aa.specific %>% complete(init_aa, final_aa) %>% filter(final_aa %in% nucleophile_aa) %>% mutate_if(is.numeric, ~round(., 2)) %>% select(-hugo) %>% magrittr::set_colnames(c("Initial Amino Acid", "Acquired Nucleophile", "EG rate (%)", "NPC rate (%)", "lower bound EG rate (%)", "upper bound EG rate (%)")) %>% write_csv(., file = "../x2nuc/stable1.csv")
#


#reweight for nuc
nucidx <- data.frame(aa_vec = aa_vec) %>% mutate(idx = 1:n()) %>% filter(aa_vec %in% nucleophile_aa) %>% pull(idx)

trans.aa.nuc <- lapply(seq_along(nucidx), function(idx) generic_reweight(x[[nucidx[idx]]], total_nuc) %>% mutate(aa = aa_vec[nucidx[idx]])) %>% do.call(rbind, .) %>% as_tibble()


trans.aa.nuc <- trans.aa.nuc %>% mutate(init_aa = str_sub(hugo, 1, 3)) %>% dplyr::rename(final_aa = aa) %>% filter(init_aa %in% standard_aa)

#Transistion Probability to a Specific Nucleophile Across All Acquired Nucleophiles and all initial amino acids

plt.aa.nuc <- ggplot(trans.aa.nuc %>% complete(init_aa, final_aa), aes(x=init_aa, y= final_aa, fill = pct_us)) + geom_tile() +  scale_fill_gradient(low = "lightgrey", high = "red", na.value = "lightgrey", breaks = scales::pretty_breaks(n = 8), limits = c(0, 15))  + labs(x = "", y="", fill = "Mutation\n Rate (%)") + theme_minimal() + theme(plot.title = element_text(size = 16, hjust = 0.5), axis.text.y = element_text(size = 9), axis.text.x = element_text(size = 10, hjust = 1)) + theme(panel.background = element_blank(), panel.grid = element_blank())   + ylab("Acquired Nucleophile") + xlab("Mutated Amino Acid") + theme(axis.title = element_text(size = 14)) + theme(legend.position = "right", legend.text = element_text(size = 8, vjust = 0.5), legend.key.size = unit(0.8, "cm"), legend.title = element_text(size = 10, hjust = 0.5)) + ggtitle("Transitions to a Specific Nucleophile Across All Acquired Nucleophiles") 

#Supplemental Figure 3
ggsave(plot = plt.aa.nuc, filename = "../x2nuc/sfig4.pdf", device = cairo_pdf, units = "in", width = 14, height = 7)

#reweight for all aa
#27 transistions with a prevalence greater than 1%
#show as bar plot
#diff between tcga and US rate are small, the biggest is 0.4%
trans.aa.all <- lapply(seq_along(x), function(idx) generic_reweight(x[[idx]], total_aa) %>% mutate(aa = aa_vec[idx])) %>% do.call(rbind, .) %>% as_tibble()

trans.aa.all <- trans.aa.all %>% mutate(init_aa = str_sub(hugo, 1, 3)) %>% dplyr::rename(final_aa = aa)

plt.aa.all <- ggplot(trans.aa.all %>% complete(init_aa, final_aa), aes(x=init_aa, y= final_aa, fill = pct_us)) + geom_tile() +  scale_fill_gradient(low = "lightgrey", high = "red", na.value = "lightgrey", breaks = scales::pretty_breaks(n = 8), limits = c(0, 5))  + labs(x = "", y="", fill = "Acquired\nMutation (%)") + theme_minimal() + theme(plot.title = element_text(size = 9, hjust = 0.5), axis.text.y = element_text(size = 9), axis.text.x = element_text(size = 10, angle = 45, hjust = 1)) + theme(panel.background = element_blank(), panel.grid = element_blank()) + theme(axis.ticks = element_line(color = "grey70"))  + ylab("Acquired Amino Acid") + xlab("Mutated Amino Acid") + theme(axis.title = element_text(size = 14)) + theme(legend.position = "right", legend.text = element_text(size = 8, vjust = 0.5), legend.key.size = unit(0.8, "cm"), legend.title = element_text(size = 10, hjust = 0.5)) 


#### Fig 3

cancercontr <- list.files("../analyses", pattern = "*cancer-contr", full.names = TRUE, recursive = TRUE)

mutrates <- list.files("../analyses", full.names = TRUE, recursive = TRUE, pattern = "*mut-rates.rds")
generates <- list.files("../analyses", full.names = TRUE, recursive = TRUE, pattern = "*gene-rates.rds")

dfmut <- do.call(rbind, lapply(mutrates, function(x) readRDS(x) %>% mutate(aa = gsub(".*p.", "", hugo))))

dfmut <- dfmut %>% mutate(diff = pct_us - pct_tcga) %>% arrange(desc(pct_us))


aa_converter <- AMINO_ACID_CODE %>% as.data.frame() %>% magrittr::set_colnames("three_letter") %>% mutate(one_letter = rownames(.)) %>% as_tibble()

#Top 25 mutations with the largest difference between US estimate and TCGA estimate - the data show that we can just reference this in the paper and need not specifically have a figure on this

nuc <- dfmut %>% mutate(gene_label = str_sub(gene, 1, -2), final_aa = str_sub(hugo, -3, -1)) %>% filter(final_aa %in% nucleophile_aa) %>% arrange(desc(pct_us)) %>% dplyr::slice(1:25)

nuc <- nuc %>% mutate(gene_label = str_sub(gene, 1, -2), init_aa = str_sub(gsub("(.*)\\.", "", hugo), 1, 3), number = (str_extract_all(gsub("(.*)\\.", "", hugo), "\\d+")), final_aa = str_sub(hugo, -3, -1))

nuclong <- nuc

nuc$number <- lapply(nuc$number, function(x) x[[1]]) %>% unlist()

nuc <- nuc %>% left_join(aa_converter, by=c("init_aa" = "three_letter")) %>% left_join(aa_converter, by = c("final_aa"="three_letter")) %>% mutate(label_text = paste0(gene_label, " (", one_letter.x, number, one_letter.y, ")")) %>% select(-c(gene_label, init_aa, number, final_aa, one_letter.x, one_letter.y))

nuc_bar <- nuc 

mutbar_second_y_axis_label <- paste0("Estimated Number of New\nCases in the US, 2023")

label_order <- nuc_bar %>% arrange(desc(pct_us)) %>% pull(label_text)

nuc_bar <- nuc_bar %>% mutate(cases = round(pct_us * 0.01 * 1958310))


#mut.bar %>% ggplot(aes(x=label_text, y=pct_us)) + geom_bar(position = position_dodge(width = 5), stat="identity",color="grey50", fill="grey50", alpha = 0.75) + geom_errorbar(aes(x= label_text, ymin = pct_lb, ymax = pct_ub), color = "grey55", size = 0.2, width = 0.4) + theme_minimal() + theme(panel.grid = element_blank(), panel.border = element_rect(color="black", fill="transparent")) + geom_hline(yintercept = c(0.5, 1, 1.5), linetype="dashed", color="grey60") + scale_y_continuous(name = "Estimated Mutation Proportion\nfor U.S. Population (%)", limits=c(0,1.8), breaks=seq(0,2,0.5), sec.axis = sec_axis(~ . * 0.01 * 1958310, name = mutbar_second_y_axis_label, breaks = seq(0, 30000, 30000/6))) + theme(axis.ticks = element_line(color="black")) + xlab("Top 50 Point Mutations")  + theme(axis.title = element_text(size=18), axis.text.y = element_text(size = 16), axis.text.x = element_text(size = 12, angle = 70, hjust = 1)) 

fig3a <- nuc_bar %>% mutate(label_text = factor(label_text, levels = label_order)) %>% ggplot(aes(x=pct_us, y=label_text)) + geom_errorbar(aes(y = label_text, xmin = pct_lb, xmax = pct_ub), color = "grey55", size = 0.2, width = 0.4) + geom_bar(position = position_dodge(width = 5), stat="identity",color="grey70", fill="grey70", alpha = 0.75) + theme_minimal() + theme(panel.grid = element_blank(), panel.border = element_rect(color="black", fill="transparent")) + geom_vline(xintercept = c(1), linetype="dashed", color="grey60", alpha = 0.8) + scale_x_continuous(name = "Estimated Mutation Proportion\nfor U.S. Population (%)", limits=c(0,5), breaks=seq(0,5,0.5), sec.axis = sec_axis(~ . * 0.01 * 1958310, name = mutbar_second_y_axis_label, breaks = seq(0, 90000, 90000/9))) + theme(axis.ticks = element_line(color="black")) + ylab("Top 25 Acquired Nucleophiles")  + theme(axis.title = element_text(size=18), axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12, hjust = 1))  + geom_text(aes(x = pct_us + 0.27, y = label_text, label = cases), size = 3.5)



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

nuccontr <- cancercontr[data.frame(idx = cancercontr %>% dirname() %>% basename()) %>% mutate(i = 1:n()) %>% filter(idx %in% nucleophile_aa) %>% pull(i)]

cancercontr_df <- lapply(nuccontr, function(p) parse_cancer_contribution_file(p, nuclong, incidence_df))

t25cancercontr_df <- lapply(cancercontr_df, function(p) p$t25_df) %>% do.call(rbind, .) 

df_25_cc <- t25cancercontr_df %>% select(-label_text) %>% relocate(hugo) %>% left_join(.,nuc %>% select(hugo, label_text, pct_us)) %>% relocate(label_text, pct_us) %>% select(-hugo) 

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



short_levels <- data.frame(cancer_types = cancer_types, cancer_types_short = cancer_types_short) %>% mutate(cancer_types = factor(cancer_types, levels = (incidence_df %>% select(2,3) %>% arrange(desc(incidence)) %>% pull(name)) ) ) %>% as_tibble() %>% pull(cancer_types_short) 

df_25_cc <- df_25_cc %>% left_join(., data.frame(cancer_types = cancer_types, cancer_types_short = cancer_types_short) %>% mutate(cancer_types_short  = factor(cancer_types_short, levels = short_levels)), by=c("name.y" = "cancer_types")) 

#11 hieght by 8 width
fig3c <- ggplot(df_25_cc, aes(x= cancer_types_short, y= label_text, fill = value)) + geom_tile() +  scale_fill_gradient(low = "lightgrey", high = "red", na.value = "lightgrey", breaks = scales::pretty_breaks(n = 6))  + labs(x = "", y="", fill = "% of samples") + theme_minimal() + theme(plot.title = element_text(size = 9, hjust = 1), axis.text.y = element_text(size = 9), axis.text.x = element_text(size = 10, angle = 45, hjust = 1, vjust = 0.96),legend.title = element_text(size = 10), legend.text = element_text(size = 8), legend.key.size = unit(0.25, "cm")) + theme(panel.background = element_blank(), panel.grid = element_blank()) + theme(axis.ticks.x = element_line(color = "grey70"))  + ylab("Top 25 Acquired Nucleophiles") + xlab("Cancers with >1% Incidence") + theme(axis.title = element_text(size = 14)) 

##fig1b

nuc_all <- dfmut %>% mutate(gene_label = str_sub(gene, 1, -2), final_aa = str_sub(hugo, -3, -1)) %>% filter(final_aa %in% nucleophile_aa) %>% arrange(desc(pct_us))


t1p <- nuc_all %>% arrange(desc(pct_us)) %>% mutate(idx = 1:n()) %>% filter(idx < (nrow(nuc_all)*0.01))

points_mut <- t1p %>% arrange(desc(pct_us)) %>% mutate(us_rank = 1:n()) %>% filter(abs(diff) > 0.5 | us_rank < 7) %>% mutate(gene_label = str_sub(gene, 1, -2), init_aa = str_sub(gsub("(.*)\\.", "", hugo), 1, 3), number = str_extract_all(gsub("(.*)\\.", "", hugo), "\\d+") %>% unlist(), final_aa = str_sub(hugo, -3, -1)) %>% left_join(aa_converter, by=c("init_aa" = "three_letter")) %>% left_join(aa_converter, by = c("final_aa"="three_letter")) %>% mutate(label_text = paste0(gene_label, " (", one_letter.x, number, one_letter.y, ")")) %>% select(-c(gene_label, init_aa, number, final_aa, one_letter.x, one_letter.y)) %>% mutate(nx = -0.04) %>% mutate(ny = 0.02) %>% mutate(position = "right") %>% mutate(ny = ifelse(us_rank == 6, 0.06, ny)) %>% mutate(nx = ifelse(us_rank == 6, 0.03, nx)) %>% mutate(position = ifelse(us_rank == 6, "left", position))

fig3b <- t1p %>% ggplot(aes(x=pct_tcga, y=pct_us)) + geom_point(alpha=0.8, size=1.75, color="black") + theme_minimal() + geom_abline(color="grey70", size=1, alpha=0.8, linetype="dashed", slope = 1, intercept = 0) + theme(panel.grid = element_blank(), panel.border = element_rect(color="black", fill="transparent")) + scale_x_continuous(limits=c(0,5), breaks=seq(0, 5, 0.5)) + scale_y_continuous(limits=c(0,5), breaks=seq(0, 5, 0.5)) + ylab("Estimed Epidemio-Genomic (EG)\nMutation Rate (%)") + xlab("Estimated Naive Pan-Cancer (NPC) Mutation Rate (%)") + theme(axis.title=element_text(size=16), axis.text=element_text(size=14), axis.ticks=element_line(color="black")) + geom_text(data=points_mut, mapping=aes(x=pct_tcga+nx, y=pct_us+ny, label=label_text, hjust = position), lineheight = 0.5) + geom_errorbar(aes(ymin = pct_lb, ymax = pct_ub), size = 0.03, width = 0, color = "grey70")


#Figure 3 and Figure 4 X2Nuc
fig3 <- cowplot::plot_grid(fig3b, fig3a, nrow = 1, labels = c("A", "B"), label_size = 18, rel_widths = c(1, 1.3))

fig4 <- fig3c

ggsave(plot = fig3, filename = "../x2nuc/fig3.pdf", units = "in", height = 6, width = 15, device = cairo_pdf)

ggsave(plot = fig4, filename = "../x2nuc/fig4.pdf", units = "in", height = 8, width = 11, device = cairo_pdf)


#Supplemental Table 3 x2nuc
nuc_all %>% mutate_if(is.numeric, ~signif(., 3)) %>% magrittr::set_colnames(c("Mutation", "gene_tmp", "EG rate (%)", "NPC rate (%)", "lower bound EG rate (%)", "upper bound EG rate (%)", "aa", "diff", "Gene", "Acquired Nucleophile")) %>% select(-c(gene_tmp, aa)) %>% relocate(Mutation, Gene) %>% write_csv(., "../x2nuc/stable3.csv")

#Supplemental Table 5 x2nuc
df_25_cc %>% select(label_text, value, name.y) %>% magrittr::set_colnames(c("Mutation", "value", "name.y")) %>% pivot_wider(names_from = name.y, values_from = value) %>% write_csv(., "../x2nuc/stable5.csv")

cor.test(t1p$pct_us, t1p$pct_tcga) #Pearson correlation = 0.85 and p-val < 2.2e-16

#Correlation between tcga and us rate for top x% of US cases
df <- nuc_all %>% arrange(desc(pct_us)) %>% mutate(idx = 1:n())

pct_keep <- c(1, 10, 25, 50, 100)


#ranked by pct_us
get_corr_df_pearson <- function(t, df2){
  
  df2 <- df %>% filter(idx < (nrow(df) * t * (1/100)) )
  v <- cor.test(df2$pct_tcga, df2$pct_us)
  corr <- v$estimate %>% round(., 2) %>% unname()
  lb <- v$conf.int[1] %>% round(., 2) %>% unname()
  ub <- v$conf.int[2] %>% round(., 2) %>% unname()
  n_pts <- nrow(df2)
  data.frame(top_pct_of_data = t, corr = corr, lb = lb, ub = ub, p = "2.2e-16", n_points = n_pts, method = "pearson")

}

get_corr_df_spearman <- function(t, df2){
  
  df2 <- df %>% filter(idx < (nrow(df) * t * (1/100)) )
  v <- cor.test(df2$pct_tcga, df2$pct_us, method = "spearman")
  corr <- v$estimate %>% round(., 2) %>% unname()
  n_pts <- nrow(df2)
  data.frame(top_pct_of_data = t, corr = corr, p = "2.2e-16", n_points = n_pts, method = "spearman")

}

pearson <- do.call(rbind, lapply(pct_keep, function(q) get_corr_df_pearson(q, df2))) %>% as_tibble()
spearman <- do.call(rbind, lapply(pct_keep, function(q) get_corr_df_spearman(q, df2))) %>% as_tibble()

write_csv(pearson, "../x2nuc/stable4.csv")

# i manually tested these values and the p was < 2.2e-16
#top 1% p < 2.2e-16 and cor = 0.85
#top 10% p < 2.2e-16 and cor = 0.84
#top 25% p < 2.2e-16 and cor = 0.92 
#top 50% p < 2.2e-16 and cor = 0.79
#top 100 p < 2.2e-16 and cor = 0.73









#