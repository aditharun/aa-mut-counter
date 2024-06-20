library(tidyverse)
library(agvgd)

datadir <- "results"

allcounts <- file.path(datadir, "all-muts", "all-counts-matrix.txt") %>% read_tsv()

total <- allcounts %>% filter(Hugo_Symbol == "Total")

files <- list.files("results", pattern = "*.txt", full.names = TRUE)
aa_files <- files[!grepl("mutations", files)] 

extract_transition_matrix <- function(filepath, total_cols){
	final_aa <- filepath %>% basename() %>% str_sub(., 1, -5) %>% str_sub(., 8, -1)
	test <- filepath %>% read_tsv()
	test <- test %>% mutate(init_aa = sub(".*p\\.", "", Hugo_Symbol)) %>% mutate(init_aa = str_sub(init_aa, 1, 3)) %>% mutate(final_aa = final_aa) %>% relocate(init_aa, final_aa) %>% select(-Hugo_Symbol)
	test <- test %>% group_by(init_aa, final_aa) %>% summarise_all(sum) 

	test <- test %>% ungroup()

	cols_to_add <- setdiff(total_cols, colnames(test)[-c(1:2)])

	if (!identical(cols_to_add, character(0))){
		test[,cols_to_add] <- 0
	}

	test

}

data <- lapply(aa_files, function(x) extract_transition_matrix(x, colnames(total)[-1]))

df <- do.call(rbind,  data) %>% as_tibble()

df <- df %>% select(-All)

df.l <- df %>% pivot_longer(-c(init_aa, final_aa)) %>% magrittr::set_colnames(c("init_aa", "final_aa", "rosetta", "count"))

total.l <- total %>% select(-All) %>% pivot_longer(-Hugo_Symbol) %>% select(-Hugo_Symbol) %>% magrittr::set_colnames(c("rosetta", "total"))

rate.df <- df.l %>% left_join(total.l) %>% mutate(rate = count / total) 



#70317 is breast cancer

#4 columns by 2 rows 
# need to think about design here, do they all have their own color scheme! YES

gi_tumor_codes <- c("70367", "70387", "70187", "70347", "70337", "70377", "70217", "81603")
gi_tumor_names <- c("Colorectal Adenocarcinoma", "Pancreatic Adenocarcinoma", "Hepatocellular Carcinoma", "Gastric Adenocarcinoma", "Esophageal Adenocarcinoma", "Gallbladder Cancer", "Esophageal Squamous Cell Carcinoma", "Cholangiocarcinoma")

levels_tumor_names <- c("Esophageal Adenocarcinoma", "Esophageal Squamous Cell Carcinoma", "Gastric Adenocarcinoma", "Gallbladder Cancer", "Cholangiocarcinoma", "Hepatocellular Carcinoma", "Pancreatic Adenocarcinoma", "Colorectal Adenocarcinoma")

gi_df <- data.frame(rosetta = gi_tumor_codes, names = gi_tumor_names) %>% as_tibble()

gi_df <- gi_df[match(levels_tumor_names, gi_df$names), ]


create_heatmap <- function(rate.df, gi_df, index){

	code <- gi_df$rosetta[index]
	cancer_name <- gi_df$names[index]

	total_grid <- expand.grid(init_aa = agvgd::amino_acids("three_letter"), final_aa = agvgd::amino_acids("three_letter")) %>% as_tibble() 

	cancer.df <- rate.df %>% filter(rosetta == code) 

	mat <- total_grid %>% left_join(cancer.df) %>% mutate(rate = ifelse(is.na(rate), 0, rate))

	mat$init_aa <- factor(mat$init_aa, levels = agvgd::amino_acids("three_letter") %>% sort())
	mat$final_aa <- factor(mat$final_aa, levels = agvgd::amino_acids("three_letter") %>% sort() %>% rev())

	plt <- mat %>% ggplot(., aes(x = init_aa, y = final_aa, fill = rate)) +  geom_tile() +  scale_fill_gradient(low = "lightgrey", high = "red", na.value = "lightgrey", breaks = scales::pretty_breaks(n = 6)) + labs(x = "Initial Amino Acid", y = "Final Amino Acid", fill = "Average #\nof Mutations") +  theme_minimal() + theme(panel.grid = element_blank()) + ggtitle(cancer_name) + theme(plot.title = element_text(size = 12, hjust = 0.5), axis.text = element_text(size = 6), axis.title = element_text(size = 9), legend.title = element_text(size = 10), legend.text = element_text(size = 8), legend.key.size = unit(0.25, "cm"))

	if (!index %in% c(1, 5)){
		plt <- plt + theme(legend.title = element_blank())
	}

	plt

}


plt.list <- lapply(1:nrow(gi_df), function(x) create_heatmap(rate.df, gi_df, x))

cowplot::plot_grid(plotlist = plt.list, nrow = 2, ncol = 4, rel_widths = c(1.08, 1, 1, 1)) %>% ggsave(filename = "~/Downloads/mutation-landscape-of-gi-tumors.pdf", plot = ., units = "in", height = 6, width = 20, device = cairo_pdf)















#