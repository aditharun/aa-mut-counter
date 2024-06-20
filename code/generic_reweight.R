#A generic reweighting structure



#df: list of matrices with each row being a mutation/gene/transition we care about, and each column being a rosetta code except for the first two which are Hugo_Symbol and All
#total: matrix with 1xn with first two columns Hugo_Symbol and All where Hugo Symbol is "Total" and the values for each rosetta code is the total (however, we'd like to define)

generic_reweight <- function(df, total){

	seer.df <- "../data/seer-abundances.xlsx" %>% read_excel() %>% magrittr::set_colnames(c("rosetta", "cancer", "incidence")) %>% mutate(incidence_frac = incidence / 100)

	all_cases_cols <- "../results/all-muts/all-counts-matrix.txt" %>% read_table() %>% colnames()


	align_columns <- function(all_cases_cols, df){

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

	#if only 1 element in df then do not send through lapply()
	if (class(df)[1] == "list"){
		df2 <- lapply(df, function(x) align_columns(all_cases_cols, x)) %>% do.call(rbind, .) %>% as_tibble()
		print("s")
	} else{
		df2 <- align_columns(all_cases_cols, df) %>% as_tibble()
	}

	if ( (df2$All %>% unique() %>% sum()) == 0){
		alltmp <- df2 %>% select(-c(Hugo_Symbol, All)) %>% rowSums() %>% unname() %>% unlist()
		df2$All <- alltmp
	}


	if ( (total$All %>% unique() %>% sum()) == 0){
		alltmp <- total %>% select(-c(Hugo_Symbol, All)) %>% rowSums() %>% unname() %>% unlist()
		total$All <- alltmp
	}

	total <- total %>% pivot_longer(-c(Hugo_Symbol, All)) %>% left_join(seer.df, by=c("name"="rosetta")) 

	total <- total %>% magrittr::set_colnames(c("hugo", "all", "rosetta", "count", "cancer", "incidence", "incidence_frac"))

	tcga.aa <- df2 %>% select(-All) %>% pivot_longer(-Hugo_Symbol) %>% left_join(total %>% select(-c(count, hugo, all)), by=c("name"="rosetta"))

	tcga.aa <- tcga.aa %>% magrittr::set_colnames(c("hugo", "rosetta", "count", "cancer", "incidence", "incidence_frac"))

	total <- total %>% mutate(weight.tot = incidence_frac * count)
	tcga.aa <- tcga.aa %>% mutate(weight.aa=incidence_frac * count) 
	tcga.aa <- tcga.aa %>% mutate(idx = 1:n())
	tcga.aa.confint <- tcga.aa %>% filter(count > 0)  %>% mutate(result = purrr::map(count, conf_int), lb = purrr::map_dbl(result, "lb"), ub = purrr::map_dbl(result, "ub")) %>% select(-result)
	tcga.aa <- tcga.aa %>% left_join(tcga.aa.confint %>% select(idx, lb, ub), by=c("idx"="idx"))
	tcga.aa <- tcga.aa %>% mutate(weight.lb = lb * incidence_frac, weight.ub = ub * incidence_frac)
	tcga.aa <- tcga.aa %>% mutate(weight.lb = ifelse(is.na(weight.lb), 0, weight.lb), weight.ub = ifelse(is.na(weight.ub), 0, weight.ub))
	aa <- tcga.aa %>% left_join(total %>% select(rosetta, weight.tot, all), by=c("rosetta"="rosetta")) %>% group_by(hugo) %>% summarize(pct_us = sum(weight.aa) / sum(weight.tot), pct_tcga = sum(count) / unique(all) , pct_lb = sum(weight.lb) / sum(weight.tot), pct_ub = sum(weight.ub) / sum(weight.tot)) %>% mutate(across(where(is.double), ~ . * 100)) 

	aa
}