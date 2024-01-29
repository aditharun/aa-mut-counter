#stitch figures together from count matrices

library(optparse)

option_list <- list(
  make_option(c("-a", "--aminoacid"), type = "character", default = "Cys", 
              help = "3 letter amino acid"),
  make_option(c("-l", "--fullname"), type = "character", default = "Cysteine", 
              help = "Full amino acid name")
)


opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

amino_acid <- opt$aminoacid
amino_acid_long <- opt$fullname

system(paste0("Rscript reweighting.R ", amino_acid))
system(paste0("Rscript create-tables.R ", amino_acid))

system(paste0("Rscript fig2.R ", amino_acid, " ", amino_acid_long))
system(paste0("Rscript fig3.R ", amino_acid))
system(paste0("Rscript genes-fig.R ", amino_acid))


if (amino_acid == "Cys"){
	system(paste0("Rscript fig1.R ", amino_acid))
	system(paste0("Rscript data-points-compute.R ", amino_acid))
}

