library(tidyverse)
#library(karyoploteR)
library(genio)
library(optparse) 


# terminal inputs
option_list = list(
  make_option(c( "-p", "--pval"), type = "character", default = '5e-8', 
              help = "pvalue threshold for identifying significant snps", metavar = "character"),
  make_option(c( "-r", "--rep"), type = "character", default = '1', 
              help = "rep number for simulations", metavar = "character"),
  make_option(c("-s", "--sim"), type = "character", default = '1',
              help = "simulation case: 1,2,3,..", metavar = "numeric")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
# get values

pval <- opt$p
rep <- opt$r
case <- opt$sim 
print("input values: pvalue threshold, rep, sim case")
print(pval)
print(rep)
print(case)
setwd(paste0('/hpc/group/ochoalab/tt207/platform_bias/sim', case, '/rep', rep, '/saige/', pval))
pval = as.numeric(pval)
###########################################
# identify significant SNPs in current iteration

# identify SNPs that can be flipped
remove_snps = read.table("platform_snps_phase1.txt", header = TRUE)  
if (case == 1) {
  data <- read_plink('../../1G_6000n_500000m_controls')
} else {
  data <- read_plink('../../30G_6000n_500000m_controls_bias')
}
bim <- data$bim
comp <- function(x) chartr("ATGC","TACG",x)
bim_reverse_comp = bim[which(bim$ref == comp(bim$alt)),]
sig_snps_flip = remove_snps %>% filter(x %in% bim_reverse_comp$id) %>% pull(x)
sig_snps_remove = remove_snps %>% filter(!x %in% bim_reverse_comp$id) %>% pull(x)

# write output and rerun saige
write.table(sig_snps_flip, "platform_phase2_flip.txt", row.names = FALSE, quote = FALSE)
write.table(sig_snps_remove, "platform_phase2_remove.txt", row.names = FALSE, quote = FALSE)

print("phase 1 remove SNPs:")
print(length(remove_snps %>% pull(x)))
print("phase 2 flip SNPs:")
print(length(sig_snps_flip))
print("phase 2 remove SNPs:")
print(length(sig_snps_remove))


