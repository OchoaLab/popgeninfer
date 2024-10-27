library(tidyverse)
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
case <- opt$sim # '1'
print("input values: pvalue threshold, rep, and simulation")
print('pval:')
print(pval)
print('rep:')
print(rep)
print('sim:')
print(case)
setwd(paste0('/hpc/group/ochoalab/tt207/platform_bias/sim', case, '/rep', rep, '/hwe-test/', pval))
pval = as.numeric(pval)

###############################################
# for flips we need to test forward and reverse:
flip_rev = read.table("hwe_flip_reverse.hardy")
colnames(flip_rev) = c("CHR", "SNP", "N_OBS", "O(HOM1)", "O(HET)", "O(HOM2)", 
                       "E(HOM1)", "E(HET)", "Chi^2", "P")

flip_fwd = read.table('hwe_flip_fwd_pval.txt') %>% pull()

flip_snp = cbind(flip_rev %>% dplyr::rename(flip_rev = P), flip_fwd) %>% 
  mutate(category = ifelse(flip_rev > flip_fwd, "Flipped", ifelse(flip_rev < pval & flip_fwd < pval, "Removed", "Unchanged")))

table(flip_snp$category)
# write output and perform flip and removal
final_flip = flip_snp %>% filter(category == "Flipped") %>% pull(SNP)
write.table(final_flip, "hwe_flip_final.txt", row.names = FALSE, quote = FALSE, col.names = FALSE)
sig_snps_remove = read.table("hwe_remove.txt") %>% pull(V1)
final_remove = c(sig_snps_remove, flip_snp %>% filter(category == "Removed") %>% pull(SNP))
write.table(final_remove, "hwe_remove_final.txt", row.names = FALSE, quote = FALSE, col.names = FALSE)
print('finalized remove and flip files saved')
