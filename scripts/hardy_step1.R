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
# identify significant SNPs from hwe
hardy = read.table("bias_hwe.hardy")
colnames(hardy) = c("CHR", "SNP", "N_OBS", "O(HOM1)", "O(HET)", "O(HOM2)", 
                    "E(HOM1)", "E(HET)", "Chi^2", "P")
sig_snps = hardy %>% filter(P < pval)

# identify SNPs that can be flipped
if (case == 1){
  data <- read_plink('../../1G_6000n_500000m_controls')
} else {
  data <- read_plink('../../30G_6000n_500000m_controls_bias')
}
bim <- data$bim
comp <- function(x) chartr("ATGC","TACG",x)
bim_reverse_comp = bim[which(bim$ref == comp(bim$alt)),]
sig_snps_flip = sig_snps %>% filter(SNP %in% bim_reverse_comp$id) %>% pull(SNP)
sig_snps_flip_pval = sig_snps %>% filter(SNP %in% bim_reverse_comp$id) %>% pull(P)
sig_snps_remove = sig_snps %>% filter(!SNP %in% bim_reverse_comp$id) %>% pull(SNP)

# write output and test flips in reverse orientation
write.table(sig_snps_flip, "hwe_flip_fwd.txt", row.names = FALSE, quote = FALSE, col.names = FALSE)
write.table(sig_snps_flip_pval, "hwe_flip_fwd_pval.txt", row.names = FALSE, quote = FALSE, col.names = FALSE)
write.table(sig_snps_remove, "hwe_remove.txt", row.names = FALSE, quote = FALSE, col.names = FALSE)
print('remove and flip files saved')