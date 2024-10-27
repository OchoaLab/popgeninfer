library(tidyverse)
#library(karyoploteR)
library(genio)
library(optparse) 


# terminal inputs
option_list = list(
  make_option(c( "-i", "--iteration"), type = "character", default = '1', 
              help = "numeric number for iteration of saige control vs control", metavar = "character"),
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
iter <- opt$i
pval <- opt$p
rep <- opt$r
case <- opt$sim # '1'
print("input values: iteration, pvalue threshold, rep, and sim case")
print(iter)
print(pval)
print(rep)
print(case)
setwd(paste0('/hpc/group/ochoalab/tt207/platform_bias/sim', case, '/rep', rep, '/saige/', pval))
pval = as.numeric(pval)
###########################################
# identify significant SNPs in current iteration

# Determine the input file based on iteration
input_file <- ifelse(iter == "phase2", 
                     "platform_phase2_output.txt", 
                     paste0('platform_phase3_', iter, '_output.txt'))

next_iter <- ifelse(iter == "phase2", "0", as.numeric(iter)+1)
# Read the significant SNPs and filter them
sig_snps <- read.table(input_file, header = TRUE) %>% 
  filter(p.value < pval) %>% 
  pull(MarkerID)

# Determine the output file for exclusion list
#next_iter <- as.numeric(iter) + 1
output_file <- paste0('platform_sig_snps_remove_phase3_iter', next_iter, '.txt')

# Write the filtered SNPs to the exclusion file
write.table(sig_snps, output_file, row.names = FALSE, quote = FALSE)

# print("phase 3 remove SNPs:")
# print(length(sig_snps))

# sig_snps_p3 = read.table("platform_phase3_output.txt", header = TRUE) %>% filter(p.value < pval) #%>% pull(MarkerID)
# # write output and rerun saige
# data <- read_plink('1G_3000n_500000m_controls_phase3')
# bim <- data$bim
# comp <- function(x) chartr("ATGC","TACG",x)
# bim_reverse_comp = bim[which(bim$ref == comp(bim$alt)),]
# sig_snps_flip_p3 = sig_snps_p3 %>% filter(MarkerID %in% bim_reverse_comp$id) %>% pull(MarkerID)
# sig_snps_remove_p3 = sig_snps_p3 %>% filter(!MarkerID %in% bim_reverse_comp$id) %>% pull(MarkerID)
# 
# length(intersect(sig_snps_flip_p3, phase1_flip$x))
# length(intersect(sig_snps_remove_p3, phase1_flip$x))
# write.table(sig_snps, "platform_sig_snps_remove_phase3_iter1.txt", row.names = FALSE, quote = FALSE)
# 
# sig_snps_p3_1 = read.table("platform_phase3_1_output.txt", header = TRUE) %>% filter(p.value < pval) #%>% pull(MarkerID)

