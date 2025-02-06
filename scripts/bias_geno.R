library(tidyverse)
library(simgenphen)
library(bnpsd)
library(RColorBrewer)
library(popkin)
library(genio)
library(optparse) 

# terminal inputs
option_list = list(
  make_option(c("-n", "--num"), type = "character", default = '1', 
              help = "numeric number that indicates number of rep", metavar = "character"),
  make_option(c("-s", "--sim"), type = "character", default = '1',
              help = "simulation case: 1,2,3,..", metavar = "numeric"),
  make_option(c("-f", "--file"), type = "character", default = '1G_6000n_500000m_controls',
              help = "file name of plink file", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
# get values
rep_num <- opt$num # '1'
case <- opt$sim # '1'
file_name <- opt$file

setwd(paste0('/hpc/group/ochoalab/tt207/platform_bias/sim', case, '/rep', rep_num))
# read plink files
data = read_plink( file_name )
X = data$X
# read platform 2 individual id's, filter X genotype matrix with p2 id's
p2_id = read.table('id_p2.txt', header = TRUE) %>% pull(IID)
p2_X = X[, colnames(X) %in% p2_id]
p1_X = X[, !colnames(X) %in% p2_id]

# 6% of SNPs biased
percent_bias_snps = 0.06
n_select <- round(percent_bias_snps * nrow(p2_X))
selected_rows <- sample.int(nrow(p2_X), n_select)
write.table(selected_rows, "id_selected_rows.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)

# 4 % of SNPs flipped
# identify reverse complements in p2 for 'flip'
id_reverse_comp = read.table('id_reverse_comp.txt') %>% pull()
flip_percentage <- 0.04
selected_flip = sample(id_reverse_comp, round(flip_percentage * length(id_reverse_comp)))
write.table(selected_flip, "id_selected_flip.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)

# copy original matrix, since most data isn't edited
p2_biased_matrix <- p2_X

# apply bias
for (snp in selected_rows) {
    # w to select probability of an individual having the biased SNP
    w = rbeta(1, shape1 = 2, shape2 = 5)
    # draw once the biased genotype, the same value is used for all biased cases
    x_biased = 2 * rbinom(1, 1, 0.5)
    # bias these only
    selected_col = sample.int(ncol(p2_X), ncol(p2_X) * w)
    # apply the bias to those individuals in this row
    p2_biased_matrix[ snp, selected_col ] <- x_biased
}

# apply all flips together now
p2_biased_matrix[ selected_flip, ] <- 2 - p2_biased_matrix[ selected_flip, ]


# combine two platforms again and sort of individual id's
p_merge = cbind(p1_X, p2_biased_matrix)
write_bed(file_name, p_merge, verbose = TRUE)
