library(tidyverse)
library(optparse) 
library(gdata)

# terminal inputs
option_list = list(
  make_option(c( "-r", "--rep"), type = "character", default = '1', 
              help = "rep number for simulations", metavar = "character"),
  make_option(c("-p", "--pval"), type = "character", default = '1',
              help = "Comma-separated list of pvals for filter threshold", metavar = "character"),
  make_option(c("-s", "--sim"), type = "character", default = '1',
              help = "simulation case: 1,2,3,..", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
# get values
rep <- opt$r
print("input values: rep")
print(rep)
case <- opt$sim # '1'
print("sim case:")
print(case)
pval_vector = as.numeric(unlist(strsplit(opt$p, ",")))
print("input values: list of pvals")
print(pval_vector)
setwd(paste0('/hpc/group/ochoalab/tt207/platform_bias/sim', case, '/rep', rep, '/af-test'))

# read af outputs
p1_S1 = read.table("p1_S1.acount") %>% dplyr::rename(CHR = V1, POS = V2, A1 = V3, A2 = V4, AC = V5, N = V6)
p1_S2 = read.table("p1_S2.acount") %>% dplyr::rename(CHR = V1, POS = V2, A1 = V3, A2 = V4, AC = V5, N = V6)
p1_S3 = read.table("p1_S3.acount") %>% dplyr::rename(CHR = V1, POS = V2, A1 = V3, A2 = V4, AC = V5, N = V6)
#p1_af = gdata::combine(p1_S1, p1_S2, p1_S3) %>% dplyr::rename(CHR = V1, POS = V2, A1 = V3, A2 = V4, AC = V5, N = V6)
p2_S1 = read.table("p2_S1.acount") %>% dplyr::rename(CHR = V1, POS = V2, A1 = V3, A2 = V4, AC = V5, N = V6)
p2_S2 = read.table("p2_S2.acount") %>% dplyr::rename(CHR = V1, POS = V2, A1 = V3, A2 = V4, AC = V5, N = V6)
p2_S3 = read.table("p2_S3.acount") %>% dplyr::rename(CHR = V1, POS = V2, A1 = V3, A2 = V4, AC = V5, N = V6)
#p2_af = gdata::combine(p2_S1, p2_S2, p2_S3) %>% dplyr::rename(CHR = V1, POS = V2, A1 = V3, A2 = V4, AC = V5, N = V6)

p2_alleles = p2_S1 %>% select(CHR, POS, A1, A2)
# library(devtools)
# install_github('OchoaLab/popgeninfer')
library(popgeninfer)
x1 <- cbind(p1_S1$AC, p1_S2$AC, p1_S3$AC)
n1 <- cbind(p1_S1$N, p1_S2$N, p1_S3$N)
x2 <- cbind(p2_S1$AC, p2_S2$AC, p2_S3$AC)
n2 <- cbind(p2_S1$N, p2_S2$N, p2_S3$N)
pvals_af_test <- af_test( x1, n1, x2, n2 )$pval
pvals_af = cbind(p2_alleles, pvals_af_test)
#save(pvals_af, file = "AF_all.RData")

#load("AF_OR.RData")
# reverse comp
comp <- function(x) chartr("ATGC","TACG",x)

indexes_comp <- pvals_af$A1 == comp( pvals_af$A2 )
m_comp = sum(indexes_comp)
m_noncomp <- sum( !indexes_comp )
ids_comp <- pvals_af$POS[ indexes_comp ]
# get subset of "forward" p-values for "comp" cases
pvals_comp_fwd <- pvals_af_test[ indexes_comp ]
p2_alleles_sub = p2_alleles[indexes_comp, ]
pvals_comp_fwd_af = cbind(p2_alleles_sub, pvals_comp_fwd) 

########## if using original AF for plotting
# ids_comp = read.table("../id_reverse_comp.txt") %>% pull()
# pvals_comp_fwd <- pvals_af_test[ ids_comp ]
# p2_alleles_sub = p2_alleles[ids_comp, ]
# pvals_comp_fwd_af = cbind(p2_alleles_sub, pvals_comp_fwd) 
# indexes_comp <- rep(FALSE, length(pvals_af_test))

# Set positions in ids_comp to TRUE
indexes_comp[ids_comp] <- TRUE


#save(pvals_comp_fwd_af, file = "AF_forward.RData")
# now calculate "reverse" orientation p-values, for "comp" cases only
# subset to comp SNPs
x1c <- x1[ indexes_comp,]
n1c <- n1[ indexes_comp,]
x2c <- x2[ indexes_comp,]
n2c <- n2[ indexes_comp,]

# perform reversed test, reverse second platform only 
pvals_comp_rev <- af_test( x1c, n1c, n2c - x2c, n2c )$pval
pvals_comp_rev_af = cbind(p2_alleles_sub, pvals_comp_rev)
#save(pvals_comp_rev_af, file = "AF_reversecomp.RData")

comp_fwd_rev = merge(pvals_comp_fwd_af, pvals_comp_rev_af, by = c("CHR", "POS", "A1", "A2")) #%>% distinct(POS, .keep_all = TRUE)
####################

pcut_values = pval_vector
summary_df <- data.frame(
  pcut = numeric(),
  remove = integer(),
  flip = integer(),
  stringsAsFactors = FALSE
)
#pcut_values = 1e-04
for (pcut in pcut_values) {
  print(pcut)
  ############ reverse comp
  # classify SNPs based in their p-values (the same for all ancestries!), using this decision tree
  category_comp <- ifelse( comp_fwd_rev$pvals_comp_fwd < pcut & comp_fwd_rev$pvals_comp_rev < pcut, "Removed",
                           ifelse( comp_fwd_rev$pvals_comp_fwd < comp_fwd_rev$pvals_comp_rev, "Flipped", "Unchanged" ) )

  category_comp <- factor(category_comp, levels = c("Unchanged", "Flipped", "Removed"))
  
  # comp df
  category_comp_df = cbind(comp_fwd_rev, category_comp)
  print(table(category_comp_df$category_comp))
  # Extract POS (indices) for "flip"
  pos_flip <- subset(category_comp_df, category_comp == "Flipped")$POS

  # Extract POS (indices) for "remove"
  pos_remove <- subset(category_comp_df, category_comp == "Removed")$POS

  ############ non-reverse comp
  pvals_noncomp <- pvals_af_test[ !indexes_comp ]
  
  # classify SNPs (same for all ancestries)
  category_noncomp <- ifelse( pvals_noncomp < pcut , "Removed", "Unchanged")
  # set order for plotting
  category_noncomp <- factor(category_noncomp, levels = c("Unchanged", "Removed"))
  category_noncomp_df = cbind(p2_alleles[!indexes_comp,], category_noncomp)
  table(category_noncomp_df$category_noncomp)
  # Extract POS (indices) for "remove"
  pos_remove_noncomp <- unique(subset(category_noncomp_df, category_noncomp == "Removed")$POS)
  
  # write SNP id's to remove/flip into text file
  write_lines( unique(c( pos_remove, pos_remove_noncomp )), paste0("remove_", pcut, ".txt"))
  write_lines( unique(pos_flip), paste0("flip_", pcut, ".txt"))
  
  # Add results to dataframe
  summary_df <- rbind(summary_df, data.frame(
    pcut = pcut,
    remove = length(c(pos_remove, pos_remove_noncomp)),
    flip = length(pos_flip),
    keep = length(c(subset(category_comp_df, category_comp == "Unchanged")$POS, unique(subset(category_noncomp_df, category_noncomp == "Unchanged")$POS))),
    stringsAsFactors = FALSE
  ))
}

write.table(summary_df, "summary.txt", row.names = FALSE, quote = FALSE)
