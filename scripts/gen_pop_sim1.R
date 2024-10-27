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
              help = "simulation case: 1,2,3,..", metavar = "numeric")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
# get values
rep_num <- opt$num # '1'
case <- opt$sim # '1'

source("/hpc/group/ochoalab/tt207/platform_bias/make_corr_psi_tree.R")
setwd(paste0('/hpc/group/ochoalab/tt207/platform_bias/sim', case))
new_dir <- paste0('rep', rep_num)
dir.create(new_dir)
setwd(new_dir)

# define population structure
# 3 subpopulations
n1 <- 2000
n2 <- 2000
n3 <- 2000

# here's the labels (for simplicity, list all individuals of S1 first, then S2, then S3)
labs <- c(
  rep.int('S1', n1),
  rep.int('S2', n2),
  rep.int('S3', n3)
)

n_ind = length(labs)
# total number of individuals
n_ind

# create random 1 male/2 female for sex 
sex = rbinom(length(labs), 1, 0.5) +1
table(sex)

# number of subpopulations "k_subpops"
k_subpops <- length(unique(labs))
k_subpops

# get its current coancestry, again to be rescaled in the end
tree_subpops <- make_corr_psi_tree( k_subpops )
coanc_subpops <- coanc_tree( tree_subpops )

# construct Q matrix and scale things appropriately
Fst <- 0.3

Q <- admix_prop_indep_subpops(labs)
Fst_temp <- fst_admix(Q, coanc_subpops)
coanc_factor = Fst/Fst_temp
coanc_subpops = coanc_subpops*coanc_factor

# desired admixture matrix
tree_subpops <- scale_tree( tree_subpops, coanc_factor )

# this reports the final edge lengths shown in the paper
tree_subpops$edge.length
# [1] 0.1800000 0.1800000 0.2195122 0.2195122

# and the total covariances/FSTs
coanc_subpops
##      S1   S2   S3
## S1 0.18 0.00 0.00
## S2 0.00 0.36 0.18
## S3 0.00 0.18 0.36

# save parameters
write_grm( 'Psi', coanc_subpops )
# Q = admix prop
write_matrix( 'Q', Q , ext = 'txt.gz' )

# the true coancestry matrix
coancestry <- coanc_admix( Q, coanc_subpops )
write_grm( 'Theta', coancestry )

####### draw genotypes
m_loci <- 500000
#m_loci <- 50
reverse_comp_percentage <- 0.1
# flip_percentage <- 0.04
# percent_bias_snps <- 0.06
###### create bim file
# create a dummy tibble with the right columns
library(tibble)
temp_bim <- tibble(
  chr = rep(1, m_loci),
  id = 1:m_loci,
  posg = rep(0, m_loci),
  pos = 1:m_loci,
  alt = rep(0, m_loci),
  ref = rep(0, m_loci)
)

source('/hpc/group/ochoalab/tt207/platform_bias/generate_bim.R')
bim_new = generate_bim(temp_bim, reverse_comp_percentage)
# check number of reverse comp
bim_reverse_comp = bim_new[which(bim_new$ref == comp(bim_new$alt)),]
# write file with reverse comp id's
write.table(bim_reverse_comp$id, "id_reverse_comp.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)

########################
# ancestral AFs
p_anc <- draw_p_anc( m_loci )
####################### platform 1 (unbiased)
# independent subpops (intermediate) AFs
p_subpops <- draw_p_subpops_tree( p_anc, tree_subpops )
# genotypes
X <- draw_genotypes_admix( p_subpops, admix_proportions = Q )

# write bim/bed/fam
fam <- tibble(fam = labs)
fam <- genio::make_fam(fam) %>% select(-sex, -pheno)
# s3 = fam_new %>% filter(fam == "S3") %>% pull(id)
# write.table(s3, "/hpc/group/ochoalab/tt207/platform_bias/S3_2000.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)

# create 'trait' column 0/1 to indicate the two platforms
n <- c(n1, n2, n3)
trait = unlist(lapply(n, function(n) rep(0:1, each = n / 2)))
fam_new = cbind(fam, sex, trait) %>% dplyr::rename(pheno = trait)
fam_order = fam_new %>% arrange(pheno)
write_fam("1G_6000n_500000m_controls_bias", fam_order)
write_bim("1G_6000n_500000m_controls_bias", bim_new)
write_plink("1G_6000n_500000m_controls", X, fam = fam_new, bim = bim_new)

# write platform individual id's into text file:
id_p1 = fam_new %>% filter(pheno == 0) %>% select(FID = fam, IID = id)
id_p2 = fam_new %>% filter(pheno == 1) %>% select(FID = fam, IID = id)
write.table(id_p1, "id_p1.txt", col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(id_p2, "id_p2.txt", col.names = TRUE, row.names = FALSE, quote = FALSE)

# write covar
covar_file = fam_new %>% select(-pat, -mat) %>% dplyr::rename(famid = fam, iid = id)
write.table(covar_file, "1G_covar_6000n_500000m_controls.txt", col.names = TRUE, row.names = FALSE, quote = FALSE)
# p_anc 
write.table(p_anc, "1G_p_anc_6000n_500000m_controls.txt", sep = " ", col.names = TRUE, row.names = FALSE)
# and true ancestral allele frequencies
write_matrix( 'P', p_subpops, ext = 'txt.gz' )


##### run analysis with saige approach and af-test appraoch
dir.create('saige')
dir.create('af-test')
