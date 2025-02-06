library(bnpsd)
library(simfam)
library(genio)
library(simgenphen)
library(popkin)
library(tidyverse)
library(optparse) 

# terminal inputs
option_list = list(make_option(c( "-n", "--num"), type = "character", default = '1', 
                               help = "numeric number that indicates number of rep", metavar = "character"),
                   make_option(c("-s", "--sim"), type = "character", default = '1',
                               help = "simulation case: 1,2,3,..", metavar = "numeric"))

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
# get values
rep_num <- opt$num # '1'
case <- opt$sim # '1'
setwd(paste0('/hpc/group/ochoalab/tt207/platform_bias/sim', case, '/rep', rep_num))


# draw trait
plink_file = paste0('30G_6000n_500000m_controls')
plink = read_plink(plink_file)
print('finish reading plink file')
dim(plink$X)

#X_G = plink$X
fam_G = plink$fam %>% group_by(fam) %>%
  mutate(pheno = case_when(
    fam %in% c("S1", "S2", "S3") & row_number() <= 1000 ~ 0,
    fam %in% c("S1", "S2", "S3") & row_number() > 1000 & row_number() <= 2000 ~ 1,
    TRUE ~ NA_real_  
  )) %>%
  ungroup()

covar_file = cbind(fam_G$fam, fam_G$id, fam_G$sex, fam_G$pheno) %>% as.data.frame() %>% dplyr::rename(famid = V1, iid = V2, sex = V3, pheno = V4)
write.table(covar_file, paste0('30G_covar_6000n_500000m.txt'), col.names = TRUE, row.names = FALSE, quote = FALSE)
# write platform 2 id's (biased platform) - for easier downstream processing
p2_id = covar_file %>% filter(pheno == 1) %>% dplyr::rename(FID = famid, IID = iid) 
write.table(p2_id, paste0('id_p2.txt'), col.names = TRUE, row.names = FALSE, quote = FALSE)

# rewrite fam file so that it's ordered by platform 1 followed by platform 2
fam_platform = fam_G %>% arrange(pheno)
write_fam('30G_6000n_500000m_controls_bias', fam_platform)

# rewrite bim file with sampled percentage of reverse complements
bim_G = plink$bim
reverse_comp_percentage <- 0.1
source('/hpc/group/ochoalab/tt207/platform_bias/generate_bim.R')
bim_new = generate_bim(bim_G, reverse_comp_percentage)
# check number of reverse comp
bim_reverse_comp = bim_new[which(bim_new$ref == comp(bim_new$alt)),]
# write file with reverse comp id's
write.table(bim_reverse_comp$id, "id_reverse_comp.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)
write_bim('30G_6000n_500000m_controls_bias', bim_new)

