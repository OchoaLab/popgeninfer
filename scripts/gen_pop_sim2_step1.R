library(simfam)
library(simfam)
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

setwd(paste0('/hpc/group/ochoalab/tt207/platform_bias/sim', case))
new_dir <- paste0('rep', rep_num)
dir.create(new_dir)
setwd(new_dir)
dir.create('subpop')
setwd('subpop')

# sim 2 : 3 subpopulations, admixture within subpopulations
# simulate pedigree for each subpop separately

subpop = 3
n = 2000
# number of generations
G = 30

data <- sim_pedigree( n, G )
ids = data$ids

#head(data$ids, n = 5)
#head(ids)


for (x in 1:subpop) {
  data <- sim_pedigree( n, G )
  ids = data$ids
  print(paste('write fam file for subpop', x))
  write.table(data$fam, paste0("pedigree_fam_30G_S", x, ".txt"),
              col.names = TRUE, row.names = FALSE, quote = FALSE)
  
  print(paste('write fam ids file for subpop', x))
  save( ids, file = paste0("pedigree_ids_30G_S", x, ".RData")) 
}

