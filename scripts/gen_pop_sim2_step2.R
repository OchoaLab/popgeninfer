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
setwd(paste0('/hpc/group/ochoalab/tt207/platform_bias/sim', case, '/rep', rep_num, '/subpop'))

# run loop for each pop, simulate separately
subpop = 3
# 2000 individuals in each subpop
n = 2000
# number of generations
G = 30


for (x in 1:subpop){
  ## pedigree data
  fam = read.table(paste0('pedigree_fam_30G_S', x, ".txt"), header = TRUE)
  load(paste0('pedigree_ids_30G_S', x, ".RData"))

  # draw genotypes X through pedigree
  # read X genotype file from sim1 (first generation) !!
  plink_file = paste0('/hpc/group/ochoalab/tt207/platform_bias/sim1/rep', rep_num, '/subpop/1G_2000n_500000m_S', x)
  plink = read_plink(plink_file)
  X_1 = plink$X
  X_famid = plink$fam$fam
  #rownames(X_1)
  # need to rewrite famid's to indicate subpop
  famid = rep(X_famid, G)
  fam$fam <- famid
  
  # update fam id to reflect subpop
  new_fam = fam %>% separate(id, c("group", "id"), sep = "-") %>% 
    separate(pat, c("pat_group", "pat"), sep = "-") %>% 
    separate(mat, c("mat_group", "mat"), sep = "-") %>% 
    mutate(id = as.numeric(id) + ifelse(x == 1, 0, 2000 * (x - 1)),
           pat = ifelse(is.na(pat) == TRUE, NA, as.numeric(pat) + 2000* (x - 1)),
           mat = ifelse(is.na(mat) == TRUE, NA, as.numeric(mat) + 2000* (x - 1))) %>% 
    mutate(id = paste0(group, "-", id),
           pat = ifelse(is.na(pat) == TRUE, NA, paste0(pat_group, "-", pat)),
           mat = ifelse(is.na(mat) == TRUE, NA, paste0(mat_group, "-", mat))) %>% select(-group, -mat_group, -pat_group)
  
  print(paste('write fam file for subpop', x))
  write.table(new_fam, paste0('pedigree_fam_30G_S', x, "_new.txt"),
              col.names = TRUE, row.names = FALSE, quote = FALSE)
  
  # update id's to reflect subpop
  increment_suffixes <- function(sublist, x) {
    if (x == 1) {
      return(sublist)  # Skip the function if x is 1
    } else {
      increment <- 2000 * (x - 1)
    }
    
    split_elements <- strsplit(sublist, "-")
    prefixes <- sapply(split_elements, `[`, 1)
    suffixes <- sapply(split_elements, `[`, 2)
    new_suffixes <- as.numeric(suffixes) + increment
    new_elements <- paste(prefixes, new_suffixes, sep = "-")
    return(new_elements)
  }
  
  new_ids <- lapply(ids, increment_suffixes, x)
  
  print(paste('write fam ids file for subpop', x))
  save( new_ids, file = paste0('pedigree_ids_30G_S', x, "_new.RData")) 
  
  fam <- prune_fam( new_fam, new_ids[[G]] )
  
  X_1_id = colnames(X_1) # 2001, 2002...
  
  colnames(X_1) = paste0("1-", X_1_id)
  print(X_1[1:5, 1:5])
  # generate fam file for last generation
  fam_G = tail(fam, n)
  # generate last generation
  print("generate last generation X_G")
  
  name_out <- paste0('30G_2000n_500000m_S', x)
  m_loci = 500000
  
  ## NEW version that is a bit more memory efficient
  # a global variable updated as we go
  m_last <- 0
  # simulator function
  sim_chunk <- function( m_chunk ) {
    # these are the locus indexes to process right now:
    indexes <- m_last + ( 1 : m_chunk )
    
    X_G <- geno_last_gen( X_1[ indexes, , drop = FALSE ], fam, new_ids )
    rownames(X_G) <- NULL
    
    # create dummy bim to go with this
    bim <- make_bim( n = m_chunk )
    # this creates IDs (and positions) that don't repeat/clash
    bim$id <- bim$id + m_last
    bim$pos <- bim$id
    #print(bim)
    # update global value (use <<-) for next round
    m_last <<- m_last + m_chunk
    return( list( X = X_G, bim = bim ) )
  }
  # simulate data and write it on the go!
  print('sim_and_plink')
  sim_and_write_plink( sim_chunk, m_loci, fam_G, name_out )
  

  }
# 
# 
# 
