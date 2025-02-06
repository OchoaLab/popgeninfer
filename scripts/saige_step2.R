library(SAIGE)
library(optparse) 
# terminal inputs
option_list = list(
  make_option(c( "-i", "--iteration"), type = "character", default = '1', 
              help = "numeric number for iteration of saige control vs control", metavar = "character"),
  make_option(c( "-p", "--path"), type = "character",
              help = "current file path", metavar = "character"),
  make_option(c( "-s", "--sim"), type = "character",
              help = "simulation case", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
# get values
iter <- opt$i
path <- opt$p
case <- opt$s
setwd(path)
#dir = "/datacommons/ochoalab/ssns_gwas/array/saige_platform/flip/"
#### only controls
if (case == 1) {
  if (iter == 0) {
    plinkFile <- '../../1G_6000n_500000m_controls' # rep1 folder
    outputPrefix= paste0('../platform_', iter) # saige folder
    GMMATmodelFile = paste0('../platform_', iter, '.rda') 
    varianceRatioFile = paste0('../platform_', iter, '.varianceRatio.txt') 
    SAIGEOutputFile = paste0('../platform_', iter, '_output.txt') 
  } else {
    # in pval folder
    plinkFile <- paste0('1G_6000n_500000m_controls_', iter)
    GMMATmodelFile = paste0('platform_', iter, '.rda') 
    varianceRatioFile = paste0('platform_', iter, '.varianceRatio.txt') 
    SAIGEOutputFile = paste0('platform_', iter, '_output.txt') 
  }
} else {
  if (iter == 0) {
    plinkFile <- '../../30G_6000n_500000m_controls_bias' # rep1 folder
    outputPrefix= paste0('../platform_', iter) # saige folder
    GMMATmodelFile = paste0('../platform_', iter, '.rda') 
    varianceRatioFile = paste0('../platform_', iter, '.varianceRatio.txt') 
    SAIGEOutputFile = paste0('../platform_', iter, '_output.txt') 
  } else {
    # in pval folder
    plinkFile <- paste0('30G_6000n_500000m_controls_', iter)
    GMMATmodelFile = paste0('platform_', iter, '.rda') 
    varianceRatioFile = paste0('platform_', iter, '.varianceRatio.txt') 
    SAIGEOutputFile = paste0('platform_', iter, '_output.txt') 
  }
  
}

print( 'saige step 2')
print(plinkFile)
print(GMMATmodelFile)
print(varianceRatioFile)
print(SAIGEOutputFile)

SPAGMMATtest(bedFile=paste0(plinkFile, ".bed"),
             bimFile=paste0(plinkFile, ".bim"),
             famFile=paste0(plinkFile, ".fam"),
             AlleleOrder= 'alt-first',
             is_imputed_data=TRUE,
             #impute_method = opt$impute_method,
             GMMATmodelFile=GMMATmodelFile,
             varianceRatioFile=varianceRatioFile,
             SAIGEOutputFile=SAIGEOutputFile,
             is_output_moreDetails =TRUE,
             is_overwrite_output = TRUE,
             #SPAcutoff = opt$SPAcutoff, default 2
             is_Firth_beta = TRUE, # for binary traits
             #pCutoffforFirth = opt$pCutoffforFirth, # default 0.01
             LOCO = FALSE,
             min_MAF=0,
             min_MAC=0.5,
             max_missing = 1,
             dosage_zerod_cutoff = 0,
             dosage_zerod_MAC_cutoff = 0
)