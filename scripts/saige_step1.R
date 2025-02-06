library(SAIGE)
library(optparse) 
library(genio)

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


#dir = "/datacommons/ochoalab/ssns_gwas/array/saige_platform/flip/"
setwd(path)
#### only controls
#phenoFile <- '../../1G_covar_6000n_500000m_controls.txt'

if (case == 1) {
  phenoFile <- '../../1G_covar_6000n_500000m_controls.txt'
  if (iter == 0) {
    plinkFile <- '../../1G_6000n_500000m_controls' # rep1 folder
    outputPrefix= paste0('../platform_', iter) # saige folder
  } else {
    # in pval folder
    plinkFile <- paste0('1G_6000n_500000m_controls_', iter)
    outputPrefix= paste0('platform_', iter) 
  }
} else {
  phenoFile <- '../../30G_covar_6000n_500000m.txt'
  if (iter == 0) {
    plinkFile <- '../../30G_6000n_500000m_controls_bias' # rep1 folder
    outputPrefix= paste0('../platform_', iter) # saige folder
  } else {
    # in pval folder
    plinkFile <- paste0('30G_6000n_500000m_controls_', iter)
    outputPrefix= paste0('platform_', iter) 
  }
}





print('start saige step1')
print(plinkFile)
print(phenoFile)
print(outputPrefix)

phenoCol= 'pheno'
sampleIDColinphenoFile='iid' 
traitType='binary'        
IsOverwriteVarianceRatioFile=TRUE

fitNULLGLMM(plinkFile = plinkFile,
            phenoFile = phenoFile,
            phenoCol = phenoCol,
            sampleIDColinphenoFile = sampleIDColinphenoFile,
            traitType = traitType,
            outputPrefix = outputPrefix,
            #covarColList = covars,
            #qCovarCol = qcovars,
            IsOverwriteVarianceRatioFile = TRUE,
            LOCO = FALSE,
            #sexCol = "sex",
            minMAFforGRM = 0,
            maxMissingRateforGRM = 1
)
