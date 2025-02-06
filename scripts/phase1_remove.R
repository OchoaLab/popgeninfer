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
print("input values: iteration, pvalue threshold, and rep")
print('iter:')
print(iter)
print(pval)
print('rep:')
print(rep)
print('sim:')
print(case)
setwd(paste0('/hpc/group/ochoalab/tt207/platform_bias/sim', case, '/rep', rep, '/saige/', pval))
pval = as.numeric(pval)
###########################################
# identify significant SNPs in current iteration

# identify SNPs that can be flipped
if (iter == 0) {
  sig_snps = read.table(paste0("../platform_", iter, "_output.txt"), header = TRUE) %>% filter(p.value < pval) %>% pull(MarkerID)
  print("remove SNPs:")
  print(length(sig_snps))
  # master file of SNPs identified as flip or remove
  write.table(sig_snps, "platform_snps_phase1.txt", row.names = FALSE, quote = FALSE)
 
} else {
  sig_snps = read.table(paste0("platform_", iter, "_output.txt"), header = TRUE) %>% 
    filter(p.value < pval) %>% filter(p.value != 0) %>% pull(MarkerID)
  
  if (length(sig_snps) == 0) {
    message("saige p-values are less than the defined p-value threshold. Stopping the script.")
    stop("p-values are less than the defined p-value threshold.")
  } else {
    # load sig snps from original control vs control saige output:
    snps_master = read.table("platform_snps_phase1.txt", header = TRUE) %>% select(1) %>% pull()
    print("snps_master")
    print(length(snps_master))
    print("sig_snp")
    print(length(sig_snps))
    # write master file for removal
    write.table(unique(c(snps_master, sig_snps)), "platform_snps_phase1.txt", row.names = FALSE, quote = FALSE)
    write.table(unique(sig_snps), paste0("platform_snps_phase1_", iter,".txt"), row.names = FALSE, quote = FALSE)
    
    print("remove SNPs accumulated count:")
    print(length(c(snps_master, sig_snps))) 
  }
  
}




##############################
#plot 
##############################
# GR_input = saige_output %>% mutate(chr = paste0("chr", CHR), start = POS, end = POS) %>% 
#   select(chr, start, end, A1 = Allele1, A2 = Allele2, P = p.value)
# platform_gr = makeGRangesFromDataFrame(data.frame(GR_input), keep.extra.columns = TRUE)
# 
# y_max = round(-log10(min(saige_output$p.value))) # 295
# 
# png( paste0('./flip/saige_platform_', iter, '.png'), width=18.5, height=7, res=1000, units = 'in')
# #par(mfrow=c(5,3), mar=c(2,2,2,2))
# chr_plot = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", 
#              "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22")
# 
# pp <- getDefaultPlotParams(plot.type=4)
# pp$ideogramlateralmargin = 0 
# pp$ideogramheight = 0
# pp$data1inmargin = -7
# pp$leftmargin <- 0.1
# plot.new()
# kp <- plotKaryotype("hg38", plot.type=4, ideogram.plotter = NULL, labels.plotter = NULL,
#                     chromosomes = chr_plot, main = "array vs WGS", cex=1.5, plot.params = pp)
# # bottom plot
# kpPlotManhattan(kp, data = platform_gr, points.col="2blues", r0=autotrack(1,1, margin = 0.15),genomewideline = -log10(5e-8),
#                 suggestiveline = -log10(1e-5), , suggestive.col = "blue",
#                 genomewide.col = "red", ymax=y_max)
# kpAddLabels(kp, labels = " ", srt=90, pos=3, r0=autotrack(1,1, margin = 0.15), cex=2, label.margin = 0.055)
# kpAddLabels(kp, labels = "-log10( p-value )", srt=90, pos=3, r0=autotrack(1,1, margin = 0.15), cex=1.5, label.margin = 0.039)
# kpAxis(kp, ymin=0, ymax=y_max, r0=autotrack(1,1 , margin = 0.15))
# 
# kpAddChromosomeNames(kp, srt=45, cex = 1)
# invisible( dev.off() )
