
# popgeninfer

The goal of popgeninfer is to perform simple statistical tests for population genetics data.
The package currently implements tests for allele frequency differences between datasets. 

The `/scripts` folder contains R and bash scripts for implementing Allele Frequency Filter (AF-filter) test and Logistic Mixed Model Filter (LMM-filter) test for harmonizing external controls in multiethnic case-control association studies with platform-specific genotyping biases. Current scripts are written for our simulated data for two scenarios (1) single generation with no family relatedness, and (2) 30 generations with relatedness within subpopulation. Scripts for simulating genotypes with platform biases and allele flips are included as well. R-package version for AF-filter and LMM-filter are in development. 

## Installation

The current development version can be installed from the GitHub repository using `devtools`:
```R
install.packages("devtools") # if needed
library(devtools)
install_github('OchoaLab/popgeninfer', build_vignettes = TRUE)
```


## Example

``` r
library(popgeninfer)

# test that the allele frequency in one dataset (counts of one allele `x1` out of `n1` total alleles) equals that in a second dataset (counts `x2` out of `n2`).
data <- af_test( x1, n1, x2, n2 )
data$pval
```

