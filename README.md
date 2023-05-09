
# popgeninfer

The goal of popgeninfer is to perform simple statistical tests for population genetics data.
The package currently implements tests for allele frequency differences between datasets.

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

