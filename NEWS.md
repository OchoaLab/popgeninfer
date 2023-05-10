# popgeninfer 0.0.0.9000 (2023-05-08)

- Added function `af_test`, which implements a clean and well-tested version of the allele frequency test first implemented by Tiffany Tu for the NS GWAS project.  This tests if the data from two datasets have the same allele frequency (null hypothesis) or not (alternative hypothesis).  Supports combining data from different ancestries into a single p-value.

# popgeninfer 0.0.1.9000 (2023-05-09)

- Function `af_test` changed default behavior for vector and matrix inputs:
  - Previously, vector inputs were treated as data to combine, and matrices were flattened to vectors
  - Now, vectors are treated as column matrices, and for matrices in general rows are treated as separate tests and columns as data to combine
  - The new setup aligns better with the actual uses of this function in various lab projects.
- Added package vignette, confirming with simple simulations that `af_test` p-values are uniform under the null and concentrated near zero under the alternative hypothesis.
