# popgeninfer 0.0.0.9000 (2023-05-08)

- Added function `af_test`, which implements a clean and well-tested version of the allele frequency test first implemented by Tiffany Tu for the NS GWAS project.  This tests if the data from two datasets have the same allele frequency (null hypothesis) or not (alternative hypothesis).  Supports combining data from different ancestries into a single p-value.

# popgeninfer 0.0.1.9000 (2023-05-09)

- Function `af_test` changed default behavior for vector and matrix inputs:
  - Previously, vector inputs were treated as data to combine, and matrices were flattened to vectors
  - Now, vectors are treated as column matrices, and for matrices in general rows are treated as separate tests and columns as data to combine
  - The new setup aligns better with the actual uses of this function in various lab projects.
- Added package vignette, confirming with simple simulations that `af_test` p-values are uniform under the null and concentrated near zero under the alternative hypothesis.

# popgeninfer 0.0.2.9000 (2024-05-15)

- Added function `af_test_or`, which complements `af_test` to calculate odds ratios, together with confidence intervals, and p-values.  However, this test has a different alternative model (more constrained because there is a single OR shared across ancestries) so the p-value are different when there is more than one column/ancestry.

# popgeninfer 0.0.3.9000 (2024-05-16)

- Function `af_test_or` corrected a minor bug that swapped the roles of the allele counts (which are response) and study (which is covariate).  Also internally optimized how this Binomial data is passed to function `glm`, and made minor clarifications to documentation
- Updated vignette to include `af_test_or`, and extended simulations to include new cases relevant to that function.  The vignette validates the accuracy of the new function and contrasts it to `af_test`.
