library(tidyverse)

# input: reverse_comp_percentage, bim (must have alt and ref columns), 
# output: bim_new (updated bim file), bim_reverse_comp$id (id's of reverse complement in the bim file)

# Reverse complement function
comp <- function(x) {
  chartr("ATGC", "TACG", x)
}

# Define alleles and generate all possible pairs once
alleles <- c("A", "T", "C", "G")
all_pairs <- expand.grid(ref = alleles, alt = alleles)

# Remove reverse complement pairs
non_reverse_comp_pairs <- all_pairs[!(all_pairs$ref == all_pairs$alt | 
                                        all_pairs$ref == comp(all_pairs$alt)), ]

# Identify reverse complement pairs
reverse_comp_pairs <- all_pairs[all_pairs$ref == comp(all_pairs$alt), ]

# Function to generate new bim with reverse complement
generate_bim <- function(bim, reverse_comp_percentage) {
  
  # Number of loci (SNPs) in the bim file
  m_loci <- nrow(bim)
  
  # Sample reverse complement and non-reverse complement pairs
  sample_reverse_comp <- reverse_comp_pairs %>%
    sample_n(round(m_loci * reverse_comp_percentage), replace = TRUE)
  sample_non_reverse_comp <- non_reverse_comp_pairs %>%
    sample_n(m_loci - round(m_loci * reverse_comp_percentage), replace = TRUE)
  
  # Combine reverse comp and non-reverse comp, shuffle rows
  ref_alt <- rbind(sample_reverse_comp, sample_non_reverse_comp) %>%
    sample_frac(1)
  
  # Rewrite the bim file with new ref and alt columns
  bim_new <- bim %>% select(-alt, -ref)
  bim_new$alt <- ref_alt$alt
  bim_new$ref <- ref_alt$ref
  
  return(bim_new)
}
