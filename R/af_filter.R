#' Allele frequency filter
#'
#' This applies the AF-filter, which is a wrapper around [af_test()] that tests whether allele frequencies match in the given orientation between two datasets, then identifies flippable (reverse complement) loci and also tests them in the reverse orientation in the second dataset.
#' This function returns the input locus info table with the two sets of p-values, as applicable.
#' However, no thresholds or classification results are applied here, for that see another function.
#'
#' @inheritParams af_test
#' @param bim Locus annotation table, containing minimally two columns named `ref` and `alt` with the two alleles at that locus (such as "A" and "G", or "TGT" and "T").
#'
#' @return The input `bim` table with the following three columns added:
#' - `revcomp`: Logical, whether locus is identified as having reverse complement (flippale) alleles.
#' - `pval_fwd`: p-values from [af_test()] in given orientation.
#' - `pval_rev`: p-values from [af_test()] in reversed orientation for the second dataset, for loci with `revcomp == TRUE`.  Rest of loci have `NA` p-values.
#'
#' @examples
#' # data with three loci, two ancestries, two datasets
#' # make our minimal BIM table of locus annotations
#' # note last three are reverse complement ref/alt pairs
#' bim <- data.frame(
#'   ref = c('A', 'G', 'C', 'T', 'G'),
#'   alt = c('G', 'T', 'G', 'A', 'C')
#' )
#' m_loci <- nrow( bim )
#' k_subpops <- 2
#' # all same sample size with no missingness for simplicity
#' n <- 10
#' # and same allele frequency in this case
#' p <- 0.5
#' # the actual data we need as input
#' n1 <- matrix( n, nrow = m_loci, ncol = k_subpops )
#' n2 <- n1
#' x1 <- matrix(
#'   rbinom( m_loci * k_subpops, n, p ),
#'   nrow = m_loci, ncol = k_subpops
#' )
#' x2 <- matrix(
#'   rbinom( m_loci * k_subpops, n, p ),
#'   nrow = m_loci, ncol = k_subpops
#' )
#'
#' # apply function, get desired p-values!
#' data <- af_filter( x1, n1, x2, n2, bim )
#'
#' @seealso
#' [af_test()] for the association test based on allele frequencies/counts.
#' 
#' [revcomp()] for determining reverse complements.
#' 
#' @export
af_filter <- function( x1, n1, x2, n2, bim ) {
    # identify reverse complement cases, which are "flippable"
    bim$revcomp <- bim$ref == revcomp( bim$alt )

    # forward test
    bim$pval_fwd <- af_test( x1, n1, x2, n2 )$pval
    
    # perform reversed test, reverse second platform only
    x1c <- x1[ bim$revcomp, ]
    n1c <- n1[ bim$revcomp, ]
    x2c <- x2[ bim$revcomp, ]
    n2c <- n2[ bim$revcomp, ]
    bim$pval_rev <- NA
    bim$pval_rev[ bim$revcomp ] <- af_test( x1c, n1c, n2c - x2c, n2c )$pval

    # this is the BIM table with the addition of columns: `revcomp` (boolean), `pval_fwd`, and `pval_rev`
    return( bim )
}
