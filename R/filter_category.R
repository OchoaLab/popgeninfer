#' Calculate locus categories depending on threshold
#'
#' This takes the output of [af_filter()] or an equivalent function and classifies loci depending on the pattern of significance, which is threshold-dependent.
#' For loci that only have forward p-values, the are classified as "remove" if they are significant (below the threshold) or "keep" otherwise.
#' For loci with both forward and reversed p-values (which must have reversed complement ref/alt alleles), they are classified as "remove" if both p-values are significant, "flip" if the forward p-value is the smallest, and "keep" otherwise.
#'
#' @param data A tibble/data.frame with at least three columns:
#' - `revcomp`: Logical, whether locus is identified as having reverse complement (flippale) alleles.
#' - `pval_fwd`: p-values in given orientation.
#' - `pval_rev`: p-values in reversed orientation for the second dataset, for loci with `revcomp == TRUE`.  Values for rest of loci are ignored.
#' @param pcut The p-value threshold
#'
#' @return A vector with length equal to the number of rows of `data`, with one of 3 possible values depending on classification outcome per locus: "remove", "flip", or "keep".
#'
#' @examples
#' library(tibble)
#' # an example with all cases
#' data <- tibble(
#'     revcomp  = c( FALSE, FALSE,  TRUE,  TRUE,  TRUE,  TRUE,  TRUE ),
#'     pval_fwd = c(   0.9, 1e-10,   0.1,   0.9,   0.9, 1e-10, 1e-10 ),
#'     pval_rev = c(    NA,    NA,   0.9,   0.1, 1e-10,   0.9, 1e-10 )
#' )
#' pcut <- 0.01
#'
#' # get categories!
#' cat <- filter_category( data, pcut )
#'
#' @seealso
#' [af_filter()]
#' 
#' @export
filter_category <- function( data, pcut ) {
    # apply thresholds, categorize loci
    # copy down values
    pval_rest <- data$pval_fwd[ !data$revcomp ]
    pval_fwd <- data$pval_fwd[ data$revcomp ]
    pval_rev <- data$pval_rev[ data$revcomp ]

    # make output
    x <- vector( 'character', nrow( data ) )
    
    x[ data$revcomp ] <- ifelse( pval_fwd < pcut & pval_rev < pcut, "remove",
                         ifelse( pval_fwd < pval_rev, "flip", "keep" ) )
    x[ !data$revcomp ] <- ifelse( pval_rest < pcut , "remove", "keep")

    # return category vector
    return( x )
}
