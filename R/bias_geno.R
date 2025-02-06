#' Simulate genotyping platform bias
#'
#' This function accepts a true genotype matrix as the input and either the proportion of biased loci or their exact indexes, and simulates biased genotypes determined by those and some additional parameters.
#' The goal is to simulate genotyping platform-specific errors where there is a greater tendency to read a certain allele over the other.
#' In particular, once biased loci are determined, each of those loci is drawn from a Beta distribution a random probability that individuals at that locus will have an erroneous genotype (`error_prob`), and a value for the erroneous genotype (`error_geno`, either 0 or 2, picked with equal probability).
#' Then, at the same biased locus, a proportion `error_prob` of random individuals have their genotypes replaced by the `error_geno` value, which therefore results in biased genotypes towards the value of `error_geno`.
#'
#' @param X True genotype matrix, with values in 0, 1, or 2 counting reference alleles.
#' This matrix must have loci along the rows and individuals along columns.
#' @param p Proportion of biased loci.
#' Either this or `biased_loci` must be provided.
#' If provided and `biased_loci` is `NULL`, then biased loci will be picked randomly from among all loci.
#' @param biased_loci The indexes of loci desired to be biased.
#' Either this or `p` must be provided.
#' If provided, this is the set of loci that will be biased, and `p` is ignored.
#' @param shape1 The first shape parameter of the Beta distribution used to generate error probabilities per locus.
#' @param shape2 The second shape parameter of the Beta distribution used to generate error probabilities per locus.
#'
#' @return A named list with the following elements
#' - `X`: The modified genotype matrix, with the same dimensions as the input `X`.  A proportion `p` of its loci will have platform biases, the rest of the loci are unaltered.
#' - `biased_loci`: The length `m_biased` vector of indexes of biased loci.
#' - `error_probs`: The length `m_biased` vector of randomly drawn error probabilities from the Beta distribution, one per locus.
#' - `error_genos`: The length `m_biased` vector of randomly drawn erroneous genotypes, one per locus.
#'
#' @examples
#' # simulate a toy genotype matrix with 10 loci and 3 individuals
#' X <- matrix(
#'     rbinom( 30, 2, 0.5 ),
#'     nrow = 10,
#'     ncol = 3
#' )
#'
#' # simulate platform bias!
#' obj <- bias_geno( X, p = 0.1 )
#' 
#' # the genotype matrix with platform-specific biases is here
#' X_biased <- obj$X
#' 
#' # the other elements in the list `obj` provide further details
#' # regarding how the platform biases were generated
#' 
#' @export
bias_geno <- function(
                      X,
                      p = NA,
                      biased_loci = NULL,
                      shape1 = 1,
                      shape2 = 5
                      ) {
    # check inputs
    if ( missing( X ) )
        stop( '`X` is required!' )
    if ( is.null( biased_loci ) ) {
        if ( is.na( p ) )
            stop( 'Either `p` or `biased_loci` must be provided!' )
        if ( p < 0 )
            stop( '`p` must be non-negative' )
        if ( p > 1 )
            stop( '`p` must be <= 1' )
        # not sure if this ever makes sense, but no biased rows means we return the input matrix unchanged
        if ( p == 0 )
            return ( X )
        
        # else p is defined and in range, in which case the following succeeds
        m_loci <- nrow( X )
        m_biased <- round( p * m_loci )
        biased_loci <- sample.int( m_loci, m_biased )
    } else {
        m_biased <- length( biased_loci )
    }

    # draw probabilities per biased locus of an individual having the biased SNP (an "error")
    error_probs <- stats::rbeta( m_biased, shape1 = shape1, shape2 = shape2 )
    # draw once the biased genotypes, the same value is used for all biased cases per locus
    error_genos <- 2L * stats::rbinom( m_biased, 1L, 0.5 )
    
    # apply bias in a loop per locus
    n_ind <- ncol( X )
    for ( i in 1 : m_biased ) {
        # select individuals randomly to bias at this locus
        selected_inds <- sample.int( n_ind, n_ind * error_probs[i] )
        # replace the original genotypes of those biased individuals at this locus with the shared biased value
        X[ biased_loci[ i ], selected_inds ] <- error_genos[i]
    }

    # done!
    return( list(
        X = X,
        biased_loci = biased_loci,
        error_probs = error_probs,
        error_genos = error_genos
    ) )
}

