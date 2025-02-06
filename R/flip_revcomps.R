#' Simulate strand flipping for ambiguous ref/alt pairs
#'
#' This function transforms a genotype matrix by "flipping" some of its loci so that the opposite allele is counted.
#' In other words, its numeric value `x`, which should be in 0, 1, or 2, is replaced with `2-x` but only at some proportion `p` of loci that are flippable.
#' Flippable loci are those where the alt allele is the reverse complement of the ref allele, so that there is strand ambiguity and alleles cannot be aligned trivially to external datasets.
#'
#' @inheritParams bias_geno
#' @param p Proportion of flippable loci to flip.
#' Either this and `bim` or `flipped_loci` must be provided.
#' If provided and `flipped_loci` is `NULL`, then flipped loci will be picked randomly from among all flippable loci.
#' @param bim A data frame or tibble with `ref` and `alt` allele columns for every locus in `X`, used only to determine which loci are flippable.
#' Either this and `p` or `flipped_loci` must be provided.
#' @param flipped_loci The indexes of loci desired to be flipped.
#' Either this or `p` and `bim` must be provided.
#' If provided, this is the set of loci that will be flipped, and `p` and `bim` are ignored.
#' In this case `flipped_loci` is not checked against the list of flippable loci.
#'
#' @return A named list with the following elements
#' - `X`: The modified genotype matrix, with the same dimensions as the input `X`.  A proportion `p` of its flippable loci will be flipped, the rest of the loci are unaltered.
#' - `flipped_loci`: The vector of indexes of flipped loci.
#'
#' @examples
#' # simulate a toy genotype matrix with 5 loci and 3 individuals
#' X <- matrix(
#'     rbinom( 15, 2, 0.5 ),
#'     nrow = 5,
#'     ncol = 3
#' )
#'
#' # for our fake BIM table, only the first three are flippable (reverse complement ref/alt)
#' library(tibble)
#' bim <- tibble(
#'     alt = c('A', 'G', 'C', 'T', 'GATTACA'),
#'     ref = c('T', 'C', 'G', 'G', 'G')
#' ) 
#'
#' # randomly pick two of the three flippable loci to flip!
#' obj <- flip_revcomps( X, p = 2/3, bim )
#' 
#' # this is the modified genotype matrix
#' obj$X
#'
#' # and the indexes of the flipped loci
#' obj$flipped_loci
#'
#' @seealso
#' [revcomp()] for reverse complements.
#' 
#' @export
flip_revcomps <- function(
                          X,
                          p = NA,
                          bim = NULL,
                          flipped_loci = NULL
                          ) {
    # check inputs
    if ( missing( X ) )
        stop( '`X` is required!' )
    if ( is.null( flipped_loci ) ) {
        if ( is.na( p ) || is.null( bim ) )
            stop( 'Either `flipped_loci` or both `p` and `bim` must be provided!' )
        if ( nrow( bim ) != nrow( X ) )
            stop( 'The number of rows of `X` and `bim` must match!' )
        
        # if p and bim are provided, lets determine which loci to flip
        
        # first find loci that are flippable.  If they're not reverse complements, they can be aligned trivially, so there wouldn't be flip ambiguity
        flippable_loci <- which( bim$ref == revcomp( bim$alt ) )

        # now select loci to actually flip from among flippable loci
        flipped_loci <- sample( flippable_loci, round( p * length( flippable_loci ) ) )
    }

    # apply all flips together now
    X[ flipped_loci, ] <- 2 - X[ flipped_loci, ]

    # return data
    return( list(
        X = X,
        flipped_loci = flipped_loci
    ) )
}
