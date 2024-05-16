#' Estimate odds ratios (ORs) between two conditions
#'
#' This function tests if the ORs between allele frequencies at each locus is non-zero between two datasets based on their counts and sample sizes, returning ORs (odds of first condition over odds of second), confidence intervals, and p-values.
#' Under the null hypothesis, the OR is 1 in both datasets, while under the alternative it is non-zero.
#' Counts are assumed to be Binomially distributed.
#' If inputs are matrices of equal dimensions, rows are tested separately but columns are combined, assigned a single OR assuming that it is shared across columns (useful for combining data from different ancestries, or different loci).
#' Vector inputs are treated as column matrices.
#' Parameter fitting and significance is calculated using [stats::glm()], with dataset as the effect variable and column (when more than one) as a categorical covariate conditioned upon.
#'
#' This function accepts exactly the same inputs as [af_test()], which it complements by providing odds ratios and 95% confidence intervals, which the former does not return at all.
#' However, a shared OR across columns of the same row/locus requires a more restricted alternative hypothesis with one alternative parameter, the OR, and thus a 1 degree of freedom test, so these p-values are different from those of [af_test()].
#' (In contrast, [af_test()] implicitly allows ORs to be different for different columns, and they can even go in opposite directions, so the test has degrees of freedom equal to the number of columns.)
#' 
#' @inheritParams af_test
#' @param level The level for the confidence interval.
#'
#' @return A tibble with rows matching inputs and with these columns:
#' - `OR`: the estimated odds ratio.
#' - `CIL`: lower confidence interval.
#' - `CIU`: upper confidence interval.
#' - `pval`: the p-value.
#'
#' @examples
#' ### Test single locus ###
#' 
#' # Simulate data where the alternative hypothesis is true
#' # true allele frequencies are different
#' p1 <- 0.2
#' p2 <- 0.3
#' # sample sizes can be very different
#' n1 <- 100
#' n2 <- 50
#' # simulate counts from Binomial distribution
#' x1 <- rbinom( 1, n1, p1 )
#' x2 <- rbinom( 1, n2, p2 )
#'
#' # this performs desired statistical test
#' data <- af_test_or( x1, n1, x2, n2 )
#' # inspect stats!
#' # though null hypothesis is true, there is low power because sample sizes are small,
#' # so example is most likely insignificant
#' data
#'
#' ### Multiple loci and ancestries ###
#'
#' # it's useful to think of each column as the same locus but a different ancestry,
#' # so the allele frequencies can be different.
#' # Each row is a different locus to be tested separately.
#' # here the null hypothesis is true so both datasets will share the same allele frequency.
#' # will simulate 2 loci, 3 ancestries
#' p <- matrix(
#'     c(
#'         0.2, 0.7, 0.5,
#'         0.1, 0.3, 0.5
#'     ),
#'     nrow = 2,
#'     ncol = 3,
#'     byrow = TRUE
#' )
#' # all sample sizes can be different
#' n1 <- matrix( c( 100, 80, 200, 33, 66, 99 ), nrow = 2, byrow = TRUE )
#' n2 <- matrix( c( 50, 15, 99, 44, 55, 66 ), nrow = 2, byrow = TRUE )
#' # simulate counts from Binomial distribution (vectorized across all 3 ancestries,
#' # reshaped as matrix of desired dimensions)
#' x1 <- matrix( rbinom( 6, n1, p ), nrow = 2 )
#' x2 <- matrix( rbinom( 6, n2, p ), nrow = 2 )
#'
#' # perform desired statistical test, a single p-value testing whether the shared OR across 
#' # ancestries differs from 1.
#' data <- af_test_or( x1, n1, x2, n2 )
#' data
#'
#' @seealso
#' [af_test()], which unlike this function calculates p-values (but not ORs) under a more general alternative hypothesis that the OR is non-zero with different values in each column.
#' 
#' @export
af_test_or <- function( x1, n1, x2, n2, level = 0.95 ) {
    # coerce vector inputs to column matrices
    if ( !is.matrix( x1 ) )
        x1 <- as.matrix( x1 )
    if ( !is.matrix( n1 ) )
        n1 <- as.matrix( n1 )
    if ( !is.matrix( x2 ) )
        x2 <- as.matrix( x2 )
    if ( !is.matrix( n2 ) )
        n2 <- as.matrix( n2 )
    
    # confirm that all matrices have same dimensions
    mk <- dim( x1 )
    if ( any( dim( n1 ) != mk ) )
        stop( 'Dimensions of `x1` (', toString( mk ), ') and `n1` (', toString( dim( n1 ) ), ') differ!' )
    if ( any( dim( x2 ) != mk ) )
        stop( 'Dimensions of `x1` (', toString( mk ), ') and `x2` (', toString( dim( x2 ) ), ') differ!' )
    if ( any( dim( n2 ) != mk ) )
        stop( 'Dimensions of `x1` (', toString( mk, ), ') and `n2` (', toString( dim( n2 ) ), ') differ!' )

    # separate dimensions now
    m <- mk[1] # snps
    
    # test each locus separately
    # output is, for now, a matrix with a row for every locus, and 4 values that are the desired stats
    data <- matrix( NA, m, 4 )
    for ( i in 1 : m ) {
        # calculate and store values in i'th row
        data[ i, ] <- af_test_or_single( x1[ i, ], n1[ i, ], x2[ i, ], n2[ i, ] )
    }

    # now reorganize, order defined in af_test_or_single, into tibble
    # also, better to report transformed data (originals are log-OR)
    data <- tibble::tibble(
        OR = exp( data[ , 1 ] ),
        CIL = exp( data[ , 2 ] ),
        CIU = exp( data[ , 3 ] ),
        pval = data[ , 4 ]
    )
    
    # done, return!
    return( data )
}
