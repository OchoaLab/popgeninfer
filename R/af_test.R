#' Tests if allele frequencies match between datasets
#'
#' This function tests if allele frequencies at each locus match between two datasets based on their counts and sample sizes.
#' Under the null hypothesis, the frequency is the same in both datasets, while under the alternative they are different.
#' Counts are assumed to be Binomially distributed.
#' Significance is calculated using the likelihood ratio test.
#' If inputs are matrices of equal dimensions, rows are tested separately but columns are combined, assigned a single p-value calculated with degrees of freedom equal to the number of columns, which tests the joint null hypothesis (useful for combining data from different ancestries, or different loci).
#' Vector inputs are treated as column matrices.
#'
#' @param x1 The number of alleles of one type in the first dataset.
#' @param n1 The total numbers of alleles (of all types) in the first dataset.
#' @param x2 The number of alleles of one type in the second dataset.
#' @param n2 The total numbers of alleles (of all types) in the second dataset.
#'
#' @return A list with these values:
#' - `stat`: the test statistics (length equal to number of rows of inputs).
#' - `df`: degrees of freedom used for chi squared test.
#' - `pval`: the p-values (length equal to number of rows of inputs).
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
#' data <- af_test( x1, n1, x2, n2 )
#' # the p-value is most desired quantity.
#' # though null hypothesis is true, there is low power because sample sizes are small,
#' # so example is most likely insignificant
#' data$pval
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
#' # perform desired statistical test, a single p-value testing whether all of the three ancestries
#' # have the same allele frequency in both datasets (within each ancestry) or not, which is a
#' # 3-degree-of-freedom test.
#' data <- af_test( x1, n1, x2, n2 )
#' data$pval
#' 
#' @export
af_test <- function( x1, n1, x2, n2 ) {
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
    m <- mk[1]
    k <- mk[2]

    # let vectorization work its magic!
    # works with input matrices but output is flattened to vector, so return to matrix
    #stat <- sum( af_test_single( x1, n1, x2, n2 ) )
    stat <- matrix(
        af_test_single( x1, n1, x2, n2 ),
        nrow = m,
        ncol = k
    )
    # sum across ancestries but not loci
    stat <- rowSums( stat )
    
    # under the null, distribution is asymptotically chi-squared distributed
    # it is the same number of degrees of freedom in all rows
    pval <- stats::pchisq( q = stat, df = k, lower.tail = FALSE )
    return(
        list(
            stat = stat,
            df = k,
            pval = pval
        )
    )
}

# original version from Tiffany for k=3 only
# kept for internal tests
likelihoodratio_3 <- function(x1, n1, x2, n2, x3, n3, x4, n4, x5, n5, x6, n6){
    ## alt pmf
    p1 = x1/n1
    p2 = x2/n2
    pmf1 = stats::dbinom(x1, n1, p1, log = TRUE)
    pmf2 = stats::dbinom(x2, n2, p2, log = TRUE)
    
    p3 = x3/n3
    p4 = x4/n4
    pmf3 = stats::dbinom(x3, n3, p3, log = TRUE)
    pmf4 = stats::dbinom(x4, n4, p4, log = TRUE)
    
    p5 = x5/n5
    p6 = x6/n6
    pmf5 = stats::dbinom(x5, n5, p5, log = TRUE)
    pmf6 = stats::dbinom(x6, n6, p6, log = TRUE)
    
    altpmf_b = pmf1 + pmf2 
    altpmf_a = pmf3 + pmf4 
    altpmf_w = pmf5 + pmf6
    
    ## null pmf
    x_b = x1+x2
    x_a = x3+x4
    x_w = x5+x6
    n_b = n1+n2
    n_a = n3+n4
    n_w = n5+n6
    p_b = x_b/n_b
    p_a = x_a/n_a
    p_w = x_w/n_w
    
    nullpmf_b = stats::dbinom(x1, n1, p_b, log = TRUE) + stats::dbinom(x2, n2, p_b, log = TRUE)
    nullpmf_a = stats::dbinom(x3, n3, p_a, log = TRUE) + stats::dbinom(x4, n4, p_a, log = TRUE)
    nullpmf_w = stats::dbinom(x5, n5, p_w, log = TRUE) + stats::dbinom(x6, n6, p_w, log = TRUE)
    
    #likelihood ratio
    stat = -2*((nullpmf_b - altpmf_b) + (nullpmf_a - altpmf_a) + (nullpmf_w - altpmf_w)) 
    # converge to chi-square distribution
    pval = stats::pchisq(q = stat, df = 3, lower.tail = FALSE)
    return(
        list(
            stat = stat,
            df = 3,
            pval = pval
        )
    )
}
