# returns statistic only for max flexibility
af_test_single <- function( x1, n1, x2, n2 ){
    # confirm that vectors have same length
    k <- length( x1 )
    if ( length( n1 ) != k )
        stop( 'Length of `x1` (', k, ') and `n1` (', length( n1 ), ') differ!' )
    if ( length( x2 ) != k )
        stop( 'Length of `x1` (', k, ') and `x2` (', length( x2 ), ') differ!' )
    if ( length( n2 ) != k )
        stop( 'Length of `x1` (', k, ') and `n2` (', length( n2 ), ') differ!' )
    
    # under the alternative hypothesis, allele frequencies are different
    # these are the respective MLEs
    p1 <- x1 / n1
    p2 <- x2 / n2
    # total log probability mass functions (PMF), treating datasets as independent
    pmf1 <- stats::dbinom( x1, n1, p1, log = TRUE ) + stats::dbinom( x2, n2, p2, log = TRUE )
  
    # under the null hypothesis, the allele frequency is the same in both cases
    # this is the MLE in this case
    p0 <- ( x1 + x2 ) / ( n1 + n2 )
    # and the total log PMF
    pmf0 <- stats::dbinom( x1, n1, p0, log = TRUE ) + stats::dbinom( x2, n2, p0, log = TRUE )
    
    # test statistic
    stat <- 2 * ( pmf1 - pmf0 )
    return( stat )
    ## # converge to chi-square distribution
    ## pvalue <- pchisq( q = stat, df = 1, lower.tail = FALSE )
    ## return( pvalue )
}
