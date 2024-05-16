# this is testing a single SNP, so the inputs are each vectors with length equal the number of ancestries
af_test_or_single <- function( x1, n1, x2, n2, level = 0.95 ) {
    # in small simulations we encounter edge cases frequently, which produce errors for the regular code; handle those now!
    # in all of these, all sample probabilities are 1 or all 0, that's the problem!
    # since the probabilities are equal, makes sense to assume that the OR=1, but otherwise assign the worst CIs and p-values
    if ( all( c(x1, x2) == 0 ) || all( c(n1-x1, n2-x2) == 0 ) )
        return( c(0, -Inf, Inf, 1) )

    # inputs were tested outside, here assume they are all the same length
    k <- length( x1 )
    # reorganized data for glm
    # observed successes and failures (alleles), one row per dataset
    # goes through all ancestries first for first dataset, then all ancestries for second dataset
    y <- rbind( cbind( x1, n1 - x1 ), cbind( x2, n2 - x2 ) )
    # dataset 1 indicator (second is dataset 2)
    x <- c( rep.int( 1, k ), rep.int( 0, k ) )
    # ancestry indicator (force categorical)
    a <- as.factor( c( 1:k, 1:k ) )
    # fit model!
    # in these cases glm can produce warnings, such as "glm.fit: fitted probabilities numerically 0 or 1 occurred", always silence it
    suppressWarnings(
        obj <- if ( k == 1 ) {
                   stats::glm( y ~ x, family = stats::binomial( link = 'logit' ) )
               } else
                   stats::glm( y ~ x + a, family = stats::binomial( link = 'logit' ) )
    )
    # get summary
    objsum <- summary( obj )
    # and get coefficients from there
    objcoef <- stats::coef( objsum )
    # this is desired coefficient, the others we can ignore
    beta <- objcoef[ 'x', 'Estimate' ]
    # and p-value using this model (different from af_test because alternative is restricted, so it's a 1-df test)
    pval <- objcoef[ 'x', 'Pr(>|z|)' ]
    # use special function to get CIs
    # but silence it, there's no other way to silence it unfortunately
    # can also produce additional warnings such as "glm.fit: fitted probabilities numerically 0 or 1 occurred"
    obj2 <- suppressWarnings( suppressMessages( stats::confint( obj, 'x', level = level ) ) )
    # in extreme cases a bound can be NA, correct that now with an appropriate infinity instead
    if ( is.na( obj2[1] ) ) obj2[1] <- -Inf
    if ( is.na( obj2[2] ) ) obj2[2] <- Inf
    # return all values in a vector for now
    return( c(beta, obj2, pval) )
}
