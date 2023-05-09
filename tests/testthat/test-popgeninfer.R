test_that("af_test_single, af_test work with k=1", {
    # simulate some null data
    p <- 0.4
    n1 <- 10
    n2 <- 13
    x1 <- rbinom( 1, n1, p )
    x2 <- rbinom( 1, n2, p )

    # test it!
    expect_silent(
        stat <- af_test_single( x1, n1, x2, n2 )
    )
    expect_true( stat >= 0 )
    # repeat with more general wrapper
    expect_silent(
        obj <- af_test( x1, n1, x2, n2 )
    )
    expect_equal( obj$stat, stat )
    pval <- obj$pval
    #expect_true( pval >= 0 )
    expect_true( pval <= 1 )
    # in this case null holds, so it should be very unlikely that the result is extremely significant
    expect_true( pval > 1e-10 )
    
    # simulate some alternative data
    # big sample sizes too to best ensure a significant result
    p1 <- 0.4
    p2 <- 0.1
    n1 <- 1000
    n2 <- 1390
    x1 <- rbinom( 1, n1, p1 )
    x2 <- rbinom( 1, n2, p2 )

    # test it!
    expect_silent(
        stat <- af_test_single( x1, n1, x2, n2 )
    )
    expect_true( stat >= 0 )
    # repeat with more general wrapper
    expect_silent(
        obj <- af_test( x1, n1, x2, n2 )
    )
    expect_equal( obj$stat, stat )
    pval <- obj$pval
    expect_true( pval >= 0 )
    #expect_true( pval <= 1 )
    # in this case alt holds, and sample sizes are large, so a significant result is most likely
    expect_true( pval < 1e-5 )
})

test_that( "af_test works with k=3", {
    k <- 3
    # simulate some null data
    # totally random data!
    p1 <- runif( k )
    p2 <- p1
    # same for sample sizes
    n1 <- rpois( k, 10 )
    n2 <- rpois( k, 10 )
    # don't allow zeroes of course
    n1[ n1 == 0 ] <- 1
    n2[ n2 == 0 ] <- 1
    # simulate data in batch
    x1 <- rbinom( k, n1, p1 )
    x2 <- rbinom( k, n2, p2 )
    
    # this is main version
    expect_silent(
        obj <- af_test( x1, n1, x2, n2 )
    )
    stat <- obj$stat
    df <- obj$df
    pval <- obj$pval
    expect_true( stat >= 0 )
    expect_equal( df, k )
    #expect_true( pval >= 0 )
    expect_true( pval <= 1 )
    # in this case null holds, so it should be very unlikely that the result is extremely significant
    expect_true( pval > 1e-10 )

    # another version
    expect_silent(
        obj2 <- likelihoodratio_3(
            x1[1], n1[1], x2[1], n2[1],
            x1[2], n1[2], x2[2], n2[2],
            x1[3], n1[3], x2[3], n2[3]
        )
    )
    expect_equal( obj2, obj )

    # simulate some alternative data
    # totally random data!
    p1 <- runif( k )
    p2 <- runif( k )
    # same for sample sizes
    n1 <- rpois( k, 1000 )
    n2 <- rpois( k, 1000 )
    # don't allow zeroes of course
    n1[ n1 == 0 ] <- 1
    n2[ n2 == 0 ] <- 1
    # simulate data in batch
    x1 <- rbinom( k, n1, p1 )
    x2 <- rbinom( k, n2, p2 )
    
    # this is main version
    expect_silent(
        obj <- af_test( x1, n1, x2, n2 )
    )
    stat <- obj$stat
    df <- obj$df
    pval <- obj$pval
    expect_true( stat >= 0 )
    expect_equal( df, k )
    expect_true( pval >= 0 )
    #expect_true( pval <= 1 )
    # in this case alt holds, and sample sizes are large, so a significant result is most likely
    expect_true( pval < 1e-5 )
    
    # another version
    expect_silent(
        obj2 <- likelihoodratio_3(
            x1[1], n1[1], x2[1], n2[1],
            x1[2], n1[2], x2[2], n2[2],
            x1[3], n1[3], x2[3], n2[3]
        )
    )
    expect_equal( obj2, obj )
})
