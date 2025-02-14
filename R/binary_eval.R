#' Calculate binary classification statistics
#'
#' This function calculates precision, recall, and the underlying true positive (TP), false positive (FP), and false negative (FN) counts from two binary vectors.
#' This general function works well, but it is potentially inefficient for generating curves from large data with a monotonic threshold (such as a single p-value), so it is intended for more complex thresholding scenarios.
#'
#' @param true The vector of true logical classifications.
#' @param pred The vector of predicted logical classifications.
#' @param names Logical for whether the return vector should have names or not.
#'
#' @return A named vector (unless `names = FALSE`) with five values in this order:
#' - `precision`: `TP / ( TP + FP )`
#' - `recall`: `TP / ( TP + FN )`
#' - `TP`: the number of rows where both `true` and `pred` are `TRUE`
#' - `FP`: the number of rows where `true` is `FALSE` and `pred` is `TRUE`
#' - `FN`: the number of rows where `true` is `TRUE` and `pred` is `FALSE`
#'
#' @examples
#' # input classifications, true and prediction
#' true <- c(TRUE, TRUE, TRUE,  TRUE, FALSE, FALSE, FALSE)
#' pred <- c(TRUE, TRUE, TRUE, FALSE,  TRUE,  TRUE, FALSE)
#'
#' # get stats!
#' data <- binary_eval( true, pred )
#' data
#'
#' @seealso
#' [filter_eval()] for a wrapper applied to 3-category classifications (remove, flip, keep).
#'
#' @export
# inputs in these cases are booleans already!  These cases are easy, though they're potentially inefficient in cases where thresholds give nested answers (here it's probably not so, not generally)
binary_eval <- function( true, pred, names = TRUE ) {
    TP <- sum( true & pred )
    FP <- sum( !true & pred )
    FN <- sum( true & !pred )
    precision <- TP / ( TP + FP )
    recall <- TP / ( TP + FN )
    x <- c( precision, recall, TP, FP, FN )
    if ( names )
        names( x ) <- c( 'precision', 'recall', 'TP', 'FP', 'FN' )
    return( x )
}
