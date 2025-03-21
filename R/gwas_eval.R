#' Calculate GWAS-inspired classification statistics for flip/remove/keep categories
#'
#' This function calculates precision, recall, and the underlying true positive (TP), false positive (FP), and false negative (FN) counts from two vectors whose values are "flip", "remove" or "keep", using a special criterion that is specific to what occurs in these cases when biased or flipped loci are tested in a GWAS, differentiating between removals and the rest and whether non-removals are correct or not.
#' In particular, here `TP` counts the number of keep and flip that are correctly predicted; `FP` is the number of predicted keep or flip that are incorrect; `FN` is the number of predicted remove that are incorrect; lastly, the implied `TN` (not explicitly calculated) is the number of correct removals.
#'
#' Categories must match "keep", "flip" and "remove" in spelling and capitalization, otherwise an error is thrown.
#'
#' @param true Vector of true string categories, with values in "flip", "remove", or "keep".
#' @param pred Vector of predicted string categories, with values in "flip", "remove", or "keep".
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
#' # the input data should look something like this:
#' true <- c('keep', 'keep',   'keep', 'flip',   'flip', 'remove', 'remove', 'remove')
#' pred <- c('keep', 'flip', 'remove', 'flip', 'remove', 'remove',   'flip',   'keep')
#'
#' # get stats!
#' data <- gwas_eval( true, pred )
#' data
#'
#' @seealso
#' [binary_eval()], which is used to calculate the desired statistics for logical (binary) classification.
#' 
#' [filter_eval()] for a wrapper applying [gwas_eval()] and variants of [binary_eval()] to 3-category classifications (remove, flip, keep).
#'
#' @export
gwas_eval <- function( true, pred, names = TRUE ) {
    # validate inputs
    check_cats( true )
    check_cats( pred )
    
    # the cases are more complicated here, this counts them as desired
    TP <- sum( true == pred & true != 'remove' )
    FN <- sum( pred == 'remove' & true != 'remove' )
    FP <- sum( pred != 'remove' & true != pred )
    # the rest is a copy of [binary_eval()]
    precision <- TP / ( TP + FP )
    recall <- TP / ( TP + FN )
    x <- c( precision, recall, TP, FP, FN )
    if ( names )
        names( x ) <- c( 'precision', 'recall', 'TP', 'FP', 'FN' )
    return( x )
}

check_cats <- function( vec ) {
    if ( anyNA( vec ) )
        stop( 'Input vector must not have missing values!' )
    if ( !all( vec %in% c('keep', 'flip', 'remove') ) )
        stop( 'Input vector must have all values equal to either "keep", "flip", or "remove"!' )
}
