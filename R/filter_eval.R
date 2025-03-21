#' Calculate precision and recall for flip/remove/keep categories using three criteria
#'
#' This function calculates precision and recall, which are ordinarily binary classification statistics calculated with [binary_eval()], on our 3-category data (flip, remove, keep) by using a GWAS-inspired criterion described in [gwas_eval()] and also treating the minor categories as their own in turn and the rest are merged using [binary_eval()].
#' In particular, for "flip" we consider the binary class `TRUE` if it was a flip and `FALSE` if it was anything else (in this case either remove or keep).
#' Similarly, for "remove" we consider the binary class `TRUE` if it was a remove and `FALSE for anything else (flip or keep).
#' No statistics are calculated for "keep", which are ordinarily the most common category, though that will be reflected in the other categories.
#' 
#' Categories must match "keep", "flip" and "remove" in spelling and capitalization, otherwise an error is thrown.
#'
#' @param true Vector of true string categories, with values in "flip", "remove", or "keep".
#' @param pred Vector of predicted string categories, with values in "flip", "remove", or "keep".
#'
#' @return A tibble with two rows (for flip and remove), and three columns:
#' - `type`: either "flip" or "remove"
#' - `precision`: the precision for flips or removes, as indicated by `type`
#' - `recall`: the recall for flips or removes, as indicated by `type`
#'
#' @examples
#' # the input data should look something like this:
#' true <- c('keep', 'keep',   'keep', 'flip',   'flip', 'remove', 'remove', 'remove')
#' pred <- c('keep', 'flip', 'remove', 'flip', 'remove', 'remove',   'flip',   'keep')
#'
#' # calculate statistics!
#' data <- filter_eval( true, pred )
#' data
#' 
#' @seealso
#' [gwas_eval()], which is used internally to calculate statistics using a custom, GWAS-inspired criterion.
#' 
#' [binary_eval()], which is used internally to calculate statistics once our multi-category classes are binarized.
#'
#' @export
filter_eval <- function( true, pred ) {
    # calculate gwas-inspired metric, which also validates inputs
    gwas <- gwas_eval( true, pred, names = FALSE )
    # calculate data for flips, then removes
    # names get in the way of testing, skip here
    flip <- binary_eval( true == 'flip', pred == 'flip', names = FALSE )
    remove <- binary_eval( true == 'remove', pred == 'remove', names = FALSE )
    
    # this organization facilitates plotting
    return(
        tibble::tibble(
                    type = c('gwas', 'flip', 'remove'),
                    precision = c( gwas[1], flip[1], remove[1] ),
                    recall = c( gwas[2], flip[2], remove[2] )
                )
    )
}

