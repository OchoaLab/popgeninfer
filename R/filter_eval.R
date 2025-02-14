#' Calculate precision and recall for flip/remove/keep categories
#'
#' This function calculates precision and recall, which are ordinarily binary classification statistics calculated with [binary_eval()], on our 3-category data (flip, remove, keep) by treating the minor categories as their own in turn and the rest are merged.
#' In particular, for "flip" we consider the binary class `TRUE` if it was a flip and `FALSE` if it was anything else (in this case either remove or keep).
#' Similarly, for "remove" we consider the binary class `TRUE` if it was a remove and `FALSE for anything else (flip or keep).
#' No statistics are calculated for "keep", which are ordinarily the most common category.
#' Note that categories must match "flip" and "remove" in spelling and capitalization to be counted correctly, everything else is implicitly treated as the third category "keep", even if there are additional categories.
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
#' # the input data should look like this:
#' true <- c('keep', 'keep',   'keep', 'flip',   'flip', 'remove', 'remove', 'remove')
#' pred <- c('keep', 'flip', 'remove', 'flip', 'remove', 'remove',   'flip',   'keep')
#'
#' # calculate statistics!
#' data <- filter_eval( true, pred )
#' data
#' 
#' @seealso
#' [binary_eval()], which is used internally to calculate the desired statistics once our multi-category classes are binarized.
#'
#' @export
filter_eval <- function( true, pred ) {
    # calculate data for flips, then removes
    # names get in the way of testing, skip here
    flip <- binary_eval( true == 'flip', pred == 'flip', names = FALSE )
    remove <- binary_eval( true == 'remove', pred == 'remove', names = FALSE )

    # this organization facilitates plotting
    return(
        tibble::tibble(
                    type = c('flip', 'remove'),
                    precision = c( flip[1], remove[1] ),
                    recall = c( flip[2], remove[2] )
                )
    )
}

