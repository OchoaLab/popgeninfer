% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filter_eval.R
\name{filter_eval}
\alias{filter_eval}
\title{Calculate precision and recall for flip/remove/keep categories}
\usage{
filter_eval(true, pred)
}
\arguments{
\item{true}{Vector of true string categories, with values in "flip", "remove", or "keep".}

\item{pred}{Vector of predicted string categories, with values in "flip", "remove", or "keep".}
}
\value{
A tibble with two rows (for flip and remove), and three columns:
\itemize{
\item \code{type}: either "flip" or "remove"
\item \code{precision}: the precision for flips or removes, as indicated by \code{type}
\item \code{recall}: the recall for flips or removes, as indicated by \code{type}
}
}
\description{
This function calculates precision and recall, which are ordinarily binary classification statistics calculated with \code{\link[=binary_eval]{binary_eval()}}, on our 3-category data (flip, remove, keep) by treating the minor categories as their own in turn and the rest are merged.
In particular, for "flip" we consider the binary class \code{TRUE} if it was a flip and \code{FALSE} if it was anything else (in this case either remove or keep).
Similarly, for "remove" we consider the binary class \code{TRUE} if it was a remove and `FALSE for anything else (flip or keep).
No statistics are calculated for "keep", which are ordinarily the most common category.
Note that categories must match "flip" and "remove" in spelling and capitalization to be counted correctly, everything else is implicitly treated as the third category "keep", even if there are additional categories.
}
\examples{
# the input data should look like this:
true <- c('keep', 'keep',   'keep', 'flip',   'flip', 'remove', 'remove', 'remove')
pred <- c('keep', 'flip', 'remove', 'flip', 'remove', 'remove',   'flip',   'keep')

# calculate statistics!
data <- filter_eval( true, pred )
data

}
\seealso{
\code{\link[=binary_eval]{binary_eval()}}, which is used internally to calculate the desired statistics once our multi-category classes are binarized.
}
