#' Returns reverse complements of DNA variants
#'
#' This is a very simple, standalone implementation that covers DNA cases: ACGT lowercase and uppercase, one or multiple letters.
#' Letters that do not fall in those cases are not complemented, but all strings are reversed.
#'
#' @param x A vector of DNA variants
#'
#' @return A vector of the same length where every string has been reversed-complemented.
#'
#' @examples
#' revcomp( c( 'GATTACA', 'TATA' ) )
#'
#' @export
# this is the solution proposed by Megatron here:
# https://stackoverflow.com/questions/20371854/complement-a-dna-sequence
revcomp <- function( x )
    stringi::stri_reverse( chartr( "acgtACGT", "tgcaTGCA", x ) )
