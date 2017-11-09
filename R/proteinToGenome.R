## Code related to mapping of coordinates within proteins to genomic
## coordinates.


#' @description Convert within-protein coordinates to transcript (CDS) relative
#'    coordinates.
#'
#' @param x `IRanges` object with coordinates within a protein sequence.
#'
#' @return `IRanges` with the coordinates within the coding region (!) of the
#'     encoding transcript.
#'
#' @author Johannes Rainer
#'
#' @md
#'
#' @noRd
#'
#' @examples
#'
#' prt <- IRanges(start = 23, end = 27)
#' .proteinCoordsToTx(prt)
.proteinCoordsToTx <- function(x) {
    end(x) <- end(x) * 3
    start(x) <- 1 + (start(x) - 1) * 3
    x
}


#' @description Checks for each provided range whether its width matches the
#'     expected length. This is to check which transcripts' CDS length matches
#'     the expected length base on the protein's sequence.
#'
#' @param x `IRanges` object.
#'
#' @param length `integer(1)`
#'
#' @return `logical` with length equal to length `x`.
#'
#' @author Johannes Rainer
#'
#' @md
#'
#' @noRd
.expectedWidth <- function(x, length) {
}

## What has to be done?
## 1) for a given protein, fetch CDS of all annotated tx
## 2) check if the lengths of the CDS match, select all tx with perfect match,
##    or those with a longer CDS.
## 3) Transform protein-relative coords to tx relative coords.
## 4) Match tx cds-relative coords to individual exons - CAVE strand
