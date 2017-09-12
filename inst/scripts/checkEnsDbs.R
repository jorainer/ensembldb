#' @description Check EnsDb sqlite files found in the specified folder.
#'
#' @param x \code{character(1)} with the folder in which we're looking for EnsDb
#'     objects.
#'
#' @author Johannes Rainer
#' 
#' @noRd
#'
#' @examples
#' dir <- "/Users/jo/tmp/ensdb_20"
checkEnsDbs <- function(x) {
    edbs <- dir(x, pattern = ".sqlite$", full.names = TRUE)
    for (i in 1:length(edbs)) {
        message("\nChecking EnsDb: ", basename(edbs[i]))
        edb <- EnsDb(edbs[i])
        ensembldb:::validateEnsDb(edb)
        ensembldb:::checkValidEnsDb(edb)
        ## Now check also some query calls:
        gns <- genes(edb)
        message(" version: ", ensembldb:::dbSchemaVersion(edb))
        message(" OK")
    }
}
