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

#' @param x `character(1)` with the directory containing EnsDb SQLite files.
#'
#' @return the species/organism name of the databases that have to be re-created
check_gc_content <- function(x) {
    edbs <- dir(x, pattern = ".sqlite$", full.names = TRUE)
    failed <- lapply(edbs, function(edb) {
        edb <- EnsDb(edb)
        if (!any(colnames(mcols(transcripts(edb))) == "gc_content"))
            paste0(tolower(organism(edb)), collapse = "_")
        else NULL
    })
    failed[lengths(failed)]
}
