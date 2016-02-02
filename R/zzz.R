
.onLoad <- function(libname, pkgname){
    op <- options()
    ## What should be returned by default if the seqnames can not be mapped based
    ## on the style set by seqlevelsStyle<-
    ## Options:
    ## + NA or any other value: return this value for such cases.
    ## + MISSING: stop and throw an error.
    ## + ORIGINAL: return the original seqnames.
    opts.ens <- list(useFancyQuotes=FALSE,
                     ensembldb.seqnameNotFound="ORIGINAL")
    options(opts.ens)
    invisible()
}

