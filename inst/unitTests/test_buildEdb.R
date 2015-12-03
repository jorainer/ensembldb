test_ensDbFromGRanges <- function(){
    load(system.file("YGRanges.RData", package="ensembldb"))
    DB <- ensDbFromGRanges(Y, path=tempdir(), version=75,
                           organism="Homo_sapiens", verbose=TRUE)
    edb <- EnsDb(DB)
    checkEquals(unname(genome(edb)), "GRCh37")
}

