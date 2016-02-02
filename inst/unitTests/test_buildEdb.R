test_ensDbFromGRanges <- function(){
    load(system.file("YGRanges.RData", package="ensembldb"))
    DB <- ensDbFromGRanges(Y, path=tempdir(), version=75,
                           organism="Homo_sapiens")
    edb <- EnsDb(DB)
    checkEquals(unname(genome(edb)), "GRCh37")
}


## Test some internal functions...
test_processEnsemblFileNames <- function(){
    Test <- "Homo_sapiens.GRCh38.83.gtf.gz"
    checkTrue(ensembldb:::isEnsemblFileName(Test))
    checkEquals(ensembldb:::organismFromGtfFileName(Test), "Homo_sapiens")
    checkEquals(ensembldb:::genomeVersionFromGtfFileName(Test), "GRCh38")
    checkEquals(ensembldb:::ensemblVersionFromGtfFileName(Test), "83")

    Test <- "Homo_sapiens.GRCh38.83.chr.gff3.gz"
    checkTrue(ensembldb:::isEnsemblFileName(Test))
    checkEquals(ensembldb:::organismFromGtfFileName(Test), "Homo_sapiens")
    checkEquals(ensembldb:::genomeVersionFromGtfFileName(Test), "GRCh38")
    checkEquals(ensembldb:::ensemblVersionFromGtfFileName(Test), "83")

    Test <- "Gadus_morhua.gadMor1.83.gff3.gz"
    checkTrue(ensembldb:::isEnsemblFileName(Test))
    checkEquals(ensembldb:::organismFromGtfFileName(Test), "Gadus_morhua")
    checkEquals(ensembldb:::genomeVersionFromGtfFileName(Test), "gadMor1")
    checkEquals(ensembldb:::ensemblVersionFromGtfFileName(Test), "83")

    Test <- "Solanum_lycopersicum.GCA_000188115.2.30.chr.gtf.gz"
    checkTrue(ensembldb:::isEnsemblFileName(Test))
    checkEquals(ensembldb:::organismFromGtfFileName(Test), "Solanum_lycopersicum")
    checkEquals(ensembldb:::genomeVersionFromGtfFileName(Test), "GCA_000188115.2")
    checkEquals(ensembldb:::ensemblVersionFromGtfFileName(Test), "30")

    Test <- "ref_GRCh38.p2_top_level.gff3.gz"
    checkEquals(ensembldb:::isEnsemblFileName(Test), FALSE)
    ensembldb:::organismFromGtfFileName(Test)
    checkException(ensembldb:::genomeVersionFromGtfFileName(Test))
    ##checkException(ensembldb:::ensemblVersionFromGtfFileName(Test))
}




