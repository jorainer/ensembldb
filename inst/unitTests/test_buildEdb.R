library(EnsDb.Hsapiens.v75)
edb <- EnsDb.Hsapiens.v75

test_ensDbFromGRanges <- function(){
    load(system.file("YGRanges.RData", package="ensembldb"))
    DB <- ensDbFromGRanges(Y, path=tempdir(), version=75,
                           organism="Homo_sapiens")
    edb <- EnsDb(DB)
    checkEquals(unname(genome(edb)), "GRCh37")

    checkTrue(ensembldb:::checkValidEnsDb(edb))
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

test_organismName <- function() {
    res <- ensembldb:::.organismName("homo_sapiens")
    checkEquals(res, "Homo_sapiens")
    res <- ensembldb:::.abbrevOrganismName("homo_sapiens")
    checkEquals(res, "hsapiens")
    res <- ensembldb:::.makePackageName(dbconn(edb))
    checkEquals(res, "EnsDb.Hsapiens.v75")
}

test_isEnsemblFileName <- function() {
    res <- ensembldb:::isEnsemblFileName("Caenorhabditis_elegans.WS210.60.gtf.gz")
    checkTrue(res)
    res <- ensembldb:::isEnsemblFileName("Caenorhabditis_elegans_fdf.60.dfd.gtf.gz")
    checkTrue(!res)

    fn <- "Caenorhabditis_elegans.WS210.60.gtf.gz"
    res <- ensembldb:::ensemblVersionFromGtfFileName(fn)
    checkEquals(res, "60")
    res <- ensembldb:::organismFromGtfFileName(fn)
    checkEquals(res, "Caenorhabditis_elegans")
    res <- ensembldb:::genomeVersionFromGtfFileName(fn)
    checkEquals(res, "WS210")

    res <- ensembldb:::elementFromEnsemblFilename(fn, which = 1)
    checkEquals(res, "Caenorhabditis_elegans")
    res <- ensembldb:::elementFromEnsemblFilename(fn, which = 2)
    checkEquals(res, "WS210")
    res <- ensembldb:::elementFromEnsemblFilename(fn, which = 3)
    checkEquals(res, "60")
    res <- ensembldb:::elementFromEnsemblFilename(fn, which = 4)
    checkEquals(res, "gtf")
}

test_ensDbFromGtf_Gff <- function() {
    ## db_1 <- ensDbFromGff(system.file("gff/Devosia_geojensis.ASM96941v1.32.gff3.gz",
    ##                                   package="ensembldb"),
    ##                       outfile = tempfile())
    gff <- system.file("gff/Devosia_geojensis.ASM96941v1.32.gff3.gz",
                       package="ensembldb")
    ##Test <- import()
    db_2 <- ensDbFromGtf(system.file("gtf/Devosia_geojensis.ASM96941v1.32.gtf.gz",
                                     package="ensembldb"),
                         outfile = tempfile())
    db_2 <- EnsDb(db_2)
    checkEquals(ensemblVersion(db_2), "32")

}
