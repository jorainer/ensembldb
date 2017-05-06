
test_that(".organismName, .abbrevOrganismName and .makePackageName works", {
    res <- ensembldb:::.organismName("homo_sapiens")
    expect_equal(res, "Homo_sapiens")
    res <- ensembldb:::.abbrevOrganismName("homo_sapiens")
    expect_equal(res, "hsapiens")
    res <- ensembldb:::.makePackageName(dbconn(edb))
    expect_equal(res, "EnsDb.Hsapiens.v75")
})

test_that("ensDbFromGRanges works", {
    load(system.file("YGRanges.RData", package="ensembldb"))
    suppressWarnings(
        DB <- ensDbFromGRanges(Y, path=tempdir(), version=75,
                               organism="Homo_sapiens", skip = TRUE)
    )
    db <- EnsDb(DB)
    expect_equal(unname(genome(db)), "GRCh37")

    Test <- makeEnsembldbPackage(DB, destDir = tempdir(),
                                 version = "0.0.1", author = "J Rainer",
                                 maintainer = "")
    expect_true(ensembldb:::checkValidEnsDb(db))
})

test_that("ensDbFromGtf and Gff works", {
    gff <- system.file("gff/Devosia_geojensis.ASM96941v1.32.gff3.gz",
                       package="ensembldb")
    gtf <- system.file("gtf/Devosia_geojensis.ASM96941v1.32.gtf.gz",
                       package="ensembldb")
    suppressWarnings(
        db_gff <- EnsDb(ensDbFromGff(gff, outfile = tempfile(), skip = TRUE))
    )
    suppressWarnings(
        db_gtf <- EnsDb(ensDbFromGtf(gtf, outfile = tempfile(), skip = TRUE))
    )
    expect_equal(ensemblVersion(db_gtf), "32")
    expect_equal(ensemblVersion(db_gff), "32")

    res <- ensembldb:::compareChromosomes(db_gtf, db_gff)
    expect_equal(res, "OK")
    res <- ensembldb:::compareGenes(db_gtf, db_gff)
    expect_equal(res, "WARN")  ## differences in gene names and Entrezid.
    res <- ensembldb:::compareTx(db_gtf, db_gff)
    expect_equal(res, "OK")
    res <- ensembldb:::compareExons(db_gtf, db_gff)
    expect_equal(res, "OK")
    ## Compare them all in one call
    res <- ensembldb:::compareEnsDbs(db_gtf, db_gff)
    expect_equal(unname(res["metadata"]), "NOTE")
    expect_equal(unname(res["chromosome"]), "OK")
    expect_equal(unname(res["transcript"]), "OK")
    expect_equal(unname(res["exon"]), "OK")
})

test_that("isEnsemblFileName", {
    res <- ensembldb:::isEnsemblFileName("Caenorhabditis_elegans.WS210.60.gtf.gz")
    expect_true(res)
    res <- ensembldb:::isEnsemblFileName("Caenorhabditis_elegans_fdf.60.dfd.gtf.gz")
    expect_true(!res)

    fn <- "Caenorhabditis_elegans.WS210.60.gtf.gz"
    res <- ensembldb:::ensemblVersionFromGtfFileName(fn)
    expect_equal(res, "60")
    res <- ensembldb:::organismFromGtfFileName(fn)
    expect_equal(res, "Caenorhabditis_elegans")
    res <- ensembldb:::genomeVersionFromGtfFileName(fn)
    expect_equal(res, "WS210")

    res <- ensembldb:::elementFromEnsemblFilename(fn, which = 1)
    expect_equal(res, "Caenorhabditis_elegans")
    res <- ensembldb:::elementFromEnsemblFilename(fn, which = 2)
    expect_equal(res, "WS210")
    res <- ensembldb:::elementFromEnsemblFilename(fn, which = 3)
    expect_equal(res, "60")
    res <- ensembldb:::elementFromEnsemblFilename(fn, which = 4)
    expect_equal(res, "gtf")
})

test_that("processEnsemblFileNames works", {
    Test <- "Homo_sapiens.GRCh38.83.gtf.gz"
    expect_true(ensembldb:::isEnsemblFileName(Test))
    expect_equal(ensembldb:::organismFromGtfFileName(Test), "Homo_sapiens")
    expect_equal(ensembldb:::genomeVersionFromGtfFileName(Test), "GRCh38")
    expect_equal(ensembldb:::ensemblVersionFromGtfFileName(Test), "83")

    Test <- "Homo_sapiens.GRCh38.83.chr.gff3.gz"
    expect_true(ensembldb:::isEnsemblFileName(Test))
    expect_equal(ensembldb:::organismFromGtfFileName(Test), "Homo_sapiens")
    expect_equal(ensembldb:::genomeVersionFromGtfFileName(Test), "GRCh38")
    expect_equal(ensembldb:::ensemblVersionFromGtfFileName(Test), "83")

    Test <- "Gadus_morhua.gadMor1.83.gff3.gz"
    expect_true(ensembldb:::isEnsemblFileName(Test))
    expect_equal(ensembldb:::organismFromGtfFileName(Test), "Gadus_morhua")
    expect_equal(ensembldb:::genomeVersionFromGtfFileName(Test), "gadMor1")
    expect_equal(ensembldb:::ensemblVersionFromGtfFileName(Test), "83")

    Test <- "Solanum_lycopersicum.GCA_000188115.2.30.chr.gtf.gz"
    expect_true(ensembldb:::isEnsemblFileName(Test))
    expect_equal(ensembldb:::organismFromGtfFileName(Test), "Solanum_lycopersicum")
    expect_equal(ensembldb:::genomeVersionFromGtfFileName(Test), "GCA_000188115.2")
    expect_equal(ensembldb:::ensemblVersionFromGtfFileName(Test), "30")

    Test <- "ref_GRCh38.p2_top_level.gff3.gz"
    expect_equal(ensembldb:::isEnsemblFileName(Test), FALSE)
    ensembldb:::organismFromGtfFileName(Test)
    expect_error(ensembldb:::genomeVersionFromGtfFileName(Test))
    ##checkException(ensembldb:::ensemblVersionFromGtfFileName(Test))
})

test_that("checkExtractVersions works", {
    fn <- "Devosia_geojensis.ASM96941v1.32.gff3.gz"
    res <- ensembldb:::.checkExtractVersions(fn)
    expect_equal(unname(res["organism"]), "Devosia_geojensis")
    expect_equal(unname(res["genomeVersion"]), "ASM96941v1")
    expect_equal(unname(res["version"]), "32")
    suppressWarnings(
        res <- ensembldb:::.checkExtractVersions(fn, organism = "Homo_sapiens")
        )
    expect_equal(unname(res["organism"]), "Homo_sapiens")
    expect_error(ensembldb:::.checkExtractVersions("afdfhjd"))
})

test_that("buildMetadata works", {
    res <- ensembldb:::buildMetadata(organism = "Mus_musculus",
                                     ensemblVersion = "88",
                                     genomeVersion = "38")
    expect_equal(colnames(res), c("name", "value"))
    expect_equal(res[res$name == "Organism", "value"], "Mus_musculus")
})

test_that("guessDatabaseName works", {
    ## Testing real case examples.
    genome <- "Rnor_5.0"
    organism <- "Rattus_norvegicus"
    ensembl <- "75"
    res <- ensembldb:::.guessDatabaseName(organism, ensembl)
    expect <- "rattus_norvegicus_core_75"
    expect_equal(res, expect)

    genome <- "GRCm38"
    organism <- "Mus_musculus"
    expect <- "mus_musculus_core_75_38"
    res <- ensembldb:::.guessDatabaseName(organism, ensembl,
                                          genome = genome)
    expect_equal(expect, res)
})

test_that("getEnsemblMysqlUrl works", {
    check_getReadMysqlTable <- function(url) {
        res <- ensembldb:::.getReadMysqlTable(url, "coord_system.txt.gz",
                                              colnames = c("coord_system_id",
                                                           "species_id",
                                                           "name", "version",
                                                           "rank", "attrib"))
        expect_true(nrow(res) > 0)
    }

    ## Only run this if we have access to Ensembl.
    tmp <- try(
        RCurl::getURL(ensembldb:::.ENSEMBL_URL, dirlistonly = TRUE,
                      .opts = list(timeout = 5, maxredirs = 2))
    )
    if (!is(tmp, "try-error")) {
        res <- ensembldb:::.getEnsemblMysqlUrl(type = "ensembl",
                                               organism = "macaca mulatta",
                                               ensembl = 85)
        expect_equal(res, paste0(ensembldb:::.ENSEMBL_URL, "release-85/",
                                "mysql/macaca_mulatta_core_85_10"))
        check_getReadMysqlTable(res)
        ## Next.
        res <- ensembldb:::.getEnsemblMysqlUrl(type = "ensembl",
                                               organism = "Bos taurus",
                                               ensembl = 61)
        expect_equal(res, paste0(ensembldb:::.ENSEMBL_URL, "release-61/",
                                "mysql/bos_taurus_core_61_4j"))
        check_getReadMysqlTable(res)
        ## Next
        res <- ensembldb:::.getEnsemblMysqlUrl(type = "ensembl",
                                               organism = "Ficedula albicollis",
                                               ensembl = 77)
        expect_equal(res, paste0(ensembldb:::.ENSEMBL_URL, "release-77/",
                                "mysql/ficedula_albicollis_core_77_1"))
    }
    ## ensemblgenomes
    tmp <- try(
        RCurl::getURL(ensembldb:::.ENSEMBLGENOMES_URL, dirlistonly = TRUE,
                      .opts = list(timeout = 5, maxredirs = 2))
    )
    if (!is(tmp, "try-error")) {
        ## check fungi
        res <- ensembldb:::.getEnsemblMysqlUrl(type = "ensemblgenomes",
                                               organism = "fusarium_oxysporum",
                                               ensembl = 21)
        db_name <- "fusarium_oxysporum_core_21_74_2"
        expect_equal(res, paste0(ensembldb:::.ENSEMBLGENOMES_URL, "release-21/",
                                "fungi/mysql/", db_name))
        check_getReadMysqlTable(res)
        ## Next one
        db_name <- "solanum_lycopersicum_core_28_81_250"
        res <- ensembldb:::.getEnsemblMysqlUrl(type = "ensemblgenomes",
                                               organism = "solanum_lycopersicum",
                                               ensembl = 28)
        expect_equal(res, paste0(ensembldb:::.ENSEMBLGENOMES_URL, "release-28/",
                                "plants/mysql/", db_name))
        check_getReadMysqlTable(res)
    }
})

test_that("getSeqlengthsFromMysqlFolder works", {
    library(curl)
    ch <- new_handle(timeout = 5)
    handle_setopt(ch, timeout = 5)
    tmp <- try(
        ## RCurl::getURL(ensembldb:::.ENSEMBL_URL, dirlistonly = TRUE,
        ##               .opts = list(timeout = 5, maxredirs = 2))
        readLines(curl(ensembldb:::.ENSEMBL_URL, handle = ch))
    )
    if (!is(tmp, "try-error")) {
        ## Compare seqlengths we've in EnsDb.Hsapiens.v75 with the expected
        ## ones.
        seq_info <- seqinfo(edb)
        seq_lengths <- ensembldb:::.getSeqlengthsFromMysqlFolder(
            organism = "Homo sapiens", ensembl = 75,
            seqnames = seqlevels(seq_info))
        sl <- seqlengths(seq_info)
        sl_2 <- seq_lengths$length
        names(sl_2) <- rownames(seq_lengths)
        expect_true(all(names(sl) %in% names(sl_2)))
        expect_equal(sl, sl_2[names(sl)])
    }
})

