library(EnsDb.Hsapiens.v75)
edb <- EnsDb.Hsapiens.v75

test_ensDbFromGRanges <- function(){
    load(system.file("YGRanges.RData", package="ensembldb"))
    suppressWarnings(
        DB <- ensDbFromGRanges(Y, path=tempdir(), version=75,
                               organism="Homo_sapiens")
    )
    edb <- EnsDb(DB)
    checkEquals(unname(genome(edb)), "GRCh37")

    Test <- makeEnsembldbPackage(DB, destDir = tempdir(),
                                 version = "0.0.1", author = "J Rainer",
                                 maintainer = "")
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
    gff <- system.file("gff/Devosia_geojensis.ASM96941v1.32.gff3.gz",
                       package="ensembldb")
    gtf <- system.file("gtf/Devosia_geojensis.ASM96941v1.32.gtf.gz",
                       package="ensembldb")
    suppressWarnings(
        db_gff <- EnsDb(ensDbFromGff(gff, outfile = tempfile()))
    )
    suppressWarnings(
        db_gtf <- EnsDb(ensDbFromGtf(gtf, outfile = tempfile()))
    )
    checkEquals(ensemblVersion(db_gtf), "32")
    checkEquals(ensemblVersion(db_gff), "32")

    res <- ensembldb:::compareChromosomes(db_gtf, db_gff)
    checkEquals(res, "OK")
    res <- ensembldb:::compareGenes(db_gtf, db_gff)
    checkEquals(res, "WARN")  ## differences in gene names and Entrezid.
    res <- ensembldb:::compareTx(db_gtf, db_gff)
    checkEquals(res, "OK")
    res <- ensembldb:::compareExons(db_gtf, db_gff)
    checkEquals(res, "OK")
    ## Compare them all in one call
    res <- ensembldb:::compareEnsDbs(db_gtf, db_gff)
    checkEquals(unname(res["metadata"]), "NOTE")
    checkEquals(unname(res["chromosome"]), "OK")
    checkEquals(unname(res["transcript"]), "OK")
    checkEquals(unname(res["exon"]), "OK")
}

test_checkExtractVersions <- function() {
    fn <- "Devosia_geojensis.ASM96941v1.32.gff3.gz"
    res <- ensembldb:::.checkExtractVersions(fn)
    checkEquals(unname(res["organism"]), "Devosia_geojensis")
    checkEquals(unname(res["genomeVersion"]), "ASM96941v1")
    checkEquals(unname(res["version"]), "32")
    suppressWarnings(
        res <- ensembldb:::.checkExtractVersions(fn, organism = "Homo_sapiens")
        )
    checkEquals(unname(res["organism"]), "Homo_sapiens")
    checkException(ensembldb:::.checkExtractVersions("afdfhjd"))
}

test_buildMetadata <- function() {
    res <- ensembldb:::buildMetadata(organism = "Mus_musculus",
                                     ensemblVersion = "88",
                                     genomeVersion = "38")
    checkEquals(colnames(res), c("name", "value"))
    checkEquals(res[res$name == "Organism", "value"], "Mus_musculus")
}

test_guessDatabaseName <- function() {
    ## Testing real case examples.
    genome <- "Rnor_5.0"
    organism <- "Rattus_norvegicus"
    ensembl <- "75"
    res <- ensembldb:::.guessDatabaseName(organism, ensembl)
    expect <- "rattus_norvegicus_core_75"
    checkEquals(res, expect)

    genome <- "GRCm38"
    organism <- "Mus_musculus"
    expect <- "mus_musculus_core_75_38"
    res <- ensembldb:::.guessDatabaseName(organism, ensembl,
                                          genome = genome)
    checkEquals(expect, res)
}

test_getEnsemblMysqlUrl <- function() {
    ## Only run this if we have access to Ensembl.
    tmp <- try(
        RCurl::getURL(ensembldb:::.ENSEMBL_URL, dirlistonly = TRUE)
    )
    if (!is(tmp, "try-error")) {
        res <- ensembldb:::.getEnsemblMysqlUrl(type = "ensembl",
                                               organism = "macaca mulatta",
                                               ensembl = 85)
        checkEquals(res, paste0(ensembldb:::.ENSEMBL_URL, "release-85/",
                                "mysql/macaca_mulatta_core_85_10"))
        check_getReadMysqlTable(res)
        ## Next.
        res <- ensembldb:::.getEnsemblMysqlUrl(type = "ensembl",
                                               organism = "Bos taurus",
                                               ensembl = 61)
        checkEquals(res, paste0(ensembldb:::.ENSEMBL_URL, "release-61/",
                                "mysql/bos_taurus_core_61_4j"))
        check_getReadMysqlTable(res)
        ## Next
        res <- ensembldb:::.getEnsemblMysqlUrl(type = "ensembl",
                                               organism = "Ficedula albicollis",
                                               ensembl = 77)
        checkEquals(res, paste0(ensembldb:::.ENSEMBL_URL, "release-77/",
                                "mysql/ficedula_albicollis_core_77_1"))
    }
    ## ensemblgenomes
    tmp <- try(
        RCurl::getURL(ensembldb:::.ENSEMBLGENOMES_URL, dirlistonly = TRUE)
    )
    if (!is(tmp, "try-error")) {
        ## check fungi
        res <- ensembldb:::.getEnsemblMysqlUrl(type = "ensemblgenomes",
                                               organism = "fusarium_oxysporum",
                                               ensembl = 21)
        db_name <- "fusarium_oxysporum_core_21_74_2"
        checkEquals(res, paste0(ensembldb:::.ENSEMBLGENOMES_URL, "release-21/",
                                "fungi/mysql/", db_name))
        check_getReadMysqlTable(res)
        ## Next one
        db_name <- "solanum_lycopersicum_core_28_81_250"
        res <- ensembldb:::.getEnsemblMysqlUrl(type = "ensemblgenomes",
                                               organism = "solanum_lycopersicum",
                                               ensembl = 28)
        checkEquals(res, paste0(ensembldb:::.ENSEMBLGENOMES_URL, "release-28/",
                                "plants/mysql/", db_name))
        check_getReadMysqlTable(res)
    }
}

check_getReadMysqlTable <- function(url) {
    res <- ensembldb:::.getReadMysqlTable(url, "coord_system.txt.gz",
                                          colnames = c("coord_system_id",
                                                       "species_id",
                                                       "name", "version",
                                                       "rank", "attrib"))
    checkTrue(nrow(res) > 0)
}

test_getSeqlengthsFromMysqlFolder <- function(){
    tmp <- try(
        RCurl::getURL(ensembldb:::.ENSEMBL_URL, dirlistonly = TRUE)
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
        checkTrue(all(names(sl) %in% names(sl_2)))
        checkEquals(sl, sl_2[names(sl)])
    }
}

notrun_test_getSeqlengthsFromMysqlFolder <- function() {
    ## Test this for some more seqlengths.
    library(EnsDb.Rnorvegicus.v79)
    db <- EnsDb.Rnorvegicus.v79
    seq_info <- seqinfo(db)
    seq_lengths <- ensembldb:::.getSeqlengthsFromMysqlFolder(
        organism = "Rattus norvegicus", ensembl = 79,
        seqnames = seqlevels(seq_info))
    sl <- seqlengths(seq_info)
    sl_2 <- seq_lengths$length
    names(sl_2) <- rownames(seq_lengths)
    checkEquals(sl, sl_2)
    ## Mus musculus
}

notrun_test_ensDbFromGtf_Gff_AH <- function() {
    gtf <- paste0("/Users/jo/Projects/EnsDbs/80/caenorhabditis_elegans/",
                  "Caenorhabditis_elegans.WBcel235.80.gtf.gz")
    outf <- tempfile()
    db <- ensDbFromGtf(gtf = gtf, outfile = outf)
    ## use Gff
    gff <- paste0("/Users/jo/Projects/EnsDbs/84/canis_familiaris/gff3/",
                  "Canis_familiaris.CanFam3.1.84.gff3.gz")
    outf <- tempfile()
    db <- ensDbFromGff(gff, outfile = outf)

    ## Checking one from ensemblgenomes:
    gtf <- paste0("/Users/jo/Projects/EnsDbs/ensemblgenomes/30/",
                  "solanum_lycopersicum/",
                  "Solanum_lycopersicum.GCA_000188115.2.30.chr.gtf.gz"
                  )
    outf <- tempfile()
    db <- ensDbFromGtf(gtf = gtf, outfile = outf)
    gtf <- paste0("/Users/jo/Projects/EnsDbs/ensemblgenomes/30/",
                  "solanum_lycopersicum/",
                  "Solanum_lycopersicum.GCA_000188115.2.30.gtf.gz"
                  )
    outf <- tempfile()
    db <- ensDbFromGtf(gtf = gtf, outfile = outf)

    ## AH
    library(AnnotationHub)
    ah <- AnnotationHub()
    query(ah, c("release-83", "gtf"))
    ah_1 <- ah["AH50418"]
    db <- ensDbFromAH(ah_1, outfile = outf)
    ah_2 <- ah["AH50352"]
    db <- ensDbFromAH(ah_2, outfile = outf)
}
