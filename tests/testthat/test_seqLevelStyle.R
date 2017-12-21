## Tests related to setting seqLevelsStyle

test_that("seqlevelsStyle works", {
    orig <- getOption("ensembldb.seqnameNotFound")
    options(ensembldb.seqnameNotFound = NA)
    edb <- EnsDb.Hsapiens.v86
    SL <- seqlevels(edb)
    ucscs <- paste0("chr", c(1:22, "X", "Y", "M"))
    seqlevelsStyle(edb) <- "UCSC"
    suppressWarnings(
        SL2 <- seqlevels(edb)
    )
    expect_equal(sort(ucscs), sort(SL2[!is.na(SL2)]))
    ## Check if we throw an error message
    options(ensembldb.seqnameNotFound = "MISSING")
    expect_error(seqlevels(edb))
    ## Check if returning original names works.
    options(ensembldb.seqnameNotFound = "ORIGINAL")
    suppressWarnings(
        SL3 <- seqlevels(edb)
    )
    idx <- which(SL3 %in% ucscs)
    expect_equal(sort(SL[-idx]), sort(SL3[-idx]))
    options(ensembldb.seqnameNotFound=orig)
})

test_that("seqinfo works with seqlevelsStyle", {
    edb <- EnsDb.Hsapiens.v86
    orig <- getOption("ensembldb.seqnameNotFound")
    options(ensembldb.seqnameNotFound="MISSING")
    seqlevelsStyle(edb) <- "UCSC"
    expect_error(seqinfo(edb))
    options(ensembldb.seqnameNotFound="ORIGINAL")
    suppressWarnings(
        si <- seqinfo(edb)
    )
    options(ensembldb.seqnameNotFound=orig)
})

test_that("getWhat works with seqlevelsStyle", {
    orig <- getOption("ensembldb.seqnameNotFound")
    edb <- EnsDb.Hsapiens.v86
    seqlevelsStyle(edb) <- "Ensembl"
    ensRes <- ensembldb:::getWhat(edb, columns=c("seq_name", "seq_strand"))
    seqlevelsStyle(edb) <- "UCSC"
    suppressWarnings(
        ucscRes <- ensembldb:::getWhat(edb, columns=c("seq_name", "seq_strand"))
    )
    seqlevelsStyle(edb) <- "NCBI"
    suppressWarnings(
        ncbiRes <- ensembldb:::getWhat(edb, columns=c("seq_name", "seq_strand"))
    )
    options(ensembldb.seqnameNotFound=orig)
})

test_that("SeqNameFilter works with seqlevelsStyle", {
    orig <- getOption("ensembldb.seqnameNotFound")
    options(ensembldb.seqnameNotFound="MISSING")
    edb <- EnsDb.Hsapiens.v86
    seqlevelsStyle(edb) <- "Ensembl"
    snf <- SeqNameFilter("chrX")
    snfEns <- SeqNameFilter(c("X", "Y"))
    snfNo <- SeqNameFilter(c("bla", "blu"))
    snfSomeNo <- SeqNameFilter(c("bla", "X"))

    seqlevelsStyle(edb) <- "Ensembl"
    expect_equal(value(snf), "chrX")
    ## That makes no sense for a query though.
    expect_equal(value(snf), "chrX")
    expect_equal(value(snfEns), c("X", "Y"))
    seqlevelsStyle(edb) <- "UCSC"
    expect_equal(ensembldb:::ensDbQuery(snf, edb), "gene.seq_name = 'X'")
    expect_error(ensembldb:::ensDbQuery(snfEns, edb))
    expect_error(ensembldb:::ensDbQuery(snfNo, edb))
    expect_error(ensembldb:::ensDbQuery(snfSomeNo, edb))

    ## Setting the options to "ORIGINAL"
    options(ensembldb.seqnameNotFound="ORIGINAL")
    expect_equal(ensembldb:::ensDbQuery(snf, edb), "gene.seq_name = 'X'")
    suppressWarnings(
        expect_equal(ensembldb:::ensDbQuery(snfEns, edb),
                    "gene.seq_name in ('X','Y')")
    )
    suppressWarnings(
        expect_equal(ensembldb:::ensDbQuery(snfNo, edb),
                    "gene.seq_name in ('bla','blu')")
    )
    suppressWarnings(
        expect_equal(ensembldb:::ensDbQuery(snfSomeNo, edb),
                    "gene.seq_name in ('bla','X')")
    )
    ##
    snf <- SeqNameFilter(c("chrX", "Y"))
    suppressWarnings(
        expect_equal(ensembldb:::ensDbQuery(snf, edb),
                    "gene.seq_name in ('X','Y')")
    )
    options(ensembldb.seqnameNotFound=orig)
})

test_that("genes works with seqlevelsStyles", {
    orig <- getOption("ensembldb.seqnameNotFound")
    edb <- EnsDb.Hsapiens.v86
    ## Here we want to test whether the result returned by the function does really
    ## work when changing the seqnames.
    seqlevelsStyle(edb) <- "Ensembl"
    ensAll <- genes(edb)
    ens21Y <- genes(edb, filter=SeqNameFilter(c("Y", "21")))
    expect_equal(sort(as.character(unique(seqnames(ens21Y)))), c("21", "Y"))
    gr <- GRanges(seqnames="Y", ranges=IRanges(start=1, end=59373566), strand="+")
    ensY <- genes(edb, filter=GRangesFilter(gr))
    expect_equal(seqlevels(ensY), "Y")
    expect_equal(unique(as.character(strand(ensY))), "+")

    ## Check UCSC stuff
    seqlevelsStyle(edb) <- "UCSC"
    options(ensembldb.seqnameNotFound="ORIGINAL")
    ## Just visually inspect the seqinfo and seqnames for the "all" query.
    ucscAll <- genes(edb)
    as.character(unique(seqnames(ucscAll)))
    ucsc21Y <- genes(edb, filter=SeqNameFilter(c("chrY", "chr21")))
    expect_equal(sort(seqlevels(ucsc21Y)), c("chr21", "chrY"))
    expect_equal(sort(names(ens21Y)), sort(names(ucsc21Y)))
    ## GRangesFilter.
    gr <- GRanges(seqnames="chrY", ranges=IRanges(start=1, end=59373566), strand="+")
    ucscY <- genes(edb, filter=GRangesFilter(gr))
    expect_equal(seqlevels(ucscY), "chrY")
    expect_equal(unique(as.character(strand(ucscY))), "+")
    expect_equal(sort(names(ensY)), sort(names(ucscY)))
    options(ensembldb.seqnameNotFound=orig)
})

test_that("transcripts works with seqlevelsStyle", {
    orig <- getOption("ensembldb.seqnameNotFound")
    edb <- EnsDb.Hsapiens.v86
    seqlevelsStyle(edb) <- "Ensembl"
    ens21Y <- transcripts(edb, filter = SeqNameFilter(c("Y", "21")))
    expect_equal(sort(seqlevels(ens21Y)), c("21", "Y"))
    gr <- GRanges(seqnames="Y", ranges=IRanges(start=1, end=59373566), strand="+")
    ensY <- transcripts(edb, filter=GRangesFilter(gr))
    expect_equal(seqlevels(ensY), "Y")
    expect_equal(unique(as.character(strand(ensY))), "+")

    ## Check UCSC stuff
    seqlevelsStyle(edb) <- "UCSC"
    options(ensembldb.seqnameNotFound="ORIGINAL")
    ucsc21Y <- transcripts(edb, filter=SeqNameFilter(c("chrY", "chr21")))
    expect_equal(sort(seqlevels(ucsc21Y)), c("chr21", "chrY"))
    expect_equal(sort(names(ens21Y)), sort(names(ucsc21Y)))
    ## GRangesFilter.
    gr <- GRanges(seqnames="chrY", ranges=IRanges(start=1, end=59373566),
                  strand="+")
    ucscY <- transcripts(edb, filter=GRangesFilter(gr))
    expect_equal(seqlevels(ucscY), "chrY")
    expect_equal(unique(as.character(strand(ucscY))), "+")
    expect_equal(sort(names(ensY)), sort(names(ucscY)))
    options(ensembldb.seqnameNotFound=orig)
})

test_that("transcriptsBy works with seqlevelsStyle", {
    orig <- getOption("ensembldb.seqnameNotFound")
    edb <- EnsDb.Hsapiens.v86
    seqlevelsStyle(edb) <- "Ensembl"
    ens21Y <- transcriptsBy(edb, filter=SeqNameFilter(c("Y", "21")))
    expect_equal(sort(seqlevels(ens21Y)), c("21", "Y"))
    gr <- GRanges(seqnames="Y", ranges=IRanges(start=1, end=59373566), strand="+")
    ensY <- transcriptsBy(edb, filter=GRangesFilter(gr))
    expect_equal(seqlevels(ensY), "Y")
    suppressWarnings(
        expect_equal(unique(as.character(unlist(strand(ensY)))), "+")
    )
    
    ## Check UCSC stuff
    seqlevelsStyle(edb) <- "UCSC"
    options(ensembldb.seqnameNotFound="ORIGINAL")
    ucsc21Y <- transcriptsBy(edb, filter=SeqNameFilter(c("chrY", "chr21")))
    expect_equal(sort(seqlevels(ucsc21Y)), c("chr21", "chrY"))
    expect_equal(sort(names(ens21Y)), sort(names(ucsc21Y)))
    ## GRangesFilter.
    gr <- GRanges(seqnames="chrY", ranges=IRanges(start=1, end=59373566), strand="+")
    ucscY <- transcriptsBy(edb, filter=GRangesFilter(gr))
    expect_equal(seqlevels(ucscY), "chrY")
    suppressWarnings(
        expect_equal(unique(as.character(unlist(strand(ucscY)))), "+")
    )
    expect_equal(sort(names(ensY)), sort(names(ucscY)))
    options(ensembldb.seqnameNotFound=orig)
})

test_that("exons works with seqlevelsStyle", {
    orig <- getOption("ensembldb.seqnameNotFound")
    edb <- EnsDb.Hsapiens.v86
    seqlevelsStyle(edb) <- "Ensembl"
    ens21Y <- exons(edb, filter=SeqNameFilter(c("Y", "21")))
    expect_equal(sort(seqlevels(ens21Y)), c("21", "Y"))
    gr <- GRanges(seqnames="Y", ranges=IRanges(start=1, end=59373566), strand="+")
    ensY <- exons(edb, filter=GRangesFilter(gr))
    expect_equal(seqlevels(ensY), "Y")
    expect_equal(unique(as.character(strand(ensY))), "+")

    ## Check UCSC stuff
    seqlevelsStyle(edb) <- "UCSC"
    options(ensembldb.seqnameNotFound="ORIGINAL")
    ucsc21Y <- exons(edb, filter=SeqNameFilter(c("chrY", "chr21")))
    expect_equal(sort(seqlevels(ucsc21Y)), c("chr21", "chrY"))
    expect_equal(sort(names(ens21Y)), sort(names(ucsc21Y)))
    ## GRangesFilter.
    gr <- GRanges(seqnames="chrY", ranges=IRanges(start=1, end=59373566), strand="+")
    ucscY <- exons(edb, filter=GRangesFilter(gr))
    expect_equal(seqlevels(ucscY), "chrY")
    expect_equal(unique(as.character(strand(ucscY))), "+")
    expect_equal(sort(names(ensY)), sort(names(ucscY)))
    options(ensembldb.seqnameNotFound=orig)
})

test_that("exonsBy works with seqlevelsStyle", {
    orig <- getOption("ensembldb.seqnameNotFound")
    edb <- EnsDb.Hsapiens.v86
    seqlevelsStyle(edb) <- "Ensembl"
    ens21Y <- exonsBy(edb, filter=SeqNameFilter(c("Y", "21")))
    expect_equal(sort(seqlevels(ens21Y)), c("21", "Y"))
    gr <- GRanges(seqnames="Y", ranges=IRanges(start=1, end=59373566), strand="+")
    ensY <- exonsBy(edb, filter=GRangesFilter(gr))
    expect_equal(seqlevels(ensY), "Y")
    suppressWarnings(
        expect_equal(unique(as.character(unlist(strand(ensY)))), "+")
    )
    ## Check UCSC stuff
    seqlevelsStyle(edb) <- "UCSC"
    options(ensembldb.seqnameNotFound="ORIGINAL")
    ucsc21Y <- exonsBy(edb, filter=SeqNameFilter(c("chrY", "chr21")))
    expect_equal(sort(seqlevels(ucsc21Y)), c("chr21", "chrY"))
    expect_equal(sort(names(ens21Y)), sort(names(ucsc21Y)))
    ## GRangesFilter.
    gr <- GRanges(seqnames="chrY", ranges=IRanges(start=1, end=59373566), strand="+")
    ucscY <- exonsBy(edb, filter=GRangesFilter(gr))
    expect_equal(seqlevels(ucscY), "chrY")
    suppressWarnings(
        expect_equal(unique(as.character(unlist(strand(ucscY)))), "+")
    )
    expect_equal(sort(names(ensY)), sort(names(ucscY)))
    options(ensembldb.seqnameNotFound=orig)
})

test_that("cdsBy works with seqlevelsStyle", {
    orig <- getOption("ensembldb.seqnameNotFound")
    edb <- EnsDb.Hsapiens.v86
    seqlevelsStyle(edb) <- "Ensembl"
    ens21Y <- cdsBy(edb, filter=SeqNameFilter(c("Y", "21")))
    expect_equal(sort(seqlevels(ens21Y)), c("21", "Y"))
    gr <- GRanges(seqnames="Y", ranges=IRanges(start=1, end=59373566), strand="+")
    ensY <- cdsBy(edb, filter=GRangesFilter(gr))
    expect_equal(seqlevels(ensY), "Y")
    suppressWarnings(
        expect_equal(unique(as.character(unlist(strand(ensY)))), "+")
    )
    ## Check UCSC stuff
    seqlevelsStyle(edb) <- "UCSC"
    options(ensembldb.seqnameNotFound="ORIGINAL")
    ucsc21Y <- cdsBy(edb, filter=SeqNameFilter(c("chrY", "chr21")))
    expect_equal(sort(seqlevels(ucsc21Y)), c("chr21", "chrY"))
    expect_equal(sort(names(ens21Y)), sort(names(ucsc21Y)))
    ## GRangesFilter.
    gr <- GRanges(seqnames="chrY", ranges=IRanges(start=1, end=59373566), strand="+")
    ucscY <- cdsBy(edb, filter=GRangesFilter(gr))
    expect_equal(seqlevels(ucscY), "chrY")
    suppressWarnings(
        expect_equal(unique(as.character(unlist(strand(ucscY)))), "+")
    )
    expect_equal(sort(names(ensY)), sort(names(ucscY)))
    options(ensembldb.seqnameNotFound=orig)
})

test_that("threeUTRsByTranscript works with seqlevelsStyle", {
    orig <- getOption("ensembldb.seqnameNotFound")
    edb <- EnsDb.Hsapiens.v86
    seqlevelsStyle(edb) <- "Ensembl"
    ens21Y <- threeUTRsByTranscript(edb, filter=SeqNameFilter(c("Y", "21")))
    expect_equal(sort(seqlevels(ens21Y)), c("21", "Y"))
    gr <- GRanges(seqnames="Y", ranges=IRanges(start=1, end=59373566), strand="+")
    ensY <- threeUTRsByTranscript(edb, filter=GRangesFilter(gr))
    expect_equal(seqlevels(ensY), "Y")
    suppressWarnings(
        expect_equal(unique(as.character(unlist(strand(ensY)))), "+")
    )
    ## Check UCSC stuff
    seqlevelsStyle(edb) <- "UCSC"
    options(ensembldb.seqnameNotFound="ORIGINAL")
    ucsc21Y <- threeUTRsByTranscript(edb, filter=SeqNameFilter(c("chrY", "chr21")))
    expect_equal(sort(seqlevels(ucsc21Y)), c("chr21", "chrY"))
    expect_equal(sort(names(ens21Y)), sort(names(ucsc21Y)))
    ## GRangesFilter.
    gr <- GRanges(seqnames="chrY", ranges=IRanges(start=1, end=59373566), strand="+")
    ucscY <- threeUTRsByTranscript(edb, filter=GRangesFilter(gr))
    expect_equal(seqlevels(ucscY), "chrY")
    suppressWarnings(
        expect_equal(unique(as.character(unlist(strand(ucscY)))), "+")
    )
    expect_equal(sort(names(ensY)), sort(names(ucscY)))
    options(ensembldb.seqnameNotFound=orig)
})

test_that("fiveUTRsByTranscript works with seqlevelsStyle", {
    orig <- getOption("ensembldb.seqnameNotFound")
    edb <- EnsDb.Hsapiens.v86
    seqlevelsStyle(edb) <- "Ensembl"
    ens21Y <- fiveUTRsByTranscript(edb, filter=SeqNameFilter(c("Y", "21")))
    expect_equal(sort(seqlevels(ens21Y)), c("21", "Y"))
    gr <- GRanges(seqnames="Y", ranges=IRanges(start=1, end=59373566), strand="+")
    ensY <- fiveUTRsByTranscript(edb, filter=GRangesFilter(gr))
    expect_equal(seqlevels(ensY), "Y")
    suppressWarnings(
        expect_equal(unique(as.character(unlist(strand(ensY)))), "+")
    )
    ## Check UCSC stuff
    seqlevelsStyle(edb) <- "UCSC"
    options(ensembldb.seqnameNotFound="ORIGINAL")
    ucsc21Y <- fiveUTRsByTranscript(edb, filter=SeqNameFilter(c("chrY", "chr21")))
    expect_equal(sort(seqlevels(ucsc21Y)), c("chr21", "chrY"))
    expect_equal(sort(names(ens21Y)), sort(names(ucsc21Y)))
    ## GRangesFilter.
    gr <- GRanges(seqnames="chrY", ranges=IRanges(start=1, end=59373566), strand="+")
    ucscY <- fiveUTRsByTranscript(edb, filter=GRangesFilter(gr))
    expect_equal(seqlevels(ucscY), "chrY")
    suppressWarnings(
        expect_equal(unique(as.character(unlist(strand(ucscY)))), "+")
    )
    expect_equal(sort(names(ensY)), sort(names(ucscY)))
    options(ensembldb.seqnameNotFound=orig)
})

test_that("seting and getting seqlevelsStyle works", {
    edb <- EnsDb.Hsapiens.v86
    ## Testing the getter/setter for the seqlevelsStyle.
    expect_equal(seqlevelsStyle(edb), "Ensembl")
    expect_equal(NA, ensembldb:::getProperty(edb, "seqlevelsStyle"))

    seqlevelsStyle(edb) <- "Ensembl"
    expect_equal(seqlevelsStyle(edb), "Ensembl")
    expect_equal("Ensembl", ensembldb:::getProperty(edb, "seqlevelsStyle"))

    ## Try NCBI.
    seqlevelsStyle(edb) <- "NCBI"
    expect_equal(seqlevelsStyle(edb), "NCBI")

    ## Try UCSC.
    seqlevelsStyle(edb) <- "UCSC"
    expect_equal(seqlevelsStyle(edb), "UCSC")

    ## Error checking:
    expect_error(seqlevelsStyle(edb) <- "bla")
})

test_that("formatting seqnames for query works with seqlevelsStyle", {
    ## Testing if the formating/mapping between seqnames works as expected
    ## We want to map anything TO Ensembl.
    ## Check also the warning messages!
    ucscs <- c("chr1", "chr3", "chr1", "chr9", "chrM", "chr1", "chrX")
    enses <- c("1", "3", "1", "9", "MT", "1", "X")
    ## reset
    edb <- EnsDb.Hsapiens.v86
    ## Shouldn't do anything here.
    seqlevelsStyle(edb)
    ensembldb:::dbSeqlevelsStyle(edb)
    got <- ensembldb:::formatSeqnamesForQuery(edb, enses)
    expect_equal(got, enses)
    ## Change the seqlevels to UCSC
    seqlevelsStyle(edb) <- "UCSC"
    ## If ifNotFound is not specified we suppose to get an error.
    options(ensembldb.seqnameNotFound="MISSING")
    expect_error(ensembldb:::formatSeqnamesForQuery(edb, enses))
    ## With specifying ifNotFound
    suppressWarnings(
        got <- ensembldb:::formatSeqnamesForQuery(edb, enses, ifNotFound=NA)
    )
    expect_equal(all(is.na(got)), TRUE)
    ## Same by setting the option
    options(ensembldb.seqnameNotFound=NA)
    suppressWarnings(
        got <- ensembldb:::formatSeqnamesForQuery(edb, enses)
    )
    expect_equal(all(is.na(got)), TRUE)

    ## Now the working example:
    got <- ensembldb:::formatSeqnamesForQuery(edb, ucscs)
    expect_equal(got, enses)
    ## What if one is not mappable:
    suppressWarnings(
        got <- ensembldb:::formatSeqnamesForQuery(edb, c(ucscs, "asdfd"),
                                                  ifNotFound=NA)
    )
    expect_equal(got, c(enses, NA))
})

test_that("formating seqnames from query works with seqlevelsStyle", {
    ucscs <- c("chr1", "chr3", "chr1", "chr9", "chrM", "chr1", "chrX")
    enses <- c("1", "3", "1", "9", "MT", "1", "X")
    edb <- EnsDb.Hsapiens.v86
    ## Shouldn't do anything here.
    seqlevelsStyle(edb)
    ensembldb:::dbSeqlevelsStyle(edb)
    got <- ensembldb:::formatSeqnamesFromQuery(edb, enses)
    expect_equal(got, enses)
    ## Change the seqlevels to UCSC
    seqlevelsStyle(edb) <- "UCSC"
    ## If ifNotFound is not specified we suppose to get an error.
    options(ensembldb.seqnameNotFound="MISSING")
    expect_error(ensembldb:::formatSeqnamesFromQuery(edb, ucsc))
    ## With specifying ifNotFound
    suppressWarnings(
        got <- ensembldb:::formatSeqnamesFromQuery(edb, ucscs, ifNotFound=NA)
    )
    expect_equal(all(is.na(got)), TRUE)
    ## Same using options
    options(ensembldb.seqnameNotFound=NA)
    suppressWarnings(
        got <- ensembldb:::formatSeqnamesFromQuery(edb, ucscs, ifNotFound=NA)
    )
    expect_equal(all(is.na(got)), TRUE)
    ## Now the working example:
    got <- ensembldb:::formatSeqnamesFromQuery(edb, enses)
    expect_equal(got, ucscs)
    ## What if one is not mappable:
    suppressWarnings(
        got <- ensembldb:::formatSeqnamesFromQuery(edb, c(enses, "asdfd"),
                                                   ifNotFound=NA)
    )
    expect_equal(got, c(ucscs, NA))
    suppressWarnings(
        got <- ensembldb:::formatSeqnamesFromQuery(edb, c(enses, "asdfd"))
    )
    expect_equal(got, c(ucscs, NA))
})

test_that("prefixChromName works", {
    res <- ensembldb:::ucscToEns("chrY")
    expect_equal(res, "Y")
    res <- ensembldb:::prefixChromName("Y")
    expect_equal(res, "Y")
    useU <- getOption("ucscChromosomeNames", default = FALSE)
    options(ucscChromosomeNames = TRUE)
    res <- ensembldb:::prefixChromName("Y")
    expect_equal(res, "chrY")
    options(ucscChromosomeNames = useU)
})
