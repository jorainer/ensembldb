test_that(".genome_to_tx and genomeToTranscript works", {
    library(testthat)
    library(EnsDb.Hsapiens.v75)
    edbx <- filter(EnsDb.Hsapiens.v75, filter = ~ seq_name == "X")

    ## 1) ENST00000381578, + strand last four nt of CDS.
    ##    Position in ENST00000381578: 259 + 709 + 209 + 58 + 89 + 245 = 1569
    gnm <- GRanges("X", IRanges(start = 605370, end = 605374))
    res <- .genome_to_tx(gnm, db = edbx)
    expect_true(length(res) == 2)
    expect_true(is(res, "IRanges"))
    expect_equal(start(res["ENST00000381578"]), 1569)
    expect_equal(width(res["ENST00000381578"]), 5)
    ##    Restrict to negative strand - nothing there.
    strand(gnm) <- "-"
    res <- .genome_to_tx(gnm, db = edbx)
    expect_equal(length(res), 1)
    expect_equal(start(res), -1)
    expect_true(is(res, "IRanges"))

    ## 2) ENST00000486554, - strand, nt 1-3 of tx
    gnm <- GRanges("X:106959629-106959631")
    res <- ensembldb:::.genome_to_tx(gnm, db = edbx)
    expect_equal(names(res), c("ENST00000372390", "ENST00000486554"))
    expect_equal(start(res["ENST00000486554"]), 1)
    expect_equal(width(res["ENST00000486554"]), 3)

    ## 2) ENST00000486554, - strand, nt 1-3 of CDS, nt 501-503 of tx
    gnm <- GRanges("X:106959129-106959131")
    res <- ensembldb:::.genome_to_tx(gnm, db = edbx)
    expect_equal(length(res), 11)
    expect_equal(start(res["ENST00000486554"]), 501)
    expect_equal(width(res["ENST00000486554"]), 3)
    
    ## 3) As above, just not all within the exon.
    gnm <- GRanges("X:106959629-106959635")
    res <- .genome_to_tx(gnm, db = edbx)
    expect_true(length(res) == 1)
    expect_equal(names(res), "ENST00000372390")

    ## 4) ENST00000486554, - strand, nts 5-8 in exon 2
    ##    exon 1: 503nt: region is expected at 508-511.
    gnm <- GRanges("X", IRanges(end = 106957975, width = 4))
    res <- .genome_to_tx(gnm, db = edbx)
    expect_true(length(res) == 9)       # overlapping many tx...
    expect_equal(start(res["ENST00000486554"]), 508)
    expect_equal(end(res["ENST00000486554"]), 511)

    ## 5) ENST00000372397, - strand, last 3nt
    ##    exon 3: 446 + 52 + 1529 total length: 2025
    gnm <- GRanges("X", IRanges(start = 106956451, width = 3))
    res <- .genome_to_tx(gnm, db = edbx)
    expect_true(length(res) == 3)
    expect_equal(start(res["ENST00000372397"]), 2025)
    expect_equal(end(res["ENST00000372397"]), 2027)

    ## 6) A region with a tx/exon on both strands.
    gnm <- GRanges("7", IRanges(start = 127231550, width = 5))
    edb7 <- filter(EnsDb.Hsapiens.v75, filter = ~ seq_name == "7")
    res <- .genome_to_tx(gnm, edb7)
    tx <- transcripts(edb7, filter = TxIdFilter(names(res)))
    expect_equal(as.character(strand(tx)), c("+", "+", "+", "-"))
    expect_true(all(names(tx) %in% c("ENST00000473728", "ENST00000463733",
                                     "ENST00000000233", "ENST00000478328")))

    ## 7) Check with input being a GRanges of length > 1
    gnm <- GRanges("X", IRanges(start = c(605370, 106959629),
                                end = c(605374, 106959631) ))
    res <- .genome_to_tx(gnm, edbx)
    expect_true(is(res, "IRangesList"))
    expect_true(length(res) == 2)
    expect_equal(start(res[[1]]["ENST00000381578"]), 1569)
    expect_equal(width(res[[1]]["ENST00000381578"]), 5)
    expect_equal(names(res[[2]]), c("ENST00000372390", "ENST00000486554"))
    expect_equal(start(res[[2]]["ENST00000486554"]), 1)
    expect_equal(width(res[[2]]["ENST00000486554"]), 3)

    ## Check errors etc.
    expect_error(genomeToTranscript())
    expect_error(genomeToTranscript(x = 4, db = edbx))
    expect_error(genomeToTranscript(x = gnm))
    expect_error(genomeToTranscript(x = gnm, db = 5))

    gnm <- GRanges("X", IRanges(start = 605370, end = 605374))
    res <- genomeToTranscript(gnm, edbx)
    expect_equal(length(gnm), length(res))
    expect_true(is(res, "IRangesList"))

    gnm <- GRanges("X", IRanges(start = c(605370, 106959629, 50),
                                end = c(605374, 106959631, 50) ))
    expect_warning(res <- genomeToTranscript(gnm, edbx))
    expect_equal(length(gnm), length(res))
    expect_true(is(res, "IRangesList"))
    expect_true(start(res[[3]]) < 0)
})

test_that("genomeToProtein works", {
    edbx <- filter(EnsDb.Hsapiens.v75, filter = ~ seq_name == "X")
    ## 1) ENST00000381578, + strand first 5nt of the CDS:
    ## 591633: 5'UTR: 259 + 432 (691); 591,201 + 432 = 591633 is first nt of CDS
    ## 605371: last nt of CDS: 605126 + 246 - 1
    ## 605368: last nt prior to the stop codon
    ## 595564: range that is within the intron
    gnm <- GRanges("X", IRanges(start = c(591633, 605371, 605368, 595564),
                                width = c(5, 1, 1, 3)))
    expect_warning(res <- genomeToProtein(gnm, edbx))
    expect_true(is(res, "IRangesList"))
    expect_equal(length(res), 4)
    ## Elements 2 and 4 can not be mapped
    expect_equal(start(res[[2]]), c(-1, -1))
    expect_equal(start(res[[4]]), c(-1))
    ## Check that order is correct
    expect_true(mcols(res[[1]])$seq_start[1] == 591633)
    expect_true(mcols(res[[2]])$seq_start[1] == 605371)
    expect_true(mcols(res[[3]])$seq_start[1] == 605368)
    expect_true(mcols(res[[4]])$seq_start[1] == 595564)
    ## Check coords within the protein:
    expect_true(start(res[[1]])[1] == 1)
    expect_true(end(res[[1]])[1] == 2)
    prt <- proteins(edbx, filter = ProteinIdFilter(names(res[[3]])[1]))
    expect_equal(start(res[[3]])[1], nchar(prt$protein_sequence))

    ## 2) ENST00000486554, - strand, nt 1-3
    gnm <- GRanges("X:106959629-106959631")
    expect_warning(res <- genomeToProtein(gnm, db = edbx))
    expect_true(is(res, "IRangesList"))
    expect_true(length(res) == length(gnm))

    gnm <- GRanges("X:106959629-106959931")
    expect_warning(res <- genomeToProtein(gnm, db = edbx))
    expect_true(is(res, "IRangesList"))
    expect_true(length(res) == length(gnm))
    expect_true(start(res[[1]]) < 0)

    ## ENST00000486554
    ## 106959629:3: first 3 nt of transcript
    ## 106959629:700: spans exon 1 and part of intron 1
    ## 106959130:2: first two nt of CDS.
    ## 106957979:1 first nt in exon 2, second amino acid.
    gnm <- GRanges("X", IRanges(start = c(106959629, 106959629, 106959130,
                                          106957979),
                                width = c(3, 700, 2, 1)))
    res <- genomeToProtein(gnm, db = edbx)
    expect_true(is(res, "IRangesList"))
    expect_true(all(start(res[[1]]) < 0))
    expect_true(all(start(res[[2]]) < 0))
    expect_equal(start(res[[3]]["ENSP00000425414"]), 1)
    expect_equal(width(res[[3]]["ENSP00000425414"]), 1)
    expect_equal(start(res[[4]]["ENSP00000425414"]), 2)
    expect_equal(width(res[[4]]["ENSP00000425414"]), 1)

    gnm <- GRanges("X", IRanges(start = 106959130, width = 2))
    expect_warning(res <- genomeToProtein(gnm, edbx))
    expect_true(is(res, "IRangesList"))
    expect_true(length(res) == length(gnm))
    expect_equal(start(res[[1]]["ENSP00000425414"]), 1)
    expect_false(mcols(res[[1]]["ENSP00000425414"])$cds_ok)
})
