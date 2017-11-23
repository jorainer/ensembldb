test_that(".genome_to_tx and genomeToTranscript works", {
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

    ## 2) ENST00000486554, - strand, nt 1-3
    gnm <- GRanges("X:106959629-106959631")
    res <- .genome_to_tx(gnm, db = edbx)
    expect_equal(names(res), c("ENST00000372390", "ENST00000486554"))
    expect_equal(start(res["ENST00000486554"]), 1)
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
})
