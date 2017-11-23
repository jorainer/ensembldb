test_that("transcriptToProtein works",  {
    edbx <- filter(EnsDb.Hsapiens.v75, filter = ~ seq_name == "X")

    ## SHOX2, ENST00000381578
    ## exon 1: 259
    ## exon 2: 709 UTR: 432nt long
    ## exon 3: 209
    ## exon 4: 58
    ## exon 5: 89
    ## exon 6: 2,433 CDS: 246, then UTR.
    x <- IRanges(start = 2, width = 4)
    expect_error(.tx_to_protein(x, edbx))
    expect_error(.txs_to_proteins(x, edbx))
    
    ## Non existant transcript.
    x <- IRanges(start = 2, width = 5, names = "some")
    expect_warning(res_a <- .tx_to_protein(x, edbx))
    expect_true(is(res_a, "IRanges"))
    expect_equal(colnames(mcols(res_a)), c("tx_id", "tx_start", "tx_end",
                                           "cds_ok"))
    expect_true(start(res_a) < 0)
    expect_warning(res_b <- .txs_to_proteins(x, edbx))
    expect_equal(res_a, res_b)
    
    ## ENST00000431238: non-coding transcript
    x <- IRanges(start = 2, width = 5, names = "ENST00000431238")
    expect_warning(res_a <- .tx_to_protein(x, edbx))
    expect_true(is(res_a, "IRanges"))
    expect_equal(colnames(mcols(res_a)), c("tx_id", "tx_start", "tx_end",
                                           "cds_ok"))
    expect_true(start(res_a) < 0)
    expect_warning(res_b <- .txs_to_proteins(x, edbx))
    expect_equal(res_a, res_b)

    ## just outside the CDS
    x <- IRanges(start = 691, width = 1, names = "ENST00000381578")
    expect_warning(res_a <- .tx_to_protein(x, edbx))
    expect_true(is(res_a, "IRanges"))
    expect_equal(colnames(mcols(res_a)), c("tx_id", "tx_start", "tx_end",
                                           "cds_ok"))
    expect_true(start(res_a) < 0)
    expect_warning(res_b <- .txs_to_proteins(x, edbx))
    expect_equal(res_a, res_b)

    ## first 2nt of the CDS
    x <- IRanges(start = 692, width = 2, names = "ENST00000381578")
    res_a <- .tx_to_protein(x, edbx)
    expect_equal(start(res_a), 1)
    expect_equal(end(res_a), 1)
    expect_true(mcols(res_a)$cds_ok)
    res_b <- .txs_to_proteins(x, edbx)
    expect_equal(res_a, res_b)

    ## nt 3-4 of the CDS, should map to 1-2 of the prot seq.
    x <- IRanges(start = 694, width = 2, names = "ENST00000381578")
    res_a <- .tx_to_protein(x, edbx)
    expect_equal(start(res_a), 1)
    expect_equal(end(res_a), 2)

    ## nts of the stop codon.
    x <- IRanges(start = 1569, width = 4, names = "ENST00000381578")
    expect_warning(res_a <- .tx_to_protein(x, edbx))

    ## last nt before the stop codon
    x <- IRanges(start = 1567, width = 1, names = "ENST00000381578")
    res_a <- .tx_to_protein(x, edbx)
    prt <- proteins(edbx, filter = ProteinIdFilter(names(res_a)))
    expect_equal(start(res_a), nchar(prt$protein_sequence))

    x <- IRanges(start = 692, width = 2)
    mcols(x)$id <- "ENST00000381578"
    res_a <- .tx_to_protein(x, edbx, id = "id")
    res_b <- .txs_to_proteins(x, edbx, id = "id")
    expect_equal(res_a, res_b)
    
    ## Multiple input.
    ## just outside the CDS (last nt of the 5' UTR).
    x <- IRanges(start = c(691, 692, 32, 5000, 1, 1565),
                 width = c(1, 1, 3, 3, 2, 2),
                 names = c("ENST00000381578", "ENST00000381578", "some",
                           "ENST00000381578", "ENST00000431238",
                           "ENST00000381578"))
    expect_warning(res_a <- .tx_to_protein(x, edbx))
    expect_warning(res_b <- .txs_to_proteins(x, edbx))
    expect_equal(res_a, res_b)    
})

test_that("transcriptToProtein works", {
    edbx <- filter(EnsDb.Hsapiens.v75, filter = ~ seq_name == "X")
    
    ## Errors.
    expect_error(transcriptToProtein())
    expect_error(transcriptToProtein(db = edbx))
    expect_error(transcriptToProtein(txpos))
    
    ## 
    x <- IRanges(start = 692, width = 2)
    mcols(x)$id <- "ENST00000381578"
    res <- transcriptToProtein(x, edbx, id = "id")
    expect_equal(start(res), 1)
    expect_equal(end(res), 1)
})
