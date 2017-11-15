test_that(".proteinCoordsToTx works", {
    prts <- IRanges(start = c(1, 2), end = c(1, 2))
    res <- .proteinCoordsToTx(prts)
    expect_equal(start(res), c(1, 4))
    expect_equal(end(res), c(3, 6))

    res <- .proteinCoordsToTx(IRanges(start = 3, end = 5))
    expect_equal(start(res), 7)
    expect_equal(end(res), 15)
})

test_that(".filter_for_idType works", {
    expect_equal(.filter_for_idType("a", idType = "something"), NULL)
    expect_true(is(.filter_for_idType("a", idType = "protein_id"),
                   "ProteinIdFilter"))
    expect_true(is(.filter_for_idType("a", idType = "uniprot_id"),
                   "UniprotFilter"))
    expect_true(is(.filter_for_idType("a", idType = "tx_id"),
                   "TxIdFilter"))
})


test_that(".cds_for_id and .cds_matching_protein work", {
    prts <- c("ENSP00000466417", "don't exist", "ENSP00000437716")
    txs <- c("ENST00000589955", "ENST00000544220")

    expect_warning(res <- .cds_for_id(edb, prts))
    expect_equal(names(res), prts)
    expect_equal(names(res[[1]]), txs[1])
    expect_equal(names(res[[3]]), txs[2])
    expect_true(is.null(res[[2]]))
    ## Use uniprot...
    ## want a uniprot that is mapped to two tx.
    unis <- c("E7EMD7_HUMAN", "ALAT2_HUMAN", "H0YIP2_HUMAN")
    unis_counts <- c(7, 2, 1)
    res2 <- .cds_for_id(edb, unis, idType = "uniprot_id")
    expect_equal(unname(lengths(res2)), unis_counts)
    expect_equal(names(res2), unis)

    ## Run .cds_matching_protein on the results to avoid re-querying.
    ## H0YIP2_HUMAN has incomplete 5' and 3' CDS.
    expect_warning(clnd <- .cds_matching_protein(edb, res2))
    expect_true(all(lengths(clnd) == 1))

    ## ENSP00000437716 encoding transcript has an incomplete 3' CDS
    expect_warning(clnd <- .cds_matching_protein(edb, res))
    expect_equal(unname(lengths(clnd)), c(1, 0, 1))
})

test_that(".cds_for_id_range works", {
    rng <- IRanges(start = 2, end = 207)
    mcols(rng)$id <- "ENSP00000466417"
    res <- .cds_for_id_range(edb, rng, id = "id")
    expect_equal(length(res[[1]]), 0)
})

test_that("proteinToGenome works", {
    rng <- IRanges(start = 1, end = 1)
    names(rng) <- "H0YIP2_HUMAN"

    protein <- rng
    id <- "name"
    idType <- "uniprot_id"
    genome <- edb
})

test_that(".splice works", {
    ir <- IRanges(start = c(5, 10, 15), end = c(7, 13, 16))
    res <- .splice(ir)
    expect_equal(width(ir), width(res))
    expect_equal(start(res), c(1, 4, 8))

    res <- .splice(IRangesList(ir, ir))
})

test_that(".to_genome works", {
    ## Forward strand
    g_coords <- GRanges("2", IRanges(start = c(3, 8, 15, 19),
                                     end = c(5, 12, 16, 21)),
                        strand = "+")
    c_coords <- IRanges(start = 5, end = 7)
    res <- .to_genome(g_coords, c_coords)
    expect_equal(seqlevels(res), "2")
    expect_equal(start(res), 9)
    expect_equal(end(res), 11)
    ## 2nd example
    c_coords <- IRanges(start = 2, end = 13)
    res <- .to_genome(g_coords, c_coords)
    expect_equal(start(res), c(4, 8, 15, 19))
    expect_equal(end(res), c(5, 12, 16, 21))
    ## 3rd example
    c_coords <- IRanges(start = 5, end = 12)
    res <- .to_genome(g_coords, c_coords)
    expect_equal(start(res), c(9, 15, 19))
    expect_equal(end(res), c(12, 16, 20))
    
    ## Reverse strand
    g_coords <- GRanges("1", IRanges(start = c(28, 21, 16, 10, 3),
                                     end = c(30, 25, 17, 12, 6)),
                        strand = "-")
    c_coords <- IRanges(start = 2, end = 2)
    res <- .to_genome(g_coords, c_coords)
    expect_true(is(res, "GRanges"))
    expect_equal(start(res), 29)
    expect_equal(end(res), 29)
    expect_equal(seqlevels(res), "1")
    expect_equal(as.character(strand(res)), "-")
    ## 2nd example
    c_coords <- IRanges(start = 8, end = 16)
    res <- .to_genome(g_coords, c_coords)
    expect_equal(start(res), c(4, 10, 16, 21))
    expect_equal(end(res), c(6, 12, 17, 21))
    expect_equal(seqlevels(res), "1")
    ## passing a GRangesList
    g_coords <- GRangesList(g_coords, g_coords[2:5])
    c_coords <- IRanges(start = 4, end = 9)
    res <- .to_genome(g_coords, c_coords)
    expect_true(is.list(res))
    expect_equal(start(res[[1]]), c(17, 21))
    expect_equal(end(res[[1]]), c(17, 25))
    expect_equal(start(res[[2]]), c(11, 16, 21))
    expect_equal(end(res[[2]]), c(12, 17, 22))
    
    ## Check errors.
    expect_error(.to_genome(c_coords, c_coords))
    g_coords <- GRanges("2", IRanges(start = c(3, 8, 15, 19),
                                     end = c(5, 12, 16, 21)),
                        strand = "+")
    c_coords <- IRanges(start = 5, end = 40)
    expect_error(.to_genome(g_coords, c_coords))
    c_coords <- IRanges(start = 40, end = 50)
    expect_error(.to_genome(g_coords, c_coords))
})
