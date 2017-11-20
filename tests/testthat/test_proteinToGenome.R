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
    rngs <- IRanges(start = c(4, 1, 9), end = c(4, 1, 11))
    names(rngs) <- prts
    ## seq_names 11, 18
    dbsub <- filter(edb, filter = ~ seq_name %in% c("11", "18"))
    expect_warning(res <- .cds_for_id(dbsub, prts))
    expect_equal(names(res), prts)
    expect_equal(names(res[[1]]), txs[1])
    expect_equal(names(res[[3]]), txs[2])
    expect_true(is.null(res[[2]]))

    ## Use uniprot...
    ## want a uniprot that is mapped to two tx.
    unis <- c("E7EMD7_HUMAN", "ALAT2_HUMAN", "H0YIP2_HUMAN")
    unis_counts <- c(7, 2, 1)
    dbsub2 <- filter(edb, filter = ~ seq_name %in% c("11", "12", "16"))
    res2 <- .cds_for_id(dbsub2, unis, idType = "uniprot_id")
    expect_equal(unname(lengths(res2)), unis_counts)
    expect_equal(names(res2), unis)

    rngs <- IRanges(start = c(11, 6, 4), end = c(19, 6, 9))
    names(rngs) <- unis
    ## Can I use that right away?
    z <- res2[[1]]
    irng <- IRanges(start = 5, end = 8)
    tmp <- .to_genome(z, irng)
    
    ## Run .cds_matching_protein on the results to avoid re-querying.
    ## H0YIP2_HUMAN has incomplete 5' and 3' CDS.
    expect_warning(clnd <- .cds_matching_protein(dbsub2, res2))
    expect_equal(unname(lengths(clnd)), c(7, 2, 1))
    expect_true(is.list(clnd))
    expect_true(is(clnd[[1]], "GRangesList"))
    expect_true(all(clnd[[1]][[1]]$cds_ok))
    expect_true(all(clnd[[2]][[1]]$cds_ok))
    expect_true(all(clnd[[3]][[1]]$cds_ok == FALSE))
    
    ## ENSP00000437716 encoding transcript has an incomplete 3' CDS
    expect_warning(clnd <- .cds_matching_protein(dbsub, res))
    expect_equal(unname(lengths(clnd)), c(1, 0, 1))
    expect_true(all(clnd[[1]][[1]]$cds_ok))
    expect_true(all(clnd[[3]][[1]]$cds_ok == FALSE))
})

test_that("proteinToGenome works", {
    ## Restrict to chromosome X - speeds up stuff
    edbx <- filter(edb, filter = ~ seq_name == "X")
    ## protein_id
    ids <- c("ENSP00000425155", "ENSP00000418169", "ENSP00000441743",
             "ENSP00000385415")
    prngs <- IRanges(start = c(23, 23, 141, 23), end = c(24, 24, 421, 111))
    names(prngs) <- ids
    
    expect_warning(res <- proteinToGenome(prngs, edbx))
    ## Result has to be a triplet
    expect_true(all(sapply(res, function(z) sum(width(z))) %% 3 == 0))
    ## Manually check for the ranges here.
    ## ENSP00000418169, SYP, encoded by ENST00000479808, - strand 
    ## coords cds: 23 -> 67, 24 -> 72
    ## exon 1 width 49, 5' UTR: 13nt long. 36 nt
    ## exon 2 width 66. relative pos is 31nt in exon 2
    ## exon 2 starts: 49,055,492 [nt1], [nt31] is 49055492 - 30 = 49055462
    expect_equal(end(res[[2]]), 49055462)
    expect_equal(start(res[[2]]), 49055462 - 5)
    expect_equal(res[[2]]$protein_id, "ENSP00000418169")
    expect_equal(res[[2]]$tx_id, "ENST00000479808")
    expect_equal(res[[2]]$exon_rank, 2)
    expect_true(res[[2]]$cds_ok)

    ## ENSP00000385415, GAGE10, ENST00000407599, + strand
    ## coords cds: 23 -> 67, 111 -> 333
    ## exon 1 non-coding.
    ## exon 2 width 89nt, 5' UTR: 8nt. 67 - 1 + 8nt -> 49161331 + 74 = 49161405,
    ##        89 - 8 - (67-1) = 15nt of the 333 - 67 + 1 = 267, 252 left.
    ## exon 3 width 121nt, 131 left.
    ## exon 4 width 126, 5 left
    ## exon 5 width 117, start 49176207: cds end maps at 49176211
    strts <- c(49161405, 49161883, 49173642, 49176207)
    ends <- c(49161419, 49162003, 49173767, 49176211)
    expect_equal(start(res[[4]]), strts)
    expect_equal(end(res[[4]]), ends)
    expect_equal(res[[4]]$exon_rank, c(2, 3, 4, 5))
    expect_true(all(res[[4]]$cds_ok))
        
    ## Uniprot identifier
    ids <- c("D6RDZ7_HUMAN", "SHOX_HUMAN", "TMM27_HUMAN", "unexistant")
    prngs <- IRanges(start = c(1, 13, 43, 100), end = c(2, 21, 80, 100))
    names(prngs) <- ids
    
    res <- proteinToGenome(prngs, edbx, idType = "uniprot_id")
    ## Now, expect elements 1 and 2 to be a GRangesList and not a GRanges.
    expect_true(is(res[[1]], "GRangesList"))
    expect_true(is(res[[2]], "GRangesList"))
    expect_true(is(res[[3]], "GRanges"))
    expect_true(is(res[[4]], "GRanges"))
    expect_true(length(res[[4]]) == 0)
    res <- res[1:3]
    expect_true(all(unlist(lapply(res, function(z) sum(width(z)))) %% 3 == 0))
    ## We've got multi-mapping here
    expect_equal(unname(lengths(res)), c(4, 4, 2))
    ## Although we have different proteins, all coords are the same
    strts <- unlist(start(res[[1]]))
    expect_true(all(strts == strts[1]))
    nds <- unlist(end(res[[1]]))
    expect_true(all(nds == nds[1]))
    strts <- unlist(start(res[[2]]))
    expect_true(all(strts == strts[1]))
    nds <- unlist(end(res[[2]]))
    expect_true(all(nds == nds[1]))
})

test_that(".cds_for_id_range works", {
    rng <- IRanges(start = 2, end = 207)
    mcols(rng)$id <- "ENSP00000466417"
    res <- .cds_for_id_range(edb, rng, id = "id")
    expect_equal(length(res[[1]]), 0)
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
    expect_true(is(res, "GRangesList"))
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
