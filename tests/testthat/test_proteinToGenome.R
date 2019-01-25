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
    expect_warning(res <- ensembldb:::.cds_for_id(edb, prts))
    expect_equal(names(res), prts)
    expect_equal(names(res[[1]]), txs[1])
    expect_equal(names(res[[3]]), txs[2])
    expect_true(is.null(res[[2]]))
    res_2 <- ensembldb:::.cds_for_id2(edb, prts)
    fix_seqlevels <- function(z) {
        if (length(z))
            z <- keepSeqlevels(
                z, unique(as.character(unlist(seqnames(z),
                                              use.names = FALSE))))
        z
    }
    expect_equal(res, lapply(res_2, fix_seqlevels))

    ## Use uniprot...
    ## want a uniprot that is mapped to two tx.
    ## unis <- c("E7EMD7_HUMAN", "ALAT2_HUMAN", "H0YIP2_HUMAN")
    unis <- c("Q9NPH5", "Q8TD30", "H0YIP2")
    unis_counts <- c(10, 2, 1)
    ## dbsub2 <- filter(edb, filter = ~ seq_name %in% c("11", "12", "16"))
    res <- ensembldb:::.cds_for_id(edb, unis, idType = "uniprot_id")
    expect_equal(unname(lengths(res)), unis_counts)
    expect_equal(names(res), unis)

    res_2 <- ensembldb:::.cds_for_id2(edb, unis, idType = "uniprot_id")
    expect_equal(lapply(res_2, fix_seqlevels), res)

    rngs <- IRanges(start = c(11, 6, 4), end = c(19, 6, 9))
    names(rngs) <- unis
    ## Can I use that right away?
    z <- res[[1]]
    irng <- IRanges(start = 5, end = 8)
    tmp <- ensembldb:::.to_genome(z, irng)

    ## Run .cds_matching_protein on the results to avoid re-querying.
    ## H0YIP2_HUMAN has incomplete 5' and 3' CDS.
    expect_warning(clnd <- ensembldb:::.cds_matching_protein(edb, res))
    expect_equal(unname(lengths(clnd)), c(10, 2, 1))
    expect_true(is.list(clnd))
    expect_true(is(clnd[[1]], "GRangesList"))
    expect_true(all(clnd[[1]][[1]]$cds_ok))
    expect_true(all(clnd[[2]][[1]]$cds_ok))
    expect_true(all(clnd[[3]][[1]]$cds_ok == FALSE))

    ## ENSP00000437716 encoding transcript has an incomplete 3' CDS
    expect_warning(res <- ensembldb:::.cds_for_id(edb, prts))
    expect_warning(clnd <- ensembldb:::.cds_matching_protein(edb, res))
    expect_equal(unname(lengths(clnd)), c(1, 0, 1))
    expect_true(all(clnd[[1]][[1]]$cds_ok))
    expect_true(all(clnd[[3]][[1]]$cds_ok == FALSE))
})

test_that("proteinToGenome works", {
    ## Restrict to chromosome X - speeds up stuff
    edbx <- filter(EnsDb.Hsapiens.v86, filter = ~ seq_name == "X")
    ## protein_id
    ids <- c("ENSP00000425155", "ENSP00000418169", "ENSP00000441743",
             "ENSP00000385415")
    prngs <- IRanges(start = c(23, 23, 141, 23), end = c(24, 24, 421, 111))
    names(prngs) <- ids

    ## Test errors
    expect_error(proteinToGenome())
    expect_error(proteinToGenome(5, db = edbx))
    expect_error(proteinToGenome(prngs))
    expect_error(proteinToGenome(prngs, db = 5))

    expect_warning(res <- proteinToGenome(prngs, edbx))
    ## Result has to be a triplet
    expect_true(all(sapply(res, function(z) sum(width(z))) %% 3 == 0))
    ## Manually check for the ranges here.
    ## ENSP00000418169, SYP, encoded by ENST00000479808, - strand
    ## coords cds: 23 -> 67, 24 -> 72
    ## exon 1 width 49, 5' UTR: 13nt long. 36 nt
    ## exon 2 width 66. relative pos is 31nt in exon 2
    ## exon 2 starts: 49,055,492 [nt1], [nt31] is 49055492 - 30 = 49055462
    expect_equal(end(res[[2]]), 49199003)
    expect_equal(start(res[[2]]), 49199003 - 5)
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
    ## strts <- c(49161405, 49161883, 49173642, 49176207)
    ## ends <- c(49161419, 49162003, 49173767, 49176211)
    strts <- c(49304926, 49305404, 49317163, 49319728)
    ends <- c(49304940, 49305524, 49317288, 49319732)
    expect_equal(start(res[[4]]), strts)
    expect_equal(end(res[[4]]), ends)
    expect_equal(res[[4]]$exon_rank, c(2, 3, 4, 5))
    expect_true(all(res[[4]]$cds_ok))

    ## Uniprot identifier
    ## ids <- c("D6RDZ7_HUMAN", "SHOX_HUMAN", "TMM27_HUMAN", "unexistant")
    ids <- c("D6RDZ7", "O15266", "Q9HBJ8", "unexistant")
    prngs <- IRanges(start = c(1, 13, 43, 100), end = c(2, 21, 80, 100))
    names(prngs) <- ids

    expect_warning(res <- proteinToGenome(prngs, edbx, idType = "uniprot_id"))
    ## Now, expect elements 1 and 2 to be a GRangesList and not a GRanges.
    expect_true(is(res[[1]], "GRanges"))
    expect_true(is(res[[2]], "GRangesList"))
    expect_true(is(res[[3]], "GRanges"))
    expect_true(is(res[[4]], "GRanges"))
    expect_true(length(res[[4]]) == 0)
    res <- res[1:3]
    expect_true(all(unlist(lapply(res, function(z) sum(width(z)))) %% 3 == 0))
    ## We've got multi-mapping here
    expect_equal(unname(lengths(res)), c(1, 4, 2))
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

test_that("proteinToTranscript works", {
    prng <- IRanges(start = c(1, 2, 3, 11, 12), end = c(1, 4, 9, 12, 12),
                    names = c("ENSP00000014935", "ENSP00000173898",
                              "ENSP00000217901", "ENSP00000425155", "ffff"))

    expect_error(proteinToTranscript())
    expect_error(proteinToTranscript(5, db = edb))
    expect_error(proteinToTranscript(prng))
    expect_error(proteinToTranscript(prng, db = 5))

    expect_warning(tx_rel <- proteinToTranscript(prng, edb))
    expect_true(is(tx_rel, "IRangesList"))
    expect_true(all(unlist(lapply(tx_rel, function(z) is(z, "IRanges")))))
    expect_equal(unname(lengths(tx_rel)), c(1, 1, 1, 1, 1))
    expect_true(start(tx_rel[[5]]) == -1)
    ## All have cds OK, except 4
    expect_true(mcols(tx_rel[[1]])$cds_ok)
    expect_true(mcols(tx_rel[[2]])$cds_ok)
    expect_true(mcols(tx_rel[[3]])$cds_ok)
    expect_false(mcols(tx_rel[[4]])$cds_ok)
    expect_true(is.na(mcols(tx_rel[[5]])$cds_ok))
    expect_equal(unlist(width(tx_rel), use.names = FALSE), c(3, 9, 21, 6, 1))
    expect_equal(start(tx_rel[["ENSP00000014935"]]), (622L + 85L + 87L + 1L))
    expect_equal(start(tx_rel[["ENSP00000173898"]]), (68L + 44L + 4L))
    expect_equal(start(tx_rel[["ENSP00000217901"]]), (197L + 7L))
    expect_equal(start(tx_rel[["ENSP00000425155"]]), (86L + 92L + 31L))

    ## Now for Uniprot...
    ## ids <- c("D6RDZ7_HUMAN", "SHOX_HUMAN", "TMM27_HUMAN", "unexistant")
    ids <- c("D6RDZ7", "O15266", "Q9HBJ8", "unexistant")
    prngs <- IRanges(start = c(1, 13, 43, 100), end = c(2, 21, 80, 100))
    names(prngs) <- ids
    expect_warning(tx_rel <- proteinToTranscript(prngs, edb,
                                                 idType = "uniprot_id"))
    expect_equal(length(tx_rel), length(prngs))
    expect_equal(unname(lengths(tx_rel)), c(1, 4, 1, 1))
    lens <- unlist(width(tx_rel))
    expect_equal(unique(lens[1]), width(prngs)[1] * 3)
    expect_equal(unique(lens[2:5]), width(prngs)[2] * 3)
    expect_equal(unique(lens[6]), width(prngs)[3] * 3)

    ## Mapping fails for all...
    ids <- c("a", "unexistant")
    prngs <- IRanges(start = c(1, 13), end = c(2, 21))
    names(prngs) <- ids

    expect_warning(res <- proteinToTranscript(prngs, edb))
    expect_true(is(res, "IRangesList"))
    expect_equal(unlist(unname(start(res))), c(-1, -1))
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

    ## outside throws warning not error.
    c_coords_2 <- IRanges(start = 1, end = 1200)
    res <- expect_warning(.to_genome(g_coords, c_coords_2))
    expect_equal(lengths(res), c(0, 0))

    ## Check errors.
    expect_error(.to_genome(c_coords, c_coords))
    g_coords <- GRanges("2", IRanges(start = c(3, 8, 15, 19),
                                     end = c(5, 12, 16, 21)),
                        strand = "+")
    c_coords <- IRanges(start = 5, end = 40)
    expect_warning(.to_genome(g_coords, c_coords))
    c_coords <- IRanges(start = 40, end = 50)
    expect_warning(.to_genome(g_coords, c_coords))
})
