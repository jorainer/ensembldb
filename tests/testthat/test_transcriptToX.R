test_that(".tx_to_protein works",  {
    edbx <- filter(EnsDb.Hsapiens.v86, filter = ~ seq_name %in% c("X", "Y"))

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
    mcols(res_b) <- mcols(res_b)[1:4]
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

    ## Preloaded data test
    proteins <- proteins(edbx)
    exons <- exonsBy(edbx)
    transcripts <- transcripts(edbx)
    expect_warning(res_c <- .txs_to_proteins(x, edbx, 
                                             proteins = proteins, 
                                             exons = exons, 
                                             transcripts = transcripts))
    expect_equal(res_c, res_b)
})

test_that("transcriptToProtein works", {
    edbx <- filter(EnsDb.Hsapiens.v86, filter = ~ seq_name == "X")

    ## Errors.
    expect_error(transcriptToProtein())
    expect_error(transcriptToProtein(db = edbx))
    expect_error(transcriptToProtein(txpos))

    ##
    x <- IRanges(start = 692, width = 2)
    mcols(x)$id <- "ENST00000381578"
    expect_error(transcriptToProtein(x))
    res <- transcriptToProtein(x, edbx, id = "id")
    expect_equal(start(res), 1)
    expect_equal(end(res), 1)

    ## Preloaded data test

    proteins <- proteins(edbx)
    exons <- exonsBy(edbx)
    transcripts <- transcripts(edbx)
    expect_error(transcriptToProtein(x, edbx, id = "id", 
                                     exons = exons, 
                                     transcripts = transcripts))
    expect_error(transcriptToProtein(x, edbx, 
                                     proteins = proteins, 
                                     exons = exons, 
                                     transcripts = transcripts))
    expect_error(transcriptToProtein(x, edbx, id = "id", 
                                     proteins = proteins, 
                                     exons = proteins, 
                                     transcripts = transcripts))
    res2 <- transcriptToProtein(x, edbx, id = "id", proteins = proteins, exons = exons, transcripts = transcripts)
    expect_equal(res2, res)
})

test_that(".tx_to_genome works", {
    edbx <- filter(EnsDb.Hsapiens.v86, filter = ~ seq_name %in% c("X", "Y"))
    ## ENST00000486554:501:505 106959129:106959131 106957979:106957979 -
    ## ENST00000486554:1:5 106959627:106959631 -
    ## ENST00000381578:1:5 585079:585083 +
    ## ENST00000381578:259:260 585337:585337 591201:591201 +
    ## some:1:4 NA
    ## ENST00000155093:2:3 2935478:29935479 +
    rng <- IRanges(start = c(501, 1, 1, 259, 1, 2),
                   end = c(505, 5, 5, 260, 4, 3),
                   names = c("ENST00000486554", "ENST00000486554",
                             "ENST00000381578", "ENST00000381578", "some",
                             "ENST00000155093"))
    res <- ensembldb:::.tx_to_genome(rng, edbx)
    expect_equal(length(res), length(rng))
    expect_equal(unname(lengths(res)), c(2, 1, 1, 2, 0, 1))
    ## 1
    expect_equal(start(res[[1]]), c(107715899, 107714749))
    expect_equal(end(res[[1]]), c(107715901, 107714749))
    expect_equal(as.character(strand(res[[1]])), c("-", "-"))
    ## 2
    expect_equal(start(res[[2]]), 107716397)
    expect_equal(end(res[[2]]), 107716401)
    expect_equal(as.character(strand(res[[2]])), "-")
    ## 3
    expect_equal(start(res[[3]]), 624344)
    expect_equal(end(res[[3]]), 624348)
    expect_equal(as.character(strand(res[[3]])), "+")
    ## 4
    expect_equal(start(res[[4]]), c(624602, 630466))
    expect_equal(end(res[[4]]), c(624602, 630466))
    expect_equal(as.character(strand(res[[4]])), c("+", "+"))
    ## 6
    expect_equal(start(res[[6]]), 2935478)
    expect_equal(end(res[[6]]), 2935479)
    expect_equal(as.character(strand(res[[6]])), "+")
    expect_equal(as.character(seqnames(res[[6]])), "Y")

    ## add pre-load data test 
    exons <- exonsBy(edbx)
    res_preload <- unlist(ensembldb:::.tx_to_genome(rng, exons))
    res_unlisted <- unlist(res)
    res_unlisted$tx_id <- NULL
    seqlevels(res_preload) <- seqlevels(res_unlisted)
    expect_equal(res_unlisted, res_preload)

    ## wrong ID and range outside tx
    rng <- IRanges(start = c(501, 200, 1), end = c(505, 1200, 5),
                   names = c("ENST00000486554", "ENST00000486554", "B"))
    expect_warning(res_2 <- ensembldb:::.tx_to_genome(rng, edbx))
    a <- unlist(res[1])
    b <- unlist(res_2[1])
    seqlevels(a) <- seqlevels(b)
    expect_equal(a, b)
    expect_equal(lengths(res_2), c(ENST00000486554 = 2, ENST00000486554 = 0,
                                   B = 0))
    ## add pre-load data test 
    expect_warning(res_3 <- ensembldb:::.tx_to_genome(rng, exons))
    c <- unlist(res_3[1])
    seqlevels(res_preload) <- seqlevels(c)
    expect_equal(res_preload, res_preload)
    expect_equal(lengths(res_3), c(ENST00000486554 = 2, ENST00000486554 = 0,
                                   B = 0))
})

test_that("transcriptToGenome works", {
    edbx <- filter(EnsDb.Hsapiens.v86, filter = ~ seq_name == "X")

    x <- IRanges(start = c(259, 1, 259), end = c(260, 4, 261),
                 names = c("ENST00000381578", "some", "ENST00000381578"))
    ## Errors.
    expect_error(transcriptToGenome())
    expect_error(transcriptToGenome(db = edbx))
    expect_error(transcriptToGenome(x))

    expect_warning(res <- transcriptToGenome(x, edbx))
    expect_true(is(res, "GRangesList"))
    expect_true(length(res) == length(x))
    expect_true(length(res[[2]]) == 0)
    expect_equal(names(res), names(x))
    expect_equal(start(res[[1]]), start(res[[3]]))
    expect_equal(end(res[[1]])[1], end(res[[3]])[1])
    expect_equal(end(res[[1]])[2] + 1, end(res[[3]])[2])

    expect_warning(res <- transcriptToGenome(x[2], edbx))
    expect_true(is(res, "GRangesList"))
    expect_true(length(res) == 1)

    expect_warning(res <- transcriptToGenome(x[c(2, 2)], edbx))
    expect_true(is(res, "GRangesList"))
    expect_true(length(res) == 2)

    res <- transcriptToGenome(x[1], edbx)
    expect_true(is(res, "GRangesList"))
    expect_true(length(res) == 1)
    expect_true(length(res[[1]]) == 2)

    x <- IRanges(start = c(256, 2), end = c(265, 12000),
                 names = c("ENST00000381578", "ENST00000381578"))
    expect_warning(res <- transcriptToGenome(x, edbx)) # indeed warn expected
    expect_equal(lengths(res), c(ENST00000381578 = 2, ENST00000381578 = 0))
    
    ## add pre-load data test 
    exons <- exonsBy(edbx)
    x <- IRanges(start = c(259, 1, 259), end = c(260, 4, 261),
                names = c("ENST00000381578", "some", "ENST00000381578"))
    expect_warning(res <- transcriptToGenome(x, exons))
    expect_true(is(res, "GRangesList"))
    expect_true(length(res) == length(x))
    expect_true(length(res[[2]]) == 0)
    expect_equal(names(res), names(x))
    expect_equal(start(res[[1]]), start(res[[3]]))
    expect_equal(end(res[[1]])[1], end(res[[3]])[1])
    expect_equal(end(res[[1]])[2] + 1, end(res[[3]])[2])

    expect_warning(res <- transcriptToGenome(x[2], exons))
    expect_true(is(res, "GRangesList"))
    expect_true(length(res) == 1)

    expect_warning(res <- transcriptToGenome(x[c(2, 2)], exons))
    expect_true(is(res, "GRangesList"))
    expect_true(length(res) == 2)

    res <- transcriptToGenome(x[1], exons)
    expect_true(is(res, "GRangesList"))
    expect_true(length(res) == 1)
    expect_true(length(res[[1]]) == 2)

    x <- IRanges(start = c(256, 2), end = c(265, 12000),
                 names = c("ENST00000381578", "ENST00000381578"),
                 info = 'ORFs')
    expect_warning(res <- transcriptToGenome(x, exons)) 
    expect_true(res$ENST00000381578$info[1] == 'ORFs')
    expect_equal(lengths(res), c(ENST00000381578 = 2, ENST00000381578 = 0))
})

test_that("transcriptToCds works", {
    expect_error(transcriptToCds())

    edb18 <- filter(EnsDb.Hsapiens.v86, filter = ~ seq_name == "18")
    expect_error(transcriptToCds(db = edb18))
    ## 1) unknown tx ids
    txcoords <- IRanges(start = c(4, 3), width = c(1, 1), names = c("a", "b"))
    expect_error(transcriptToCds(x = txcoords))
    expect_warning(res <- transcriptToCds(txcoords, edb18))
    expect_true(all(start(res) == -1))
    expect_true(all(end(res) == -1))
    ## 2) all tx not coding
    txcoords <- IRanges(start = c(132, 133), end = c(323, 323),
                        names = rep("ENST00000590515", 2))
    expect_warning(res <- transcriptToCds(txcoords, edb18))
    expect_true(all(start(res) == -1))
    expect_true(all(end(res) == -1))
    ## 3) some tx not coding
    ## 4) coordinate not within coding
    txcoords <- IRanges(start = c(1463, 3, 143, 147, 1463), width = 1,
                        names = c("ENST00000398117", "ENST00000333681",
                                  "ENST00000590515", "ENST00000589955",
                                  "ENST00000398117"))
    expect_warning(res <- transcriptToCds(txcoords, edb18))
    expect_equal(start(res), c(1, -1, -1, 1, 1))
    expect_equal(end(res), c(1, -1, -1, 1, 1))
    expect_equal(res[1], res[5])
    ## End position outside of CDS
    txcoords <- IRanges(start = c(1463, 3, 143, 147), width = c(4, 1, 1, 765),
                        names = c("ENST00000398117", "ENST00000333681",
                                  "ENST00000590515", "ENST00000589955"))
    expect_warning(res <- transcriptToCds(txcoords, edb18))
    expect_equal(start(res), c(1, -1, -1, -1))
    expect_equal(end(res), c(4, -1, -1, -1))
    txcoords <- IRanges(start = c(3, 143, 147), width = c(1, 1, 765),
                        names = c("ENST00000333681",
                                  "ENST00000590515", "ENST00000589955"))
    expect_warning(res <- transcriptToCds(txcoords, edb18))
    expect_true(all(start(res) < 0))
    expect_true(all(end(res) < 0))

    ## Preloaded data test
    exons <- exonsBy(edb18)
    transcripts <- transcripts(edb18)
    expect_error(transcriptToCds(db = edb18, exons = exons,transcripts = transcripts))
    ## 1) unknown tx ids
    txcoords <- IRanges(start = c(4, 3), width = c(1, 1), names = c("a", "b"))
    expect_error(transcriptToCds(x = txcoords))
    expect_warning(res <- transcriptToCds(txcoords, edb18, exons = exons,transcripts = transcripts))
    expect_true(all(start(res) == -1))
    expect_true(all(end(res) == -1))
    ## 2) all tx not coding
    txcoords <- IRanges(start = c(132, 133), end = c(323, 323),
                        names = rep("ENST00000590515", 2))
    expect_warning(res <- transcriptToCds(txcoords, edb18, exons = exons,transcripts = transcripts))
    expect_true(all(start(res) == -1))
    expect_true(all(end(res) == -1))
    ## 3) some tx not coding
    ## 4) coordinate not within coding
    txcoords <- IRanges(start = c(1463, 3, 143, 147, 1463), width = 1,
                        names = c("ENST00000398117", "ENST00000333681",
                                  "ENST00000590515", "ENST00000589955",
                                  "ENST00000398117"))
    expect_warning(res <- transcriptToCds(txcoords, edb18, exons = exons,transcripts = transcripts))
    expect_equal(start(res), c(1, -1, -1, 1, 1))
    expect_equal(end(res), c(1, -1, -1, 1, 1))
    expect_equal(res[1], res[5])
    ## End position outside of CDS
    txcoords <- IRanges(start = c(1463, 3, 143, 147), width = c(4, 1, 1, 765),
                        names = c("ENST00000398117", "ENST00000333681",
                                  "ENST00000590515", "ENST00000589955"))
    expect_warning(res <- transcriptToCds(txcoords, edb18, exons = exons,transcripts = transcripts))
    expect_equal(start(res), c(1, -1, -1, -1))
    expect_equal(end(res), c(4, -1, -1, -1))
    txcoords <- IRanges(start = c(3, 143, 147), width = c(1, 1, 765),
                        names = c("ENST00000333681",
                                  "ENST00000590515", "ENST00000589955"))
    expect_warning(res <- transcriptToCds(txcoords, edb18, exons = exons,transcripts = transcripts))
    expect_true(all(start(res) < 0))
    expect_true(all(end(res) < 0))
})

test_that("cdsToTranscript works", {
    edb18 <- filter(EnsDb.Hsapiens.v86, filter = ~ seq_name == "18")
    expect_error(cdsToTranscript())
    expect_error(cdsToTranscript(db = edb18))
    txcoords <- IRanges(start = c(4, 3, 143, 147), width = 1,
                        names = c("ENST00000398117", "ENST00000333681",
                                  "ENST00000590515", "ENST00000589955"))
    expect_error(cdsToTranscript(x = txcoords))
    expect_warning(res <- cdsToTranscript(txcoords, edb18))
    expect_equal(start(res), c(1466, 902, -1, 293))

    txcoords <- IRanges(start = c(4, 3, 50000, 147), width = 1,
                        names = c("ENST00000398117", "ENST00000398117",
                                  "ENST00000398117", "b"))
    expect_warning(res <- cdsToTranscript(txcoords, edb18))
    expect_equal(start(res), c(1466, 1465, -1, -1))

    ## Map variants:
    ## ENST00000070846:c.1643DelG
    ## ENST00000070846:c.1881DelC
    ## ENST00000379802:c.6995C>A
    ## ENST00000261590:c.1088C>T
    ## ENST00000261590:c.561T>G
    rngs <- IRanges(start = c(1643, 1881, 6995, 1088, 561), width = 1,
                    names = c("ENST00000070846", "ENST00000070846",
                              "ENST00000379802", "ENST00000261590",
                              "ENST00000261590"))
    rngs_tx <- cdsToTranscript(rngs, EnsDb.Hsapiens.v86)
    gnm <- transcriptToGenome(rngs_tx, EnsDb.Hsapiens.v86)

    library(BSgenome.Hsapiens.NCBI.GRCh38)
    res <- getSeq(BSgenome.Hsapiens.NCBI.GRCh38, unlist(gnm))
    exp <- c("G", "C", "C", "C", "T")
    expect_equal(exp, unname(as.character(res)))
    
    ## Preloaded data test
    exons <- exonsBy(EnsDb.Hsapiens.v86)
    transcripts <- transcripts(EnsDb.Hsapiens.v86)

    edb18 <- filter(EnsDb.Hsapiens.v86, filter = ~ seq_name == "18")
    expect_error(cdsToTranscript())
    expect_error(cdsToTranscript(db = edb18, exons = exons,transcripts = transcripts))
    txcoords <- IRanges(start = c(4, 3, 143, 147), width = 1,
                        names = c("ENST00000398117", "ENST00000333681",
                                  "ENST00000590515", "ENST00000589955"))
    expect_error(cdsToTranscript(x = txcoords))
    expect_warning(res <- cdsToTranscript(txcoords, edb18, exons = exons,transcripts = transcripts))
    expect_equal(start(res), c(1466, 902, -1, 293))

    txcoords <- IRanges(start = c(4, 3, 50000, 147), width = 1,
                        names = c("ENST00000398117", "ENST00000398117",
                                  "ENST00000398117", "b"))
    expect_warning(res <- cdsToTranscript(txcoords, edb18, exons = exons,transcripts = transcripts))
    expect_equal(start(res), c(1466, 1465, -1, -1))

    ## Map variants:
    ## ENST00000070846:c.1643DelG
    ## ENST00000070846:c.1881DelC
    ## ENST00000379802:c.6995C>A
    ## ENST00000261590:c.1088C>T
    ## ENST00000261590:c.561T>G
    rngs <- IRanges(start = c(1643, 1881, 6995, 1088, 561), width = 1,
                    names = c("ENST00000070846", "ENST00000070846",
                              "ENST00000379802", "ENST00000261590",
                              "ENST00000261590"))
    rngs_tx <- cdsToTranscript(rngs, EnsDb.Hsapiens.v86, exons = exons,transcripts = transcripts)
    gnm <- transcriptToGenome(rngs_tx, exons)

    library(BSgenome.Hsapiens.NCBI.GRCh38)
    res <- getSeq(BSgenome.Hsapiens.NCBI.GRCh38, unlist(gnm))
    exp <- c("G", "C", "C", "C", "T")
    expect_equal(exp, unname(as.character(res)))
})

test_that(".ids_message works", {
    res <- .ids_message(c("a", "b"))
    expect_equal(res, paste(c("a", "b"), collapse = ", "))
    res <- .ids_message(c("a", "b", "c", "d", "e", "f"))
    expect_equal(res, "a, b, c ... (3 more)")
    expect_equal(.ids_message(c("c", "b", "c", "d")), "c, b, c ... (1 more)")
})
