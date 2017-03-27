test_that("set returnFilterColumns works", {
    orig <- returnFilterColumns(edb)
    returnFilterColumns(edb) <- TRUE
    expect_equal(TRUE, returnFilterColumns(edb))
    returnFilterColumns(edb) <- FALSE
    expect_equal(FALSE, returnFilterColumns(edb))
    expect_error(returnFilterColumns(edb) <- "d")
    expect_error(returnFilterColumns(edb) <- c(TRUE, FALSE))
    ## Restore the "original" setting
    returnFilterColumns(edb) <- orig
})

test_that("returnFilterColumns works with_genes", {
    orig <- returnFilterColumns(edb)

    returnFilterColumns(edb) <- FALSE
    ## What happens if we use a GRangesFilter with return filter cols FALSE?
    grf <- GRangesFilter(GRanges(17, IRanges(57180000, 57233000)),
                         type = "within")
    res <- genes(edb, filter = grf)
    expect_equal(res$gene_id,
                c("ENSG00000224738", "ENSG00000182628", "ENSG00000252212",
                  "ENSG00000211514", "ENSG00000207996"))
    cols <- c("gene_id", "gene_name")
    res <- genes(edb, filter = grf, return.type = "data.frame",
                 columns = cols)
    ## Expect only the columns
    expect_equal(colnames(res), cols)
    returnFilterColumns(edb) <- TRUE
    res <- genes(edb, filter = grf, return.type = "data.frame",
                 columns = cols)
    ## Now I expect also the gene coords.
    expect_equal(colnames(res), c(cols, "gene_seq_start", "gene_seq_end",
                                 "seq_name", "seq_strand"))

    ## Use a gene biotype filter
    gbt <- GeneBiotypeFilter("protein_coding")

    returnFilterColumns(edb) <- TRUE
    res <- genes(edb, filter = list(gbt, grf), return.type = "data.frame",
                 columns = cols)
    expect_equal(res$gene_name, "SKA2")
    expect_equal(colnames(res), c(cols, "gene_biotype", "gene_seq_start",
                                 "gene_seq_end", "seq_name", "seq_strand"))
    returnFilterColumns(edb) <- FALSE
    res <- genes(edb, filter = list(gbt, grf), return.type = "data.frame",
                 columns = cols)
    expect_equal(colnames(res), cols)
    returnFilterColumns(edb) <- orig
})

test_that("returnFilterColumns works with_tx", {
    orig <- returnFilterColumns(edb)
    returnFilterColumns(edb) <- FALSE
    ## What happens if we use a GRangesFilter with return filter cols FALSE?
    grf <- GRangesFilter(GRanges(17, IRanges(57180000, 57233000)),
                         type = "within")
    res <- transcripts(edb, filter = grf)
    cols <- c("tx_id", "gene_name")
    res <- transcripts(edb, filter = grf, return.type = "data.frame",
                       columns = cols)
    ## Expect only the columns
    expect_equal(colnames(res), cols)
    returnFilterColumns(edb) <- TRUE
    res <- transcripts(edb, filter = grf, return.type = "data.frame",
                       columns = cols)
    ## Now I expect also the gene coords.
    expect_equal(colnames(res), c(cols, "tx_seq_start", "tx_seq_end", "seq_name",
                                 "seq_strand"))
    ## Use a gene biotype filter
    gbt <- GeneBiotypeFilter("protein_coding")

    returnFilterColumns(edb) <- TRUE
    res <- transcripts(edb, filter = list(gbt, grf), return.type = "data.frame",
                       columns = cols)
    expect_equal(unique(res$gene_name), "SKA2")
    expect_equal(colnames(res), c(cols, "gene_biotype", "tx_seq_start",
                                 "tx_seq_end", "seq_name", "seq_strand"))
    returnFilterColumns(edb) <- FALSE
    res <- transcripts(edb, filter = list(gbt, grf), return.type = "data.frame",
                       columns = cols)
    expect_equal(colnames(res), cols)
    returnFilterColumns(edb) <- orig
})

test_that("returnFilterColumns works with exons", {
    orig <- returnFilterColumns(edb)
    returnFilterColumns(edb) <- FALSE
    ## What happens if we use a GRangesFilter with return filter cols FALSE?
    grf <- GRangesFilter(GRanges(17, IRanges(57180000, 57233000)),
                         type = "within")
    res <- exons(edb, filter = grf)
    cols <- c("exon_id", "gene_name")
    res <- exons(edb, filter = grf, return.type = "data.frame",
                 columns = cols)
    ## Expect only the columns
    expect_equal(colnames(res), cols)
    returnFilterColumns(edb) <- TRUE
    res <- exons(edb, filter = grf, return.type = "data.frame",
                 columns = cols)
    ## Now I expect also the gene coords.
    expect_equal(colnames(res), c(cols, "exon_seq_start", "exon_seq_end",
                                 "seq_name", "seq_strand"))
    ## Use a gene biotype filter
    gbt <- GeneBiotypeFilter("protein_coding")
    returnFilterColumns(edb) <- TRUE
    res <- exons(edb, filter = list(gbt, grf), return.type = "data.frame",
                 columns = cols)
    expect_equal(unique(res$gene_name), c("TRIM37", "SKA2"))
    expect_equal(colnames(res), c(cols, "gene_biotype", "exon_seq_start",
                                 "exon_seq_end", "seq_name", "seq_strand"))
    returnFilterColumns(edb) <- FALSE
    res <- exons(edb, filter = list(gbt, grf), return.type = "data.frame",
                 columns = cols)
    expect_equal(colnames(res), cols)
    returnFilterColumns(edb) <- orig
})

test_that("returnFilterColumns works with exonsBy", {
    orig <- returnFilterColumns(edb)
    returnFilterColumns(edb) <- FALSE
    ## What happens if we use a GRangesFilter with return filter cols FALSE?
    grf <- GRangesFilter(GRanges(17, IRanges(57180000, 57233000)),
                         type = "within")
    ## By genes
    cols <- c("exon_id", "gene_name")
    res <- exonsBy(edb, by = "gene", filter = grf, columns = cols)
    res <- unlist(res)
    ## Expect only the columns
    expect_equal(colnames(mcols(res)), cols)
    returnFilterColumns(edb) <- TRUE
    res <- exonsBy(edb, by = "gene", filter = grf, columns = cols)
    res <- unlist(res)
    ## Now I expect also the gene coords, but not the seq_name and seq_strand,
    ## as these are redundant with data which is in the GRanges!
    expect_equal(colnames(mcols(res)), c(cols, "gene_seq_start", "gene_seq_end"))
    ## Use a gene biotype filter
    gbt <- GeneBiotypeFilter("protein_coding")
    returnFilterColumns(edb) <- TRUE
    res <- unlist(exonsBy(edb, by = "gene", filter = list(gbt, grf), columns = cols))
    expect_equal(unique(res$gene_name), c("SKA2"))
    expect_equal(colnames(mcols(res)), c(cols, "gene_biotype", "gene_seq_start", "gene_seq_end"))
    returnFilterColumns(edb) <- FALSE
    res <- unlist(exonsBy(edb, by = "gene", filter = list(gbt, grf), columns = cols))
    expect_equal(colnames(mcols(res)), cols)
    ## By tx
    returnFilterColumns(edb) <- FALSE
    cols <- c("exon_id", "gene_name")
    res <- exonsBy(edb, by = "tx", filter = grf, columns = cols)
    res <- unlist(res)
    ## Expect only the columns
    expect_equal(colnames(mcols(res)), c(cols, "exon_rank"))
    returnFilterColumns(edb) <- TRUE
    res <- exonsBy(edb, by = "tx", filter = grf, columns = cols)
    res <- unlist(res)
    ## Now I expect also the gene coords.
    expect_equal(colnames(mcols(res)), c(cols, "tx_seq_start", "tx_seq_end",
                                        "exon_rank"))
    ## Use a gene biotype filter
    gbt <- GeneBiotypeFilter("protein_coding")
    returnFilterColumns(edb) <- TRUE
    res <- unlist(exonsBy(edb, by = "tx", filter = list(gbt, grf),
                          columns = cols))
    expect_equal(unique(res$gene_name), c("SKA2"))
    expect_equal(colnames(mcols(res)), c(cols, "gene_biotype", "tx_seq_start",
                                        "tx_seq_end", "exon_rank"))
    returnFilterColumns(edb) <- FALSE
    res <- unlist(exonsBy(edb, by = "tx", filter = list(gbt, grf),
                          columns = cols))
    expect_equal(colnames(mcols(res)), c(cols, "exon_rank"))
    returnFilterColumns(edb) <- orig
})

test_that("returnFilterColumns works with transcriptsBy", {
    orig <- returnFilterColumns(edb)
    returnFilterColumns(edb) <- FALSE
    ## What happens if we use a GRangesFilter with return filter cols FALSE?
    grf <- GRangesFilter(GRanges(17, IRanges(57180000, 57233000)),
                         type = "within")
    ## By genes
    cols <- c("tx_id", "gene_name")
    res <- transcriptsBy(edb, by = "gene", filter = grf, columns = cols)
    res <- unlist(res)
    ## Expect only the columns
    expect_equal(colnames(mcols(res)), cols)
    returnFilterColumns(edb) <- TRUE
    res <- transcriptsBy(edb, by = "gene", filter = grf, columns = cols)
    res <- unlist(res)
    ## Now I expect also the gene coords.
    expect_equal(colnames(mcols(res)), c(cols, "gene_seq_start", "gene_seq_end"))
    ## Use a gene biotype filter
    gbt <- GeneBiotypeFilter("protein_coding")
    returnFilterColumns(edb) <- TRUE
    res <- unlist(transcriptsBy(edb, by = "gene", filter = list(gbt, grf),
                                columns = cols))
    expect_equal(unique(res$gene_name), c("SKA2"))
    expect_equal(colnames(mcols(res)),
                c(cols, "gene_biotype", "gene_seq_start", "gene_seq_end"))
    returnFilterColumns(edb) <- FALSE
    res <- unlist(transcriptsBy(edb, by = "gene", filter = list(gbt, grf),
                                columns = cols))
    expect_equal(colnames(mcols(res)), cols)
    returnFilterColumns(edb) <- orig
})

test_that("returnFilterColumns works with_cdsBy", {
    orig <- returnFilterColumns(edb)
    grf <- GRangesFilter(GRanges(17, IRanges(57180000, 57233000)),
                         type = "within")
    ## By tx
    returnFilterColumns(edb) <- FALSE
    cols <- c("gene_id", "gene_name")
    res <- cdsBy(edb, by = "tx", filter = grf, columns = cols)
    res <- unlist(res)
    ## Expect only the columns
    expect_equal(colnames(mcols(res)), c(cols, "exon_id", "exon_rank"))
    returnFilterColumns(edb) <- TRUE
    res <- cdsBy(edb, by = "tx", filter = grf, columns = cols)
    res <- unlist(res)
    ## Now I expect also the gene coords.
    expect_equal(colnames(mcols(res)), c(cols, "tx_seq_start", "tx_seq_end",
                                        "seq_name", "seq_strand", "exon_id",
                                        "exon_rank"))
    ## Use a gene biotype filter
    gbt <- GeneBiotypeFilter("protein_coding")
    returnFilterColumns(edb) <- TRUE
    res <- unlist(cdsBy(edb, by = "tx", filter = list(gbt, grf), columns = cols))
    expect_equal(unique(res$gene_name), c("SKA2"))
    expect_equal(colnames(mcols(res)), c(cols, "gene_biotype", "tx_seq_start",
                                        "tx_seq_end", "seq_name", "seq_strand",
                                        "exon_id", "exon_rank"))
    returnFilterColumns(edb) <- FALSE
    res <- unlist(cdsBy(edb, by = "tx", filter = list(gbt, grf), columns = cols))
    expect_equal(colnames(mcols(res)), c(cols, "exon_id", "exon_rank"))
    returnFilterColumns(edb) <- orig
})

test_that("returnFilterColumns works with threeUTRsByTranscript", {
    orig <- returnFilterColumns(edb)
    grf <- GRangesFilter(GRanges(17, IRanges(57180000, 57233000)),
                         type = "within")
    ## By tx
    returnFilterColumns(edb) <- FALSE
    cols <- c("gene_id", "gene_name")
    res <- threeUTRsByTranscript(edb, filter = grf, columns = cols)
    res <- unlist(res)
    ## Expect only the columns
    expect_equal(colnames(mcols(res)), c(cols, "exon_id", "exon_rank"))
    returnFilterColumns(edb) <- TRUE
    res <- threeUTRsByTranscript(edb, filter = grf, columns = cols)
    res <- unlist(res)
    ## Now I expect also the gene coords.
    expect_equal(colnames(mcols(res)), c(cols, "tx_seq_start", "tx_seq_end",
                                        "seq_name", "seq_strand", "exon_id",
                                        "exon_rank"))
    ## Use a gene biotype filter
    gbt <- GeneBiotypeFilter("protein_coding")
    returnFilterColumns(edb) <- TRUE
    res <- unlist(threeUTRsByTranscript(edb, filter = list(gbt, grf),
                                        columns = cols))
    expect_equal(unique(res$gene_name), c("SKA2"))
    expect_equal(colnames(mcols(res)), c(cols, "gene_biotype", "tx_seq_start",
                                        "tx_seq_end", "seq_name", "seq_strand",
                                        "exon_id", "exon_rank"))
    returnFilterColumns(edb) <- FALSE
    res <- unlist(threeUTRsByTranscript(edb, filter = list(gbt, grf),
                                        columns = cols))
    expect_equal(colnames(mcols(res)), c(cols, "exon_id", "exon_rank"))
    returnFilterColumns(edb) <- orig
})

