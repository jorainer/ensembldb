test_that("OnlyCodingTxFilter constructor works", {
    fl <- OnlyCodingTxFilter()
    expect_true(is(fl, "OnlyCodingTxFilter"))
})

test_that("ProtDomIdFilter constructor works", {
    fl <- ProtDomIdFilter("a")
    expect_true(is(fl, "ProtDomIdFilter"))
    fl <- AnnotationFilter(~ prot_dom_id %in% 1:4)
    expect_true(is(fl, "ProtDomIdFilter"))
    expect_equal(value(fl), c("1", "2", "3", "4"))
    expect_equal(condition(fl), "==")
})

test_that("ProteinDomainIdFilter constructor works", {
    fl <- ProteinDomainIdFilter("a")
    expect_true(is(fl, "ProteinDomainIdFilter"))
    fl <- AnnotationFilter(~ protein_domain_id %in% 1:4)
    expect_true(is(fl, "ProteinDomainIdFilter"))
    expect_equal(value(fl), c("1", "2", "3", "4"))
    expect_equal(condition(fl), "==")
})

test_that("ProteinDomainSourceFilter constructor works", {
    fl <- ProteinDomainSourceFilter("a")
    expect_true(is(fl, "ProteinDomainSourceFilter"))
    fl <- AnnotationFilter(~ protein_domain_source %in% 1:4)
    expect_true(is(fl, "ProteinDomainSourceFilter"))
    expect_equal(value(fl), c("1", "2", "3", "4"))
    expect_equal(condition(fl), "==")
})

test_that("UniprotDbFilter constructor works", {
    fl <- UniprotDbFilter("a")
    expect_true(is(fl, "UniprotDbFilter"))
    fl <- AnnotationFilter(~ uniprot_db != "4")
    expect_true(is(fl, "UniprotDbFilter"))
    expect_equal(value(fl), "4")
    expect_equal(condition(fl), "!=")
})

test_that("UniprotMappingTypeFilter constructor works", {
    fl <- UniprotMappingTypeFilter("a")
    expect_true(is(fl, "UniprotMappingTypeFilter"))
    fl <- AnnotationFilter(~ uniprot_mapping_type != "4")
    expect_true(is(fl, "UniprotMappingTypeFilter"))
    expect_equal(value(fl), "4")
    expect_equal(condition(fl), "!=")
})

test_that("GRangesFilter works for EnsDb", {
    ## Testing slots
    gr <- GRanges("X", ranges = IRanges(123, 234), strand = "-")
    grf <- GRangesFilter(gr, type = "within")
    ## Now check some stuff
    expect_equal(start(grf), start(gr))
    expect_equal(end(grf), end(gr))
    expect_equal(as.character(strand(gr)), strand(grf))
    expect_equal(as.character(seqnames(gr)), seqnames(grf))

    ## Test column:
    ## filter alone.
    exp <- c(start = "gene_seq_start", end = "gene_seq_end",
             seqname = "seq_name", strand = "seq_strand")
    expect_equal(ensembldb:::ensDbColumn(grf), exp)
    grf@feature <- "tx"
    exp <- c(start = "tx_seq_start", end = "tx_seq_end",
             seqname = "seq_name", strand = "seq_strand")
    expect_equal(ensembldb:::ensDbColumn(grf), exp)
    grf@feature <- "exon"
    exp <- c(start = "exon_seq_start", end = "exon_seq_end",
             seqname = "seq_name", strand = "seq_strand")
    expect_equal(ensembldb:::ensDbColumn(grf), exp)
    ## filter and ensdb.
    exp <- c(start = "exon.exon_seq_start", end = "exon.exon_seq_end",
             seqname = "gene.seq_name", strand = "gene.seq_strand")
    expect_equal(ensembldb:::ensDbColumn(grf, edb), exp)
    grf@feature <- "tx"
    exp <- c(start = "tx.tx_seq_start", end = "tx.tx_seq_end",
             seqname = "gene.seq_name", strand = "gene.seq_strand")
    expect_equal(ensembldb:::ensDbColumn(grf, edb), exp)
    grf@feature <- "gene"
    exp <- c(start = "gene.gene_seq_start", end = "gene.gene_seq_end",
             seqname = "gene.seq_name", strand = "gene.seq_strand")
    expect_equal(ensembldb:::ensDbColumn(grf, edb), exp)

    exp <- paste0("(gene_seq_start>=123 and gene_seq_end<=234 and",
                  " seq_name='X' and seq_strand = -1)")
    expect_equal(ensembldb:::ensDbQuery(grf), exp)
    ## what if we set strand to *
    grf2 <- GRangesFilter(GRanges("1", IRanges(123, 234)), type = "within")
    exp <- paste0("(gene.gene_seq_start>=123 and gene.gene_seq_end<=234",
                  " and gene.seq_name='1')")
    expect_equal(ensembldb:::ensDbQuery(grf2, edb), exp)

    ## Now, using overlapping.
    grf2 <- GRangesFilter(GRanges("X", IRanges(123, 234), strand = "-"),
                          type = "any", feature = "transcript")
    exp <- paste0("(tx.tx_seq_start<=234 and tx.tx_seq_end>=123 and",
                  " gene.seq_name='X' and gene.seq_strand = -1)")
    expect_equal(ensembldb:::ensDbQuery(grf2, edb), exp)
})

test_that("TxSupportLevelFilter works for EnsDb", {
    fl <- TxSupportLevelFilter(3)
    expect_true(is(fl, "TxSupportLevelFilter"))
    fl <- AnnotationFilter(~ tx_support_level == 3)
    expect_true(is(fl, "TxSupportLevelFilter"))
})
