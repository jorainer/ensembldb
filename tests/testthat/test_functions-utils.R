
test_that("orderDataFrameBy works", {
    res <- exons(edb, filter = GenenameFilter("ZBTB16"),
                 return.type = "DataFrame")
    ## Order by end
    res_2 <- ensembldb:::orderDataFrameBy(res, by = "exon_seq_end")
    idx <- order(res_2$exon_seq_end)
    expect_equal(idx, 1:nrow(res_2))
})

test_that("addFilterColumns works for AnnotationFilterList", {
    afl <- AnnotationFilterList(GenenameFilter(2), SymbolFilter(23))
    afl2 <- AnnotationFilterList(SeqNameFilter(4), afl)
    res <- ensembldb:::addFilterColumns(cols = "gene_biotype", filter = afl, edb)
    expect_equal(res, c("gene_biotype", "gene_name", "symbol"))
    res <- ensembldb:::addFilterColumns(cols = "gene_biotype", filter = afl2,
                                        edb)
    expect_equal(res, c("gene_biotype", "seq_name", "gene_name", "symbol"))
})

test_that("functions work for encapsuled AnnotationFilterLists", {
    fl <- AnnotationFilterList(GenenameFilter("a"),
                               AnnotationFilterList(TxIdFilter("b")))
    res <- ensembldb:::.processFilterParam(fl, edb)
    expect_equal(res, fl)
    res <- ensembldb:::setFeatureInGRangesFilter(fl, "gene")
    expect_equal(res, fl)
    res <- ensembldb:::addFilterColumns("z", fl, edb)
    expect_equal(res, c("z", "gene_name", "tx_id"))
    res <- ensembldb:::getWhat(edb, filter = fl)
    ## Check if content is the same
    flts1 <- ~ genename == "BCL2" & tx_biotype == "protein_coding"
    res1 <- transcripts(edb, filter = flts1)
    flts2 <- AnnotationFilterList(
        AnnotationFilterList(GenenameFilter("BCL2"),
                             AnnotationFilterList(
                                 TxBiotypeFilter("protein_coding"))
                             ))
    res2 <- transcripts(edb, filter = flts2)
    expect_equal(res1, res2)
})

## Here we want to test if we get always also the filter columns back.
test_that("multiFilterReturnCols works also with symbolic filters", {
    cols <- ensembldb:::addFilterColumns(edb, cols = c("exon_id"),
                                         filter = SymbolFilter("SKA2"))
    expect_equal(cols, c("exon_id", "symbol"))
    ## Two filter
    cols <- ensembldb:::addFilterColumns(edb, cols = c("exon_id"),
                                         filter = list(SymbolFilter("SKA2"),
                                                       GenenameFilter("SKA2")))
    expect_equal(cols, c("exon_id", "symbol", "gene_name"))
    cols <- ensembldb:::addFilterColumns(edb, cols = c("exon_id"),
                                         filter = list(SymbolFilter("SKA2"),
                                                       GenenameFilter("SKA2"),
                                                       GRangesFilter(
                                                           GRanges("3",
                                                                   IRanges(3, 5)
                                                                   ))))
    expect_equal(cols, c("exon_id", "symbol", "gene_name", "gene_seq_start",
                         "gene_seq_end", "seq_name", "seq_strand"))
    cols <- ensembldb:::addFilterColumns(edb, cols = c("exon_id"),
                                         filter = list(SymbolFilter("SKA2"),
                                                       GenenameFilter("SKA2"),
                                                       GRangesFilter(
                                                           GRanges("3",
                                                                   IRanges(3, 5)
                                                                   ),
                                                           feature = "exon")))
    expect_equal(cols, c("exon_id", "symbol", "gene_name", "exon_seq_start",
                         "exon_seq_end", "seq_name", "seq_strand"))
    ## SeqStartFilter and GRangesFilter
    ssf <- TxStartFilter(123)
    cols <- ensembldb:::addFilterColumns(edb, cols = c("exon_id"),
                                         filter = list(SymbolFilter("SKA2"),
                                                       GenenameFilter("SKA2"),
                                                       GRangesFilter(
                                                           GRanges("3",
                                                                   IRanges(3, 5)
                                                                   ),
                                                           feature = "exon"),
                                                       ssf))
    expect_equal(cols, c("exon_id", "symbol", "gene_name", "exon_seq_start",
                         "exon_seq_end", "seq_name", "seq_strand",
                         "tx_seq_start"))
})

test_that("SQLiteName2MySQL works", {
    have <- "EnsDb.Hsapiens.v75"
    want <- "ensdb_hsapiens_v75"
    expect_equal(ensembldb:::SQLiteName2MySQL(have), want)
})

test_that("anyProteinColumns works", {
    expect_true(ensembldb:::anyProteinColumns(c("gene_id", "protein_id")))
    expect_true(!ensembldb:::anyProteinColumns(c("gene_id", "exon_id")))
})

test_that("listProteinColumns works", {
    if (hasProteinData(edb)) {
        res <- listProteinColumns(edb)
        expect_true(any(res == "protein_id"))
        expect_true(any(res == "uniprot_id"))
        expect_true(any(res == "protein_domain_id"))
        ## That's new columns fetched for Uniprot:
        expect_true(any(res == "uniprot_db"))
        expect_true(any(res == "uniprot_mapping_type"))
    } else {
        expect_error(listProteinColumns(edb))
    }
})

test_that("strand2num works", {
    expect_equal(ensembldb:::strand2num("+"), 1)
    expect_equal(ensembldb:::strand2num("+1"), 1)
    expect_equal(ensembldb:::strand2num("-"), -1)
    expect_equal(ensembldb:::strand2num("-1"), -1)
    expect_equal(ensembldb:::strand2num(1), 1)
    expect_equal(ensembldb:::strand2num(5), 1)
    expect_equal(ensembldb:::strand2num(-1), -1)
    expect_equal(ensembldb:::strand2num(-5), -1)
    expect_error(ensembldb:::strand2num("a"))
})

test_that("num2strand works", {
    expect_equal(ensembldb:::num2strand(1), "+")
    expect_equal(ensembldb:::num2strand(-1), "-")
})
