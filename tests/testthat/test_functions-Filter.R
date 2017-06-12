
test_that(".fieldInEnsDb works", {
    expect_equal(unname(ensembldb:::.fieldInEnsDb("symbol")), "gene_name")
    expect_equal(unname(ensembldb:::.fieldInEnsDb("gene_biotype")), "gene_biotype")
    expect_equal(unname(ensembldb:::.fieldInEnsDb("entrez")), "entrezid")
    expect_equal(unname(ensembldb:::.fieldInEnsDb("gene_id")), "gene_id")
    expect_equal(unname(ensembldb:::.fieldInEnsDb("genename")), "gene_name")
    expect_equal(unname(ensembldb:::.fieldInEnsDb("seq_name")), "seq_name")
    expect_equal(unname(ensembldb:::.fieldInEnsDb("tx_id")), "tx_id")
    expect_equal(unname(ensembldb:::.fieldInEnsDb("tx_biotype")), "tx_biotype")
    expect_equal(unname(ensembldb:::.fieldInEnsDb("tx_name")), "tx_id")
    expect_equal(unname(ensembldb:::.fieldInEnsDb("exon_id")), "exon_id")
    expect_equal(unname(ensembldb:::.fieldInEnsDb("exon_rank")), "exon_idx")
    expect_equal(unname(ensembldb:::.fieldInEnsDb("protein_id")), "protein_id")
    expect_equal(unname(ensembldb:::.fieldInEnsDb("uniprot")), "uniprot_id")
    expect_equal(unname(ensembldb:::.fieldInEnsDb("uniprot_db")), "uniprot_db")
    expect_equal(unname(ensembldb:::.fieldInEnsDb("uniprot_mapping_type")),
                 "uniprot_mapping_type")
    expect_equal(unname(ensembldb:::.fieldInEnsDb("prot_dom_id")),
                 "protein_domain_id")
    expect_equal(unname(ensembldb:::.fieldInEnsDb("description")),
                 "description")
    expect_equal(unname(ensembldb:::.fieldInEnsDb("tx_support_level")),
                 "tx_support_level")
    expect_error(ensembldb:::.fieldInEnsDb("aaa"))
})

test_that(".conditionForEnsDb works", {
    smb <- SymbolFilter("a")
    expect_equal(condition(smb), "==")
    expect_equal(ensembldb:::.conditionForEnsDb(smb), "=")
    smb <- SymbolFilter(c("a", "b", "c"))
    expect_equal(ensembldb:::.conditionForEnsDb(smb), "in")
    smb <- SymbolFilter(c("a", "b", "c"), condition = "!=")
    expect_equal(ensembldb:::.conditionForEnsDb(smb), "not in")
    smb <- SymbolFilter(c("a"), condition = "!=")
    expect_equal(ensembldb:::.conditionForEnsDb(smb), "!=")
    smb <- SymbolFilter(c("a"), condition = "startsWith")
    expect_equal(ensembldb:::.conditionForEnsDb(smb), "like")
    smb <- SymbolFilter(c("a"), condition = "endsWith")
    expect_equal(ensembldb:::.conditionForEnsDb(smb), "like")
    ## Tests for numeric filters
    fl <- GeneStartFilter(4)
    expect_equal(ensembldb:::.conditionForEnsDb(fl), "=")
    fl <- GeneStartFilter(4, condition = ">")
    expect_equal(ensembldb:::.conditionForEnsDb(fl), ">")
    fl <- GeneStartFilter(4, condition = ">=")
    expect_equal(ensembldb:::.conditionForEnsDb(fl), ">=")
    fl <- GeneStartFilter(4, condition = "<")
    expect_equal(ensembldb:::.conditionForEnsDb(fl), "<")
    fl <- GeneStartFilter(4, condition = "<=")
    expect_equal(ensembldb:::.conditionForEnsDb(fl), "<=")
    fl <- TxSupportLevelFilter(4, condition = "<=")
    expect_equal(ensembldb:::.conditionForEnsDb(fl), "<=")
})

test_that(".valueForEnsDb works", {
    smb <- SymbolFilter("a")
    expect_equal(ensembldb:::.valueForEnsDb(smb), "'a'")
    smb <- SymbolFilter(c("a", "b", "b", "c"))
    expect_equal(ensembldb:::.valueForEnsDb(smb), "('a','b','c')")
    smb <- SymbolFilter("a", condition = "startsWith")
    expect_equal(ensembldb:::.valueForEnsDb(smb), "'a%'")
    smb <- SymbolFilter("a", condition = "endsWith")
    expect_equal(ensembldb:::.valueForEnsDb(smb), "'%a'")
    ## Tests for numeric filters
    fl <- GeneStartFilter(4)
    expect_equal(ensembldb:::.valueForEnsDb(fl), 4)
})

test_that(".queryForEnsDb works", {
    smb <- SymbolFilter("a")
    expect_equal(ensembldb:::.queryForEnsDb(smb), "gene_name = 'a'")
    smb <- SymbolFilter(c("a", "x"), condition = "!=")
    expect_equal(ensembldb:::.queryForEnsDb(smb), "gene_name not in ('a','x')")
    ## Tests for numeric filters
    fl <- GeneStartFilter(5, condition = "<=")
    expect_equal(ensembldb:::.queryForEnsDb(fl), "gene_seq_start <= 5")
    fl <- TxSupportLevelFilter(5, condition = "<=")
    expect_equal(ensembldb:::.queryForEnsDb(fl), "tx_support_level <= 5")
})

test_that(".queryForEnsDbWithTables works", {
    smb <- SymbolFilter("a")
    expect_equal(ensembldb:::.queryForEnsDbWithTables(smb), "gene_name = 'a'")
    smb <- SymbolFilter(c("a", "x"), condition = "!=")
    expect_equal(ensembldb:::.queryForEnsDbWithTables(smb),
                 "gene_name not in ('a','x')")
    ## With edb
    smb <- SymbolFilter("a")
    expect_equal(ensembldb:::.queryForEnsDbWithTables(smb, edb),
                 "gene.gene_name = 'a'")
    smb <- SymbolFilter(c("a", "x"), condition = "!=")
    expect_equal(ensembldb:::.queryForEnsDbWithTables(smb, edb),
                 "gene.gene_name not in ('a','x')")
    ## With edb, tables
    smb <- SymbolFilter("a")
    expect_equal(ensembldb:::.queryForEnsDbWithTables(smb, edb, c("gene", "tx")),
                 "gene.gene_name = 'a'")
    fl <- GeneIdFilter("b")
    expect_equal(ensembldb:::.queryForEnsDbWithTables(fl, edb),
                 "gene.gene_id = 'b'")
    expect_equal(ensembldb:::.queryForEnsDbWithTables(fl, edb, c("tx", "gene")),
                 "tx.gene_id = 'b'")
    fl <- GeneIdFilter("b", condition = "contains")
    expect_equal(ensembldb:::.queryForEnsDbWithTables(fl, edb, c("tx", "gene")),
                 "tx.gene_id like '%b%'")
    ## Entrez
    if (as.numeric(ensembldb:::dbSchemaVersion(edb)) > 1) {
        fl <- EntrezFilter("g")
        expect_equal(ensembldb:::.queryForEnsDbWithTables(fl, edb),
                     "entrezgene.entrezid = 'g'")
        fl <- EntrezFilter("g", condition = "endsWith")
        expect_equal(ensembldb:::.queryForEnsDbWithTables(fl, edb),
                     "entrezgene.entrezid like '%g'")
    } else {
        fl <- EntrezFilter("g")
        expect_equal(ensembldb:::.queryForEnsDbWithTables(fl, edb),
                     "gene.entrezid = 'g'")
        fl <- EntrezFilter("g", condition = "endsWith")
        expect_equal(ensembldb:::.queryForEnsDbWithTables(fl, edb),
                     "gene.entrezid like '%g'")
    }
    ## Numeric filters
    fl <- TxStartFilter(123)
    expect_equal(ensembldb:::.queryForEnsDbWithTables(fl, edb),
                 "tx.tx_seq_start = 123")
    expect_error(ensembldb:::.queryForEnsDbWithTables(fl, edb, "gene"))
    if (any(listColumns(edb) == "tx_support_level")) {
        fl <- TxSupportLevelFilter(3)
        expect_equal(ensembldb:::.queryForEnsDbWithTables(fl, edb),
                     "tx.tx_support_level = 3")
    }
})

test_that(".processFilterParam works", {
    ## Check that the processFilterParam does what we expect. Check input and
    ## return ALWAYS an AnnotationFilterList object.
    snf <- SeqNameFilter(c("Y", 9))
    res <- ensembldb:::.processFilterParam(snf, db = edb)
    expect_true(is(res, "AnnotationFilterList"))
    expect_equal(res[[1]], snf)

    ## - single filter
    gif <- GeneIdFilter("BCL2", condition = "!=")
    res <- ensembldb:::.processFilterParam(gif, db = edb)
    expect_true(is(res, "AnnotationFilterList"))
    expect_equal(res[[1]], gif)
    
    ## - list of filters
    snf <- SeqNameFilter("X")
    res <- ensembldb:::.processFilterParam(list(gif, snf), edb)
    expect_true(is(res, "AnnotationFilterList"))
    expect_true(length(res) == 2)
    expect_equal(res[[1]], gif)
    expect_equal(res[[2]], snf)
    expect_equal(res@logOp, "&")
    
    ## - AnnotationFilterList
    afl <- AnnotationFilterList(gif, snf, logOp = "|")
    res <- ensembldb:::.processFilterParam(afl, edb)
    expect_true(is(res, "AnnotationFilterList"))
    expect_equal(afl, res)
    afl <- AnnotationFilterList(gif, snf, logOp = "&")
    res <- ensembldb:::.processFilterParam(afl, edb)
    expect_true(is(res, "AnnotationFilterList"))
    expect_equal(afl, res)
    
    ## - filter expression
    res <- ensembldb:::.processFilterParam(~ gene_id != "BCL2" |
                                               seq_name == "X", edb)
    expect_true(is(res, "AnnotationFilterList"))
    expect_equal(res, AnnotationFilterList(gif, snf, logOp = "|"))
    flt <- ~ gene_id != "BCL2" | seq_name == "X"
    res <- ensembldb:::.processFilterParam(flt, edb)
    expect_true(is(res, "AnnotationFilterList"))
    expect_equal(res, AnnotationFilterList(gif, snf, logOp = "|"))
    res <- ensembldb:::.processFilterParam(~ gene_id != "BCL2", edb)
    expect_true(is(res, "AnnotationFilterList"))
    expect_equal(res, AnnotationFilterList(gif))

    ## - Errors
    expect_error(ensembldb:::.processFilterParam(db = edb))
    expect_error(ensembldb:::.processFilterParam(4, edb))
    expect_error(ensembldb:::.processFilterParam(list(afl, "a"), edb))
    expect_error(ensembldb:::.processFilterParam("a", edb))
    expect_error(ensembldb:::.processFilterParam(~ gene_bla == "14", edb))
    ## Errors for filters that are not supported.
    expect_error(ensembldb:::.processFilterParam(CdsEndFilter(123), edb))
    
    ## Same with calls from within a function.
    testFun <- function(filter = AnnotationFilterList()) {
        ensembldb:::.processFilterParam(filter, db = edb)
    }

    res <- testFun()
    expect_true(is(res, "AnnotationFilterList"))
    expect_true(length(res) == 0)
    res <- testFun(filter = ~ gene_id == 4)
    expect_true(is(res, "AnnotationFilterList"))
    expect_equal(res[[1]], GeneIdFilter(4))
    res <- testFun(filter = GenenameFilter("BCL2"))
    expect_true(is(res, "AnnotationFilterList"))
    expect_equal(res[[1]], GenenameFilter("BCL2"))
    res <- testFun(filter = AnnotationFilterList(GenenameFilter("BCL2")))
    expect_true(is(res, "AnnotationFilterList"))
    expect_equal(res[[1]], GenenameFilter("BCL2"))
    
    gene <- "ZBTB16"
    otherFun <- function(gn) {
        testFun(filter = GenenameFilter(gn))
    }
    res <- otherFun(gene)
})

test_that("setFeatureInGRangesFilter works", {
    afl <- AnnotationFilterList(GeneIdFilter(123), SeqNameFilter(3),
                                GRangesFilter(GRanges()))
    afl2 <- AnnotationFilterList(afl, GRangesFilter(GRanges()))

    res <- ensembldb:::setFeatureInGRangesFilter(afl, feature = "tx")
    expect_equal(res[[3]]@feature, "tx")
    res <- ensembldb:::setFeatureInGRangesFilter(afl2, feature = "tx")
    expect_equal(res[[2]]@feature, "tx")
    expect_equal(res[[1]][[3]]@feature, "tx")    
})

test_that(".AnnottionFilterClassNames works", {
    afl1 <- AnnotationFilter(~ genename == 3 & seq_name != 5)
    expect_equal(.AnnotationFilterClassNames(afl1),
                 c("GenenameFilter", "SeqNameFilter"))
    afl2 <- AnnotationFilter(~ gene_start > 13 | seq_strand == "+")
    expect_equal(.AnnotationFilterClassNames(afl2),
                 c("GeneStartFilter", "SeqStrandFilter"))
    afl3 <- AnnotationFilterList(afl1, SymbolFilter(4))
    expect_equal(.AnnotationFilterClassNames(afl3),
                 c("GenenameFilter", "SeqNameFilter", "SymbolFilter"))
    afl4 <- AnnotationFilterList(afl2, afl3)
    expect_equal(.AnnotationFilterClassNames(afl4),
                 c("GeneStartFilter", "SeqStrandFilter", "GenenameFilter",
                   "SeqNameFilter", "SymbolFilter"))
})

test_that(".anyIs works", {
    sf <- SymbolFilter("d")
    expect_true(ensembldb:::.anyIs(sf, "SymbolFilter"))
    expect_false(ensembldb:::.anyIs(GenenameFilter(3), "SymbolFilter"))

    flts <- AnnotationFilterList(sf, TxIdFilter("b"))
    expect_true(any(ensembldb:::.anyIs(flts, "SymbolFilter")))
    expect_false(any(ensembldb:::.anyIs(flts, "BLa")))

    ## Additional nesting.
    flts <- AnnotationFilterList(flts, GenenameFilter("2"))
    expect_true(any(ensembldb:::.anyIs(flts, "SymbolFilter")))
    expect_true(any(ensembldb:::.anyIs(flts, "GenenameFilter")))
    expect_true(any(ensembldb:::.anyIs(flts, "TxIdFilter")))    
})
