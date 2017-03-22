library(RUnit)
library(ensembldb)

test_fieldInEnsDb <- function() {
    checkEquals(unname(ensembldb:::.fieldInEnsDb("symbol")), "gene_name")
    checkEquals(unname(ensembldb:::.fieldInEnsDb("gene_biotype")), "gene_biotype")
    checkEquals(unname(ensembldb:::.fieldInEnsDb("entrez")), "entrezid")
    checkEquals(unname(ensembldb:::.fieldInEnsDb("gene_id")), "gene_id")
    checkEquals(unname(ensembldb:::.fieldInEnsDb("genename")), "gene_name")
    checkEquals(unname(ensembldb:::.fieldInEnsDb("seq_name")), "seq_name")
    checkEquals(unname(ensembldb:::.fieldInEnsDb("tx_id")), "tx_id")
    checkEquals(unname(ensembldb:::.fieldInEnsDb("tx_biotype")), "tx_biotype")
    checkEquals(unname(ensembldb:::.fieldInEnsDb("tx_name")), "tx_id")
    checkEquals(unname(ensembldb:::.fieldInEnsDb("exon_id")), "exon_id")
    checkEquals(unname(ensembldb:::.fieldInEnsDb("exon_rank")), "exon_idx")
    checkEquals(unname(ensembldb:::.fieldInEnsDb("protein_id")), "protein_id")
    checkEquals(unname(ensembldb:::.fieldInEnsDb("uniprot")), "uniprot_id")
    checkEquals(unname(ensembldb:::.fieldInEnsDb("uniprot_db")), "uniprot_db")
    checkEquals(unname(ensembldb:::.fieldInEnsDb("uniprot_mapping_type")),
                "uniprot_mapping_type")
    checkEquals(unname(ensembldb:::.fieldInEnsDb("prot_dom_id")),
                "protein_domain_id")
    checkException(ensembldb:::.fieldInEnsDb("aaa"))
}

test_conditionForEnsDb <- function() {
    smb <- SymbolFilter("a")
    checkEquals(condition(smb), "==")
    checkEquals(ensembldb:::.conditionForEnsDb(smb), "=")
    smb <- SymbolFilter(c("a", "b", "c"))
    checkEquals(ensembldb:::.conditionForEnsDb(smb), "in")
    smb <- SymbolFilter(c("a", "b", "c"), condition = "!=")
    checkEquals(ensembldb:::.conditionForEnsDb(smb), "not in")
    smb <- SymbolFilter(c("a"), condition = "!=")
    checkEquals(ensembldb:::.conditionForEnsDb(smb), "!=")
    smb <- SymbolFilter(c("a"), condition = "startsWith")
    checkEquals(ensembldb:::.conditionForEnsDb(smb), "like")
    smb <- SymbolFilter(c("a"), condition = "endsWith")
    checkEquals(ensembldb:::.conditionForEnsDb(smb), "like")
    ## Tests for numeric filters
    fl <- GeneStartFilter(4)
    checkEquals(ensembldb:::.conditionForEnsDb(fl), "=")
    fl <- GeneStartFilter(4, condition = ">")
    checkEquals(ensembldb:::.conditionForEnsDb(fl), ">")
    fl <- GeneStartFilter(4, condition = ">=")
    checkEquals(ensembldb:::.conditionForEnsDb(fl), ">=")
    fl <- GeneStartFilter(4, condition = "<")
    checkEquals(ensembldb:::.conditionForEnsDb(fl), "<")
    fl <- GeneStartFilter(4, condition = "<=")
    checkEquals(ensembldb:::.conditionForEnsDb(fl), "<=")
}

test_valueForEnsDb <- function() {
    smb <- SymbolFilter("a")
    checkEquals(ensembldb:::.valueForEnsDb(smb), "'a'")
    smb <- SymbolFilter(c("a", "b", "b", "c"))
    checkEquals(ensembldb:::.valueForEnsDb(smb), "('a','b','c')")
    smb <- SymbolFilter("a", condition = "startsWith")
    checkEquals(ensembldb:::.valueForEnsDb(smb), "'a%'")
    smb <- SymbolFilter("a", condition = "endsWith")
    checkEquals(ensembldb:::.valueForEnsDb(smb), "'%a'")
    ## Tests for numeric filters
    fl <- GeneStartFilter(4)
    checkEquals(ensembldb:::.valueForEnsDb(fl), 4)
}

test_queryForEnsDb <- function() {
    smb <- SymbolFilter("a")
    checkEquals(ensembldb:::.queryForEnsDb(smb), "gene_name = 'a'")
    smb <- SymbolFilter(c("a", "x"), condition = "!=")
    checkEquals(ensembldb:::.queryForEnsDb(smb), "gene_name not in ('a','x')")
    ## Tests for numeric filters
    fl <- GeneStartFilter(5, condition = "<=")
    checkEquals(ensembldb:::.queryForEnsDb(fl), "gene_seq_start <= 5")
}

test_queryForEnsDbWithTables <- function() {
    smb <- SymbolFilter("a")
    checkEquals(ensembldb:::.queryForEnsDbWithTables(smb), "gene_name = 'a'")
    smb <- SymbolFilter(c("a", "x"), condition = "!=")
    checkEquals(ensembldb:::.queryForEnsDbWithTables(smb),
                "gene_name not in ('a','x')")
    ## With edb
    smb <- SymbolFilter("a")
    checkEquals(ensembldb:::.queryForEnsDbWithTables(smb, edb),
                "gene.gene_name = 'a'")
    smb <- SymbolFilter(c("a", "x"), condition = "!=")
    checkEquals(ensembldb:::.queryForEnsDbWithTables(smb, edb),
                "gene.gene_name not in ('a','x')")
    ## With edb, tables
    smb <- SymbolFilter("a")
    checkEquals(ensembldb:::.queryForEnsDbWithTables(smb, edb, c("gene", "tx")),
                "gene.gene_name = 'a'")
    fl <- GeneIdFilter("b")
    checkEquals(ensembldb:::.queryForEnsDbWithTables(fl, edb),
                "gene.gene_id = 'b'")
    checkEquals(ensembldb:::.queryForEnsDbWithTables(fl, edb, c("tx", "gene")),
                "tx.gene_id = 'b'")
    ## Entrez
    fl <- EntrezFilter("g")
    checkEquals(ensembldb:::.queryForEnsDbWithTables(fl, edb),
                "gene.entrezid = 'g'")
    fl <- EntrezFilter("g", condition = "endsWith")
    checkEquals(ensembldb:::.queryForEnsDbWithTables(fl, edb),
                "gene.entrezid like '%g'")
    ## Numeric filters
    fl <- TxStartFilter(123)
    checkEquals(ensembldb:::.queryForEnsDbWithTables(fl, edb),
                "tx.tx_seq_start = 123")
    checkException(ensembldb:::.queryForEnsDbWithTables(fl, edb, "gene"),
                    "tx.tx_seq_start = 123")
}

test_processFilterParam <- function() {
    ## Check that the processFilterParam does what we expect. Check input and
    ## return ALWAYS an AnnotationFilterList object.

    ## - No input
    ## res <- ensembldb:::.processFilterParam()
    ## checkTrue(is(res, "AnnotationFilterList"))
    ## checkTrue(length(res) == 0)
    snf <- SeqNameFilter(c("Y", 9))
    res <- ensembldb:::.processFilterParam(snf)
    checkTrue(is(res, "AnnotationFilterList"))
    checkEquals(res[[1]], snf)

    ## - single filter
    gif <- GeneIdFilter("BCL2", condition = "!=")
    res <- ensembldb:::.processFilterParam(gif)
    checkTrue(is(res, "AnnotationFilterList"))
    checkEquals(res[[1]], gif)
    
    ## - list of filters
    snf <- SeqNameFilter("X")
    res <- ensembldb:::.processFilterParam(list(gif, snf))
    checkTrue(is(res, "AnnotationFilterList"))
    checkTrue(length(res) == 2)
    checkEquals(res[[1]], gif)
    checkEquals(res[[2]], snf)
    checkEquals(res@logOp, "&")
    
    ## - AnnotationFilterList
    afl <- AnnotationFilterList(gif, snf, logOp = "|")
    res <- ensembldb:::.processFilterParam(afl)
    checkTrue(is(res, "AnnotationFilterList"))
    checkEquals(afl, res)
    
    ## - filter expression
    res <- ensembldb:::.processFilterParam(gene_id != "BCL2" | seq_name == "X")
    checkTrue(is(res, "AnnotationFilterList"))
    checkEquals(res, AnnotationFilterList(gif, snf, logOp = "|"))
    res <- ensembldb:::.processFilterParam(gene_id != "BCL2")
    checkTrue(is(res, "AnnotationFilterList"))
    checkEquals(res, AnnotationFilterList(gif))

    ## - Errors
    checkException(ensembldb:::.processFilterParam())
    checkException(ensembldb:::.processFilterParam(4))
    checkException(ensembldb:::.processFilterParam(list(afl, "a")))
    checkException(ensembldb:::.processFilterParam("a"))
    checkException(ensembldb:::.processFilterParam(gene_bla == "14"))
}
