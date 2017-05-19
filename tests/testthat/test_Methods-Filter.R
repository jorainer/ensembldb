
test_that("ensDbColumn works", {
    smb <- SymbolFilter("a")
    expect_equal(unname(ensembldb:::ensDbColumn(smb)), "gene_name")
    expect_equal(unname(ensembldb:::ensDbColumn(smb, edb)), "gene.gene_name")
    expect_error(unname(ensembldb:::ensDbColumn(smb, edb, "tx")))
    expect_equal(unname(ensembldb:::ensDbColumn(smb, edb, "gene")),
                 "gene.gene_name")
    ##
    fl <- OnlyCodingTxFilter()
    expect_equal(ensembldb:::ensDbColumn(fl), "tx.tx_cds_seq_start")
    ## gene filters:
    fl <- GeneIdFilter("a")
    expect_equal(unname(ensembldb:::ensDbColumn(fl)), "gene_id")
    expect_equal(unname(ensembldb:::ensDbColumn(fl, edb)), "gene.gene_id")
    expect_equal(unname(ensembldb:::ensDbColumn(fl, edb, "tx")), "tx.gene_id")
    fl <- GenenameFilter("a")
    expect_equal(unname(ensembldb:::ensDbColumn(fl)), "gene_name")
    expect_equal(unname(ensembldb:::ensDbColumn(fl, edb)), "gene.gene_name")
    expect_error(unname(ensembldb:::ensDbColumn(fl, edb, "tx")))
    fl <- GeneStartFilter(123)
    expect_equal(unname(ensembldb:::ensDbColumn(fl)), "gene_seq_start")
    expect_equal(unname(ensembldb:::ensDbColumn(fl, edb)), "gene.gene_seq_start")
    expect_error(unname(ensembldb:::ensDbColumn(fl, edb, "tx")))
    fl <- GeneEndFilter(123)
    expect_equal(unname(ensembldb:::ensDbColumn(fl)), "gene_seq_end")
    expect_equal(unname(ensembldb:::ensDbColumn(fl, edb)), "gene.gene_seq_end")
    expect_error(unname(ensembldb:::ensDbColumn(fl, edb, "tx")))
    fl <- EntrezFilter(123)
    expect_equal(unname(ensembldb:::ensDbColumn(fl)), "entrezid")
    if (as.numeric(ensembldb:::dbSchemaVersion(edb)) > 1) {
        expect_equal(unname(ensembldb:::ensDbColumn(fl, edb)),
                     "entrezgene.entrezid")
    } else {
        expect_equal(unname(ensembldb:::ensDbColumn(fl, edb)), "gene.entrezid")
    }
    expect_error(unname(ensembldb:::ensDbColumn(fl, edb, "tx")))
    fl <- SeqNameFilter("a")
    expect_equal(unname(ensembldb:::ensDbColumn(fl)), "seq_name")
    expect_equal(unname(ensembldb:::ensDbColumn(fl, edb)), "gene.seq_name")
    expect_error(unname(ensembldb:::ensDbColumn(fl, edb, "tx")))
    fl <- SeqStrandFilter("a")
    expect_equal(unname(ensembldb:::ensDbColumn(fl)), "seq_strand")
    expect_equal(unname(ensembldb:::ensDbColumn(fl, edb)), "gene.seq_strand")
    expect_error(unname(ensembldb:::ensDbColumn(fl, edb, "tx")))
    ## tx filters:
    fl <- TxIdFilter("a")
    expect_equal(unname(ensembldb:::ensDbColumn(fl)), "tx_id")
    expect_equal(unname(ensembldb:::ensDbColumn(fl, edb)), "tx.tx_id")
    expect_error(unname(ensembldb:::ensDbColumn(fl, edb, "gene")))
    expect_equal(unname(ensembldb:::ensDbColumn(fl, edb, "protein")),
                 "protein.tx_id")
    fl <- TxNameFilter("a")
    expect_equal(unname(ensembldb:::ensDbColumn(fl)), "tx_id")
    expect_equal(unname(ensembldb:::ensDbColumn(fl, edb)), "tx.tx_id")
    expect_error(unname(ensembldb:::ensDbColumn(fl, edb, "gene")))
    fl <- TxBiotypeFilter("a")
    expect_equal(unname(ensembldb:::ensDbColumn(fl)), "tx_biotype")
    expect_equal(unname(ensembldb:::ensDbColumn(fl, edb)), "tx.tx_biotype")
    expect_error(unname(ensembldb:::ensDbColumn(fl, edb, "gene")))
    fl <- TxStartFilter(123)
    expect_equal(unname(ensembldb:::ensDbColumn(fl)), "tx_seq_start")
    expect_equal(unname(ensembldb:::ensDbColumn(fl, edb)), "tx.tx_seq_start")
    expect_error(unname(ensembldb:::ensDbColumn(fl, edb, "gene")))
    fl <- TxEndFilter(123)
    expect_equal(unname(ensembldb:::ensDbColumn(fl)), "tx_seq_end")
    expect_equal(unname(ensembldb:::ensDbColumn(fl, edb)), "tx.tx_seq_end")
    expect_error(unname(ensembldb:::ensDbColumn(fl, edb, "gene")))
    ## exon filters:
    fl <- ExonIdFilter(123)
    expect_equal(unname(ensembldb:::ensDbColumn(fl)), "exon_id")
    expect_equal(unname(ensembldb:::ensDbColumn(fl, edb)), "tx2exon.exon_id")
    expect_equal(unname(ensembldb:::ensDbColumn(fl, edb, "tx2exon")),
                 "tx2exon.exon_id")
    expect_equal(unname(ensembldb:::ensDbColumn(fl, edb, "exon")),
                 "exon.exon_id")
    expect_error(unname(ensembldb:::ensDbColumn(fl, edb, "gene")))
    fl <- ExonEndFilter(123)
    expect_equal(unname(ensembldb:::ensDbColumn(fl)), "exon_seq_end")
    expect_equal(unname(ensembldb:::ensDbColumn(fl, edb)), "exon.exon_seq_end")
    expect_error(unname(ensembldb:::ensDbColumn(fl, edb, "gene")))
    fl <- ExonStartFilter(123)
    expect_equal(unname(ensembldb:::ensDbColumn(fl)), "exon_seq_start")
    expect_equal(unname(ensembldb:::ensDbColumn(fl, edb)), "exon.exon_seq_start")
    expect_error(unname(ensembldb:::ensDbColumn(fl, edb, "gene")))
    fl <- ExonRankFilter(123)
    expect_equal(unname(ensembldb:::ensDbColumn(fl)), "exon_idx")
    expect_equal(unname(ensembldb:::ensDbColumn(fl, edb)), "tx2exon.exon_idx")
    expect_error(unname(ensembldb:::ensDbColumn(fl, edb, "gene")))
    ## protein filters:
    fl <- ProteinIdFilter("a")
    expect_equal(unname(ensembldb:::ensDbColumn(fl)), "protein_id")
    expect_equal(unname(ensembldb:::ensDbColumn(fl, edb)), "protein.protein_id")
    expect_equal(unname(ensembldb:::ensDbColumn(fl, edb, "uniprot")),
                 "uniprot.protein_id")
    expect_error(unname(ensembldb:::ensDbColumn(fl, edb, "gene")))
    fl <- UniprotFilter("a")
    expect_equal(unname(ensembldb:::ensDbColumn(fl)), "uniprot_id")
    expect_equal(unname(ensembldb:::ensDbColumn(fl, edb)), "uniprot.uniprot_id")
    expect_error(unname(ensembldb:::ensDbColumn(fl, edb, "gene")))
    fl <- ProtDomIdFilter("a")
    expect_equal(unname(ensembldb:::ensDbColumn(fl)), "protein_domain_id")
    expect_equal(unname(ensembldb:::ensDbColumn(fl, edb)),
                 "protein_domain.protein_domain_id")
    expect_error(unname(ensembldb:::ensDbColumn(fl, edb, "gene")))
    fl <- UniprotDbFilter("a")
    expect_equal(unname(ensembldb:::ensDbColumn(fl)), "uniprot_db")
    expect_equal(unname(ensembldb:::ensDbColumn(fl, edb)),
                 "uniprot.uniprot_db")
    expect_error(unname(ensembldb:::ensDbColumn(fl, edb, "gene")))
    fl <- UniprotMappingTypeFilter("a")
    expect_equal(unname(ensembldb:::ensDbColumn(fl)), "uniprot_mapping_type")
    expect_equal(unname(ensembldb:::ensDbColumn(fl, edb)),
                 "uniprot.uniprot_mapping_type")
    expect_error(unname(ensembldb:::ensDbColumn(fl, edb, "gene")))    
})

test_that("ensDbQuery works", {
    smb <- SymbolFilter("I'm a gene")
    expect_equal(ensembldb:::ensDbQuery(smb), "gene_name = 'I''m a gene'")
    expect_equal(ensembldb:::ensDbQuery(smb, edb), "gene.gene_name = 'I''m a gene'")
    expect_equal(ensembldb:::ensDbQuery(smb, edb, c("gene", "tx")),
                 "gene.gene_name = 'I''m a gene'")
    expect_error(ensembldb:::ensDbQuery(smb, edb, "tx"))
    smb <- SymbolFilter(c("a", "x"), condition = "!=")
    expect_equal(ensembldb:::ensDbQuery(smb), "gene_name not in ('a','x')")
    smb <- SymbolFilter(c("a", "x"), condition = "!=")
    expect_equal(ensembldb:::ensDbQuery(smb, edb),
                 "gene.gene_name not in ('a','x')")
    ## gene_id with tx table.
    fl <- GeneIdFilter("a")
    expect_equal(ensembldb:::ensDbQuery(fl), "gene_id = 'a'")
    expect_equal(ensembldb:::ensDbQuery(fl, edb), "gene.gene_id = 'a'")
    expect_equal(ensembldb:::ensDbQuery(fl, edb, "tx"), "tx.gene_id = 'a'")
    ## numeric filter(s)
    fl <- ExonRankFilter(21)
    expect_equal(ensembldb:::ensDbQuery(fl), "exon_idx = 21")
    expect_equal(ensembldb:::ensDbQuery(fl, edb), "tx2exon.exon_idx = 21")
    ##
    fl <- OnlyCodingTxFilter()
    expect_equal(ensembldb:::ensDbQuery(fl), "tx.tx_cds_seq_start is not null")
    ## 
    fL <- AnnotationFilterList(GeneIdFilter("a"),
                               TxBiotypeFilter("coding", condition = "!="),
                               GeneStartFilter(123, condition = "<"))
    res <- ensembldb:::ensDbQuery(fL)
    expect_equal(res, paste0("(gene_id = 'a' and tx_biotype != 'coding' ",
                             "and gene_seq_start < 123)"))
    res <- ensembldb:::ensDbQuery(fL, edb)
    expect_equal(res, paste0("(gene.gene_id = 'a' and tx.tx_biotype != ",
                             "'coding' and gene.gene_seq_start < 123)"))
    fL <- AnnotationFilterList(GeneIdFilter("a"),
                               TxBiotypeFilter("coding", condition = "!="),
                               TxStartFilter(123, condition = "<"))
    res <- ensembldb:::ensDbQuery(fL, edb)
    expect_equal(res, paste0("(gene.gene_id = 'a' and tx.tx_biotype != ",
                             "'coding' and tx.tx_seq_start < 123)"))
    res <- ensembldb:::ensDbQuery(fL, edb, c("tx", "gene", "exon"))
    expect_equal(res, paste0("(tx.gene_id = 'a' and tx.tx_biotype != ",
                             "'coding' and tx.tx_seq_start < 123)"))
    fl <- ProteinIdFilter("a")
    expect_equal(unname(ensembldb:::ensDbQuery(fl)), "protein_id = 'a'")
    expect_equal(unname(ensembldb:::ensDbQuery(fl, edb)),
                 "protein.protein_id = 'a'")
    expect_equal(unname(ensembldb:::ensDbQuery(fl, edb, "uniprot")),
                 "uniprot.protein_id = 'a'")    
    fl <- ProtDomIdFilter("a")
    expect_equal(ensembldb:::ensDbQuery(fl), "protein_domain_id = 'a'")
    expect_equal(ensembldb:::ensDbQuery(fl, edb),
                 "protein_domain.protein_domain_id = 'a'")
    fl <- UniprotDbFilter("a")
    expect_equal(ensembldb:::ensDbQuery(fl), "uniprot_db = 'a'")
    expect_equal(ensembldb:::ensDbQuery(fl, edb),
                 "uniprot.uniprot_db = 'a'")
    fl <- UniprotMappingTypeFilter("a")
    expect_equal(ensembldb:::ensDbQuery(fl), "uniprot_mapping_type = 'a'")
    expect_equal(ensembldb:::ensDbQuery(fl, edb),
                 "uniprot.uniprot_mapping_type = 'a'")
    ## Seq name and seq strand.
    fl <- SeqNameFilter("3")
    expect_equal(ensembldb:::ensDbQuery(fl), "seq_name = '3'")
    expect_equal(ensembldb:::ensDbQuery(fl, edb), "gene.seq_name = '3'")
    fl <- SeqNameFilter("chr3")
    expect_equal(ensembldb:::ensDbQuery(fl), "seq_name = 'chr3'")
    expect_equal(ensembldb:::ensDbQuery(fl, edb), "gene.seq_name = 'chr3'")
    seqlevelsStyle(edb) <- "UCSC"
    expect_equal(ensembldb:::ensDbQuery(fl), "seq_name = 'chr3'")
    expect_equal(ensembldb:::ensDbQuery(fl, edb), "gene.seq_name = '3'")
    seqlevelsStyle(edb) <- "Ensembl"
    fl <- SeqStrandFilter("+")
    expect_equal(ensembldb:::ensDbQuery(fl), "seq_strand = 1")
    expect_equal(ensembldb:::ensDbQuery(fl, edb), "gene.seq_strand = 1")
    fl <- SeqStrandFilter("+1")
    expect_equal(ensembldb:::ensDbQuery(fl), "seq_strand = 1")
    expect_equal(ensembldb:::ensDbQuery(fl, edb), "gene.seq_strand = 1")
    fl <- SeqStrandFilter("-")
    expect_equal(ensembldb:::ensDbQuery(fl), "seq_strand = -1")
    expect_equal(ensembldb:::ensDbQuery(fl, edb), "gene.seq_strand = -1")
    fl <- SeqStrandFilter("-1")
    expect_equal(ensembldb:::ensDbQuery(fl), "seq_strand = -1")
    expect_equal(ensembldb:::ensDbQuery(fl, edb), "gene.seq_strand = -1")
    ## GRangesFilter: see test_ensDb_for_GRangesFilter
})

test_that("ensDbQuery works for AnnotationFilterList", {
    gnf <- GenenameFilter("BCL2", condition = "!=")
    snf <- SeqNameFilter(4)
    ssf <- SeqStrandFilter("+")
    afl <- AnnotationFilterList(gnf, snf, ssf, logOp = c("|", "&"))
    Q <- ensembldb:::ensDbQuery(afl)
    expect_equal(Q, "(gene_name != 'BCL2' or seq_name = '4' and seq_strand = 1)")
    
    ## Nested AnnotationFilterLists.
    afl1 <- AnnotationFilterList(GenenameFilter("BCL2"),
                                 GenenameFilter("BCL2L11"), logOp = "|")
    afl2 <- AnnotationFilterList(afl1, SeqNameFilter(18))
    Q <- ensembldb:::ensDbQuery(afl2, db = edb)
    expect_equal(Q, paste0("((gene.gene_name = 'BCL2' or gene.gene_name = ",
                           "'BCL2L11') and gene.seq_name = '18')"))
    library(RSQLite)
    res <- dbGetQuery(dbconn(edb), paste0("select distinct gene_name from gene",
                                          " where ", Q))
    expect_equal(res$gene_name, "BCL2")
    res2 <- genes(edb,
                  filter = AnnotationFilterList(GenenameFilter(c("BCL2L11",
                                                                 "BCL2")),
                                                SeqNameFilter(18)))
    expect_equal(res$gene_name, res2$gene_name)
    ## Same with a GRangesFilter.
    grf <- GRangesFilter(GRanges(18, IRanges(60790600, 60790700)))
    afl2 <- AnnotationFilterList(afl1, grf)
    Q <- ensembldb:::ensDbQuery(afl2, db = edb)
    expect_equal(Q, paste0("((gene.gene_name = 'BCL2' or gene.gene_name = ",
                           "'BCL2L11') and (gene.gene_seq_start<=60790700",
                           " and gene.gene_seq_end>=60790600 and gene.seq_name",
                           "='18'))"))
    res <- dbGetQuery(dbconn(edb), paste0("select distinct gene_name from gene",
                                          " where ", Q))
    expect_equal(res$gene_name, "BCL2")    
})

test_that("ensDbQuery works for SeqNameFilter", {
    fl <- SeqNameFilter("3")
    res <- ensembldb:::ensDbQuery(fl)
    expect_equal(res, "seq_name = '3'")
    res <- ensembldb:::ensDbQuery(fl, edb)
    expect_equal(res, "gene.seq_name = '3'")
    fl <- SeqNameFilter("chr3")
    res <- ensembldb:::ensDbQuery(fl)
    expect_equal(res, "seq_name = 'chr3'")
    res <- ensembldb:::ensDbQuery(fl, edb)
    expect_equal(res, "gene.seq_name = 'chr3'")
    seqlevelsStyle(edb) <- "UCSC"
    res <- ensembldb:::ensDbQuery(fl)
    expect_equal(res, "seq_name = 'chr3'")
    res <- ensembldb:::ensDbQuery(fl, edb)
    expect_equal(res, "gene.seq_name = '3'")
    seqlevelsStyle(edb) <- "Ensembl"
})

test_that("ensDbQuery works for GRangesFilter", {
    gr <- GRanges(seqnames = "a",
                  ranges = IRanges(start = 1, end = 5))
    F <- GRangesFilter(value = gr, type = "within")
    expect_equal(unname(ensembldb:::ensDbColumn(F)), c("gene_seq_start",
                                                       "gene_seq_end",
                                                       "seq_name",
                                                       "seq_strand"))
    expect_equal(unname(ensembldb:::ensDbColumn(F, edb)),
                 c("gene.gene_seq_start", "gene.gene_seq_end",
                   "gene.seq_name", "gene.seq_strand"))

    expect_equal(ensembldb:::ensDbQuery(F),
                 paste0("(gene_seq_start>=1 and gene_seq_end",
                        "<=5 and seq_name='a')"))
    expect_equal(ensembldb:::ensDbQuery(F, edb),
                 paste0("(gene.gene_seq_start>=1 and gene.gene_seq_end",
                        "<=5 and gene.seq_name='a')"))
    F <- GRangesFilter(value = gr, type = "any")
    expect_equal(ensembldb:::ensDbQuery(F),
                 paste0("(gene_seq_start<=5 and gene_seq_end",
                        ">=1 and seq_name='a')"))
    expect_equal(ensembldb:::ensDbQuery(F, edb),
                 paste0("(gene.gene_seq_start<=5 and gene.gene_seq_end",
                        ">=1 and gene.seq_name='a')"))
    
    ## tx
    F <- GRangesFilter(value = gr, feature = "tx", type = "within")
    expect_equal(unname(ensembldb:::ensDbColumn(F)), c("tx_seq_start",
                                                       "tx_seq_end",
                                                       "seq_name",
                                                       "seq_strand"))
    expect_equal(unname(ensembldb:::ensDbColumn(F, edb)), c("tx.tx_seq_start",
                                                            "tx.tx_seq_end",
                                                            "gene.seq_name",
                                                            "gene.seq_strand"))
    ## exon
    F <- GRangesFilter(value = gr, feature = "exon", type = "within")
    expect_equal(unname(ensembldb:::ensDbColumn(F)), c("exon_seq_start",
                                                       "exon_seq_end",
                                                       "seq_name",
                                                       "seq_strand"))
    expect_equal(unname(ensembldb:::ensDbColumn(F, edb)), c("exon.exon_seq_start",
                                                            "exon.exon_seq_end",
                                                            "gene.seq_name",
                                                            "gene.seq_strand"))
    ## Check the buildWhere
    res <- ensembldb:::buildWhereForGRanges(F, columns = ensembldb:::ensDbColumn(F))
    expect_equal(res,"(exon_seq_start>=1 and exon_seq_end<=5 and seq_name='a')")
    res <- ensembldb:::ensDbQuery(F)
    expect_equal(res,"(exon_seq_start>=1 and exon_seq_end<=5 and seq_name='a')")
    res <- ensembldb:::ensDbQuery(F, edb)
    expect_equal(res, paste0("(exon.exon_seq_start>=1 and exon.exon_seq_end",
                             "<=5 and gene.seq_name='a')"))
})

test_that("buildWhereForGRanges works", {
    grng <- GRanges(seqname = "X", IRanges(start = 10, end = 100))
    ## start
    flt <- GRangesFilter(grng, type = "start")
    res <- ensembldb:::buildWhereForGRanges(
                           flt, columns = ensembldb:::ensDbColumn(flt))
    expect_equal(res, "(gene_seq_start=10 and seq_name='X')")
    res <- ensembldb:::buildWhereForGRanges(
                           flt, columns = ensembldb:::ensDbColumn(flt, db = edb))
    expect_equal(res, "(gene.gene_seq_start=10 and gene.seq_name='X')")
    ## end
    flt <- GRangesFilter(grng, type = "end")
    res <- ensembldb:::buildWhereForGRanges(
                           flt, columns = ensembldb:::ensDbColumn(flt))
    expect_equal(res, "(gene_seq_end=100 and seq_name='X')")
    ## equal
    flt <- GRangesFilter(grng, type = "equal")
    res <- ensembldb:::buildWhereForGRanges(
                           flt, columns = ensembldb:::ensDbColumn(flt))
    expect_equal(res, "(gene_seq_start=10 and gene_seq_end=100 and seq_name='X')")
    ## within
    flt <- GRangesFilter(grng, type = "within")
    res <- ensembldb:::buildWhereForGRanges(
                           flt, columns = ensembldb:::ensDbColumn(flt))
    expect_equal(res, "(gene_seq_start>=10 and gene_seq_end<=100 and seq_name='X')")
    ## any
    flt <- GRangesFilter(grng, type = "any")
    res <- ensembldb:::buildWhereForGRanges(
                           flt, columns = ensembldb:::ensDbColumn(flt))
    expect_equal(res, "(gene_seq_start<=100 and gene_seq_end>=10 and seq_name='X')")

    ## Same with a strand specified.
    grng <- GRanges(seqname = "X", IRanges(start = 10, end = 100), strand = "-")
    ## start
    flt <- GRangesFilter(grng, type = "start")
    res <- ensembldb:::buildWhereForGRanges(
                           flt, columns = ensembldb:::ensDbColumn(flt))
    expect_equal(res, "(gene_seq_start=10 and seq_name='X' and seq_strand = -1)")
    ## end
    flt <- GRangesFilter(grng, type = "end")
    res <- ensembldb:::buildWhereForGRanges(
                           flt, columns = ensembldb:::ensDbColumn(flt))
    expect_equal(res, "(gene_seq_end=100 and seq_name='X' and seq_strand = -1)")
    ## equal
    flt <- GRangesFilter(grng, type = "equal")
    res <- ensembldb:::buildWhereForGRanges(
                           flt, columns = ensembldb:::ensDbColumn(flt))
    expect_equal(res, paste0("(gene_seq_start=10 and gene_seq_end=100 and ",
                            "seq_name='X' and seq_strand = -1)"))
    ## within
    flt <- GRangesFilter(grng, type = "within")
    res <- ensembldb:::buildWhereForGRanges(
                           flt, columns = ensembldb:::ensDbColumn(flt))
    expect_equal(res, paste0("(gene_seq_start>=10 and gene_seq_end<=100 and ",
                            "seq_name='X' and seq_strand = -1)"))
    ## any
    flt <- GRangesFilter(grng, type = "any")
    res <- ensembldb:::buildWhereForGRanges(
                           flt, columns = ensembldb:::ensDbColumn(flt))
    expect_equal(res, paste0("(gene_seq_start<=100 and gene_seq_end>=10 and ",
                            "seq_name='X' and seq_strand = -1)"))
})

test_that("ensDbColumn works with AnnotationFilterList", {
    afl <- AnnotationFilterList(GeneIdFilter(123), SeqNameFilter(3))
    afl2 <- AnnotationFilterList(afl, SeqNameFilter(5))
    res <- ensembldb:::ensDbColumn(afl2)
    expect_equal(res, c("gene_id", "seq_name"))
})

############################################################
## Using protein data based filters.
test_that("ProteinIdFilter works", {
    pf <- ProteinIdFilter("ABC")
    expect_equal(value(pf), "ABC")
    expect_equal(field(pf), "protein_id")
    expect_equal(ensembldb:::ensDbQuery(pf), "protein_id = 'ABC'")
    if (hasProteinData(edb)) {
        expect_equal(ensembldb:::ensDbColumn(pf, edb), "protein.protein_id")
        expect_equal(ensembldb:::ensDbColumn(pf, edb,
                                             with.tables = "protein_domain"),
                    "protein_domain.protein_id")
        expect_equal(ensembldb:::ensDbColumn(pf, edb, with.tables = "uniprot"),
                     "uniprot.protein_id")
        expect_equal(ensembldb:::ensDbQuery(pf, edb), "protein.protein_id = 'ABC'")
        expect_equal(ensembldb:::ensDbQuery(pf, edb, with.tables = "uniprot"),
                    "uniprot.protein_id = 'ABC'")
    } else {
        expect_error(ensembldb:::ensDbColumn(pf, edb))
        expect_error(ensembldb:::ensDbQuery(pf, edb))
        expect_error(ensembldb:::ensDbColumn(pf, edb, with.tables = "uniprot"))
        expect_error(ensembldb:::ensDbQuery(pf, edb, with.tables = "uniprot"))
    }
    pf <- ProteinIdFilter(c("A", "B"))
    expect_equal(ensembldb:::ensDbQuery(pf), "protein_id in ('A','B')")
    expect_error(ProteinIdFilter("B", condition = ">"))
})

test_that("UniprotFilter works", {
    pf <- UniprotFilter("ABC")
    expect_equal(value(pf), "ABC")
    expect_equal(field(pf), "uniprot")
    expect_equal(unname(ensembldb:::ensDbColumn(pf)), "uniprot_id")
    expect_equal(ensembldb:::ensDbQuery(pf), "uniprot_id = 'ABC'")
    if (hasProteinData(edb)) {
        expect_equal(ensembldb:::ensDbColumn(pf, edb), "uniprot.uniprot_id")
        expect_equal(ensembldb:::ensDbColumn(pf, edb, with.tables = "uniprot"),
                    "uniprot.uniprot_id")
        expect_equal(ensembldb:::ensDbQuery(pf, edb), "uniprot.uniprot_id = 'ABC'")
        expect_equal(ensembldb:::ensDbQuery(pf, edb, with.tables = "uniprot"),
                    "uniprot.uniprot_id = 'ABC'")
    } else {
        expect_error(ensembldb:::ensDbColumn(pf, edb))
        expect_error(ensembldb:::ensDbQuery(pf, edb))
        expect_error(ensembldb:::ensDbColumn(pf, edb, with.tables = "uniprot"))
        expect_error(ensembldb:::ensDbQuery(pf, edb, with.tables = "uniprot"))
    }
    pf <- UniprotFilter(c("A", "B"))
    expect_equal(ensembldb:::ensDbQuery(pf), "uniprot_id in ('A','B')")
    expect_error(UniprotFilter("B", condition = ">"))
})

test_that("ProtDomIdFilter works", {
    pf <- ProtDomIdFilter("ABC")
    expect_equal(value(pf), "ABC")
    expect_equal(field(pf), "prot_dom_id")
    expect_equal(unname(ensembldb:::ensDbColumn(pf)), "protein_domain_id")
    expect_equal(ensembldb:::ensDbQuery(pf), "protein_domain_id = 'ABC'")
    if (hasProteinData(edb)) {
        expect_equal(ensembldb:::ensDbColumn(pf, edb),
                     "protein_domain.protein_domain_id")
        expect_equal(ensembldb:::ensDbColumn(pf, edb,
                                             with.tables = "protein_domain"),
                    "protein_domain.protein_domain_id")
        expect_equal(ensembldb:::ensDbQuery(pf, edb),
                     "protein_domain.protein_domain_id = 'ABC'")
        expect_equal(ensembldb:::ensDbQuery(pf, edb,
                                            with.tables = "protein_domain"),
                    "protein_domain.protein_domain_id = 'ABC'")
    } else {
        expect_error(ensembldb:::ensDbColumn(pf, edb))
        expect_error(ensembldb:::ensDbQuery(pf, edb))
        expect_error(ensembldb:::ensDbColumn(pf, edb,
                                             with.tables = "protein_domain"))
        expect_error(ensembldb:::ensDbQuery(pf, edb,
                                            with.tables = "protein_domain"))
    }
    pf <- ProtDomIdFilter(c("A", "B"))
    expect_equal(ensembldb:::ensDbQuery(pf), "protein_domain_id in ('A','B')")
    expect_error(ProtDomIdFilter("B", condition = ">"))
})

test_that("UniprotDbFilter works", {
    pf <- UniprotDbFilter("ABC")
    expect_equal(value(pf), "ABC")
    expect_equal(field(pf), "uniprot_db")
    expect_equal(ensembldb:::ensDbQuery(pf), "uniprot_db = 'ABC'")
    if (hasProteinData(edb)) {
        expect_equal(ensembldb:::ensDbColumn(pf, edb), "uniprot.uniprot_db")
        expect_equal(ensembldb:::ensDbColumn(pf, edb, with.tables = "uniprot"),
                    "uniprot.uniprot_db")
        expect_equal(ensembldb:::ensDbQuery(pf, edb), "uniprot.uniprot_db = 'ABC'")
        expect_equal(ensembldb:::ensDbQuery(pf, edb, with.tables = "uniprot"),
                    "uniprot.uniprot_db = 'ABC'")
    } else {
        expect_error(ensembldb:::ensDbColumn(pf, edb))
        expect_error(ensembldb:::ensDbQuery(pf, edb))
        expect_error(ensembldb:::ensDbColumn(pf, edb, with.tables = "uniprot"))
        expect_error(ensembldb:::ensDbQuery(pf, edb, with.tables = "uniprot"))
    }
    pf <- UniprotDbFilter(c("A", "B"))
    expect_equal(ensembldb:::ensDbQuery(pf), "uniprot_db in ('A','B')")
    expect_error(UniprotDbFilter("B", condition = ">"))
})

test_that("UniprotMappingTypeFilter works", {
    pf <- UniprotMappingTypeFilter("ABC")
    expect_equal(value(pf), "ABC")
    expect_equal(field(pf), "uniprot_mapping_type")
    expect_equal(ensembldb:::ensDbQuery(pf), "uniprot_mapping_type = 'ABC'")
    if (hasProteinData(edb)) {
        expect_equal(ensembldb:::ensDbColumn(pf, edb),
                     "uniprot.uniprot_mapping_type")
        expect_equal(ensembldb:::ensDbColumn(pf, edb, with.tables = "uniprot"),
                    "uniprot.uniprot_mapping_type")
        expect_equal(ensembldb:::ensDbQuery(pf, edb),
                     "uniprot.uniprot_mapping_type = 'ABC'")
        expect_equal(ensembldb:::ensDbQuery(pf, edb, with.tables = "uniprot"),
                    "uniprot.uniprot_mapping_type = 'ABC'")
    } else {
        expect_error(ensembldb:::ensDbColumn(pf, edb))
        expect_error(ensembldb:::ensDbQuery(pf, edb))
        expect_error(ensembldb:::ensDbColumn(pf, edb, with.tables = "uniprot"))
        expect_error(ensembldb:::ensDbQuery(pf, edb, with.tables = "uniprot"))
    }
    pf <- UniprotMappingTypeFilter(c("A", "B"))
    expect_equal(ensembldb:::ensDbQuery(pf), "uniprot_mapping_type in ('A','B')")
    expect_error(UniprotMappingTypeFilter("B", condition = ">"))
})
