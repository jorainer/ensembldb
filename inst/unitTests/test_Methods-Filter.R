## Tests for methods in Methods-Filter.R
library(RUnit)
library(EnsDb.Hsapiens.v75)
edb <- EnsDb.Hsapiens.v75

test_supportedFilters <- function() {
    res <- supportedFilters(edb)
    if (!hasProteinData(edb))
        expect_identical(length(res), 19)
    else 
        expect_identical(length(res), 24)
}

test_ensDbColumn <- function() {
    smb <- SymbolFilter("a")
    checkEquals(unname(ensembldb:::ensDbColumn(smb)), "gene_name")
    checkEquals(unname(ensembldb:::ensDbColumn(smb, edb)), "gene.gene_name")
    checkException(unname(ensembldb:::ensDbColumn(smb, edb, "tx")))
    checkEquals(unname(ensembldb:::ensDbColumn(smb, edb, "gene")),
                "gene.gene_name")
    ##
    fl <- OnlyCodingTxFilter()
    checkEquals(ensembldb:::ensDbColumn(fl), "tx.tx_cds_seq_start")
    ## gene filters:
    fl <- GeneIdFilter("a")
    checkEquals(unname(ensembldb:::ensDbColumn(fl)), "gene_id")
    checkEquals(unname(ensembldb:::ensDbColumn(fl, edb)), "gene.gene_id")
    checkEquals(unname(ensembldb:::ensDbColumn(fl, edb, "tx")), "tx.gene_id")
    fl <- GenenameFilter("a")
    checkEquals(unname(ensembldb:::ensDbColumn(fl)), "gene_name")
    checkEquals(unname(ensembldb:::ensDbColumn(fl, edb)), "gene.gene_name")
    checkException(unname(ensembldb:::ensDbColumn(fl, edb, "tx")))
    fl <- GeneStartFilter(123)
    checkEquals(unname(ensembldb:::ensDbColumn(fl)), "gene_seq_start")
    checkEquals(unname(ensembldb:::ensDbColumn(fl, edb)), "gene.gene_seq_start")
    checkException(unname(ensembldb:::ensDbColumn(fl, edb, "tx")))
    fl <- GeneEndFilter(123)
    checkEquals(unname(ensembldb:::ensDbColumn(fl)), "gene_seq_end")
    checkEquals(unname(ensembldb:::ensDbColumn(fl, edb)), "gene.gene_seq_end")
    checkException(unname(ensembldb:::ensDbColumn(fl, edb, "tx")))
    fl <- EntrezFilter(123)
    checkEquals(unname(ensembldb:::ensDbColumn(fl)), "entrezid")
    checkEquals(unname(ensembldb:::ensDbColumn(fl, edb)), "gene.entrezid")
    checkException(unname(ensembldb:::ensDbColumn(fl, edb, "tx")))
    fl <- SeqNameFilter("a")
    checkEquals(unname(ensembldb:::ensDbColumn(fl)), "seq_name")
    checkEquals(unname(ensembldb:::ensDbColumn(fl, edb)), "gene.seq_name")
    checkException(unname(ensembldb:::ensDbColumn(fl, edb, "tx")))
    fl <- SeqStrandFilter("a")
    checkEquals(unname(ensembldb:::ensDbColumn(fl)), "seq_strand")
    checkEquals(unname(ensembldb:::ensDbColumn(fl, edb)), "gene.seq_strand")
    checkException(unname(ensembldb:::ensDbColumn(fl, edb, "tx")))
    ## tx filters:
    fl <- TxIdFilter("a")
    checkEquals(unname(ensembldb:::ensDbColumn(fl)), "tx_id")
    checkEquals(unname(ensembldb:::ensDbColumn(fl, edb)), "tx.tx_id")
    checkException(unname(ensembldb:::ensDbColumn(fl, edb, "gene")))
    checkEquals(unname(ensembldb:::ensDbColumn(fl, edb, "protein")),
                "protein.tx_id")
    fl <- TxNameFilter("a")
    checkEquals(unname(ensembldb:::ensDbColumn(fl)), "tx_id")
    checkEquals(unname(ensembldb:::ensDbColumn(fl, edb)), "tx.tx_id")
    checkException(unname(ensembldb:::ensDbColumn(fl, edb, "gene")))
    fl <- TxBiotypeFilter("a")
    checkEquals(unname(ensembldb:::ensDbColumn(fl)), "tx_biotype")
    checkEquals(unname(ensembldb:::ensDbColumn(fl, edb)), "tx.tx_biotype")
    checkException(unname(ensembldb:::ensDbColumn(fl, edb, "gene")))
    fl <- TxStartFilter(123)
    checkEquals(unname(ensembldb:::ensDbColumn(fl)), "tx_seq_start")
    checkEquals(unname(ensembldb:::ensDbColumn(fl, edb)), "tx.tx_seq_start")
    checkException(unname(ensembldb:::ensDbColumn(fl, edb, "gene")))
    fl <- TxEndFilter(123)
    checkEquals(unname(ensembldb:::ensDbColumn(fl)), "tx_seq_end")
    checkEquals(unname(ensembldb:::ensDbColumn(fl, edb)), "tx.tx_seq_end")
    checkException(unname(ensembldb:::ensDbColumn(fl, edb, "gene")))
    ## exon filters:
    fl <- ExonIdFilter(123)
    checkEquals(unname(ensembldb:::ensDbColumn(fl)), "exon_id")
    checkEquals(unname(ensembldb:::ensDbColumn(fl, edb)), "exon.exon_id")
    checkEquals(unname(ensembldb:::ensDbColumn(fl, edb, "tx2exon")),
                "tx2exon.exon_id")
    checkException(unname(ensembldb:::ensDbColumn(fl, edb, "gene")))
    fl <- ExonEndFilter(123)
    checkEquals(unname(ensembldb:::ensDbColumn(fl)), "exon_seq_end")
    checkEquals(unname(ensembldb:::ensDbColumn(fl, edb)), "exon.exon_seq_end")
    checkException(unname(ensembldb:::ensDbColumn(fl, edb, "gene")))
    fl <- ExonStartFilter(123)
    checkEquals(unname(ensembldb:::ensDbColumn(fl)), "exon_seq_start")
    checkEquals(unname(ensembldb:::ensDbColumn(fl, edb)), "exon.exon_seq_start")
    checkException(unname(ensembldb:::ensDbColumn(fl, edb, "gene")))
    fl <- ExonRankFilter(123)
    checkEquals(unname(ensembldb:::ensDbColumn(fl)), "exon_idx")
    checkEquals(unname(ensembldb:::ensDbColumn(fl, edb)), "tx2exon.exon_idx")
    checkException(unname(ensembldb:::ensDbColumn(fl, edb, "gene")))
    ## protein filters:
    fl <- ProteinIdFilter("a")
    checkEquals(unname(ensembldb:::ensDbColumn(fl)), "protein_id")
    checkEquals(unname(ensembldb:::ensDbColumn(fl, edb)), "protein.protein_id")
    checkEquals(unname(ensembldb:::ensDbColumn(fl, edb, "uniprot")),
                "uniprot.protein_id")
    checkException(unname(ensembldb:::ensDbColumn(fl, edb, "gene")))
    fl <- UniprotFilter("a")
    checkEquals(unname(ensembldb:::ensDbColumn(fl)), "uniprot_id")
    checkEquals(unname(ensembldb:::ensDbColumn(fl, edb)), "uniprot.uniprot_id")
    checkException(unname(ensembldb:::ensDbColumn(fl, edb, "gene")))
    fl <- ProtDomIdFilter("a")
    checkEquals(unname(ensembldb:::ensDbColumn(fl)), "protein_domain_id")
    checkEquals(unname(ensembldb:::ensDbColumn(fl, edb)),
                "protein_domain.protein_domain_id")
    checkException(unname(ensembldb:::ensDbColumn(fl, edb, "gene")))
    fl <- UniprotDbFilter("a")
    checkEquals(unname(ensembldb:::ensDbColumn(fl)), "uniprot_db")
    checkEquals(unname(ensembldb:::ensDbColumn(fl, edb)),
                "uniprot.uniprot_db")
    checkException(unname(ensembldb:::ensDbColumn(fl, edb, "gene")))
    fl <- UniprotMappingTypeFilter("a")
    checkEquals(unname(ensembldb:::ensDbColumn(fl)), "uniprot_mapping_type")
    checkEquals(unname(ensembldb:::ensDbColumn(fl, edb)),
                "uniprot.uniprot_mapping_type")
    checkException(unname(ensembldb:::ensDbColumn(fl, edb, "gene")))    
}

test_ensDbQuery <- function() {
    smb <- SymbolFilter("I'm a gene")
    checkEquals(ensembldb:::ensDbQuery(smb), "gene_name = 'I''m a gene'")
    checkEquals(ensembldb:::ensDbQuery(smb, edb), "gene.gene_name = 'I''m a gene'")
    checkEquals(ensembldb:::ensDbQuery(smb, edb, c("gene", "tx")),
                "gene.gene_name = 'I''m a gene'")
    checkException(ensembldb:::ensDbQuery(smb, edb, "tx"))
    smb <- SymbolFilter(c("a", "x"), condition = "!=")
    checkEquals(ensembldb:::ensDbQuery(smb), "gene_name not in ('a','x')")
    smb <- SymbolFilter(c("a", "x"), condition = "!=")
    checkEquals(ensembldb:::ensDbQuery(smb, edb),
                "gene.gene_name not in ('a','x')")
    ## gene_id with tx table.
    fl <- GeneIdFilter("a")
    checkEquals(ensembldb:::ensDbQuery(fl), "gene_id = 'a'")
    checkEquals(ensembldb:::ensDbQuery(fl, edb), "gene.gene_id = 'a'")
    checkEquals(ensembldb:::ensDbQuery(fl, edb, "tx"), "tx.gene_id = 'a'")
    ## numeric filter(s)
    fl <- ExonRankFilter(21)
    checkEquals(ensembldb:::ensDbQuery(fl), "exon_idx = 21")
    checkEquals(ensembldb:::ensDbQuery(fl, edb), "tx2exon.exon_idx = 21")
    ##
    fl <- OnlyCodingTxFilter()
    checkEquals(ensembldb:::ensDbQuery(fl), "tx.tx_cds_seq_start is not null")
    ## for a list.
    fL <- list(GeneIdFilter("a"), TxBiotypeFilter("coding", condition = "!="),
               GeneStartFilter(123, condition = "<"))
    res <- ensembldb:::ensDbQuery(fL)
    checkEquals(res, paste0(" where gene_id = 'a' and tx_biotype != 'coding' ",
                            "and gene_seq_start < 123"))
    res <- ensembldb:::ensDbQuery(fL, edb)
    checkEquals(res, paste0(" where gene.gene_id = 'a' and tx.tx_biotype != ",
                            "'coding' and gene.gene_seq_start < 123"))
    fL <- list(GeneIdFilter("a"), TxBiotypeFilter("coding", condition = "!="),
               TxStartFilter(123, condition = "<"))
    res <- ensembldb:::ensDbQuery(fL, edb)
    checkEquals(res, paste0(" where gene.gene_id = 'a' and tx.tx_biotype != ",
                            "'coding' and tx.tx_seq_start < 123"))
    res <- ensembldb:::ensDbQuery(fL, edb, c("tx", "gene", "exon"))
    checkEquals(res, paste0(" where tx.gene_id = 'a' and tx.tx_biotype != ",
                            "'coding' and tx.tx_seq_start < 123"))
    fl <- ProteinIdFilter("a")
    checkEquals(unname(ensembldb:::ensDbQuery(fl)), "protein_id = 'a'")
    checkEquals(unname(ensembldb:::ensDbQuery(fl, edb)),
                "protein.protein_id = 'a'")
    checkEquals(unname(ensembldb:::ensDbQuery(fl, edb, "uniprot")),
                "uniprot.protein_id = 'a'")    
    fl <- ProtDomIdFilter("a")
    checkEquals(ensembldb:::ensDbQuery(fl), "protein_domain_id = 'a'")
    checkEquals(ensembldb:::ensDbQuery(fl, edb),
                "protein_domain.protein_domain_id = 'a'")
    fl <- UniprotDbFilter("a")
    checkEquals(ensembldb:::ensDbQuery(fl), "uniprot_db = 'a'")
    checkEquals(ensembldb:::ensDbQuery(fl, edb),
                "uniprot.uniprot_db = 'a'")
    fl <- UniprotMappingTypeFilter("a")
    checkEquals(ensembldb:::ensDbQuery(fl), "uniprot_mapping_type = 'a'")
    checkEquals(ensembldb:::ensDbQuery(fl, edb),
                "uniprot.uniprot_mapping_type = 'a'")
    ## Seq name and seq strand.
    fl <- SeqNameFilter("3")
    checkEquals(ensembldb:::ensDbQuery(fl), "seq_name = '3'")
    checkEquals(ensembldb:::ensDbQuery(fl, edb), "gene.seq_name = '3'")
    fl <- SeqNameFilter("chr3")
    checkEquals(ensembldb:::ensDbQuery(fl), "seq_name = 'chr3'")
    checkEquals(ensembldb:::ensDbQuery(fl, edb), "gene.seq_name = 'chr3'")
    seqlevelsStyle(edb) <- "UCSC"
    checkEquals(ensembldb:::ensDbQuery(fl), "seq_name = 'chr3'")
    checkEquals(ensembldb:::ensDbQuery(fl, edb), "gene.seq_name = '3'")
    seqlevelsStyle(edb) <- "Ensembl"
    fl <- SeqStrandFilter("+")
    checkEquals(ensembldb:::ensDbQuery(fl), "seq_strand = 1")
    checkEquals(ensembldb:::ensDbQuery(fl, edb), "gene.seq_strand = 1")
    fl <- SeqStrandFilter("+1")
    checkEquals(ensembldb:::ensDbQuery(fl), "seq_strand = 1")
    checkEquals(ensembldb:::ensDbQuery(fl, edb), "gene.seq_strand = 1")
    fl <- SeqStrandFilter("-")
    checkEquals(ensembldb:::ensDbQuery(fl), "seq_strand = -1")
    checkEquals(ensembldb:::ensDbQuery(fl, edb), "gene.seq_strand = -1")
    fl <- SeqStrandFilter("-1")
    checkEquals(ensembldb:::ensDbQuery(fl), "seq_strand = -1")
    checkEquals(ensembldb:::ensDbQuery(fl, edb), "gene.seq_strand = -1")
    ## GRangesFilter: see test_ensDb_for_GRangesFilter
}

test_ensDbQuery_AnnotationFilterList <- function() {
    gnf <- GenenameFilter("BCL2", condition = "!=")
    snf <- SeqNameFilter(4)
    ssf <- SeqStrandFilter("+")
    afl <- AnnotationFilterList(gnf, snf, ssf, logOp = c("|", "&"))
    Q <- ensembldb:::ensDbQuery(afl)
    checkEquals(Q, "(gene_name != 'BCL2' or seq_name = '4' and seq_strand = 1)")
    
    ## Nested AnnotationFilterLists.
    afl1 <- AnnotationFilterList(GenenameFilter("BCL2"),
                                 GenenameFilter("BCL2L11"), logOp = "|")
    afl2 <- AnnotationFilterList(afl1, SeqNameFilter(18))
    Q <- ensembldb:::ensDbQuery(afl2, db = edb)
    checkEquals(Q, paste0("((gene.gene_name = 'BCL2' or gene.gene_name = ",
                          "'BCL2L11') and gene.seq_name = '18')"))
    library(RSQLite)
    res <- dbGetQuery(dbconn(edb), paste0("select distinct gene_name from gene",
                                          " where ", Q))
    checkEquals(res$gene_name, "BCL2")
    res2 <- genes(edb,
                  filter = AnnotationFilterList(GenenameFilter(c("BCL2L11",
                                                                 "BCL2")),
                                                SeqNameFilter(18)))
    checkEquals(res$gene_name, res2$gene_name)
    ## Same with a GRangesFilter.
    grf <- GRangesFilter(GRanges(18, IRanges(60790600, 60790700)))
    afl2 <- AnnotationFilterList(afl1, grf)
    Q <- ensembldb:::ensDbQuery(afl2, db = edb)
    checkEquals(Q, paste0("((gene.gene_name = 'BCL2' or gene.gene_name = ",
                          "'BCL2L11') and (gene.gene_seq_start<=60790700",
                          " and gene.gene_seq_end>=60790600 and gene.seq_name",
                          "='18'))"))
    res <- dbGetQuery(dbconn(edb), paste0("select distinct gene_name from gene",
                                          " where ", Q))
    checkEquals(res$gene_name, "BCL2")    
}

test_value_for_SeqNameFilter <- function() {
    fl <- SeqNameFilter("3")
    checkEquals(value(fl), "3")
    checkEquals(value(fl, edb), "3")
    fl <- SeqNameFilter("chr3")
    checkEquals(value(fl), "chr3")
    checkEquals(value(fl, edb), "chr3")
    seqlevelsStyle(edb) <- "UCSC"
    checkEquals(value(fl, edb), "3")
    seqlevelsStyle(edb) <- "Ensembl"
}

test_strand2num <- function(x) {
    checkEquals(ensembldb:::strand2num("+"), 1)
    checkEquals(ensembldb:::strand2num("+1"), 1)
    checkEquals(ensembldb:::strand2num("-"), -1)
    checkEquals(ensembldb:::strand2num("-1"), -1)
    checkEquals(ensembldb:::strand2num(1), 1)
    checkEquals(ensembldb:::strand2num(5), 1)
    checkEquals(ensembldb:::strand2num(-1), -1)
    checkEquals(ensembldb:::strand2num(-5), -1)
    checkException(ensembldb:::strand2num("a"))
}

test_num2strand <- function(x) {
    checkEquals(ensembldb:::num2strand(1), "+")
    checkEquals(ensembldb:::num2strand(-1), "-")
}

test_ensDb_for_GRangesFilter <- function() {
    gr <- GRanges(seqnames = "a",
                  ranges = IRanges(start = 1, end = 5))
    F <- GRangesFilter(value = gr, type = "within")
    checkEquals(unname(ensembldb:::ensDbColumn(F)), c("gene_seq_start",
                                                      "gene_seq_end",
                                                      "seq_name",
                                                      "seq_strand"))
    checkEquals(unname(ensembldb:::ensDbColumn(F, edb)),
                c("gene.gene_seq_start", "gene.gene_seq_end",
                  "gene.seq_name", "gene.seq_strand"))

    checkEquals(ensembldb:::ensDbQuery(F),
                paste0("(gene_seq_start>=1 and gene_seq_end",
                       "<=5 and seq_name='a')"))
    checkEquals(ensembldb:::ensDbQuery(F, edb),
                paste0("(gene.gene_seq_start>=1 and gene.gene_seq_end",
                       "<=5 and gene.seq_name='a')"))
    F <- GRangesFilter(value = gr, type = "any")
    checkEquals(ensembldb:::ensDbQuery(F),
                paste0("(gene_seq_start<=5 and gene_seq_end",
                       ">=1 and seq_name='a')"))
    checkEquals(ensembldb:::ensDbQuery(F, edb),
                paste0("(gene.gene_seq_start<=5 and gene.gene_seq_end",
                       ">=1 and gene.seq_name='a')"))
    
    ## tx
    F <- GRangesFilter(value = gr, feature = "tx", type = "within")
    checkEquals(unname(ensembldb:::ensDbColumn(F)), c("tx_seq_start",
                                                      "tx_seq_end",
                                                      "seq_name",
                                                      "seq_strand"))
    checkEquals(unname(ensembldb:::ensDbColumn(F, edb)), c("tx.tx_seq_start",
                                                           "tx.tx_seq_end",
                                                           "gene.seq_name",
                                                           "gene.seq_strand"))
    ## exon
    F <- GRangesFilter(value = gr, feature = "exon", type = "within")
    checkEquals(unname(ensembldb:::ensDbColumn(F)), c("exon_seq_start",
                                                      "exon_seq_end",
                                                      "seq_name",
                                                      "seq_strand"))
    checkEquals(unname(ensembldb:::ensDbColumn(F, edb)), c("exon.exon_seq_start",
                                                           "exon.exon_seq_end",
                                                           "gene.seq_name",
                                                           "gene.seq_strand"))
    ## Check the buildWhere
    res <- ensembldb:::buildWhereForGRanges(F, columns = ensembldb:::ensDbColumn(F))
    checkEquals(res,"(exon_seq_start>=1 and exon_seq_end<=5 and seq_name='a')")
    res <- ensembldb:::ensDbQuery(F)
    checkEquals(res,"(exon_seq_start>=1 and exon_seq_end<=5 and seq_name='a')")
    res <- ensembldb:::ensDbQuery(F, edb)
    checkEquals(res, paste0("(exon.exon_seq_start>=1 and exon.exon_seq_end",
                            "<=5 and gene.seq_name='a')"))
}

test_buildWhereForGRanges <- function() {
    grng <- GRanges(seqname = "X", IRanges(start = 10, end = 100))
    ## start
    flt <- GRangesFilter(grng, type = "start")
    res <- ensembldb:::buildWhereForGRanges(
                           flt, columns = ensembldb:::ensDbColumn(flt))
    checkEquals(res, "(gene_seq_start=10 and seq_name='X')")
    res <- ensembldb:::buildWhereForGRanges(
                           flt, columns = ensembldb:::ensDbColumn(flt, db = edb))
    checkEquals(res, "(gene.gene_seq_start=10 and gene.seq_name='X')")
    ## end
    flt <- GRangesFilter(grng, type = "end")
    res <- ensembldb:::buildWhereForGRanges(
                           flt, columns = ensembldb:::ensDbColumn(flt))
    checkEquals(res, "(gene_seq_end=100 and seq_name='X')")
    ## equal
    flt <- GRangesFilter(grng, type = "equal")
    res <- ensembldb:::buildWhereForGRanges(
                           flt, columns = ensembldb:::ensDbColumn(flt))
    checkEquals(res, "(gene_seq_start=10 and gene_seq_end=100 and seq_name='X')")
    ## within
    flt <- GRangesFilter(grng, type = "within")
    res <- ensembldb:::buildWhereForGRanges(
                           flt, columns = ensembldb:::ensDbColumn(flt))
    checkEquals(res, "(gene_seq_start>=10 and gene_seq_end<=100 and seq_name='X')")
    ## any
    flt <- GRangesFilter(grng, type = "any")
    res <- ensembldb:::buildWhereForGRanges(
                           flt, columns = ensembldb:::ensDbColumn(flt))
    checkEquals(res, "(gene_seq_start<=100 and gene_seq_end>=10 and seq_name='X')")

    ## Same with a strand specified.
    grng <- GRanges(seqname = "X", IRanges(start = 10, end = 100), strand = "-")
    ## start
    flt <- GRangesFilter(grng, type = "start")
    res <- ensembldb:::buildWhereForGRanges(
                           flt, columns = ensembldb:::ensDbColumn(flt))
    checkEquals(res, "(gene_seq_start=10 and seq_name='X' and seq_strand = -1)")
    ## end
    flt <- GRangesFilter(grng, type = "end")
    res <- ensembldb:::buildWhereForGRanges(
                           flt, columns = ensembldb:::ensDbColumn(flt))
    checkEquals(res, "(gene_seq_end=100 and seq_name='X' and seq_strand = -1)")
    ## equal
    flt <- GRangesFilter(grng, type = "equal")
    res <- ensembldb:::buildWhereForGRanges(
                           flt, columns = ensembldb:::ensDbColumn(flt))
    checkEquals(res, paste0("(gene_seq_start=10 and gene_seq_end=100 and ",
                            "seq_name='X' and seq_strand = -1)"))
    ## within
    flt <- GRangesFilter(grng, type = "within")
    res <- ensembldb:::buildWhereForGRanges(
                           flt, columns = ensembldb:::ensDbColumn(flt))
    checkEquals(res, paste0("(gene_seq_start>=10 and gene_seq_end<=100 and ",
                            "seq_name='X' and seq_strand = -1)"))
    ## any
    flt <- GRangesFilter(grng, type = "any")
    res <- ensembldb:::buildWhereForGRanges(
                           flt, columns = ensembldb:::ensDbColumn(flt))
    checkEquals(res, paste0("(gene_seq_start<=100 and gene_seq_end>=10 and ",
                            "seq_name='X' and seq_strand = -1)"))
}
