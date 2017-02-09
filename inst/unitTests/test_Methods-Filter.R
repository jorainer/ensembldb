## Tests for methods in Methods-Filter.R
library(RUnit)
library(EnsDb.Hsapiens.v75)
edb <- EnsDb.Hsapiens.v75

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
    ## TODO GRangesFilter
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
}
