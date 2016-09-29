############################################################
## Tests for the dbhelpers.R functions
library(EnsDb.Hsapiens.v75)
edb <- EnsDb.Hsapiens.v75

test_addRequiredTables <- function() {
    have <- c("exon", "gene")
    need <- c("exon", "gene", "tx2exon", "tx")
    checkEquals(sort(need), sort(ensembldb:::addRequiredTables(edb, have)))

    have <- c("exon", "chromosome")
    need <- c("exon", "tx2exon", "tx", "gene", "chromosome")
    checkEquals(sort(need), sort(ensembldb:::addRequiredTables(edb, have)))

    have <- c("chromosome", "tx")
    need <- c("chromosome", "tx", "gene")
    checkEquals(sort(need), sort(ensembldb:::addRequiredTables(edb, have)))

    if (hasProteinData(edb)) {
        have <- c("uniprot", "exon")
        need <- c("uniprot", "exon", "protein", "tx", "tx2exon")
        checkEquals(sort(need),
                    sort(ensembldb:::addRequiredTables(edb, have)))

        have <- c("uniprot", "chromosome")
        need <- c("uniprot", "chromosome", "protein", "tx", "gene")
        checkEquals(sort(need),
                    sort(ensembldb:::addRequiredTables(edb, have)))

        have <- c("protein_domain", "gene")
        need <- c("protein_domain", "gene", "protein", "tx")
        checkEquals(sort(need),
                    sort(ensembldb:::addRequiredTables(edb, have)))

        have <- c("protein", "exon")
        need <- c("protein", "exon", "tx", "tx2exon")
        checkEquals(sort(need),
                    sort(ensembldb:::addRequiredTables(edb, have)))
    }
}

test_prefixColumns <- function() {
    ## gene
    res <- ensembldb:::prefixColumns(edb, "gene_id")
    checkEquals(unname(unlist(res)), "gene.gene_id")
    res <- ensembldb:::prefixColumns(edb, "gene_id",
                                     with.tables = c("tx", "exon"))
    checkEquals(unname(unlist(res)), "tx.gene_id")
    ## tx
    res <- ensembldb:::prefixColumns(edb, "tx_id")
    checkEquals(unname(unlist(res)), "tx.tx_id")
    res <- ensembldb:::prefixColumns(edb, "tx_id",
                                     with.tables = "tx2exon")
    checkEquals(unname(unlist(res)), "tx2exon.tx_id")
    ## exon
    res <- ensembldb:::prefixColumns(edb, "exon_id")
    checkEquals(unname(unlist(res)), "tx2exon.exon_id")
    res <- ensembldb:::prefixColumns(edb, "exon_id",
                                     with.tables = c("exon"))
    checkEquals(unname(unlist(res)), "exon.exon_id")

    if (hasProteinData(edb)) {
        res <- ensembldb:::prefixColumns(edb, "protein_id")
        checkEquals(unname(unlist(res)), "protein.protein_id")
        res <- ensembldb:::prefixColumns(edb, "protein_id",
                                         with.tables = c("uniprot", "protein_domain"))
        checkEquals(unname(unlist(res)), "uniprot.protein_id")
    }
}

test_joinQueryOnTables <- function() {
    res <- ensembldb:::joinQueryOnTables(edb, tab = c("exon", "gene"))
    checkEquals(res, paste0("gene join tx on (gene.gene_id=tx.gene_id) ",
                            "join tx2exon on (tx.tx_id=tx2exon.tx_id) ",
                            "join exon on (tx2exon.exon_id=exon.exon_id)"))
    res <- ensembldb:::joinQueryOnTables(edb, tab = c("gene", "chromosome"))
    checkEquals(res, "gene join chromosome on (gene.seq_name=chromosome.seq_name)")
    res <- ensembldb:::joinQueryOnTables(edb, tab = c("exon", "tx"))
    checkEquals(res, paste0("tx join tx2exon on (tx.tx_id=tx2exon.tx_id) ",
                            "join exon on (tx2exon.exon_id=exon.exon_id)"))
    res <- ensembldb:::joinQueryOnTables(edb, tab = c("chromosome", "tx"))
    checkEquals(res, paste0("gene join tx on (gene.gene_id=tx.gene_id) ",
                            "join chromosome on (gene.seq_name=chromosome.seq_name)"))
    if (hasProteinData(edb)) {
        res <- ensembldb:::joinQueryOnTables(edb, tab = c("tx", "uniprot"))
        checkEquals(res, paste0("tx join protein on (tx.tx_id=protein.tx_id) ",
                                "join uniprot on (protein.protein_id=uniprot.protein_id)"))
        res <- ensembldb:::joinQueryOnTables(edb, tab = c("protein", "gene"))
        checkEquals(res, paste0("gene join tx on (gene.gene_id=tx.gene_id) ",
                                "join protein on (tx.tx_id=protein.tx_id)"))
        res <- ensembldb:::joinQueryOnTables(edb,
                                             tab = c("uniprot", "protein_domain"))
        checkEquals(res, "uniprot join protein_domain on (uniprot.protein_id=protein_domain.protein_id)")
        res <- ensembldb:::joinQueryOnTables(edb, tab = c("uniprot", "exon"))
        checkEquals(res, paste0("tx join tx2exon on (tx.tx_id=tx2exon.tx_id) ",
                                "join exon on (tx2exon.exon_id=exon.exon_id) ",
                                "join protein on (tx.tx_id=protein.tx_id) ",
                                "join uniprot on (protein.protein_id=uniprot.protein_id)"))
    }
}
