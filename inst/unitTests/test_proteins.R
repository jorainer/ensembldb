############################################################
## Tests related to protein data.
library(EnsDb.Hsapiens.v75)
edb <- EnsDb.Hsapiens.v75

############################################################
## Getting protein data in other methods.
test_genes_with_proteins <- function() {
    suppressWarnings(
        res <- genes(edb, columns = c("gene_name", "gene_id", "protein_id",
                                      "uniprot_id"),
                     filter = GenenameFilter("ZBTB16"),
                     return.type = "data.frame")
    )
    if (hasProteinData(edb)) {
        ## have a 1:n mapping of protein_id to uniprot id:
        checkTrue(length(unique(res$protein_id)) <
                  nrow(res))
        checkEquals(colnames(res), c("gene_name", "gene_id", "protein_id",
                                     "uniprot_id"))
    } else {
        checkEquals(colnames(res), c("gene_name", "gene_id"))
        checkEquals(nrow(res), length(unique(res$gene_id)))
    }
    ## Next one fetching also protein domain data.
    suppressWarnings(
        res <- genes(edb, columns = c("gene_name", "tx_id", "protein_id",
                                      "protein_domain_id"),
                     filter = GenenameFilter("ZBTB16"),
                     return.type = "data.frame")
    )
    if (hasProteinData(edb)) {
        checkEquals(colnames(res), c("gene_name", "tx_id", "protein_id",
                                     "protein_domain_id", "gene_id"))
        checkTrue(nrow(res) > length(unique(res$protein_id)))
        checkTrue(nrow(res) > length(unique(res$tx_id)))
    } else {
        checkEquals(colnames(res), c("gene_name", "tx_id", "gene_id"))
        checkTrue(nrow(res) == length(unique(res$tx_id)))
    }
}

test_transcripts_with_proteins <- function() {
    suppressWarnings(
        res <- transcripts(edb, columns = c("tx_biotype", "protein_id",
                                            "uniprot_id"),
                           filter = GenenameFilter("ZBTB16"),
                           return.type = "data.frame")
    )
}

############################################################
## Using protein data based filters.
test_ProteinidFilter <- function() {
    pf <- ProteinidFilter("ABC")
    checkEquals(value(pf), "ABC")
    checkEquals(column(pf), "protein_id")
    checkEquals(where(pf), "protein_id = 'ABC'")
    if (hasProteinData(edb)) {
        checkEquals(column(pf, edb), "protein.protein_id")
        checkEquals(column(pf, edb, with.tables = "protein_domain"),
                    "protein_domain.protein_id")
        checkEquals(column(pf, edb, with.tables = "uniprot"),
                    "uniprot.protein_id")
        checkEquals(where(pf, edb), "protein.protein_id = 'ABC'")
        checkEquals(where(pf, edb, with.tables = "uniprot"),
                    "uniprot.protein_id = 'ABC'")
    } else {
        checkException(column(pf, edb))
        checkException(where(pf, edb))
        checkException(column(pf, edb, with.tables = "uniprot"))
        checkException(where(pf, edb, with.tables = "uniprot"))
    }
    pf <- ProteinidFilter(c("A", "B"))
    checkEquals(where(pf), "protein_id in ('A','B')")
    checkException(ProteinidFilter("B", condition = ">"))
}

test_UniprotidFilter <- function() {
    pf <- UniprotidFilter("ABC")
    checkEquals(value(pf), "ABC")
    checkEquals(column(pf), "uniprot_id")
    checkEquals(where(pf), "uniprot_id = 'ABC'")
    if (hasProteinData(edb)) {
        checkEquals(column(pf, edb), "uniprot.uniprot_id")
        checkEquals(column(pf, edb, with.tables = "protein_domain"),
                    "protein_domain.protein_id")
        checkEquals(column(pf, edb, with.tables = "uniprot"),
                    "uniprot.protein_id")
        checkEquals(where(pf, edb), "protein.protein_id = 'ABC'")
        checkEquals(where(pf, edb, with.tables = "uniprot"),
                    "uniprot.protein_id = 'ABC'")
    } else {
        checkException(column(pf, edb))
        checkException(where(pf, edb))
        checkException(column(pf, edb, with.tables = "uniprot"))
        checkException(where(pf, edb, with.tables = "uniprot"))
    }
    pf <- UniprotidFilter(c("A", "B"))
    checkEquals(where(pf), "uniprot_id in ('A','B')")
    checkException(UniprotidFilter("B", condition = ">"))
}

test_ProtdomidFilter <- function() {
    pf <- ProtdomidFilter("ABC")
    checkEquals(value(pf), "ABC")
    checkEquals(column(pf), "protein_domain_id")
    checkEquals(where(pf), "protein_domain_id = 'ABC'")
    if (hasProteinData(edb)) {
        checkEquals(column(pf, edb), "protein_domain.protein_domain_id")
        checkEquals(column(pf, edb, with.tables = "protein_domain"),
                    "protein_domain.protein_domain_id")
        checkEquals(where(pf, edb), "protein_domain.protein_domain_id = 'ABC'")
        checkEquals(where(pf, edb, with.tables = "protein_domain"),
                    "protein.protein_domain_id = 'ABC'")
    } else {
        checkException(column(pf, edb))
        checkException(where(pf, edb))
        checkException(column(pf, edb, with.tables = "protein_domain"))
        checkException(where(pf, edb, with.tables = "protein_domain"))
    }
    pf <- ProtdomidFilter(c("A", "B"))
    checkEquals(where(pf), "protein_domain_id in ('A','B')")
    checkException(ProtdomidFilter("B", condition = ">"))
}


############################################################
## The dedicated methods to fetch protein data.
