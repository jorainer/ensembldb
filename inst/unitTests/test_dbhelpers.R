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

    if (ensembldb:::hasProteinData(edb)) {
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
}

test_joinQueryOnTables <- function() {
}
