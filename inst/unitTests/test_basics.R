############################################################
## Testing basic functionality on EnsDb objects.
library(EnsDb.Hsapiens.v75)
library(RSQLite)
edb <- EnsDb.Hsapiens.v75

test_metadata <- function() {
    res <- metadata(edb)
    checkEquals(res, dbGetQuery(dbconn(edb), "select * from metadata"))
}

test_ensemblVersion <- function() {
    checkEquals(ensemblVersion(edb), "75")
}


test_getMetadataValue <- function() {
    checkException(ensembldb:::getMetadataValue(edb))
}

test_listTables <- function() {
    res <- listTables(edb)
    if (!hasProteinData(edb)) {
        checkEquals(names(res), names(ensembldb:::.ENSDB_TABLES))
    } else {
        checkEquals(names(res), c("gene", "tx", "tx2exon", "exon",
                                  "chromosome", "protein", "uniprot",
                                  "protein_domain", "metadata"))
    }
}

test_listColumns <- function() {
    res <- listColumns(edb, table = "gene")
    checkEquals(res, c(ensembldb:::.ENSDB_TABLES$gene, "symbol"))
    res <- listColumns(edb, table = "tx")
    checkEquals(res, c(ensembldb:::.ENSDB_TABLES$tx, "tx_name"))
    res <- listColumns(edb, table = "exon")
    checkEquals(res, c(ensembldb:::.ENSDB_TABLES$exon))
    res <- listColumns(edb, table = "chromosome")
    checkEquals(res, c(ensembldb:::.ENSDB_TABLES$chromosome))
    res <- listColumns(edb, table = "tx2exon")
    checkEquals(res, c(ensembldb:::.ENSDB_TABLES$tx2exon))
    if (hasProteinData(edb)) {
        res <- listColumns(edb, table = "protein")
        checkEquals(res, ensembldb:::.ENSDB_PROTEIN_TABLES$protein)
        res <- listColumns(edb, table = "uniprot")
        checkEquals(res, ensembldb:::.ENSDB_PROTEIN_TABLES$uniprot)
        res <- listColumns(edb, table = "protein_domain")
        checkEquals(res, ensembldb:::.ENSDB_PROTEIN_TABLES$protein_domain)
    }
}

test_cleanColumns <- function() {
    cols <- c("gene_id", "tx_id", "tx_name")
    res <- ensembldb:::cleanColumns(edb, cols)
    checkEquals(cols, res)
    cols <- c(cols, "not there")
    suppressWarnings(
        res <- ensembldb:::cleanColumns(edb, cols)
    )
    checkEquals(cols[1:3], res)
    cols <- c("gene_id", "protein_id", "tx_id", "protein_sequence")
    suppressWarnings(
        res <- ensembldb:::cleanColumns(edb, cols)
    )
    if (hasProteinData(edb)) {
        checkEquals(res, cols)
    } else {
        checkEquals(res, cols[c(1, 3)])
    }
}

test_tablesForColumns <- function() {
    checkException(ensembldb:::tablesForColumns(edb))
    res <- ensembldb:::tablesForColumns(edb, columns = "tx_id")
    if (hasProteinData(edb))
        checkEquals(res, c("tx", "tx2exon", "protein"))
    else
        checkEquals(res, c("tx", "tx2exon"))

    res <- ensembldb:::tablesForColumns(edb, columns = "seq_name")
    checkEquals(res, c("gene", "chromosome"))

    if (hasProteinData(edb)) {
        res <- ensembldb:::tablesForColumns(edb, columns = "protein_id")
        checkEquals(res, c("protein", "uniprot", "protein_domain"))
    }
}

test_tablesByDegree <- function() {
    res <- ensembldb:::tablesByDegree(edb,
                                      tab = c("chromosome", "gene", "tx"))
    checkEquals(res, c("gene", "tx", "chromosome"))
}

test_updateEnsDb <- function(){
    edb2 <- updateEnsDb(edb)
    checkEquals(edb2@tables, edb@tables)
    checkTrue(.hasSlot(edb2, ".properties"))
}

test_properties <- function(){
    origProps <- ensembldb:::properties(edb)
    checkEquals(ensembldb:::getProperty(edb, "foo"), NA)

    checkException(ensembldb:::setProperty(edb, "foo"))

    edb <- ensembldb:::setProperty(edb, foo="bar")
    checkEquals(ensembldb:::getProperty(edb, "foo"), "bar")
    checkEquals(length(ensembldb:::properties(edb)),
                length(origProps) + 1)
}
