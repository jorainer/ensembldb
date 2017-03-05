## Some utility functions for Filters.

## Vector to map AnnotationFilter fields to actual database columns.
## Format: field name = database column name
.ENSDB_FIELDS <- c(
    ## gene
    entrez = "entrezid",
    gene_biotype = "gene_biotype",
    gene_id = "gene_id",
    genename = "gene_name",
    symbol = "gene_name",
    seq_name = "seq_name",
    seq_strand = "seq_strand",
    gene_start = "gene_seq_start",
    gene_end = "gene_seq_end",
    ## tx
    tx_id = "tx_id",
    tx_biotype = "tx_biotype",
    tx_name = "tx_id",
    tx_start = "tx_seq_start",
    tx_end = "tx_seq_end",
    ## exon
    exon_id = "exon_id",
    exon_rank = "exon_idx",
    exon_start = "exon_seq_start",
    exon_end = "exon_seq_end",
    ## protein
    protein_id = "protein_id",
    uniprot = "uniprot_id",
    uniprot_db = "uniprot_db",
    uniprot_mapping_type = "uniprot_mapping_type",
    prot_dom_id = "protein_domain_id"
)

.supportedFilters <- function(x) {
    flts <- c(
        "EntrezFilter", "GeneBiotypeFilter", "GeneIdFilter", "GenenameFilter",
        "SymbolFilter", "SeqNameFilter", "SeqStrandFilter", "GeneStartFilter",
        "GeneEndFilter", "TxIdFilter", "TxBiotypeFilter", "TxNameFilter",
        "TxStartFilter", "TxEndFilter", "ExonIdFilter", "ExonRankFilter",
        "ExonStartFilter", "ExonEndFilter", "GRangesFilter"
    )
    if (hasProteinData(x))
        flts <- c(flts, "ProteinIdFilter", "UniprotFilter", "UniprotDbFilter",
                  "UniprotMappingTypeFilter", "ProtDomIdFilter")
    return(sort(flts))
}

#' Utility function to map from the default AnnotationFilters fields to the
#' database columns used in ensembldb.
#'
#' @param x The field name to be \emph{translated}.
#' @return The column name in the EnsDb database.
#' @noRd
.fieldInEnsDb <- function(x) {
    if (length(x) == 0 || missing(x))
        stop("Error in .fieldInEnsDb: got empty input argument!")
    cn <- .ENSDB_FIELDS[x]
    if (is.na(cn))
        stop("Unable to map field '", x, "'!")
    return(cn)
}


#' Utility function to map the condition of an AnnotationFilter to the SQL
#' condition to be used in the EnsDb database.
#'
#' @param x An \code{AnnotationFilter}.
#'
#' @return A character representing the condition for the SQL call.
#' @noRd
.conditionForEnsDb <- function(x) {
    cond <- condition(x)
    if (length(value(x)) > 1) {
        if (cond == "==")
            cond <- "in"
        if (cond == "!=")
            cond <- "not in"
    }
    if (cond == "==")
        cond <- "="
    if (cond %in% c("startsWith", "endsWith"))
        cond <- "like"
    return(cond)
}

#' Single quote character values, paste multiple values and enclose in quotes.
#'
#' @param x An \code{AnnotationFilter} object.
#' @noRd
.valueForEnsDb <- function(x) {
    vals <- unique(x@value)
    if (x@.valueIsCharacter) {
        vals <- sQuote(gsub(unique(vals), pattern = "'", replacement = "''"))
    }
    if (length(vals) > 1)
        vals <- paste0("(",  paste0(vals, collapse = ","), ")")
    ## Process the like/startsWith/endsWith
    if (condition(x) == "startsWith")
        vals <- paste0("'", unique(x@value), "%'")
    if (condition(x) == "endsWith")
        vals <- paste0("'%", unique(x@value), "'")
    return(vals)
}

#' That's to build the standard query from an AnnotationFilter for EnsDb.
#'
#' @param x An \code{AnnotationFilter}.
#' @noRd
.queryForEnsDb <- function(x) {
    paste(.fieldInEnsDb(x@field), .conditionForEnsDb(x), .valueForEnsDb(x))
}

#' This is a slightly more sophisticated function that does also prefix the
#' columns.
#' @noRd
.queryForEnsDbWithTables <- function(x, db, tables = character()) {
    clmn <- .fieldInEnsDb(x@field)
    if (!missing(db)) {
        if (length(tables) == 0)
            tables <- names(listTables(db))
        clmn <- unlist(prefixColumns(db, clmn, with.tables = tables))
    }
    res <- paste(clmn, .conditionForEnsDb(x), .valueForEnsDb(x))
    ## cat("  ", res, "\n")
    return(res)
}
