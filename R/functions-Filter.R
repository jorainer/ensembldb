## Some utility functions for Filters.

## Vector to map AnnotationFilter fields to actual database columns.
## Format: field name = database column name
.ENSDB_FIELDS <- c(
    ## gene
    entrez = "entrezid",
    gene_biotype = "gene_biotype",
    gene_id = "gene_id",
    genename = "gene_name",
    gene_name = "gene_name",
    symbol = "gene_name",
    seq_name = "seq_name",
    seq_strand = "seq_strand",
    gene_start = "gene_seq_start",
    gene_end = "gene_seq_end",
    description = "description",
    ## tx
    tx_id = "tx_id",
    tx_biotype = "tx_biotype",
    tx_name = "tx_id",
    tx_start = "tx_seq_start",
    tx_end = "tx_seq_end",
    tx_support_level = "tx_support_level",
    tx_external_name = "tx_external_name",
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
    prot_dom_id = "protein_domain_id",
    protein_domain_id = "protein_domain_id",
    protein_domain_source = "protein_domain_source"
)

.supportedFilters <- function(x) {
    flds <- .filterFields(x)
    flts <- c(.fieldToClass(flds), "GRangesFilter")
    flds <- c(flds, NA)
    idx <- order(flts)
    data.frame(filter = flts[idx], field = flds[idx], stringsAsFactors = FALSE)
}

.filterFields <- function(x) {
    flds <- c("entrez", "gene_biotype", "gene_id", "gene_name", "genename",
              "symbol", "seq_name", "seq_strand", "gene_start", "gene_end",
              "tx_id", "tx_biotype", "tx_name", "tx_start", "tx_end",
              "exon_id", "exon_rank", "exon_start", "exon_end")
    if (.has_tx_external_name(x))
        flds <- c(flds, "tx_external_name")
    if (hasProteinData(x))
        flds <- c(flds, "protein_id", "uniprot", "uniprot_db",
                  "uniprot_mapping_type", "prot_dom_id",
                  "protein_domain_id",
                  "protein_domain_source")
    if (any(listColumns(x) == "tx_support_level"))
        flds <- c(flds, "tx_support_level")
    sort(flds)
}

.fieldToClass <- function(field) {
    class <- gsub("_([[:alpha:]])", "\\U\\1", field, perl=TRUE)
    class <- sub("^([[:alpha:]])", "\\U\\1", class, perl=TRUE)
    paste0(class, if (length(class)) "Filter" else character(0))
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
    if (is.na(.ENSDB_FIELDS[x]))
        stop("Unable to map field '", x, "'!")
    else
        .ENSDB_FIELDS[x]
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
    if (length(unique(value(x))) > 1) {
        if (cond == "==")
            cond <- "in"
        if (cond == "!=")
            cond <- "not in"
    }
    if (cond == "==")
        cond <- "="
    if (cond %in% c("startsWith", "endsWith", "contains"))
        cond <- "like"
    cond
}

#' Single quote character values, paste multiple values and enclose in quotes.
#'
#' @param x An \code{AnnotationFilter} object.
#' @noRd
.valueForEnsDb <- function(x) {
    vals <- unique(value(x))
    if (is(x, "CharacterFilter")) {
        vals <- sQuote(gsub(unique(vals), pattern = "'", replacement = "''"))
    }
    if (length(vals) > 1)
        vals <- paste0("(",  paste0(vals, collapse = ","), ")")
    ## Process the like/startsWith/endsWith
    if (condition(x) == "startsWith")
        vals <- paste0("'", unique(x@value), "%'")
    if (condition(x) == "endsWith")
        vals <- paste0("'%", unique(x@value), "'")
    if (condition(x) == "contains")
        vals <- paste0("'%", unique(x@value), "%'")
    vals
}

#' That's to build the standard query from an AnnotationFilter for EnsDb.
#'
#' @param x An \code{AnnotationFilter}.
#' @noRd
.queryForEnsDb <- function(x) {
    paste(.fieldInEnsDb(field(x)), .conditionForEnsDb(x), .valueForEnsDb(x))
}

#' This is a slightly more sophisticated function that does also prefix the
#' columns.
#' @noRd
.queryForEnsDbWithTables <- function(x, db, tables = character()) {
    clmn <- .fieldInEnsDb(field(x))
    if (!missing(db)) {
        if (length(tables) == 0)
            tables <- names(listTables(db))
        clmn <- unlist(prefixColumns(db, clmn, with.tables = tables))
    }
    res <- paste(clmn, .conditionForEnsDb(x), .valueForEnsDb(x))
    res
}

#' Simple helper function to convert expressions to AnnotationFilter or
#' AnnotationFilterList.
#'
#' @param x Can be an \code{AnnotationFilter}, an \code{AnnotationFilterList},
#' a \code{list} or a filter \code{expression}. This should NOT be empty!
#'
#' @return Returns an \code{AnnotationFilterList} with all filters.
#'
#' @noRd
.processFilterParam <- function(x, db) {
    if (missing(db))
        stop("Argument 'db' missing.")
    ## Check if x is a formula and eventually translate it.
    if (is(x, "formula"))
        res <- AnnotationFilter(x)
    else res <- x
    if (is(res, "AnnotationFilter"))
        res <- AnnotationFilterList(res)
    if (!is(res, "AnnotationFilterList")) {
        ## Did not get a filter expression, thus checking what we've got.
        if (is(res, "list")) {
            if (length(res)) {
                ## Check that all elements are AnnotationFilter objects!
                if (!all(unlist(lapply(res, function(z) {
                    inherits(z, "AnnotationFilter")
                }), use.names = FALSE)))
                    stop("One of more elements in 'filter' are not ",
                         "'AnnotationFilter' objects!", call. = FALSE)
                res <- as(res, "AnnotationFilterList")
                res@logOp <- rep("&", (length(res) - 1))
            } else {
                res <- AnnotationFilterList()
            }
        } else {
            stop("'filter' has to be an 'AnnotationFilter', a list of ",
                 "'AnnotationFilter' object, an 'AnnotationFilterList' ",
                 "or a valid filter expression!", call. = FALSE)
        }
    }
    supp_filters <- supportedFilters(db)$filter
    have_filters <- unique(.AnnotationFilterClassNames(res))
    if (!all(have_filters %in% supp_filters))
        stop("AnnotationFilter classes: ",
             paste(have_filters[!(have_filters %in% supp_filters)]),
             " are not supported by this EnsDb database.", call. = FALSE)
    res
}


############################################################
## setFeatureInGRangesFilter
##
## Simple helper function to set the @feature in GRangesFilter
## depending on the calling method.
setFeatureInGRangesFilter <- function(x, feature){
    for (i in seq(along.with = x)){
        if (is(x[[i]], "GRangesFilter"))
            x[[i]]@feature <- feature
        if (is(x[[i]], "AnnotationFilterList"))
            x[[i]] <- setFeatureInGRangesFilter(x[[i]], feature = feature)
    }
    x
}

############################################################
## isProteinFilter
##' evaluates whether the filter is a protein annotation related filter.
##' @param x The object that should be evaluated. Can be an AnnotationFilter or
##'     an AnnotationFilterList.
##' @return Returns TRUE if 'x' is a filter for protein annotation tables and
##' FALSE otherwise.
##' @noRd
isProteinFilter <- function(x) {
    if (is(x, "AnnotationFilterList"))
        return(unlist(lapply(x, isProteinFilter)))
    else
        return(is(x, "ProteinIdFilter") | is(x, "UniprotFilter") |
               is(x, "ProtDomIdFilter") | is(x, "UniprotDbFilter") |
               is(x, "UniprotMappingTypeFilter"))
}

#' build the \emph{where} query for a \code{GRangedFilter}. Supported conditions
#' are: \code{"start"}, \code{"end"}, \code{"equal"}, \code{"within"},
#' \code{"any"}.
#'
#' @param grf \code{GRangesFilter}.
#'
#' @param columns named character vectors with the column names for start, end,
#'     strand and seq_name.
#'
#' @param db An optional \code{EnsDb} instance. Used to \emph{translate}
#'     seqnames depending on the specified seqlevels style.
#'
#' @return A character with the corresponding \emph{where} query.
#' @noRd
buildWhereForGRanges <- function(grf, columns, db = NULL){
    condition <- condition(grf)
    if (!(condition %in% c("start", "end", "within", "equal", "any")))
        stop("'condition' ", condition, " not supported. Condition (type) can ",
             "be one of 'any', 'start', 'end', 'equal', 'within'.")
    if( is.null(names(columns)))
        stop("The vector with the required column names for the",
             " GRangesFilter query has to have names!")
    if (!all(c("start", "end", "seqname", "strand") %in% names(columns)))
        stop("'columns' has to be a named vector with names being ",
             "'start', 'end', 'seqname', 'strand'!")
    ## Build the query to fetch all features that are located within the range
    quers <- sapply(as(value(grf), "GRangesList"), function(z) {
        if (!is.null(db)) {
            seqn <- formatSeqnamesForQuery(db, as.character(seqnames(z)))
        } else {
            seqn <- as.character(seqnames(z))
        }
        ## start: start, seqname and strand have to match.
        if (condition == "start") {
            query <- paste0(columns["start"], "=", start(z), " and ",
                            columns["seqname"], "='", seqn, "'")
        }
        ## end: end, seqname and strand have to match.
        if (condition == "end") {
            query <- paste0(columns["end"], "=", end(z), " and ",
                            columns["seqname"], "='", seqn, "'")
        }
        ## equal: start, end, seqname and strand have to match.
        if (condition == "equal") {
            query <- paste0(columns["start"], "=", start(z), " and ",
                            columns["end"], "=", end(z), " and ",
                            columns["seqname"], "='", seqn, "'")
        }
        ## within: start has to be >= start, end <= end, seqname and strand
        ##         have to match.
        if (condition == "within") {
            query <- paste0(columns["start"], ">=", start(z), " and ",
                            columns["end"], "<=", end(z), " and ",
                            columns["seqname"], "='", seqn, "'")
        }
        ## any: essentially the overlapping.
        if (condition == "any") {
            query <- paste0(columns["start"], "<=", end(z), " and ",
                            columns["end"], ">=", start(z), " and ",
                            columns["seqname"], "='", seqn, "'")
        }
        ## Include the strand, if it's not "*"
        if(as.character(strand(z)) != "*"){
            query <- paste0(query, " and ", columns["strand"], " = ",
                            strand2num(as.character(strand(z))))
        }
        return(query)
    })
    if(length(quers) > 1)
        quers <- paste0("(", quers, ")")
    ## Collapse now the queries.
    query <- paste0(quers, collapse=" or ")
    paste0("(", query, ")")
}

#' @description Helper to extract all AnnotationFilter class names from an
#'     AnnotationFilterList (recursively!)
#'
#' @param x The \code{AnnotationFilterList}.
#'
#' @return A \code{character} with the names of the classes.
#' @noRd
.AnnotationFilterClassNames <- function(x) {
    classes <- lapply(x, function(z) {
        if (is(z, "AnnotationFilterList"))
            return(.AnnotationFilterClassNames(z))
        class(z)
    })
    unlist(classes, use.names = FALSE)
}

#' @description Test if any of the filter(s) is an SymbolFilter.
#'
#' @noRd
.anyIs <- function(x, what = "SymbolFilter") {
    if (is(x, "AnnotationFilter")) {
        is(x, what)
    } else {
        unlist(lapply(x, .anyIs, what = what))
    }
}
