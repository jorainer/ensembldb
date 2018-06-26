## Methods for filter classes.

#' @description Extract the field/column name from an AnnotationFilter and
#'     ensure that it matches the correct database column name. Depending on
#'     whether argument \code{db} is present, the column names are also
#'     prefixed with the name of the corresponding table.
#'
#' @param object An \code{AnnotationFilter} object.
#'
#' @param db An \code{EnsDb} object.
#'
#' @param with.tables \code{character} specifying the tables that should be
#'     considered when prefixing the column name.
#' 
#' @noRd
setMethod("ensDbColumn", "AnnotationFilter",
          function(object, db, with.tables = character()) {
              clmn <- .fieldInEnsDb(object@field)
              if (missing(db))
                  return(clmn)
              if (length(with.tables) == 0)
                  with.tables <- names(listTables(db))
              unlist(prefixColumns(db, clmn, with.tables = with.tables),
                     use.names = FALSE)
          })

setMethod("ensDbColumn", "AnnotationFilterList",
          function(object, db, with.tables = character()) {
              if (length(object) == 0)
                  return(character())
              unique(unlist(lapply(object, ensDbColumn, db,
                                   with.tables = with.tables)))
          })

#' @title Convert an AnnotationFilter to a SQL WHERE condition for EnsDb
#'
#' @aliases convertFilter,AnnotationFilter,EnsDb-method
#'     convertFilter,AnnotationFilterList,EnsDb-method
#' 
#' @description `convertFilter` converts an `AnnotationFilter::AnnotationFilter`
#'     or `AnnotationFilter::AnnotationFilterList` to an SQL where condition
#'     for an `EnsDb` database.
#'
#' @note This function *might* be used in direct SQL queries on the SQLite
#'     database underlying an `EnsDb` but is more thought to illustrate the
#'     use of `AnnotationFilter` objects in combination with SQL databases.
#'     This method is used internally to create the SQL calls to the database.
#'
#' @param object `AnnotationFilter` or `AnnotationFilterList` objects (or
#'     objects extending these classes).
#'
#' @param db `EnsDb` object.
#'
#' @param with.tables optional `character` vector specifying the names of the
#'     database tables that are being queried.
#'
#' @return A `character(1)` with the SQL where condition.
#' 
#' @md
#'
#' @rdname convertFilter
#' 
#' @author Johannes Rainer
#'
#' @examples
#'
#' library(EnsDb.Hsapiens.v86)
#' edb <- EnsDb.Hsapiens.v86
#'
#' ## Define a filter
#' flt <- AnnotationFilter(~ gene_name == "BCL2")
#'
#' ## Use the method from the AnnotationFilter package:
#' convertFilter(flt)
#'
#' ## Create a combination of filters
#' flt_list <- AnnotationFilter(~ gene_name %in% c("BCL2", "BCL2L11") &
#'     tx_biotype == "protein_coding")
#' flt_list
#'
#' convertFilter(flt_list)
#'
#' ## Use the filters in the context of an EnsDb database:
#' convertFilter(flt, edb)
#'
#' convertFilter(flt_list, edb)
setMethod("convertFilter", signature = c(object = "AnnotationFilter",
                                         db = "EnsDb"),
          function(object, db, with.tables = character()) {
              ensDbQuery(object, db, with.tables)
          })
#' @rdname convertFilter
setMethod("convertFilter", signature = c(object = "AnnotationFilterList",
                                         db = "EnsDb"),
          function(object, db, with.tables = character()) {
              ensDbQuery(object, db, with.tables)
          })

#' @description Build the \emph{where} query for an \code{AnnotationFilter} or
#'     \code{AnnotationFilterList}.
#'
#' @noRd
setMethod("ensDbQuery", "AnnotationFilter",
          function(object, db, with.tables = character()) {
              .queryForEnsDbWithTables(object, db, with.tables)
          })

setMethod("ensDbQuery", "AnnotationFilterList",
          function(object, db, with.tables = character()) {
              wq <- NULL
              if (length(object)) {
                  wq <- ensDbQuery(object[[1]], db, with.tables = with.tables)
                  if (length(object) > 1) {
                      for (i in 2:length(object)) {
                          wq <- paste(wq, .logOp2SQL(object@logOp[(i -1)]),
                                      ensDbQuery(object[[i]], db,
                                                 with.tables = with.tables))
                      }
                  }
                  ## Encapsule all inside brackets.
                  wq <- paste0("(", wq, ")")
              }
              wq
          })

#' Need an ensDbQuery for SeqNameFilter to support different chromosome naming
#' styles
#' 
#' @noRd
setMethod("ensDbQuery", "SeqNameFilter",
          function(object, db, with.tables = character()) {
              ## val <- sQuote(value(object, db))
              ## Doing all the stuff in here:
              vals <- value(object)
              clmn <- .fieldInEnsDb(field(object))
              if (!missing(db)) {
                  ## o Eventually rename the seqname based on the seqlevelsStyle.
                  vals <- formatSeqnamesForQuery(db, vals)
                  if (length(with.tables) == 0)
                      with.tables <- names(listTables(db))
                  clmn <- unlist(prefixColumns(db, clmn,
                                               with.tables = with.tables))
              }
              ## o Quote the values.
              vals <- sQuote(vals)
              ## o Concatenate values.
              if (length(vals) > 1)
                  vals <- paste0("(", paste0(vals, collapse = ","), ")")
              paste(clmn, .conditionForEnsDb(object), vals)
          })

setMethod("ensDbQuery", "SeqStrandFilter",
          function(object, db, with.tables = character()) {
              ## We have to ensure that value is converted to +1, -1.
              val <- strand2num(value(object))
              clmn <- .fieldInEnsDb(field(object))
              if (!missing(db)) {
                  if (length(with.tables) == 0)
                      with.tables <- names(listTables(db))
                  clmn <- unlist(prefixColumns(db, clmn,
                                               with.tables = with.tables))
              }
              paste(clmn, .conditionForEnsDb(object), val)
          })


setMethod("start", signature(x="GRangesFilter"),
          function(x, ...){
              start(value(x))
          })

setMethod("end", signature(x="GRangesFilter"),
          function(x, ...){
              end(value(x))
          })

setMethod("strand", signature(x="GRangesFilter"),
          function(x, ...){
              as.character(strand(value(x)))
          })

#' @description \code{seqnames}: accessor for the sequence names of the
#' \code{GRanges} object within a \code{GRangesFilter}
#' @param x For \code{seqnames}, \code{seqlevels}: a \code{GRangesFilter} object.
#' 
#' @rdname Filter-classes
setMethod("seqnames", signature(x="GRangesFilter"),
          function(x){
              as.character(seqnames(value(x)))
          })

#' @description \code{seqnames}: accessor for the \code{seqlevels} of the
#' \code{GRanges} object within a \code{GRangesFilter}
#' 
#' @rdname Filter-classes
setMethod("seqlevels", signature(x="GRangesFilter"),
          function(x){
              seqlevels(value(x))
          })

setMethod("ensDbColumn", "GRangesFilter",
          function(object, db, with.tables = character(), ...){
              feature <- object@feature
              feature <- match.arg(feature, c("gene", "transcript", "exon",
                                              "tx"))
              if(object@feature == "transcript")
                  feature <- "tx"
              cols <- c(start = paste0(feature, "_seq_start"),
                        end = paste0(feature, "_seq_end"),
                        seqname = "seq_name",
                        strand = "seq_strand")

              if (!missing(db)) {
                  if (length(with.tables) == 0)
                      with.tables <- names(listTables(db))
                  cols <- unlist(prefixColumns(db, cols,
                                               with.tables = with.tables),
                                 use.names = FALSE)
                  ## We have to give the vector the required names!
                  names(cols) <- 1:length(cols)
                  names(cols)[grep(cols, pattern = "seq_name")] <- "seqname"
                  names(cols)[grep(cols, pattern = "seq_strand")] <- "strand"
                  names(cols)[grep(cols, pattern = "seq_start")] <- "start"
                  names(cols)[grep(cols, pattern = "seq_end")] <- "end"
              }
              cols[c("start", "end", "seqname", "strand")]
          })

setMethod("ensDbQuery", "GRangesFilter",
          function(object, db, with.tables = character()) {
              cols <- ensDbColumn(object, db, with.tables)
              if (missing(db))
                  db <- NULL
              buildWhereForGRanges(object, cols, db = db)
          })

setMethod("ensDbColumn", signature(object = "OnlyCodingTxFilter"),
          function(object, db, ...) {
              "tx.tx_cds_seq_start"
          })

setMethod("ensDbQuery", "OnlyCodingTxFilter",
          function(object, db, with.tables = character()) {
              "tx.tx_cds_seq_start is not null"
          })

setMethod("ensDbColumn", "ProteinIdFilter",
          function(object, db, with.tables = character(), ...) {
              if (missing(db)) {
                  return(callNextMethod())
              }
              if (!hasProteinData(db))
                  stop("The 'EnsDb' database used does not provide",
                       " protein annotations! A 'ProteinIdFilter' can not",
                       " be used.")
              callNextMethod()
          })

setMethod("ensDbQuery", "ProteinIdFilter",
          function(object, db, with.tables = character()) {
              if (missing(db))
                  return(callNextMethod())
              if (!hasProteinData(db))
                  stop("The 'EnsDb' database used does not provide",
                       " protein annotations! A 'ProteinIdFilter' can not",
                       " be used.")
              .queryForEnsDbWithTables(object, db, with.tables)
          })

setMethod("ensDbColumn", "UniprotFilter",
          function(object, db, with.tables = character(), ...) {
              if (missing(db))
                  return(callNextMethod())
              if (!hasProteinData(db))
                  stop("The 'EnsDb' database used does not provide",
                       " protein annotations! A 'UniprotFilter' can not",
                       " be used.")
              callNextMethod()
          })

setMethod("ensDbQuery", "UniprotFilter",
          function(object, db, with.tables = character()) {
              if (missing(db))
                  return(callNextMethod())
              if (!hasProteinData(db))
                  stop("The 'EnsDb' database used does not provide",
                       " protein annotations! A 'UniprotFilter' can not",
                       " be used.")
              .queryForEnsDbWithTables(object, db, with.tables)
          })

setMethod("ensDbColumn", "ProtDomIdFilter",
          function(object, db, with.tables = character(), ...) {
              if (missing(db))
                  return(callNextMethod())
              if (!hasProteinData(db))
                  stop("The 'EnsDb' database used does not provide",
                       " protein annotations! A 'ProtDomIdFilter' can not",
                       " be used.")
              callNextMethod()
          })

setMethod("ensDbQuery", "ProtDomIdFilter",
          function(object, db, with.tables = character()) {
              if (missing(db))
                  return(callNextMethod())
              if (!hasProteinData(db))
                  stop("The 'EnsDb' database used does not provide",
                       " protein annotations! A 'ProtDomIdFilter' can not",
                       " be used.")
              .queryForEnsDbWithTables(object, db, with.tables)
          })

setMethod("ensDbColumn", "UniprotDbFilter",
          function(object, db, with.tables = character(), ...) {
              if (missing(db))
                  return(callNextMethod())
              if (!hasProteinData(db))
                  stop("The 'EnsDb' database used does not provide",
                       " protein annotations! A 'UniprotDbFilter' can not",
                       " be used.")
              callNextMethod()
          })

setMethod("ensDbQuery", "UniprotDbFilter",
          function(object, db, with.tables = character()) {
              if (missing(db))
                  return(callNextMethod())
              if (!hasProteinData(db))
                  stop("The 'EnsDb' database used does not provide",
                       " protein annotations! A 'ProteinIdFilter' can not",
                       " be used.")
              .queryForEnsDbWithTables(object, db, with.tables)
          })

setMethod("ensDbColumn", "UniprotMappingTypeFilter",
          function(object, db, with.tables = character(), ...) {
              if (missing(db))
                  return(callNextMethod())
              if (!hasProteinData(db))
                  stop("The 'EnsDb' database used does not provide",
                       " protein annotations! A 'UniprotMappingTypeFilter' ",
                       "can not be used.")
              callNextMethod()
          })

setMethod("ensDbQuery", "UniprotMappingTypeFilter",
          function(object, db, with.tables = character()) {
              if (missing(db))
                  return(callNextMethod())
              if (!hasProteinData(db))
                  stop("The 'EnsDb' database used does not provide",
                       " protein annotations! A 'UniprotMappingTypeFilter' can not",
                       " be used.")
              .queryForEnsDbWithTables(object, db, with.tables)
          })

setMethod("ensDbColumn", "TxSupportLevelFilter",
          function(object, db, with.tables = character(), ...) {
              if (missing(db))
                  return(callNextMethod())
              if (!any(listColumns(db) %in% "tx_support_level"))
                  stop("The 'EnsDb' database used does not provide",
                       " transcript support levels! A 'TxSupportLevelFilter' ",
                       "can not be used.")
              callNextMethod()
          })

setMethod("ensDbQuery", "TxSupportLevelFilter",
          function(object, db, with.tables = character()) {
              if (missing(db))
                  return(callNextMethod())
              if (!any(listColumns(db) %in% "tx_support_level"))
                  stop("The 'EnsDb' database used does not provide",
                       " transcript support levels! A 'TxSupportLevelFilter' ",
                       "can not be used.")
              .queryForEnsDbWithTables(object, db, with.tables)
          })
