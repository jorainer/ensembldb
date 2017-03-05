
setMethod("ensDbColumn", "AnnotationFilter",
          function(object, db, with.tables = character()) {
              clmn <- .fieldInEnsDb(object@field)
              if (missing(db))
                  return(clmn)
              if (length(with.tables) == 0)
                  with.tables <- names(listTables(db))
              return(unlist(prefixColumns(db, clmn, with.tables = with.tables),
                             use.names = FALSE))
          })
setMethod("ensDbQuery", "AnnotationFilter",
          function(object, db, with.tables = character()) {
              return(.queryForEnsDbWithTables(object, db, with.tables))
          })

setMethod("ensDbQuery", "list",
          function(object, db, with.tables = character()) {
              wq <- paste0(" where ",
                           paste(unlist(lapply(object, ensDbQuery, db,
                                               with.tables = with.tables),
                                        use.names = FALSE),
                                 collapse = " and "))
              return(wq)
          })

## setMethod("value", "SeqNameFilter",
##           function(object, db){
##               if (missing(db))
##                   return(object@value)
##               val <- formatSeqnamesForQuery(db, object@value)
##               if(any(is.na(val))){
##                   stop("A value of <NA> is not allowed for a SeqNameFilter!")
##               }
##               if (length(val) > 1)
##                   val <- paste0("(",  paste0(val, collapse = ","), ")")
##               return(val)
##               ##return(ucscToEns(value(x)))
##           })
setMethod("ensDbQuery", "SeqNameFilter",
          function(object, db, with.tables = character()) {
              ## val <- sQuote(value(object, db))
              ## Doing all the stuff in here:
              vals <- object@value
              clmn <- .fieldInEnsDb(object@field)
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
              return(paste(clmn, .conditionForEnsDb(object), vals))
          })

setMethod("ensDbQuery", "SeqStrandFilter",
          function(object, db, with.tables = character()) {
              ## We have to ensure that value is converted to +1, -1.
              val <- strand2num(value(object))
              clmn <- .fieldInEnsDb(object@field)
              if (!missing(db)) {
                  if (length(with.tables) == 0)
                      with.tables <- names(listTables(db))
                  clmn <- unlist(prefixColumns(db, clmn,
                                               with.tables = with.tables))
              }
              return(paste(clmn, .conditionForEnsDb(object), val))
          })


setMethod("start", signature(x="GRangesFilter"),
          function(x, ...){
              return(start(value(x)))
          })
setMethod("end", signature(x="GRangesFilter"),
          function(x, ...){
              return(end(value(x)))
          })
setMethod("strand", signature(x="GRangesFilter"),
          function(x, ...){
              strnd <- as.character(strand(value(x)))
              return(strnd)
          })
#' @description \code{seqnames}: accessor for the sequence names of the
#' \code{GRanges} object within a \code{GRangesFilter}
#' @param x For \code{seqnames}, \code{seqlevels}: a \code{GRangesFilter} object.
#' 
#' @rdname Filter-classes
setMethod("seqnames", signature(x="GRangesFilter"),
          function(x){
              return(as.character(seqnames(value(x))))
          })
#' @description \code{seqnames}: accessor for the \code{seqlevels} of the
#' \code{GRanges} object within a \code{GRangesFilter}
#' 
#' @rdname Filter-classes
setMethod("seqlevels", signature(x="GRangesFilter"),
          function(x){
              return(seqlevels(value(x)))
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
              return(cols[c("start", "end", "seqname", "strand")])
          })
setMethod("ensDbQuery", "GRangesFilter",
          function(object, db, with.tables = character()) {
              cols <- ensDbColumn(object, db, with.tables)
              if (missing(db))
                  db <- NULL
              return(buildWhereForGRanges(object, cols, db = db))
          })


## grf: GRangesFilter
buildWhereForGRanges <- function(grf, columns, db=NULL){
    condition <- condition(grf)
    if(!any(condition == c("within", "overlapping")))
        stop(paste0("'condition' for GRangesFilter should either be ",
                    "'within' or 'overlapping', got ", condition, "."))
    if(is.null(names(columns))){
        stop(paste0("The vector with the required column names for the",
                    " GRangesFilter query has to have names!"))
    }
    if(!all(c("start", "end", "seqname", "strand") %in% names(columns)))
        stop(paste0("'columns' has to be a named vector with names being ",
                    "'start', 'end', 'seqname', 'strand'!"))
    ## Build the query to fetch all features that are located within the range
    quers <- sapply(value(grf), function(z){
        if(!is.null(db)){
            seqn <- formatSeqnamesForQuery(db, as.character(seqnames(z)))
        }else{
            seqn <- as.character(seqnames(z))
        }
        if(condition == "within"){
            query <- paste0(columns["start"], " >= ", start(z), " and ",
                            columns["end"], " <= ", end(z), " and ",
                            columns["seqname"], " = '", seqn, "'")
        }
        ## Build the query to fetch all features (partially) overlapping
        ## the range. This includes also all features (genes or transcripts)
        ## that have an intron at that position.
        if(condition == "overlapping"){
            query <- paste0(columns["start"], " <= ", end(z), " and ",
                            columns["end"], " >= ", start(z), " and ",
                            columns["seqname"], " = '", seqn, "'")
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
    query <- paste0(quers, collapse=" or ")
    ## Collapse now the queries.
    return(query)
}



## map chromosome strand...
strand2num <- function(x){
    if (is.numeric(x)) {
        if (x >= 0) return(1)
        else return(-1)
    }
    xm <- x
    if(xm == "+" | xm == "-")
        xm <- paste0(xm, 1)
    xm <- as.numeric(xm)
    if (is.na(xm))
        stop("'", x, "' can not be converted to a strand!")
    return(xm)
}
num2strand <- function(x){
    if(x < 0){
        return("-")
    }else{
        return("+")
    }
}

setMethod("ensDbColumn", signature(object = "OnlyCodingTxFilter"),
          function(object, db, ...) {
              return("tx.tx_cds_seq_start")
})
setMethod("ensDbQuery", "OnlyCodingTxFilter",
          function(object, db, with.tables = character()) {
              return("tx.tx_cds_seq_start is not null")
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
              ## return(unlist(prefixColumns(db, column(object),
              ##                             with.tables = with.tables),
              ##               use.names = FALSE))
          })
setMethod("ensDbQuery", "ProteinIdFilter",
          function(object, db, with.tables = character()) {
              if (missing(db))
                  return(callNextMethod())
              if (!hasProteinData(db))
                  stop("The 'EnsDb' database used does not provide",
                       " protein annotations! A 'ProteinIdFilter' can not",
                       " be used.")
              return(.queryForEnsDbWithTables(object, db, with.tables))
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
              ## return(unlist(prefixColumns(db, column(object),
              ##                             with.tables = with.tables),
              ##               use.names = FALSE))
          })
setMethod("ensDbQuery", "UniprotFilter",
          function(object, db, with.tables = character()) {
              if (missing(db))
                  return(callNextMethod())
              if (!hasProteinData(db))
                  stop("The 'EnsDb' database used does not provide",
                       " protein annotations! A 'UniprotFilter' can not",
                       " be used.")
              return(.queryForEnsDbWithTables(object, db, with.tables))
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
              ## return(unlist(prefixColumns(db, column(object),
              ##                             with.tables = with.tables),
              ##               use.names = FALSE))
          })
setMethod("ensDbQuery", "ProtDomIdFilter",
          function(object, db, with.tables = character()) {
              if (missing(db))
                  return(callNextMethod())
              if (!hasProteinData(db))
                  stop("The 'EnsDb' database used does not provide",
                       " protein annotations! A 'ProtDomIdFilter' can not",
                       " be used.")
              return(.queryForEnsDbWithTables(object, db, with.tables))
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
              ## return(unlist(prefixColumns(db, column(object),
              ##                             with.tables = with.tables),
              ##               use.names = FALSE))
          })
setMethod("ensDbQuery", "UniprotDbFilter",
          function(object, db, with.tables = character()) {
              if (missing(db))
                  return(callNextMethod())
              if (!hasProteinData(db))
                  stop("The 'EnsDb' database used does not provide",
                       " protein annotations! A 'ProteinIdFilter' can not",
                       " be used.")
              return(.queryForEnsDbWithTables(object, db, with.tables))
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
              ## return(unlist(prefixColumns(db, column(object),
              ##                             with.tables = with.tables),
              ##               use.names = FALSE))
          })
setMethod("ensDbQuery", "UniprotMappingTypeFilter",
          function(object, db, with.tables = character()) {
              if (missing(db))
                  return(callNextMethod())
              if (!hasProteinData(db))
                  stop("The 'EnsDb' database used does not provide",
                       " protein annotations! A 'UniprotMappingTypeFilter' can not",
                       " be used.")
              return(.queryForEnsDbWithTables(object, db, with.tables))
          })
