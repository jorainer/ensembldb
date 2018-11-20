## That's to support and interface the AnnotionDbi package.

####============================================================
##  .getColMappings
##
##  That returns a character vector of abbreviated column names
##  which can be/are used by AnnotationDbi with the names correponding
##  to the column names from ensembldb.
##  x: is supposed to be an EnsDb object.
##  all: if TRUE we return all of them, otherwise we just return those
##       that should be visible for the user.
####------------------------------------------------------------
.getColMappings <- function(x, all=FALSE){
    cols <- listColumns(x)
    if(!all){
        cols <- cols[!(cols %in% c("name", "value"))]
    }
    ret <- toupper(gsub("_", replacement="", cols))
    names(ret) <- cols
    return(ret)
}

####============================================================
##  columnForKeytype
##
##  Returns the appropriate column name in the database for the
##  given keytypes.
####------------------------------------------------------------
ensDbColumnForColumn <- function(x, column){
    maps <- .getColMappings(x)
    revmaps <- names(maps)
    names(revmaps) <- maps
    cols <- revmaps[column]
    if(any(is.na(cols))){
        warning("The following columns can not be mapped to column names in the",
                " db: ", paste(column[is.na(cols)], collapse=", "))
        cols <- cols[!is.na(cols)]
    }
    ## Fixing tx_name; tx_name should be mapped to tx_id in the database!
    ##cols[cols == "tx_name"] <- "tx_id"
    return(cols)
}


####============================================================
##  columns method
##
##  Just return the attributes, but as expected by the AnnotationDbi
##  interface (i.e. upper case, no _).
####------------------------------------------------------------
.getColumns <- function(x){
    cols <- .getColMappings(x, all=FALSE)
    names(cols) <- NULL
    return(unique(cols))
}
setMethod("columns", "EnsDb",
          function(x) .getColumns(x)
          )


####============================================================
##  keytypes method
##
##  I will essentially use all of the filters here.
####------------------------------------------------------------
setMethod("keytypes", "EnsDb",
          function(x){
              return(.filterKeytypes(withProteins = hasProteinData(x)))
          }
)
## This just returns some (eventually) usefull names for keys
.simpleKeytypes <- function(x){
    return(c("GENEID","TXID","TXNAME","EXONID","EXONNAME","CDSID","CDSNAME"))
}
.filterKeytypes <- function(withProteins = FALSE){
    return(names(.keytype2FilterMapping(withProteins = withProteins)))
}
## returns a vector mapping keytypes (names of vector) to filter names (elements).
.keytype2FilterMapping <- function(withProteins = FALSE){
    filters <- c(ENTREZID = "EntrezFilter",
                 GENEID = "GeneIdFilter",
                 GENEBIOTYPE = "GeneBiotypeFilter",
                 GENENAME = "GeneNameFilter",
                 TXID = "TxIdFilter",
                 TXBIOTYPE = "TxBiotypeFilter",
                 EXONID = "ExonIdFilter",
                 SEQNAME = "SeqNameFilter",
                 SEQSTRAND = "SeqStrandFilter",
                 TXNAME = "TxIdFilter",
                 SYMBOL = "SymbolFilter")
    if (withProteins) {
        filters <- c(filters,
                     PROTEINID = "ProteinIdFilter",
                     UNIPROTID = "UniprotFilter",
                     PROTEINDOMAINID = "ProteinDomainIdFilter",
                     PROTDOMID = "ProtDomIdFilter",
                     PROTEINDOMAINSOURCE = "ProteinDomainSourceFilter")
    }
    return(filters)
}
filterForKeytype <- function(keytype, x, vals){
    if (missing(vals))
        vals <- 1
    if (!missing(x)) {
        withProts <- hasProteinData(x)
    } else {
        withProts <- FALSE
    }
    filters <- .keytype2FilterMapping(withProts)
    if(any(names(filters) == keytype)){
        filt <- do.call(filters[keytype], args = list(value = vals))
        ## filt <- new(filters[keytype])
        return(filt)
    }else{
        stop("No filter for that keytype!")
    }
}

####============================================================
##  keys method
##
##  This keys method returns all of the keys for a specified keytype.
##  There should also be an implementation without keytypes, which
##  returns in our case the gene_ids
##
####------------------------------------------------------------
setMethod("keys", "EnsDb",
          function(x, keytype, filter, ...){
              if(missing(keytype))
                  keytype <- "GENEID"
              if(missing(filter))
                  filter <- AnnotationFilterList()
              filter <- .processFilterParam(filter, x)
              keyt <- keytypes(x)
              if (length(keytype) > 1) {
                  keytype <- keytype[1]
                  warning("Using only first provided keytype.")
              }
              if (!any(keyt == keytype))
                  stop("keytype '", keytype, "' not supported! ",
                       "Allowed choices are: ",
                       paste0("'", keyt ,"'", collapse = ", "), ".")
              keytype <- match.arg(keytype, keyt)
              ## Map the keytype to the appropriate column name.
              dbColumn <- ensDbColumnForColumn(x, keytype)
              ## Perform the query.
              res <- getWhat(x, columns = dbColumn, filter = filter)[, dbColumn]
              return(res)
          })


############################################################
## select method
##
##  We have to be carefull, if the database contains protein annotations too:
##  o If the keys are DNA/RNA related, start from a DNA/RNA related table.
##  o if keys are protein related: start from a protein column.
##  Reason is that we do have only protein annotations for protein coding genes
##  and no annotation for the remaining. Thus the type of the join (left join,
##  left outer join) is crucial, as well as the table with which we start the
##  query!
##  What if we provide more than one filter?
##  a) GeneNameFilter and ProteinidFilter: doesn't really matter from which table
##     we start, because the query will only return results with protein
##     annotions. -> if there is one DNA/RNA related filter: don't do anything.
##  b) Only protein filters: start from the highest protein table.
setMethod("select", "EnsDb",
          function(x, keys, columns, keytype, ...) {
              if (missing(keys))
                  keys <- NULL
              if (missing(columns))
                  columns <- NULL
              if (missing(keytype))
                  keytype <- NULL
              return(.select(x = x, keys = keys, columns = columns,
                             keytype = keytype, ...))
          })
.select <- function(x, keys = NULL, columns = NULL, keytype = NULL, ...) {
    extraArgs <- list(...)
    ## Perform argument checking:
    ## columns:
    if (missing(columns) | is.null(columns))
        columns <- columns(x)
    notAvailable <- !(columns %in% columns(x))
    if (all(notAvailable))
        stop("None of the specified columns are avaliable in the database!")
    if (any(notAvailable)){
        warning("The following columns are not available in the database and",
                " have thus been removed: ",
                paste(columns[notAvailable], collapse = ", "))
        columns <- columns[!notAvailable]
    }
    ## keys:
    if (is.null(keys) | missing(keys)) {
        ## Get everything from the database...
        keys <- list()
    } else {
        if (!(is(keys, "character") | is(keys, "list") | is(keys, "formula") |
              is(keys, "AnnotationFilter") | is(keys, "AnnotationFilterList")))
            stop("Argument keys should be a character vector, an object",
                 " extending AnnotationFilter, a filter expression",
                 " or an AnnotationFilterList.")
        if (is(keys, "character")) {
            if (is.null(keytype)) {
                stop("Argument keytype is mandatory if keys is a",
                     " character vector!")
            }
            ## Check also keytype:
            if (!(keytype %in% keytypes(x)))
                stop("keytype ", keytype, " not available in the database.",
                     " Use keytypes method to list all available keytypes.")
            ## Generate a filter object for the filters.
            keyFilter <- filterForKeytype(keytype, x, vals = keys)
            ## value(keyFilter) <- keys
            ## keyFilter@value <- keys
            keys <- list(keyFilter)
            ## Add also the keytype itself to the columns.
            if (!any(columns == keytype))
                columns <- c(keytype, columns)
        }
        ## Check and fix filter.
        keys <- .processFilterParam(keys, x)
    }
    ## Map the columns to column names we have in the database and
    ## add filter columns too.
    ensCols <- unique(c(ensDbColumnForColumn(x, columns),
                        addFilterColumns(character(), filter = keys, x)))
    ## TODO @jo: Do we have to check that we are allowed to have protein filters
    ##           or columns?
    ## OK, now perform the query given the filters we've got.
    ## Check if keys does only contain protein annotation columns; in that case
    ## select one of tables "protein", "uniprot", "protein_domain" in that order
    ## if (all(unlist(lapply(keys, isProteinFilter)))) {
    if (all(isProteinFilter(keys))) {
        startWith <- "protein_domain"
        if (any(unlist(lapply(keys, function(z) is(z, "UniprotFilter")))))
            startWith <- "uniprot"
        if (any(unlist(lapply(keys, function(z) is(z, "ProteinIdFilter")))))
            startWith <- "protein"
    } else {
        startWith <- NULL
    }
    ## Otherwise set startWith to NULL
    res <- getWhat(x, columns = ensCols, filter = keys, startWith = startWith)
    ## Order results if length of filters is 1.
    if (length(keys) == 1) {
        ## Define the filters on which we could sort.
        sortFilts <- c("GenenameFilter", "GeneNameFilter", "GeneIdFilter",
                       "EntrezFilter", "GeneBiotypeFilter", "SymbolFilter",
                       "TxIdFilter", "TxBiotypeFilter", "ExonIdFilter",
                       "ExonRankFilter", "SeqNameFilter")
        if (class(keys[[1]]) %in% sortFilts) {
            keyvals <- value(keys[[1]])
            ## Handle symlink Filter differently:
            if (is(keys[[1]], "SymbolFilter")) {
                ## sortCol <- ensDbColumn(keys[[1]])
                sortCol <- keys[[1]]@field
            } else {
                sortCol <- ensDbColumn(keys[[1]])
                ## sortCol <- removePrefix(ensDbColumn(keys[[1]], x))
            }
            res <- res[order(match(res[, sortCol], keyvals)), ]
        }
    } else {
        ## Show a mild warning message
        message(paste0("Note: ordering of the results might not match ordering",
                       " of keys!"))
    }
    colMap <- .getColMappings(x)
    colnames(res) <- colMap[colnames(res)]
    rownames(res) <- NULL
    if (returnFilterColumns(x))
        return(res)
    res[, columns]
}

############################################################
##  mapIds method
##
##  maps the submitted keys (names of the returned vector) to values
##  of the column specified by column.
##  x, key, column, keytype, ..., multiVals
setMethod("mapIds", "EnsDb", function(x, keys, column, keytype, ..., multiVals) {
    if(missing(keys))
        keys <- NULL
    if(missing(column))
        column <- NULL
    if(missing(keytype))
        keytype <- NULL
    if(missing(multiVals))
        multiVals <- NULL
    return(.mapIds(x = x, keys = keys, column = column, keytype = keytype,
                   multiVals = multiVals, ...))
})
## Other methods: saveDb, species, dbfile, dbconn, taxonomyId
.mapIds <- function(x, keys = NULL, column = NULL, keytype = NULL, ...,
                    multiVals = NULL) {
    if (is.null(keys))
        stop("Argument keys has to be provided!")
    ## if (!(is(keys, "character") | is(keys, "list") |
    ##       is(keys, "AnnotationFilter")))
    ##     stop("Argument keys should be a character vector, an object extending",
    ##          " AnnotationFilter or a list of objects extending AnnotationFilter.")
    if (is.null(column))
        column <- "GENEID"
    ## Have to specify the columns argument. Has to be keytype and column.
    if (is(keys, "character")){
        if (is.null(keytype))
            stop("Argument keytype is mandatory if keys is a character vector!")
        columns <- c(keytype, column)
    } else {
        ## Test if we can convert the filter. Returns ALWAYS an
        ## AnnotationFilterList
        keys <- .processFilterParam(keys, x)
        if(length(keys) > 1)
            warning("Got ", length(keys), " filter objects.",
                    " Will use the keys of the first for the mapping!")
        cn <- class(keys[[1]])[1]
        ## Use the first element to determine the keytype...
        mapping <- .keytype2FilterMapping()
        columns <- c(names(mapping)[mapping == cn], column)
        keytype <- NULL
    }
    ## if(is(keys, "list") | is(keys, "AnnotationFilter")){
    ##     if(is(keys, "list")){
    ##         if(length(keys) > 1)
    ##             warning("Got ", length(keys), " filter objects.",
    ##                     " Will use the keys of the first for the mapping!")
    ##         cn <- class(keys[[1]])[1]
    ##     }else{
    ##         cn <- class(keys)[1]
    ##     }
    ##     ## Use the first element to determine the keytype...
    ##     mapping <- .keytype2FilterMapping()
    ##     columns <- c(names(mapping)[mapping == cn], column)
    ##     keytype <- NULL
    ## }
    res <- select(x, keys = keys, columns = columns, keytype = keytype)
    if(nrow(res) == 0)
        return(character())
    ## Handling multiVals.
    if(is.null(multiVals))
        multiVals <- "first"
    if(is(multiVals, "function"))
        stop("Not yet implemented!")
    if(is.character(keys)){
        theNames <- keys
    }else{
        theNames <- unique(res[, 1])
    }
    switch(multiVals,
           first={
               vals <- res[match(theNames, res[, 1]), 2]
               names(vals) <- theNames
               return(vals)
           },
           list={
               ## vals <- split(res[, 2], f=factor(res[, 1], levels=unique(res[, 1])))
               vals <- split(res[, 2], f=factor(res[, 1], levels=unique(theNames)))
               return(vals)
           },
           filter={
               vals <- split(res[, 2], f=factor(res[, 1], levels=unique(theNames)))
               vals <- vals[unlist(lapply(vals, length)) == 1]
               return(unlist(vals))
           },
           asNA={
               ## Split the vector, set all those with multi mappings NA.
               vals <- split(res[, 2], f=factor(res[, 1], levels=unique(theNames)))
               vals[unlist(lapply(vals, length)) > 1] <- NA
               return(unlist(vals))
           },
           CharacterList={
               f <- factor(res[, 1], levels=unique(theNames))
               vals <- splitAsList(res[, 2], f=f)
               return(vals)
           })
}



