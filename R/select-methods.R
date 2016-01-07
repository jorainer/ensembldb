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
              return(.filterKeytypes())
          }
)
## This just returns some (eventually) usefull names for keys
.simpleKeytypes <- function(x){
    return(c("GENEID","TXID","TXNAME","EXONID","EXONNAME","CDSID","CDSNAME"))
}
.filterKeytypes <- function(x){
    return(names(.keytype2FilterMapping()))
}
## returns a vector mapping keytypes (names of vector) to filter names (elements).
.keytype2FilterMapping <- function(){
    filters <- c("EntrezidFilter", "GeneidFilter", "GenebiotypeFilter", "GenenameFilter",
                 "TxidFilter", "TxbiotypeFilter", "ExonidFilter", "SeqnameFilter",
                 "SeqstrandFilter")
    names(filters) <- c("ENTREZID", "GENEID", "GENEBIOTYPE", "GENENAME", "TXID",
                        "TXBIOTYPE", "EXONID", "SEQNAME", "SEQSTRAND")
    return(filters)
}
filterForKeytype <- function(keytype){
    filters <- .keytype2FilterMapping()
    if(any(names(filters) == keytype)){
        filt <- new(filters[keytype])
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
          function(x, keytype, filter,...){
              if(missing(keytype))
                  keytype <- "GENEID"
              if(missing(filter))
                  filter <- list()
              if(is(filter, "BasicFilter"))
                  filter <- list(filter)
              keyt <- keytypes(x)
              keytype <- match.arg(keytype, keyt)
              ## Map the keytype to the appropriate column name.
              dbColumn <- ensDbColumnForColumn(x, keytype)
              ## Perform the query.
              res <- getWhat(x, columns=dbColumn, filter=filter)[, dbColumn]
              return(res)
          })


####============================================================
##  select method
##
##
####------------------------------------------------------------
setMethod("select", "EnsDb",
          function(x, keys, columns, keytype, ...){
              if(missing(keys))
                  keys <- NULL
              if(missing(columns))
                  columns <- NULL
              if(missing(keytype))
                  keytype <- NULL
              return(.select(x=x, keys=keys, columns=columns, keytype=keytype, ...))
          })
.select <- function(x, keys=NULL, columns=NULL, keytype=NULL, ...){
    extraArgs <- list(...)
    ## Perform argument checking:
    ## columns:
    if(missing(columns) | is.null(columns))
        columns <- columns(x)
    notAvailable <- !(columns %in% columns(x))
    if(all(notAvailable))
        stop("None of the specified columns are avaliable in the database!")
    if(any(notAvailable)){
        warning("The following columns are not available in the database and have",
                " thus been removed: ", paste(columns[notAvailable], collapse=", "))
        columns <- columns[!notAvailable]
    }
    ## keys:
    if(is.null(keys) | missing(keys)){
        ## Get everything from the database...
        keys <- list()
    }else{
        if(!(is(keys, "character") | is(keys, "list") | is(keys, "BasicFilter")))
            stop("Argument keys should be a character vector, an object extending BasicFilter ",
                 "or a list of objects extending BasicFilter.")
        if(is(keys, "list")){
            if(!all(vapply(keys, is, logical(1L), "BasicFilter")))
                stop("If keys is a list it should be a list of objects extending BasicFilter!")
        }
        if(is(keys, "BasicFilter")){
            keys <- list(keys)
        }
        if(is(keys, "character")){
            if(is.null(keytype)){
                stop("Argument keytype is mandatory if keys is a character vector!")
            }
            ## Check also keytype:
            if(!(keytype %in% keytypes(x)))
                stop("keytype ", keytype, " not available in the database.",
                     " Use keytypes method to list all available keytypes.")
            ## Generate a filter object for the filters.
            keyFilter <- filterForKeytype(keytype)
            value(keyFilter) <- keys
            keys <- list(keyFilter)
            ## Add also the keytype itself to the columns.
            if(!any(columns == keytype))
                columns <- c(keytype, columns)
        }
    }
    ## Map the columns to column names we have in the database.
    ensCols <- ensDbColumnForColumn(x, columns)
    ## OK, now perform the query given the filters we've got.
    res <- getWhat(x, columns=ensCols, filter=keys)
    colMap <- .getColMappings(x)
    colnames(res) <- colMap[colnames(res)]
    return(res[, columns])
}


####============================================================
##  mapIds method
##
##  maps the submitted keys (names of the returned vector) to values
##  of the column specified by column.
##  x, key, column, keytype, ..., multiVals
####------------------------------------------------------------
setMethod("mapIds", "EnsDb", function(x, keys, column, keytype, ..., multiVals){
    if(missing(keys))
        keys <- NULL
    if(missing(column))
        column <- NULL
    if(missing(keytype))
        keytype <- NULL
    if(missing(multiVals))
        multiVals <- NULL
    return(.mapIds(x=x, keys=keys, column=column, keytype=keytype, multiVals=multiVals, ...))
})
## Other methods: saveDb, species, dbfile, dbconn, taxonomyId
.mapIds <- function(x, keys=NULL, column=NULL, keytype=NULL, ..., multiVals=NULL){
    if(is.null(keys))
        stop("Argument keys has to be provided!")
    if(!(is(keys, "character") | is(keys, "list") | is(keys, "BasicFilter")))
        stop("Argument keys should be a character vector, an object extending BasicFilter ",
             "or a list of objects extending BasicFilter.")
    if(is.null(column))
        column <- "GENEID"
    ## Have to specify the columns argument. Has to be keytype and column.
    if(is(keys, "character")){
        if(is.null(keytype))
            stop("Argument keytype is mandatory if keys is a character vector!")
        columns <- c(keytype, column)
    }
    if(is(keys, "list") | is(keys, "BasicFilter")){
        if(is(keys, "list")){
            if(length(keys) > 1)
                warning("Got ", length(keys), " filter objects.",
                        " Will use the keys of the first for the mapping!")
            cn <- class(keys[[1]])[1]
        }else{
            cn <- class(keys)[1]
        }
        ## Use the first element to determine the keytype...
        mapping <- .keytype2FilterMapping()
        columns <- c(names(mapping)[mapping == cn], column)
        keytype <- NULL
    }
    res <- select(x, keys=keys, columns=columns, keytype=keytype)
    if(nrow(res) == 0)
        return(character())
    ## Handling multiVals.
    if(is.null(multiVals))
        multiVals <- "first"
    if(is(multiVals, "function"))
        stop("Not yet implemented!")
    ## Eventually re-order the data.frame in the same order than the keys...
    ## That's amazingly slow!!!
    ## if(is.character(keys)){
    ##     res <- split(res, f=factor(res[, 1], levels=keys))
    ##     res <- do.call(rbind, res)
    ##     rownames(res) <- NULL
    ## }
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
               stop("Not yet implemented!")
           })
}



