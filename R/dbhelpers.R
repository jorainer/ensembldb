
## x... optional, the database file name.
## v... verbosity
## function to load the database...
EnsDb <- function(x){
    options(useFancyQuotes=FALSE)
    if(missing(x)){
        stop("No sqlite file provided!")
    }
    lite <- dbDriver("SQLite")
    con <- dbConnect(lite, dbname = x, flags=SQLITE_RO)
    tables <- dbListTables(con)
    ## read the columns for these tables.
    Tables <- vector(length=length(tables), "list")
    for(i in 1:length(Tables)){
        Tables[[ i ]] <- colnames(dbGetQuery(con, paste0("select * from ", tables[ i ], " limit 1")))
    }
    names(Tables) <- tables
    EDB <- new("EnsDb", ensdb=con, tables=Tables)
    return(EDB)
}

## x is the connection to the database, name is the name of the entry to fetch
.getMetaDataValue <- function(x, name){
    return(dbGetQuery(x, paste0("select value from metadata where name='", name, "'"))[ 1, 1])
}

####
## Note: that's the central function that checks which tables are needed for the
## least expensive join!!! The names of the tables should then also be submitted
## to any other method that calls prefixColumns (e.g. where of the Filter classes)
##
## this function checks:
## a) for multi-table columns, selects the table with the highest degree
## b) pre-pend (inverse of append ;)) the table name to the column name.
## returns a list, names being the tables and the values being the columns
## named: <table name>.<column name>
## clean: whether a cleanColumns should be called on the submitted columns.
## with.tables: force the prefix to be specifically on the submitted tables.
prefixColumns <- function(x, columns, clean=TRUE, with.tables){
    if(missing(columns))
        stop("columns is empty! No columns provided!")
    ## first get to the tables that contain these columns
    Tab <- listTables(x)   ## returns the tables by degree!
    if(!missing(with.tables)){
        with.tables <- with.tables[ with.tables %in% names(Tab) ]
        if(length(with.tables) > 0){
            Tab <- Tab[ with.tables ]
        }else{
            warning("The submitted table names are not valid in the database and were thus dropped.")
        }
        if(length(Tab)==0)
            stop("None of the tables submitted with with.tables is present in the database!")
    }
    if(clean)
        columns <- cleanColumns(x, columns)
    if(length(columns)==0){
        return(NULL)
    }
    ## group the columns by table.
    columns.bytable <- sapply(Tab, function(z){
        return(z[ z %in% columns ])
    }, simplify=FALSE, USE.NAMES=TRUE)
    ## kick out empty tables...
    columns.bytable <- columns.bytable[ unlist(lapply(columns.bytable, function(z){
        return(length(z) > 0)
    })) ]
    if(length(columns.bytable)==0)
        stop("No columns available!")
    have.columns <- NULL
    ## new approach! order the tables by number of elements, and after that, re-order them.
    columns.bytable <- columns.bytable[ order(unlist(lapply(columns.bytable, length)),
                                              decreasing=TRUE) ]
    ## has to be a for loop!!!
    ## loop throught the columns by table and sequentially kick out columns for the current table if they where already
    ## in a previous (more relevant) table
    ## however, prefer also cases were fewer tables are returned.
    for(i in 1:length(columns.bytable)){
        bm <- columns.bytable[[ i ]] %in% have.columns
        keepvals <- columns.bytable[[ i  ]][ !bm ]   ## keep those
        if(length(keepvals) > 0){
            have.columns <- c(have.columns, keepvals)
        }
        if(length(keepvals) > 0){
            columns.bytable[[ i ]] <- paste(names(columns.bytable)[ i ], keepvals, sep=".")
        }else{
            columns.bytable[[ i ]] <- keepvals
        }
    }
    ## kick out those tables with no elements left...
    columns.bytable <- columns.bytable[ unlist(lapply(columns.bytable, function(z){
        return(length(z) > 0)
    })) ]
    ## re-order by degree.
    columns.bytable <- columns.bytable[ tablesByDegree(x, names(columns.bytable)) ]
    return(columns.bytable)
}





## define a function to create a join query based on columns
## this function has to first get all tables that contain the columns,
## and then select, for columns present in more than one
## x... EnsDb
## columns... the columns
joinQueryOnColumns <- function(x, columns){
    columns.bytable <- prefixColumns(x, columns)
    ## based on that we can build the query based on the tables we've got. Note that the
    ## function internally
    ## adds tables that might be needed for the join.
    Query <- joinQueryOnTables(x, names(columns.bytable))
    return(Query)
}


## only list direct joins!!!
.JOINS <- rbind(
    c("gene", "tx", "join tx on (gene.gene_id=tx.gene_id)"),
    c("gene", "chromosome", "join chromosome on (gene.seq_name=chromosome.seq_name)"),
    c("tx", "tx2exon", "join tx2exon on (tx.tx_id=tx2exon.tx_id)"),
    c("tx2exon", "exon", "join exon on (tx2exon.exon_id=exon.exon_id)")
)
## tx is now no 1:
## .JOINS <- rbind(
##     c("tx", "gene", "join gene on (tx.gene_id=gene.gene_id)"),
##     c("gene", "chromosome", "join chromosome on (gene.seq_name=chromosome.seq_name)"),
##     c("tx", "tx2exon", "join tx2exon on (tx.tx_id=tx2exon.tx_id)"),
##     c("tx2exon", "exon", "join exon on (tx2exon.exon_id=exon.exon_id)")
##    )


joinQueryOnTables <- function(x, tab){
    ## just to be on the save side: evaluate whether we have all required tables to join;
    ## this will also ensure that the order is by degree.
    tab <- addRequiredTables(x, tab)
    Query <- tab[ 1 ]
    previous.table <- tab[ 1 ]
    for(i in 1:length(tab)){
        if(i > 1){
            Query <- paste(Query, .JOINS[ .JOINS[ , 2 ]==tab[ i ], 3 ])
        }
    }
    return(Query)
}


###
## Add additional tables in case the submitted tables are not directly connected
## and can thus not be joined. That's however not so complicated, since the database
## layout is pretty simple.
## The tables are:
##
##  exon -(exon_id=t2e_exon_id)- tx2exon -(t2e_tx_id=tx_id)- tx -(gene_id=gene_id)- gene
##                                                                                   |
##                                                   chromosome -(seq_name=seq_name)-´
addRequiredTables <- function(x, tab){
    ## dash it, as long as I can't find a way to get connected objects in a
    ## graph I'll do it manually...
    ## if we have exon and any other table, we need definitely tx2exon
    if(any(tab=="exon") & length(tab) > 1){
        tab <- unique(c(tab, "tx2exon"))
    }
    ## if we have chromosome and any other table, we'll need gene
    if(any(tab=="chromosome") & length(tab) > 1){
        tab <- unique(c(tab, "gene"))
    }
    ## if we have exon and we have gene, we'll need also tx
    if((any(tab=="exon") | (any(tab=="tx2exon"))) & any(tab=="gene")){
        tab <- unique(c(tab, "tx"))
    }
    return(tablesByDegree(x, tab))
}


.buildQuery <- function(x, columns, filter=list(), order.by="", order.type="asc",
                        group.by, skip.order.check=FALSE, return.all.columns=TRUE){
    resultcolumns <- columns    ## just to remember what we really want to give back
    ## 1) get all column names from the filters also removing the prefix.
    if(class(filter)!="list")
        stop("parameter filter has to be a list of BasicFilter classes!")
    if(length(filter) > 0){
        ## check filter!
        ## add the columns needed for the filter
        filtercolumns <- unlist(lapply(filter, column, x))
        ## remove the prefix (column name for these)
        filtercolumns <- sapply(filtercolumns, removePrefix, USE.NAMES=FALSE)
        columns <- unique(c(columns, filtercolumns))
    }
    ## 2) get all column names for the order.by:
    if(order.by!=""){
        ## if we have skip.order.check set we use the order.by as is.
        if(!skip.order.check){
            order.by <- unlist(strsplit(order.by, split=",", fixed=TRUE))
            order.by <- gsub(order.by, pattern=" ", replacement="", fixed=TRUE)
            ## allow only order.by that are also in the columns.
            order.by.nocolumns <- order.by[ !(order.by %in% columns) ]
            order.by <- order.by[ order.by %in% columns ]
            if(length(order.by.nocolumns) > 0){
                warning("columns provided in order.by (",
                        paste(order.by.nocolumns, collapse=","),
                        ") are not in columns and were thus removed.")
            }
            if(length(order.by)==0){
                order.by <- ""
            }
        }
    }else{
        order.by <- ""
    }
    ## Note: order by is now a vector!!!
    ## columns are now all columns that we want to fetch or that we need to
    ## filter or to sort.
    ## 3) check which tables we need for all of these columns:
    need.tables <- names(prefixColumns(x, columns))
    ##
    ## Now we can begin to build the query parts!
    ## a) the query part that joins all required tables.
    joinquery <- joinQueryOnColumns(x, columns=columns)
    ## b) the filter part of the query
    if(length(filter) > 0){
        filterquery <- paste(" where",
                             paste(unlist(lapply(filter, where, x,
                                                 with.tables=need.tables)),
                                   collapse=" and "))
    }else{
        filterquery <- ""
    }
    ## c) the order part of the query
    if(order.by!=""){
        if(!skip.order.check){
            order.by <- paste(unlist(prefixColumns(x=x, columns=order.by,
                                                   with.tables=need.tables),
                                     use.names=FALSE), collapse=",")
        }
        orderquery <- paste(" order by", order.by, order.type)
    }else{
        orderquery <- ""
    }
    ## And finally build the final query
    if(return.all.columns){
        resultcolumns <- columns
    }
    finalquery <- paste0("select distinct ",
                         paste(unlist(prefixColumns(x,
                                                    resultcolumns,
                                                    with.tables=need.tables),
                                      use.names=FALSE), collapse=","),
                         " from ",
                         joinquery,
                         filterquery,
                         orderquery
                         )
    return(finalquery)
}



## the buildQuery function; basically create the SQL query.
## x... EnsDb
## columns... the columns to retrieve.
## filter... list of filters
.buildQueryOld <- function(x, columns, filter=list(), order.by="", order.type="asc", group.by, skip.order.check=FALSE){
    resultcolumns <- columns    ## just to remember what we really want to give back
    if(!missing(filter)){
        if(class(filter)!="list")
            stop("parameter filter has to be a list of BasicFilter classes!")
    }
    ## first checking the filters and eventually add required columns to the columns:
    if(length(filter) > 0){
        ## check filter!
        ## add the columns needed for the filter
        filtercolumns <- unlist(lapply(filter, column, x))   ## this returns
        ## the full named columns!
        filtercolumns <- sapply(filtercolumns, removePrefix, USE.NAMES=FALSE)
        columns <- unique(c(columns, filtercolumns))
        ## next we're building the where query.
        wherequery <- paste(" where", paste(unlist(lapply(filter, where, x)),
                                            collapse=" and "))
    }else{
        wherequery <- ""
    }
    ## should we do an order.by?
    if(!missing(order.by) & order.by!=""){
        ## if we have skip.order.check set we use the order.by as is.
        if(!skip.order.check){
            order.by <- unlist(strsplit(order.by, split=",", fixed=TRUE))
            order.by <- gsub(order.by, pattern=" ", replacement="", fixed=TRUE)
            ## allow only order.by that are also in the columns.
            order.by.nocolumns <- order.by[ !(order.by %in% columns) ]
            order.by <- order.by[ order.by %in% columns ]
            if(length(order.by.nocolumns) > 0){
                warning("columns provided in order.by (",
                        paste(order.by.nocolumns, collapse=","),
                        ") are not in columns and were thus removed." )
            }
            if(length(order.by)==0){
                order.by <- ""
            }else{
                ## have to pre-pend the database table name...
                order.by <- paste(unlist(prefixColumns(x=x, columns=order.by), use.names=FALSE), collapse=",")
            }
        }
    }else{
        order.by <- ""
    }
    ## OK, build that query.
    if(order.by!=""){
        orderquery <- paste(" order by", order.by, order.type)
    }else{
        orderquery <- ""
    }
    ## now build the join query that joins all required tables.
    joinquery <- joinQueryOnColumns(x, columns=columns)
    finalquery <- paste0("select distinct ",
                         paste(unlist(prefixColumns(x, resultcolumns),
                                      use.names=FALSE), collapse=","),
                         " from ",
                         joinquery,
                         wherequery,
                         orderquery
                         )
    return(finalquery)
}

## remove the prefix again...
removePrefix <- function(x, split=".", fixed=TRUE){
    return(sapply(x, function(z){
        tmp <- unlist(strsplit(z, split=split, fixed=fixed))
        return(tmp[ length(tmp) ])
    }))
}


## just to add another layer; basically just calls buildQuery and executes the query
.getWhat <- function(x, columns, filter=list(), order.by="",
                     order.type="asc", group.by=NULL, skip.order.check=FALSE){
    ## build the query
    Q <- .buildQuery(x=x, columns=columns, filter=filter,
                     order.by=order.by, order.type=order.type, group.by=group.by,
                     skip.order.check=skip.order.check)
    ## get the data
    Res <- dbGetQuery(dbconn(x), Q)
    if(any(columns=="tx_cds_seq_start")){
        suppressWarnings(
            ## column contains "NULL" if not defined and coordinates are characters
            ## as.numeric transforms "NULL" into NA, and ensures coords are numeric.
            Res[ , "tx_cds_seq_start" ] <- as.numeric(Res[ , "tx_cds_seq_start" ])
        )
    }
    if(any(columns=="tx_cds_seq_end")){
        suppressWarnings(
            ## column contains "NULL" if not defined and coordinates are characters
            ## as.numeric transforms "NULL" into NA, and ensures coords are numeric.
            Res[ , "tx_cds_seq_end" ] <- as.numeric(Res[ , "tx_cds_seq_end" ])
        )
    }
    return(Res)
}

## Fetch data to add as a GeneTrack.
## filter ...                 Used to filter the result.
## chromosome, start, end ... Either all or none has to be specified. If specified, the function
##                            first retrieves all transcripts that have an exon in the specified
##                            range and adds them as a TranscriptidFilter to the filters. The
##                            query to fetch the "real" data is performed after.
## featureIs ...              Wheter gene_biotype or tx_biotype should be mapped to the column
##                            feature.
getGeneTrackDfForGviz <- function(x, filter=list(), chromosome=NULL, start=NULL, end=NULL,
                                  featureIs="gene_biotype"){
    featureIs <- match.arg(featureIs, c("gene_biotype", "tx_biotype"))
    if(any(!is.null(c(chromosome, start, end)))){
        if(all(!is.null(c(chromosome, start, end)))){
            ## Fetch all transcripts in that region:
            tids <- dbGetQuery(dbconn(x),
                               paste0("select distinct tx.tx_id from tx join gene on (tx.gene_id=gene.gene_id)",
                                      " where seq_name='", chromosome, "' and (",
                                      "(tx_seq_start >= ", start, " and tx_seq_start <= ", end, ") or ",
                                      "(tx_seq_end >= ", start, " and tx_seq_end <= ", end, ") or ",
                                      "(tx_seq_start <=", start, " and tx_seq_end >= ", end, ")",
                                      ")"))[, "tx_id"]
            if(length(tids) == 0)
                stop(paste0("Did not find any transcript on chromosome from ", start, " to ",
                            end, "!"))
            filter <- c(filter, TxidFilter(tids))
        }else{
            stop(paste0("Either all or none of arguments 'chromosome', 'start' and 'end' ",
                        " have to be specified!"))
        }
    }
    ## Return a data.frame with columns: chromosome, start, end, width, strand, feature,
    ## gene, exon, transcript and symbol.
    ## 1) Query the data as we usually would.
    ## 2) Perform an additional query to get cds and utr, remove all entries from the
    ##    first result for the same transcripts and rbind the data.frames.
    needCols <- c("seq_name", "exon_seq_start", "exon_seq_end", "seq_strand",
                  "featureIs", "gene_id", "exon_id",
                  "exon_idx", "tx_id", "gene_name")
    ## That's the names to which we map the original columns from the EnsDb.
    names(needCols) <- c("chromosome", "start", "end", "strand",
                         "feature", "gene", "exon", "exon_rank", "transcript",
                         "symbol")
    txs <- transcripts(x, filter=filter,
                       columns=needCols, return.type="data.frame")
    ## Rename columns
    idx <- match(needCols, colnames(txs))
    notThere <- is.na(idx)
    idx <- idx[!notThere]
    colnames(txs)[idx] <- names(needCols)[!notThere]
    ## now processing the 5utr
    fUtr <- fiveUTRsByTranscript(x, filter=filter, columns=needCols)
    if(length(fUtr) > 0){
        fUtr <- as(unlist(fUtr, use.names=FALSE), "data.frame")
        fUtr <- fUtr[, !(colnames(fUtr) %in% c("width", "seq_name", "exon_seq_start",
                                               "exon_seq_end", "strand"))]
        colnames(fUtr)[1] <- "chromosome"
        idx <- match(needCols, colnames(fUtr))
        notThere <- is.na(idx)
        idx <- idx[!notThere]
        colnames(fUtr)[idx] <- names(needCols)[!notThere]
        ## Force being in the correct ordering:
        fUtr <- fUtr[, names(needCols)]
        fUtr$feature <- "utr5"
        ## Remove transcripts from the txs data.frame
        txs <- txs[!(txs$transcript %in% fUtr$transcript), , drop=FALSE]
    }
    tUtr <- threeUTRsByTranscript(x, filter=filter, columns=needCols)
    if(length(tUtr) > 0){
        tUtr <- as(unlist(tUtr, use.names=FALSE), "data.frame")
        tUtr <- tUtr[, !(colnames(tUtr) %in% c("width", "seq_name", "exon_seq_start",
                                               "exon_seq_end", "strand"))]
        colnames(tUtr)[1] <- "chromosome"
        idx <- match(needCols, colnames(tUtr))
        notThere <- is.na(idx)
        idx <- idx[!notThere]
        colnames(tUtr)[idx] <- names(needCols)[!notThere]
        ## Force being in the correct ordering:
        tUtr <- tUtr[, names(needCols)]
        tUtr$feature <- "utr3"
        ## Remove transcripts from the txs data.frame
        if(nrow(txs) > 0){
            txs <- txs[!(txs$transcript %in% tUtr$transcript), , drop=FALSE]
        }
    }
    cds <- cdsBy(x, filter=filter, columns=needCols)
    if(length(cds) > 0){
        cds <- as(unlist(cds, use.names=FALSE), "data.frame")
        cds <- cds[, !(colnames(cds) %in% c("width", "seq_name", "exon_seq_start",
                                            "exon_seq_end", "strand"))]
        colnames(cds)[1] <- "chromosome"
        idx <- match(needCols, colnames(cds))
        notThere <- is.na(idx)
        idx <- idx[!notThere]
        colnames(cds)[idx] <- names(needCols)[!notThere]
        ## Force being in the correct ordering:
        cds <- cds[, names(needCols)]
        ## Remove transcripts from the txs data.frame
        if(nrow(txs) > 0){
            txs <- txs[!(txs$transcript %in% cds$transcript), , drop=FALSE]
        }
    }
    if(length(fUtr) > 0){
        txs <- rbind(txs, fUtr)
    }
    if(length(tUtr) > 0){
        txs <- rbind(txs, tUtr)
    }
    if(length(cds) > 0){
        txs <- rbind(txs, cds)
    }
    return(txs)
}





