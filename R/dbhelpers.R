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
        Tables[[ i ]] <- colnames(dbGetQuery(con, paste0("select * from ",
                                                         tables[ i ], " limit 1")))
    }
    names(Tables) <- tables
    EDB <- new("EnsDb", ensdb=con, tables=Tables)
    EDB <- setProperty(EDB, dbSeqlevelsStyle="Ensembl")
    ## Setting the default for the returnFilterColumns
    returnFilterColumns(EDB) <- TRUE
    ## Defining the default for the ordering
    orderResultsInR(EDB) <- FALSE
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
prefixColumns <- function(x, columns, clean = TRUE, with.tables){
    if (missing(columns))
        stop("columns is empty! No columns provided!")
    ## first get to the tables that contain these columns
    Tab <- listTables(x)   ## returns the tables by degree!
    if (!missing(with.tables)) {
        with.tables <- with.tables[ with.tables %in% names(Tab) ]
        if (length(with.tables) > 0) {
            Tab <- Tab[ with.tables ]
        } else {
            warning("The submitted table names are not valid in the database and were thus dropped.")
        }
        if (length(Tab) == 0)
            stop("None of the tables submitted with with.tables is present in the database!")
    }
    if (clean)
        columns <- cleanColumns(x, columns)
    if (length(columns) == 0) {
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

############################################################
## call the prefixColumns function and return just the column
## names, but in the same order than the provided columns.
prefixColumnsKeepOrder <- function(x, columns, clean = TRUE, with.tables) {
    res <- unlist(prefixColumns(x, columns, clean, with.tables),
                  use.names = FALSE)
    res_order <- sapply(columns, function(z) {
        idx <- grep(res, pattern = paste0("\\.", z, "$"))
        if (length(idx) == 0)
            return(NULL)
        return(res[idx[1]])
    })
    return(res_order[!is.null(res_order)])
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


############################################################
## .buildQuery
##
## The "backbone" function that builds the SQL query based on the specified
## columns, the provided filters etc.
## x an EnsDb object
.buildQuery <- function(x, columns, filter = list(), order.by = "",
                        order.type = "asc", group.by, skip.order.check=FALSE,
                        return.all.columns = TRUE) {
    resultcolumns <- columns    ## just to remember what we really want to give back
    ## 1) get all column names from the filters also removing the prefix.
    if (class(filter)!="list")
        stop("parameter filter has to be a list of BasicFilter classes!")
    if (length(filter) > 0) {
        ## check filter!
        ## add the columns needed for the filter
        filtercolumns <- unlist(lapply(filter, column, x))
        ## remove the prefix (column name for these)
        filtercolumns <- sapply(filtercolumns, removePrefix, USE.NAMES = FALSE)
        columns <- unique(c(columns, filtercolumns))
    }
    ## 2) get all column names for the order.by:
    if (any(order.by != "")) {
        ## if we have skip.order.check set we use the order.by as is.
        if (!skip.order.check) {
            order.by <- checkOrderBy(orderBy = order.by, supported = columns)
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
    if (length(filter) > 0) {
        filterquery <- paste(" where",
                             paste(unlist(lapply(filter, where, x,
                                                 with.tables = need.tables)),
                                   collapse=" and "))
    } else {
        filterquery <- ""
    }
    ## c) the order part of the query
    if (any(order.by != "")) {
        if (!skip.order.check) {
            ## order.by <- paste(unlist(prefixColumns(x=x, columns=order.by,
            ##                                        with.tables=need.tables),
            ##                          use.names=FALSE), collapse=",")
            order.by <- paste(prefixColumnsKeepOrder(x = x, columns = order.by,
                                                     with.tables = need.tables),
                              collapse=",")
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
                         ## paste(unlist(prefixColumns(x,
                         ##                            resultcolumns,
                         ##                            with.tables=need.tables),
                         ##              use.names=FALSE), collapse=","),
                         paste(prefixColumnsKeepOrder(x,
                                                      resultcolumns,
                                                      with.tables = need.tables),
                               collapse=","),
                         " from ",
                         joinquery,
                         filterquery,
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
.getWhat <- function(x, columns, filter = list(), order.by = "",
                     order.type = "asc", group.by = NULL,
                     skip.order.check = FALSE) {
    ## That's nasty stuff; for now we support the column tx_name, which we however
    ## don't have in the database. Thus, we are querying everything except that
    ## column and filling it later with the values from tx_id.
    fetchColumns <- columns
    if(any(columns == "tx_name"))
        fetchColumns <- unique(c("tx_id", fetchColumns[fetchColumns != "tx_name"]))
    if (class(filter) != "list")
        stop("parameter filter has to be a list of BasicFilter classes!")
    ## If any of the filter is a SymbolFilter, add "symbol" to the return columns.
    if (length(filter) > 0) {
        if (any(unlist(lapply(filter, function(z) {
            return(is(z, "SymbolFilter"))
        }))))
            columns <- unique(c(columns, "symbol"))  ## append a filter column.
    }
    ## Catch also a "symbol" in columns
    if(any(columns == "symbol"))
        fetchColumns <- unique(c(fetchColumns[fetchColumns != "symbol"],
                                 "gene_name"))
    ## Shall we do the ordering in R or in SQL?
    if (orderResultsInR(x) & !skip.order.check) {
        ## Build the query
        Q <- .buildQuery(x = x, columns = fetchColumns, filter = filter,
                         order.by = "", order.type = order.type,
                         group.by = group.by,
                         skip.order.check = skip.order.check)
        ## Get the data
        Res <- dbGetQuery(dbconn(x), Q)
        ## Note: we can only order by the columns that we did get back from the
        ## database; that might be different for the SQL sorting!
        Res <- orderDataFrameBy(Res, by = checkOrderBy(order.by, fetchColumns),
                                decreasing = order.type != "asc")
    } else {
        ## Build the query
        Q <- .buildQuery(x = x, columns = fetchColumns, filter = filter,
                         order.by = order.by, order.type = order.type,
                         group.by = group.by,
                         skip.order.check = skip.order.check)
        ## Get the data
        Res <- dbGetQuery(dbconn(x), Q)
    }
    ## cat("Query:\n", Q, "\n")
    if(any(columns == "tx_cds_seq_start")) {
        if (!is.integer(Res[, "tx_cds_seq_start"])) {
            suppressWarnings(
                ## column contains "NULL" if not defined and coordinates are
                ## characters as.numeric transforms "NULL" into NA, and ensures
                ## coords are numeric.
                Res[ , "tx_cds_seq_start"] <- as.integer(Res[ , "tx_cds_seq_start"])
            )
        }
    }
    if(any(columns=="tx_cds_seq_end")){
        if (!is.integer(Res[, "tx_cds_seq_end"])) {
            suppressWarnings(
                ## column contains "NULL" if not defined and coordinates are
                ## characters as.numeric transforms "NULL" into NA, and ensures
                ## coords are numeric.
                Res[ , "tx_cds_seq_end" ] <- as.integer(Res[ , "tx_cds_seq_end" ])
            )
        }
    }
    ## Fix for MySQL returning 'numeric' instead of 'integer'.
    ## THIS SHOULD BE REMOVED ONCE THE PROBLEM IS FIXED IN RMySQL!!!
    int_cols <- c("exon_seq_start", "exon_seq_end", "exon_idx", "tx_seq_start",
                  "tx_seq_end", "tx_cds_seq_start", "tx_cds_seq_end",
                  "gene_seq_start", "gene_seq_end", "seq_strand")
    for (the_col in int_cols) {
        if (any(colnames(Res) == the_col))
            if (!is.integer(Res[, the_col]))
                Res[, the_col] <- as.integer(Res[, the_col])
    }
    ## Resolving the "symlinks" again.
    if(any(columns == "tx_name")) {
        Res <- data.frame(Res, tx_name = Res$tx_id, stringsAsFactors = FALSE)
    }
    if(any(columns == "symbol")) {
        Res <- data.frame(Res, symbol = Res$gene_name, stringsAsFactors = FALSE)
    }
    ## Ensure that the ordering is as requested.
    Res <- Res[, columns, drop=FALSE]
    return(Res)
}

############################################################
## Check database validity.
.ENSDB_TABLES <- list(gene = c("gene_id", "gene_name", "entrezid",
                               "gene_biotype", "gene_seq_start",
                               "gene_seq_end", "seq_name", "seq_strand",
                               "seq_coord_system"),
                      tx = c("tx_id", "tx_biotype", "tx_seq_start",
                             "tx_seq_end", "tx_cds_seq_start",
                             "tx_cds_seq_end", "gene_id"),
                      tx2exon = c("tx_id", "exon_id", "exon_idx"),
                      exon = c("exon_id", "exon_seq_start", "exon_seq_end"),
                      chromosome = c("seq_name", "seq_length", "is_circular"),
                      metadata = c("name", "value"))
dbHasRequiredTables <- function(con, returnError = TRUE) {
    tabs <- dbListTables(con)
    if (length(tabs) == 0) {
        if (returnError)
            return("Database does not have any tables!")
        return(FALSE)
    }
    not_there <- names(.ENSDB_TABLES)[!(names(.ENSDB_TABLES) %in% tabs)]
    if (length(not_there) > 0) {
        if (returnError)
            return(paste0("Required tables ", paste(not_there, collapse = ", "),
                          " are not present in the database!"))
        return(FALSE)
    }
    return(TRUE)
}
dbHasValidTables <- function(con, returnError = TRUE) {
    for (tab in names(.ENSDB_TABLES)) {
        cols <- .ENSDB_TABLES[[tab]]
        from_db <- colnames(dbGetQuery(con, paste0("select * from ", tab,
                                                   " limit 1")))
        not_there <- cols[!(cols %in% from_db)]
        if (length(not_there) > 0) {
            if (returnError)
                return(paste0("Table ", tab, " is missing required columns ",
                              paste(not_there, collapse = ", "), "!"))
            return(FALSE)
        }
    }
    return(TRUE)
}

############################################################
## feedEnsDb2MySQL
##
##
feedEnsDb2MySQL <- function(x, mysql, verbose = TRUE) {
    if (!inherits(mysql, "MySQLConnection"))
        stop("'mysql' is supposed to be a connection to a MySQL database.")
    ## Fetch the tables and feed them to MySQL.
    sqlite_con <- dbconn(x)
    tabs <- dbListTables(sqlite_con)
    for (the_table in tabs) {
        if (verbose)
            message("Fetch table ", the_table, "...", appendLF = FALSE)
        tmp <- dbGetQuery(sqlite_con, paste0("select * from ", the_table, ";"))
        if (verbose)
            message("OK\nStoring the table in MySQL...", appendLF = FALSE)
        ## Fix tx_cds_seq_start being a character in old databases
        if (any(colnames(tmp) == "tx_cds_seq_start")) {
            suppressWarnings(
                tmp[, "tx_cds_seq_start"] <- as.integer(tmp[, "tx_cds_seq_start"])
            )
            suppressWarnings(
                tmp[, "tx_cds_seq_end"] <- as.integer(tmp[, "tx_cds_seq_end"])
            )
        }
        dbWriteTable(mysql, tmp, name = the_table, row.names = FALSE)
        if (verbose)
            message("OK")
    }
    ## Create the indices.
    if (verbose)
        message("Creating indices...", appendLF = FALSE)
    .createEnsDbIndices(mysql, indexLength = "(20)")
    if (verbose)
        message("OK")
    return(TRUE)
}
## Small helper function to cfeate all the indices.
.createEnsDbIndices <- function(con, indexLength = "") {
    dbGetQuery(con, paste0("create index seq_name_idx on chromosome (seq_name",
                           indexLength, ");"))
    dbGetQuery(con, paste0("create index gene_gene_id_idx on gene (gene_id",
                           indexLength, ");"))
    dbGetQuery(con, paste0("create index gene_gene_name_idx on gene (gene_name",
                           indexLength, ");"))
    dbGetQuery(con, paste0("create index gene_seq_name_idx on gene (seq_name",
                           indexLength, ");"))
    dbGetQuery(con, paste0("create index tx_tx_id_idx on tx (tx_id",
                           indexLength, ");"))
    dbGetQuery(con, paste0("create index tx_gene_id_idx on tx (gene_id",
                           indexLength, ");"))
    dbGetQuery(con, paste0("create index exon_exon_id_idx on exon (exon_id",
                           indexLength, ");"))
    dbGetQuery(con, paste0("create index t2e_tx_id_idx on tx2exon (tx_id",
                           indexLength, ");"))
    dbGetQuery(con, paste0("create index t2e_exon_id_idx on tx2exon (exon_id",
                           indexLength, ");"))
    dbGetQuery(con, "create index t2e_exon_idx_idx on tx2exon (exon_idx);")
}
