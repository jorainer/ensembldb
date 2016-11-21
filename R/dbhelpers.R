############################################################
## EnsDb
## Constructor function.
##' @title Connect to an EnsDb object
##'
##' @description The \code{EnsDb} constructor function connects to the database
##' specified with argument \code{x} and returns a corresponding
##' \code{\linkS4class{EnsDb}} object.
##'
##' @details By providing the connection to a MySQL database, it is possible
##' to use MySQL as the database backend and queries will be performed on that
##' database. Note however that this requires the package \code{RMySQL} to be
##' installed. In addition, the user needs to have access to a MySQL server
##' providing already an EnsDb database, or must have write privileges on a
##' MySQL server, in which case the \code{\link{useMySQL}} method can be used
##' to insert the annotations from an EnsDB package into a MySQL database.
##' @param x Either a character specifying the \emph{SQLite} database file, or
##' a \code{DBIConnection} to e.g. a MySQL database.
##' @return A \code{\linkS4class{EnsDb}} object.
##' @author Johannes Rainer
##' @examples
##' ## "Standard" way to create an EnsDb object:
##' library(EnsDb.Hsapiens.v75)
##' EnsDb.Hsapiens.v75
##'
##' ## Alternatively, provide the full file name of a SQLite database file
##' dbfile <- system.file("extdata/EnsDb.Hsapiens.v75.sqlite", package = "EnsDb.Hsapiens.v75")
##' edb <- EnsDb(dbfile)
##' edb
##'
##' ## Third way: connect to a MySQL database
##' \dontrun{
##' library(RMySQL)
##' dbcon <- dbConnect(MySQL(), user = my_user, pass = my_pass, host = my_host, dbname = "ensdb_hsapiens_v75")
##' edb <- EnsDb(dbcon)
##' }
EnsDb <- function(x){
    options(useFancyQuotes=FALSE)
    if(missing(x)){
        stop("No sqlite file provided!")
    }
    if (is.character(x)) {
        lite <- dbDriver("SQLite")
        con <- dbConnect(lite, dbname = x, flags=SQLITE_RO)
    }
    else if (is(x, "DBIConnection")) {
        con <- x
    } else {
        stop("'x' should be either a character specifying the SQLite file to",
             " be loaded, or a DBIConnection object providing the connection",
             " to the database.")
    }
    ## Check if the database is valid.
    OK <- dbHasRequiredTables(con)
    if (is.character(OK))
        stop(OK)
    OK <- dbHasValidTables(con)
    if (is.character(OK))
        stop(OK)
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
    ## Check it again...
    OK <- validateEnsDb(EDB)
    if (is.character(OK))
        stop(OK)
    return(EDB)
}

## loadEnsDb <- function(x) {
##     ## con <- ensDb( x )
##     ## EDB <- new( "EnsDb", ensdb=con )
##     return(EnsDb(x))
## }


## x is the connection to the database, name is the name of the entry to fetch
.getMetaDataValue <- function(x, name){
    return(dbGetQuery(x, paste0("select value from metadata where name='",
                                name, "'"))[ 1, 1])
}

####
## THIS HAS BEEN REPLACED WITH prefixColumns!
##
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
## prefixColumns_old <- function(x, columns, clean = TRUE, with.tables){
##     if (missing(columns))
##         stop("columns is empty! No columns provided!")
##     ## first get to the tables that contain these columns
##     Tab <- listTables(x)   ## returns the tables by degree!
##     if (!missing(with.tables)) {
##         with.tables <- with.tables[ with.tables %in% names(Tab) ]
##         if (length(with.tables) > 0) {
##             Tab <- Tab[ with.tables ]
##         } else {
##             warning("The submitted table names are not valid in the database",
##                     " and were thus dropped.")
##         }
##         if (length(Tab) == 0)
##             stop("None of the tables submitted with with.tables is present",
##                  " in the database!")
##     }
##     if (clean)
##         columns <- cleanColumns(x, columns)
##     if (length(columns) == 0) {
##         return(NULL)
##     }
##     ## group the columns by table.
##     columns.bytable <- sapply(Tab, function(z){
##         return(z[ z %in% columns ])
##     }, simplify=FALSE, USE.NAMES=TRUE)
##     ## kick out empty tables...
##     columns.bytable <- columns.bytable[ unlist(lapply(columns.bytable,
##                                                       function(z){
##                                                           return(length(z) > 0)
##                                                       })) ]
##     if(length(columns.bytable)==0)
##         stop("No columns available!")
##     have.columns <- NULL
##     ## new approach! order the tables by number of elements, and after that,
##     ## re-order them.
##     columns.bytable <- columns.bytable[ order(unlist(lapply(columns.bytable,
##                                                             length)),
##                                               decreasing=TRUE) ]
##     ## has to be a for loop!!!
##     ## loop throught the columns by table and sequentially kick out columns
##     ## for the current table if they where already
##     ## in a previous (more relevant) table
##     ## however, prefer also cases were fewer tables are returned.
##     for(i in 1:length(columns.bytable)){
##         bm <- columns.bytable[[ i ]] %in% have.columns
##         keepvals <- columns.bytable[[ i  ]][ !bm ]   ## keep those
##         if(length(keepvals) > 0){
##             have.columns <- c(have.columns, keepvals)
##         }
##         if(length(keepvals) > 0){
##             columns.bytable[[ i ]] <- paste(names(columns.bytable)[ i ],
##                                             keepvals, sep=".")
##         }else{
##             columns.bytable[[ i ]] <- keepvals
##         }
##     }
##     ## kick out those tables with no elements left...
##     columns.bytable <- columns.bytable[ unlist(lapply(columns.bytable,
##                                                       function(z){
##                                                           return(length(z) > 0)
##     })) ]
##     ## re-order by degree.
##     columns.bytable <- columns.bytable[ tablesByDegree(x,
##                                                        names(columns.bytable)) ]
##     return(columns.bytable)
## }

############################################################
## prefixColumns
##
## Determines which tables (along with the table attributes) are required for
## the join.
## Updated version of prefixColumns:
## o Uses the order of the tables returned by listTables and adds the first
##   table in which the column was found. That's different to the previous
##   default of trying to join as few tables as possible but avoids problems
##   with table joins between e.g. tx and protein in which not all tx_id are
##   present in the protein table.
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
            warning("The submitted table names are not valid in the database",
                    " and were thus dropped.")
        }
        if (length(Tab) == 0)
            stop("None of the tables submitted with with.tables is present",
                 " in the database!")
    }
    if (clean)
        columns <- cleanColumns(x, columns)
    if (length(columns) == 0) {
        return(NULL)
    }
    getCols <- columns
    result <- lapply(Tab, function(z) {
        if (length(getCols) > 0) {
            gotIt <- z[z %in% getCols]
            if (length(gotIt) > 0) {
                getCols <<- getCols[!(getCols %in% gotIt)]
                return(gotIt)
            } else {
                return(character())
            }
        }
    })
    ## If getCols length > 0 it contains columns not present in the db.
    result <- result[lengths(result) > 0]
    if (length(result) == 0)
        stop("None of the columns could be found in the database!")
    result <- mapply(result, names(result), FUN = function(z, y) {
        paste0(y, ".", z)
    }, SIMPLIFY = FALSE)
    return(result)
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
## joinQueryOnColumns <- function(x, columns, join = "join"){
##     columns.bytable <- prefixColumns(x, columns)
##     ## based on that we can build the query based on the tables we've got. Note that the
##     ## function internally
##     ## adds tables that might be needed for the join.
##     Query <- joinQueryOnTables(x, names(columns.bytable))
##     return(Query)
## }


## only list direct joins!!!
## .JOINS <- rbind(
##     c("gene", "tx", "tx on (gene.gene_id=tx.gene_id)"),
##     c("gene", "chromosome",
##       "chromosome on (gene.seq_name=chromosome.seq_name)"),
##     c("tx", "tx2exon", "tx2exon on (tx.tx_id=tx2exon.tx_id)"),
##     c("tx2exon", "exon", "exon on (tx2exon.exon_id=exon.exon_id)"),
##     c("tx", "protein", "protein on (tx.tx_id=protein.tx_id)"),
##     c("protein", "uniprot",
##       "uniprot on (protein.protein_id=uniprot.protein_id)"),
##     c("protein", "protein_domain",
##       "protein_domain on (protein.protein_id=protein_domain.protein_id)"),
##     c("uniprot", "protein_domain",
##       "protein_domain on (uniprot.protein_id=protein_domain.protein_id)")
## )

############################################################
## joinQueryOnTables
##
## Helper function to generate the join query based on the provided tables.
## joinQueryOnTables <- function(x, tab, join = "join"){
##     ## Evaluate whether we have all required tables to join;
##     ## this will also ensure that the order is by degree.
##     tab <- addRequiredTables(x, tab)
##     Query <- tab[1]
##     previous.table <- tab[1]
##     for (i in 1:length(tab)) {
##         if (i > 1) {
##             Query <- paste(Query, join, .JOINS[.JOINS[, 1] %in% previous.table &
##                                                .JOINS[, 2] == tab[i], 3])
##             ## Add the table to the previous tables.
##             previous.table <- c(previous.table, tab[i])
##         }
##     }
##     return(Query)
## }

############################################################
## ** NEW JOIN ENGINE **
##
## 1: table 1
## 2: table 2
## 3: on
## 4: suggested join
.JOINS2 <- rbind(
    c("gene", "tx", "on (gene.gene_id=tx.gene_id)", "join"),
    c("gene", "chromosome", "on (gene.seq_name=chromosome.seq_name)", "join"),
    c("tx", "tx2exon", "on (tx.tx_id=tx2exon.tx_id)", "join"),
    c("tx2exon", "exon", "on (tx2exon.exon_id=exon.exon_id)", "join"),
    c("tx", "protein", "on (tx.tx_id=protein.tx_id)", "left outer join"),
    c("protein", "protein_domain",
      "on (protein.protein_id=protein_domain.protein_id)", "left outer join"),
    c("protein", "uniprot", "on (protein.protein_id=uniprot.protein_id)",
      "left outer join"),
    c("uniprot", "protein_domain",
      "on (uniprot.protein_id=protein_domain.protein_id)", "left outer join")
)
## Takes the names of two tables, determines how to join them and returns the
## join query row, if found.
joinTwoTables <- function(a, b) {
    gotIt <- which((.JOINS2[, 1] %in% a & .JOINS2[, 2] %in% b) |
                   (.JOINS2[, 2] %in% a & .JOINS2[, 1] %in% b))
    if (length(gotIt) == 0) {
        stop("Table(s) ", paste(a, collapse = ", "), " can not be joined with ",
             paste(b, collapse = ", "), "!")
    } else {
        return(.JOINS2[gotIt[1], ])
    }
}
## x: EnsDb.
## tab: tables to join.
## join: which type of join should be used?
## startWith: optional table name from which the join should start. That's
## specifically important for a left outer join call.
joinQueryOnTables2 <- function(x, tab, join = "suggested", startWith = NULL) {
    ## join can be join, left join, left outer join or suggested in which case
    ## the join defined in the .JOINS2 table will be used.
    join <- match.arg(join, c("join", "left join", "left outer join",
                              "suggested"))
    ## Order the tables.
    ## Start with startWith, or with the first one.
    if (missing(tab))
        stop("Argument 'tab' missing! Need some tables to make a join!")
    if (!is.null(startWith)) {
        if (!any(tab == startWith))
            stop("If provided, 'startWith' has to be the name of one of the",
                 " tables that should be joined!")
    }
    ## Add eventually needed tables to link the ones provided. The tables will
    ## be ordered by degree.
    tab <- addRequiredTables(x, tab)
    if (!is.null(startWith)) {
        alreadyUsed <- startWith
        tab <- tab[tab != startWith]
    } else {
        alreadyUsed <- tab[1]
        tab <- tab[-1]
    }
    Query <- alreadyUsed
    ## Iteratively build the query.
    while (length(tab) > 0) {
        res <- joinTwoTables(a = alreadyUsed, b = tab)
        newTab <- res[1:2][!(res[1:2] %in% alreadyUsed)]
        ## Could also use the suggested join which is in element 4.
        Query <- paste(Query, ifelse(join == "suggested", res[4], join),
                       newTab, res[3])
        alreadyUsed <- c(alreadyUsed, newTab)
        tab <- tab[tab != newTab]
    }
    return(Query)
}
## this function has to first get all tables that contain the columns,
## and then select, for columns present in more than one
## x... EnsDb
## columns... the columns
## NOTE: if "startWith" is not NULL, we're adding it to the tables!!!!
joinQueryOnColumns2 <- function(x, columns, join = "suggested",
                                startWith = NULL) {
    columns.bytable <- prefixColumns(x, columns)
    ## based on that we can build the query based on the tables we've got.
    ## Note that the function internally
    ## adds tables that might be needed for the join.
    Query <- joinQueryOnTables2(x, tab = c(names(columns.bytable), startWith),
                                join = join,
                                startWith = startWith)
    return(Query)
}


###
## Add additional tables in case the submitted tables are not directly connected
## and can thus not be joined. That's however not so complicated, since the database
## layout is pretty simple.
## The tables are:
##
##      uniprot -(protein_id=protein_id)-+-(protein_id=protein_id)- protein_domain
##                                       |
##                                    protein -(tx_id=tx_id)-+
##                                                           |
##  exon -(exon_id=t2e_exon_id)- tx2exon -(t2e_tx_id=tx_id)- tx -(gene_id=gene_id)- gene
##                                                                                   |
##                                                   chromosome -(seq_name=seq_name)-+
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
    if (hasProteinData(x)) {
        ## Resolve the proteins: need tx to map between proteome and genome
        if (any(tab %in% c("uniprot", "protein_domain", "protein")) &
            any(tab %in% c("exon", "tx2exon", "gene", "chromosome")))
            tab <- unique(c(tab, "tx"))
        ## Need protein.
        if (any(tab %in% c("uniprot", "protein_domain")) &
            any(tab %in% c("exon", "tx2exon", "tx", "gene", "chromosome")))
            tab <- unique(c(tab, "protein"))
    }
    return(tablesByDegree(x, tab))
}


############################################################
## .buildQuery
##
## The "backbone" function that builds the SQL query based on the specified
## columns, the provided filters etc.
## x an EnsDb object
## startWith: optional table from which the join should start.
.buildQuery <- function(x, columns, filter = list(), order.by = "",
                        order.type = "asc", group.by, skip.order.check=FALSE,
                        return.all.columns = TRUE, join = "suggested",
                        startWith = NULL) {
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
    joinquery <- joinQueryOnColumns2(x, columns=columns, join = join,
                                     startWith = startWith)
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


## just to add another layer; basically just calls buildQuery and executes the
## query
## join: what type of join should be performed.
## startWith: the name of the table from which the query should be started.
.getWhat <- function(x, columns, filter = list(), order.by = "",
                     order.type = "asc", group.by = NULL,
                     skip.order.check = FALSE, join = "suggested",
                     startWith = NULL) {
    ## That's nasty stuff; for now we support the column tx_name, which we however
    ## don't have in the database. Thus, we are querying everything except that
    ## column and filling it later with the values from tx_id.
    fetchColumns <- columns
    if(any(columns == "tx_name"))
        fetchColumns <- unique(c("tx_id",
                                 fetchColumns[fetchColumns != "tx_name"]))
    if (class(filter) != "list")
        stop("parameter filter has to be a list of BasicFilter classes!")
    ## If any filter is a SymbolFilter, add "symbol" to the return columns.
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
                         skip.order.check = skip.order.check, join = join,
                         startWith = startWith)
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
                         skip.order.check = skip.order.check, join = join,
                         startWith = startWith)
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
.ENSDB_PROTEIN_TABLES <- list(protein = c("tx_id", "protein_id",
                                          "protein_sequence"),
                              uniprot = c("protein_id", "uniprot_id",
                                          "uniprot_db", "uniprot_mapping_type"),
                              protein_domain = c("protein_id",
                                                 "protein_domain_id",
                                                 "protein_domain_source",
                                                 "interpro_accession",
                                                 "prot_dom_start",
                                                 "prot_dom_end"))
dbHasRequiredTables <- function(con, returnError = TRUE,
                                tables = .ENSDB_TABLES) {
    tabs <- dbListTables(con)
    if (length(tabs) == 0) {
        if (returnError)
            return("Database does not have any tables!")
        return(FALSE)
    }
    not_there <- names(tables)[!(names(tables) %in% tabs)]
    if (length(not_there) > 0) {
        if (returnError)
            return(paste0("Required tables ", paste(not_there, collapse = ", "),
                          " are not present in the database!"))
        return(FALSE)
    }
    return(TRUE)
}
dbHasValidTables <- function(con, returnError = TRUE,
                             tables = .ENSDB_TABLES) {
    for (tab in names(tables)) {
        cols <- tables[[tab]]
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
    .createEnsDbIndices(mysql, indexLength = "(20)",
                        proteins = hasProteinData(x))
    if (verbose)
        message("OK")
    return(TRUE)
}
## Small helper function to cfeate all the indices.
.createEnsDbIndices <- function(con, indexLength = "", proteins = FALSE) {
    indexCols <- c(chromosome = "seq_name", gene = "gene_id", gene = "gene_name",
                   gene = "seq_name", tx = "tx_id", tx = "gene_id",
                   exon = "exon_id", tx2exon = "tx_id", tx2exon = "exon_id")
    if (proteins) {
        indexCols <- c(indexCols,
                       protein = "tx_id",
                       protein = "protein_id",
                       uniprot = "protein_id",
                       uniprot = "uniprot_id",
                       protein_domain = "protein_domain_id",
                       protein_domain = "protein_id")
    }
    for (i in 1:length(indexCols)) {
        tabname <- names(indexCols)[i]
        colname <- indexCols[i]
        dbGetQuery(con, paste0("create index ", tabname, "_", colname, "_idx ",
                                "on ", tabname, " (",colname, indexLength,")"))
    }
    ## Add the one on the numeric index:
    dbGetQuery(con, "create index tx2exon_exon_idx_idx on tx2exon (exon_idx);")
}

############################################################
## listEnsDbs
## list databases
##' @title List EnsDb databases in a MySQL server
##' @description The \code{listEnsDbs} function lists EnsDb databases in a
##' MySQL server.
##'
##' @details The use of this function requires that the \code{RMySQL} package
##' is installed and that the user has either access to a MySQL server with
##' already installed EnsDb databases, or write access to a MySQL server in
##' which case EnsDb databases could be added with the \code{\link{useMySQL}}
##' method. EnsDb databases follow the same naming conventions than the EnsDb
##' packages, with the exception that the name is all lower case and that
##' \code{"."} is replaced by \code{"_"}.
##' @param dbcon A \code{DBIConnection} object providing access to a MySQL
##' database. Either \code{dbcon} or all of the other arguments have to be
##' specified.
##' @param host Character specifying the host on which the MySQL server is
##' running.
##' @param port The port of the MySQL server (usually \code{3306}).
##' @param user The username for the MySQL server.
##' @param pass The password for the MySQL server.
##' @return A \code{data.frame} listing the database names, organism name
##' and Ensembl version of the EnsDb databases found on the server.
##' @author Johannes Rainer
##' @seealso \code{\link{useMySQL}}
##' @examples
##' \dontrun{
##' library(RMySQL)
##' dbcon <- dbConnect(MySQL(), host = "localhost", user = my_user, pass = my_pass)
##' listEnsDbs(dbcon)
##' }
listEnsDbs <- function(dbcon, host, port, user, pass) {
    if(requireNamespace("RMySQL", quietly = TRUE)) {
        if (missing(dbcon)) {
            if (missing(host) | missing(user) | missing(port) | missing(host))
                stop("Arguments 'host', 'port', 'user' and 'pass' are required",
                     " if 'dbcon' is not specified.")
            dbcon <- dbConnect(RMySQL::MySQL(), host = host, port = port, user = user,
                               pass = pass)
        }
        dbs <- dbGetQuery(dbcon, "show databases;")
        edbs <- dbs[grep(dbs$Database, pattern = "^ensdb_"), "Database"]
        edbTable <- data.frame(matrix(ncol = 3, nrow = length(edbs)))
        colnames (edbTable) <- c("dbname", "organism", "ensembl_version")
        for (i in seq_along(edbs)) {
            edbTable[i, "dbname"] <- edbs[i]
            tmp <- unlist(strsplit(edbs[i], split = "_"))
            edbTable[i, "organism"] <- tmp[2]
            edbTable[i, "ensembl_version"] <- as.numeric(gsub(pattern = "v",
                                                              replacement = "",
                                                              tmp[3]))
        }
        return(edbTable)
    } else {
        stop("Required package 'RMySQL' is not installed.")
    }
}
