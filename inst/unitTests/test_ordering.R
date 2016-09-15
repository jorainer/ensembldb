############################################################
## Some tests on the ordering/sorting of the results.
library(EnsDb.Hsapiens.v75)
edb <- EnsDb.Hsapiens.v75

## Compare the results for genes call with and without ordering in R
test_ordering_genes <- function() {
    orig <- ensembldb:::orderResultsInR(edb)
    ensembldb:::orderResultsInR(edb) <- FALSE
    res_sql <- genes(edb, return.type = "data.frame")
    ensembldb:::orderResultsInR(edb) <- TRUE
    res_r <- genes(edb, return.type = "data.frame")
    rownames(res_sql) <- NULL
    rownames(res_r) <- NULL
    checkIdentical(res_sql, res_r)
    ## Join tx table
    ensembldb:::orderResultsInR(edb) <- FALSE
    res_sql <- genes(edb, columns = c("gene_id", "tx_id"),
                     return.type = "data.frame")
    ensembldb:::orderResultsInR(edb) <- TRUE
    res_r <- genes(edb, columns = c("gene_id", "tx_id"),
                   return.type = "data.frame")
    rownames(res_sql) <- NULL
    rownames(res_r) <- NULL
    checkIdentical(res_sql, res_r)
    ## Join tx table and use an SeqnameFilter
    ensembldb:::orderResultsInR(edb) <- FALSE
    res_sql <- genes(edb, columns = c("gene_id", "tx_id"),
                     filter = SeqnameFilter("Y"))
    ensembldb:::orderResultsInR(edb) <- TRUE
    res_r <- genes(edb, columns = c("gene_id", "tx_id"),
                   filter = SeqnameFilter("Y"))
    checkIdentical(res_sql, res_r)

    ensembldb:::orderResultsInR(edb) <- orig
}

dontrun_benchmark_ordering_genes <- function() {
    .withR <- function(x, ...) {
        ensembldb:::orderResultsInR(x) <- TRUE
        genes(x, ...)
    }
    .withSQL <- function(x, ...) {
        ensembldb:::orderResultsInR(x) <- FALSE
        genes(x, ...)
    }
    library(microbenchmark)
    microbenchmark(.withR(edb), .withSQL(edb), times = 10)  ## same
    microbenchmark(.withR(edb, columns = c("gene_id", "tx_id")),
                   .withSQL(edb, columns = c("gene_id", "tx_id")),
                   times = 10)  ## R slightly faster.
    microbenchmark(.withR(edb, columns = c("gene_id", "tx_id"),
                          SeqnameFilter("Y")),
                   .withSQL(edb, columns = c("gene_id", "tx_id"),
                            SeqnameFilter("Y")),
                   times = 10)  ## same.
}

## We aim to fix issue #11 by performing the ordering in R instead
## of SQL. Thus, we don't want to run this as a "regular" test
## case.
dontrun_test_ordering_cdsBy <- function() {
    doBench <- FALSE
    if (doBench)
        library(microbenchmark)
    .withR <- function(x, ...) {
        ensembldb:::orderResultsInR(x) <- TRUE
        cdsBy(x, ...)
    }
    .withSQL <- function(x, ...) {
        ensembldb:::orderResultsInR(x) <- FALSE
        cdsBy(x, ...)
    }
    res_sql <- .withSQL(edb)
    res_r <- .withR(edb)
    checkEquals(res_sql, res_r)
    if (dobench)
        microbenchmark(.withSQL(edb), .withR(edb), times = 10)
    res_sql <- .withSQL(edb, filter = SeqnameFilter("Y"))
    res_r <- .withR(edb, filter = SeqnameFilter("Y"))
    checkEquals(res_sql, res_r)
    if (dobench)
        microbenchmark(.withSQL(edb, filter = SeqnameFilter("Y")),
                       .withR(edb, filter = SeqnameFilter("Y")),
                       times = 10)
}

dontrun_test_ordering_exonsBy <- function() {
    doBench <- FALSE
    if (doBench)
        library(microbenchmark)
    .withR <- function(x, ...) {
        ensembldb:::orderResultsInR(x) <- TRUE
        exonsBy(x, ...)
    }
    .withSQL <- function(x, ...) {
        ensembldb:::orderResultsInR(x) <- FALSE
        exonsBy(x, ...)
    }
    res_sql <- .withSQL(edb)
    res_r <- .withR(edb)
    checkEquals(res_sql, res_r)
    if (doBench)
        microbenchmark(.withSQL(edb), .withR(edb),
                       times = 3)  ## about the same; R slightly faster.
    ## with using a SeqnameFilter in addition.
    res_sql <- .withSQL(edb, filter = SeqnameFilter("Y"))
    res_r <- .withR(edb, filter = SeqnameFilter("Y")) ## query takes longer.
    checkEquals(res_sql, res_r)
    if (doBench)
        microbenchmark(.withSQL(edb, filter = SeqnameFilter("Y")),
                       .withR(edb, filter = SeqnameFilter("Y")),
                       times = 3)  ## SQL twice as fast.
    ## Now getting stuff by gene
    res_sql <- .withSQL(edb, by = "gene")
    res_r <- .withR(edb, by = "gene")
    ## checkEquals(res_sql, res_r) ## Differences due to ties
    if (doBench)
        microbenchmark(.withSQL(edb, by = "gene"),
                       .withR(edb, by = "gene"),
                       times = 3)  ## SQL faster; ???
    ## Along with a SeqnameFilter
    res_sql <- .withSQL(edb, by = "gene", filter = SeqnameFilter("Y"))
    res_r <- .withR(edb, by = "gene", filter = SeqnameFilter("Y"))
    ## Why does the query take longer for R???
    ## checkEquals(res_sql, res_r) ## Differences due to ties
    if (doBench)
        microbenchmark(.withSQL(edb, by = "gene", filter = SeqnameFilter("Y")),
                       .withR(edb, by = "gene", filter = SeqnameFilter("Y")),
                       times = 3)  ## SQL faster.
    ## Along with a GenebiotypeFilter
    if (doBench)
        microbenchmark(.withSQL(edb, by = "gene", filter = GenebiotypeFilter("protein_coding"))
                     , .withR(edb, by = "gene", filter = GenebiotypeFilter("protein_coding"))
                     , times = 3)
}

dontrun_query_tune <- function() {
    ## Query tuning:
    library(RSQLite)
    con <- dbconn(edb)

    Q <- "select distinct tx2exon.exon_id,exon.exon_seq_start,exon.exon_seq_end,gene.seq_name,tx2exon.tx_id,gene.seq_strand,tx2exon.exon_idx from gene join tx on (gene.gene_id=tx.gene_id) join tx2exon on (tx.tx_id=tx2exon.tx_id) join exon on (tx2exon.exon_id=exon.exon_id) where gene.seq_name = 'Y'"
    system.time(dbGetQuery(con, Q))

    Q2 <- "select distinct tx2exon.exon_id,exon.exon_seq_start,exon.exon_seq_end,gene.seq_name,tx2exon.tx_id,gene.seq_strand,tx2exon.exon_idx from exon join tx2exon on (tx2exon.exon_id = exon.exon_id) join tx on (tx2exon.tx_id = tx.tx_id) join gene on (gene.gene_id=tx.gene_id) where gene.seq_name = 'Y'"
    system.time(dbGetQuery(con, Q2))

    Q3 <- "select distinct tx2exon.exon_id,exon.exon_seq_start,exon.exon_seq_end,gene.seq_name,tx2exon.tx_id,gene.seq_strand,tx2exon.exon_idx from tx2exon join exon on (tx2exon.exon_id = exon.exon_id) join tx on (tx2exon.tx_id = tx.tx_id) join gene on (gene.gene_id=tx.gene_id) where gene.seq_name = 'Y'"
    system.time(dbGetQuery(con, Q3))

    Q4 <- "select distinct tx2exon.exon_id,exon.exon_seq_start,exon.exon_seq_end,gene.seq_name,tx2exon.tx_id,gene.seq_strand,tx2exon.exon_idx from tx2exon join exon on (tx2exon.exon_id = exon.exon_id) join tx on (tx2exon.tx_id = tx.tx_id) join gene on (gene.gene_id=tx.gene_id) where gene.seq_name = 'Y' order by tx.tx_id"
    system.time(dbGetQuery(con, Q4))

    Q5 <- "select distinct tx2exon.exon_id,exon.exon_seq_start,exon.exon_seq_end,gene.seq_name,tx2exon.tx_id,gene.seq_strand,tx2exon.exon_idx from tx2exon inner join exon on (tx2exon.exon_id = exon.exon_id) inner join tx on (tx2exon.tx_id = tx.tx_id) inner join gene on (gene.gene_id=tx.gene_id) where gene.seq_name = 'Y' order by tx.tx_id"
    system.time(dbGetQuery(con, Q5))

    Q6 <- "select distinct tx2exon.exon_id,exon.exon_seq_start,exon.exon_seq_end,gene.seq_name,tx2exon.tx_id,gene.seq_strand,tx2exon.exon_idx from gene inner join tx on (gene.gene_id=tx.gene_id) inner join tx2exon on (tx.tx_id=tx2exon.tx_id) inner join exon on (tx2exon.exon_id=exon.exon_id) where gene.seq_name = 'Y' order by tx.tx_id asc"
    system.time(dbGetQuery(con, Q6))

}


## Compare the performance of doing the sorting within R or
## directly in the SQL query.
dontrun_test_ordering_performance <- function() {

    library(RUnit)
    library(RSQLite)
    ## gene table: order by in SQL query vs R:
    db_con <- dbconn(edb)

    .callWithOrder <- function(con, query, orderBy = "",
                               orderSQL = TRUE) {
        if (all(orderBy == ""))
            orderBy <- NULL
        if (orderSQL & !is.null(orderBy)) {
            orderBy <- paste(orderBy, collapse = ", ")
            query <- paste0(query, " order by ", orderBy)
        }
        res <- dbGetQuery(con, query)
        if (!orderSQL & !all(is.null(orderBy))) {
            if (!all(orderBy %in% colnames(res)))
                stop("orderBy not in columns!")
            ## Do the ordering in R
            res <- res[do.call(order,
                               c(list(method = "radix"),
                                 as.list(res[, orderBy, drop = FALSE]))), ]
        }
        rownames(res) <- NULL
        return(res)
    }

    #######################
    ## gene table
    ## Simple condition
    the_q <- "select * from gene"
    system.time(res1 <- .callWithOrder(db_con, query = the_q))
    system.time(res2 <- .callWithOrder(db_con, query = the_q,
                                       orderSQL = FALSE))
    checkIdentical(res1, res2)
    ## order by gene_id
    orderBy <- "gene_id"
    system.time(res1 <- .callWithOrder(db_con, query = the_q, orderBy = orderBy))
    system.time(res2 <- .callWithOrder(db_con, query = the_q,
                                       orderBy = orderBy, orderSQL = FALSE))
    ## SQL: 0.16, R: 0.164.
    checkIdentical(res1, res2)
    ## order by gene_name
    orderBy <- "gene_name"
    system.time(res1 <- .callWithOrder(db_con, query = the_q, orderBy = orderBy))
    system.time(res2 <- .callWithOrder(db_con, query = the_q,
                                       orderBy = orderBy, orderSQL = FALSE))
    checkIdentical(res1, res2)
    ## SQL: 0.245, R: 0.185
    ## sort by gene_name and gene_seq_start
    orderBy <- c("gene_name", "gene_seq_start")
    system.time(res1 <- .callWithOrder(db_con, query = the_q, orderBy = orderBy))
    system.time(res2 <- .callWithOrder(db_con, query = the_q,
                                       orderBy = orderBy, orderSQL = FALSE))
    ## SQL: 0.26, R: 0.188
    checkEquals(res1, res2)
    ## with subsetting:
    the_q <- "select * from gene where seq_name in ('5', 'Y')"
    orderBy <- c("gene_name", "gene_seq_start")
    system.time(res1 <- .callWithOrder(db_con, query = the_q, orderBy = orderBy))
    system.time(res2 <- .callWithOrder(db_con, query = the_q,
                                       orderBy = orderBy, orderSQL = FALSE))
    ## SQL: 0.031, R: 0.024
    checkEquals(res1, res2)

    ########################
    ## joining tables.
    the_q <- paste0("select * from gene join tx on (gene.gene_id = tx.gene_id)",
                    " join tx2exon on (tx.tx_id = tx2exon.tx_id)")
    orderBy <- c("tx_id", "exon_id")
    system.time(res1 <- .callWithOrder(db_con, query = the_q, orderBy = orderBy))
    system.time(res2 <- .callWithOrder(db_con, query = the_q,
                                       orderBy = orderBy, orderSQL = FALSE))
    ## SQL: 9.6, R: 9.032
    checkEquals(res1, res2)
    ## subsetting.
    the_q <- paste0("select * from gene join tx on (gene.gene_id = tx.gene_id)",
                    " join tx2exon on (tx.tx_id = tx2exon.tx_id) where",
                    " seq_name = 'Y'")
    orderBy <- c("tx_id", "exon_id")
    system.time(res1 <- .callWithOrder(db_con, query = the_q, orderBy = orderBy))
    system.time(res2 <- .callWithOrder(db_con, query = the_q,
                                       orderBy = orderBy, orderSQL = FALSE))
    ## SQL: 0.9, R: 1.6
    checkEquals(res1, res2)
}

## implement:
## .checkOrderBy: checks order.by argument removing columns that are
## not present in the database
## orderBy columns are added to the columns.
## .orderDataFrameBy: orders the dataframe by the specified columns.
