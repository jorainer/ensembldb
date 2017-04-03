############################################################
## Compare MySQL vs SQLite backends:
## Amazing how inefficient the MySQL backend seems to be! Most
## likely it's due to RMySQL, not MySQL.
dontrun_test_MySQL_vs_SQLite <- function() {
    ## Compare the performance of the MySQL backend against
    ## the SQLite backend.
    edb_mysql <- useMySQL(edb, user = "anonuser", pass = "")

    library(microbenchmark)
    ## genes
    microbenchmark(genes(edb), genes(edb_mysql), times = 5)
    microbenchmark(genes(edb, filter = GeneBiotypeFilter("lincRNA")),
                   genes(edb_mysql, filter = GeneBiotypeFilter("lincRNA")),
                   times = 5)
    microbenchmark(genes(edb, filter = SeqNameFilter(20:23)),
                   genes(edb_mysql, filter = SeqNameFilter(20:23)),
                   times = 5)
    microbenchmark(genes(edb, columns = "tx_id"),
                   genes(edb_mysql, columns = "tx_id"),
                   times = 5)
    microbenchmark(genes(edb, filter = GenenameFilter("BCL2L11")),
                   genes(edb_mysql, filter = GenenameFilter("BCL2L11")),
                   times = 5)
    ## transcripts
    microbenchmark(transcripts(edb),
                   transcripts(edb_mysql),
                   times = 5)
    microbenchmark(transcripts(edb, filter = GenenameFilter("BCL2L11")),
                   transcripts(edb_mysql, filter = GenenameFilter("BCL2L11")),
                   times = 5)
    ## exons
    microbenchmark(exons(edb),
                   exons(edb_mysql),
                   times = 5)
    microbenchmark(exons(edb, filter = GenenameFilter("BCL2L11")),
                   exons(edb_mysql, filter = GenenameFilter("BCL2L11")),
                   times = 5)
    ## exonsBy
    microbenchmark(exonsBy(edb),
                   exonsBy(edb_mysql),
                   times = 5)
    microbenchmark(exonsBy(edb, filter = SeqNameFilter("Y")),
                   exonsBy(edb_mysql, filter = SeqNameFilter("Y")),
                   times = 5)
    ## cdsBy
    microbenchmark(cdsBy(edb), cdsBy(edb_mysql), times = 5)
    microbenchmark(cdsBy(edb, by = "gene"), cdsBy(edb_mysql, by = "gene"),
                   times = 5)
    microbenchmark(cdsBy(edb, filter = SeqStrandFilter("-")),
                   cdsBy(edb_mysql, filter = SeqStrandFilter("-")),
                   times = 5)

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

## Compare the performance of inner join with left outer join.
dontrun_test_outer_join_performance <- function() {
    Q_1 <- ensembldb:::joinQueryOnTables2(edb, tab = c("gene", "exon"))
    Q_2 <- ensembldb:::joinQueryOnTables2(edb, tab = c("gene", "exon"),
                                          startWith = "exon")
    Q_3 <- ensembldb:::joinQueryOnTables2(edb, tab = c("gene", "exon"),
                                          startWith = "exon",
                                          join = "left outer join")
    library(microbenchmark)
    library(RSQLite)
    microbenchmark(dbGetQuery(dbconn(edb), paste0("select * from ", Q_1)),
                   dbGetQuery(dbconn(edb), paste0("select * from ", Q_2)),
                   dbGetQuery(dbconn(edb), paste0("select * from ", Q_3)),
                   times = 10)
    ## Result: Q_1 is a second faster (13 instead of 14).
    ## Check performance joining tx and genes.
    Q_1 <- ensembldb:::joinQueryOnTables2(edb, tab = c("tx", "gene"))
    Q_2 <- ensembldb:::joinQueryOnTables2(edb, tab = c("tx", "gene"),
                                          startWith = "tx")
    Q_3 <- ensembldb:::joinQueryOnTables2(edb, tab = c("tx", "gene"),
                                          startWith = "tx",
                                          join = "left outer join")
    microbenchmark(dbGetQuery(dbconn(edb), paste0("select * from ", Q_1)),
                   dbGetQuery(dbconn(edb), paste0("select * from ", Q_2)),
                   dbGetQuery(dbconn(edb), paste0("select * from ", Q_3)),
                   times = 10)
    ## No difference.
}
