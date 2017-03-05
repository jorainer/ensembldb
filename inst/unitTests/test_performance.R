############################################################
## These are not test cases to be executed, but performance
## comparisons.
library(EnsDb.Hsapiens.v75)
edb <- EnsDb.Hsapiens.v75


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

