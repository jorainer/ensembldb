############################################################
## These are not test cases to be executed, but performance
## comparisons.
library(EnsDb.Hsapiens.v75)
edb <- EnsDb.Hsapiens.v75

dontrun_test_MySQL_vs_SQLite <- function() {
    ## Compare the performance of the MySQL backend against
    ## the SQLite backend.
    edb_mysql <- ensembldb:::useMySQL(edb, user = "anonuser", pass = "")

    library(microbenchmark)
    ## genes
    microbenchmark(genes(edb), genes(edb_mysql), times = 5)
    microbenchmark(genes(edb, filter = GenebiotypeFilter("lincRNA")),
                   genes(edb_mysql, filter = GenebiotypeFilter("lincRNA")),
                   times = 5)
    microbenchmark(genes(edb, filter = SeqnameFilter(20:23)),
                   genes(edb_mysql, filter = SeqnameFilter(20:23)),
                   times = 5)
    microbenchmark(genes(edb, columns = "tx_id"),
                   genes(edb_mysql, columns = "tx_id"),
                   times = 5)
    ## cdsBy
    microbenchmark(cdsBy(edb), cdsBy(edb_mysql), times = 5)
    microbenchmark(cdsBy(edb, by = "gene"), cdsBy(edb_mysql, by = "gene"),
                   times = 5)
    microbenchmark(cdsBy(edb, filter = SeqstrandFilter("-")),
                   cdsBy(edb_mysql, filter = SeqstrandFilter("-")),
                   times = 5)
}

