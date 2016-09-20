############################################################
## Can not perform these tests right away, as they require a
## working MySQL connection.
library(EnsDb.Hsapiens.v75)
edb <- EnsDb.Hsapiens.v75

dontrun_test_useMySQL <- function() {
    edb_mysql <- useMySQL(edb, user = "anonuser", host = "localhost", pass = "")
}

dontrun_test_connect_EnsDb <- function() {
    library(RMySQL)
    con <- dbConnect(MySQL(), user = "anonuser", host = "localhost", pass = "")

    ensembldb:::listEnsDbs(dbcon = con)
    ## just with user.
    ensembldb:::listEnsDbs(user = "anonuser", host = "localhost", pass = "",
                           port = 3306)

    ## Connecting directly to a EnsDb MySQL database.
    con <- dbConnect(MySQL(), user = "anonuser", host = "localhost", pass = "",
                     dbname = "ensdb_hsapiens_v75")
    edb_mysql <- EnsDb(con)
}
