library(EnsDb.Hsapiens.v75)
edb <- EnsDb.Hsapiens.v75

test_validity_functions <- function() {
    OK <- ensembldb:::dbHasRequiredTables(dbconn(edb))
    checkTrue(OK)
    ## Check the tables
    OK <- ensembldb:::dbHasValidTables(dbconn(edb))
    checkTrue(OK)
}

test_validateEnsDb <- function() {
    checkTrue(ensembldb:::validateEnsDb(edb))
}
