
test_that("validity functions work", {
    OK <- ensembldb:::dbHasRequiredTables(dbconn(edb))
    expect_true(OK)
    ## Check the tables
    OK <- ensembldb:::dbHasValidTables(dbconn(edb))
    expect_true(OK)
})

test_that("validateEnsDb works", {
    expect_true(ensembldb:::validateEnsDb(edb))
})

test_that("compareProteins works", {
    if (hasProteinData(edb)) {
        res <- ensembldb:::compareProteins(edb, edb)
        expect_equal(res, "OK")
    }
})

