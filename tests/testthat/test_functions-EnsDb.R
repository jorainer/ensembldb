test_that(".addFilter works", {
    gf <- GenenameFilter("BCL2")
    edb_2 <- ensembldb:::.addFilter(edb, filter = gf)
    expect_equal(gf, ensembldb:::getProperty(edb_2, "FILTER"))
    edb_2 <- ensembldb:::.addFilter(edb_2, filter = gf)
    expect_equal(AnnotationFilterList(gf, gf),
                 ensembldb:::getProperty(edb_2, "FILTER"))
})

test_that(".dropFilter works", {
})

test_that(".activeFilter works", {
})
