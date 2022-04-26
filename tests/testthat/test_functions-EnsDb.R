test_that(".addFilter .dropFilter and .activeFilter work", {
    gf <- GenenameFilter("BCL2")
    ## .addFilter and .activeFilter
    edb_2 <- ensembldb:::.addFilter(edb, filter = gf)
    expect_equal(AnnotationFilterList(gf),
                 ensembldb:::getProperty(edb_2, "FILTER"))
    expect_equal(AnnotationFilterList(gf),
                 ensembldb:::.activeFilter(edb_2))
    edb_2 <- ensembldb:::.addFilter(edb_2, filter = gf)
    expect_equal(AnnotationFilterList(AnnotationFilterList(gf),
                                      AnnotationFilterList(gf)),
                 ensembldb:::getProperty(edb_2, "FILTER"))
    edb_2 <- ensembldb:::.addFilter(edb, filter = ~ tx_id == 3 &
                                             tx_biotype == "protein_coding")
    flts <- ensembldb:::.activeFilter(edb_2)
    expect_equal(flts, ensembldb:::getProperty(edb_2, "FILTER"))
    expect_equal(flts, AnnotationFilter(~ tx_id == 3 &
                                            tx_biotype == "protein_coding"))

    ## Errors
    expect_error(ensembldb:::.addFilter(edb, "blabla"))
    expect_error(ensembldb:::filter(edb, "blabla"))
    expect_error(ensembldb:::.addFilter(edb))

    ## .dropFilter
    edb_2 <- ensembldb:::.dropFilter(edb_2)
    expect_equal(ensembldb:::.activeFilter(edb_2), NA)

    ## Same but with the methods.
    gf <- GenenameFilter("BCL2")
    ## .addFilter and .activeFilter
    edb_2 <- addFilter(edb, filter = gf)
    expect_equal(AnnotationFilterList(gf),
                 ensembldb:::getProperty(edb_2, "FILTER"))
    edb_2 <- filter(edb, filter = gf)
    expect_equal(AnnotationFilterList(gf),
                 ensembldb:::getProperty(edb_2, "FILTER"))
    expect_equal(AnnotationFilterList(gf),
                 activeFilter(edb_2))
    edb_2 <- addFilter(edb_2, filter = gf)
    expect_equal(AnnotationFilterList(AnnotationFilterList(gf),
                                      AnnotationFilterList(gf)),
                 ensembldb:::getProperty(edb_2, "FILTER"))
    edb_2 <- addFilter(edb, filter = ~ tx_id == 3 &
                                tx_biotype == "protein_coding")
    flts <- activeFilter(edb_2)
    expect_equal(flts, ensembldb:::getProperty(edb_2, "FILTER"))
    expect_equal(flts, AnnotationFilter(~ tx_id == 3 &
                                            tx_biotype == "protein_coding"))

    ## Errors
    expect_error(addFilter(edb, "blabla"))
    expect_error(addFilter(edb))

    ## .dropFilter
    edb_2 <- dropFilter(edb_2)
    expect_equal(activeFilter(edb_2), NA)
})

test_that(".has_tx_external_name works", {
    expect_false(.has_tx_external_name(edb))
})

test_that(".fix_is_circular works", {
    x <- data.frame(is_circular = 0L, seq_name = c("A", "B", "C", "D", "MT"))
    res <- .fix_is_circular(x)
    expect_equal(res$is_circular, c(0L, 0L, 0L, 0L, 1L))

    res <- .fix_is_circular(x, c("D", "A", "Z"))
    expect_equal(res$is_circular, c(1L, 0L, 0L, 1L, 0L))
})
