test_that("submitting seqlevel mapping custom mapping", {
    ## Setting the mapping data.frame as "genomeStyle" property should work
    mymap <- data.frame(myway = c("a", "b", "c", "e", "D"),
                        Ensembl = c(1, 5, 8, 9, "MT"))
    x <- edb
    seqlevelsStyle(x) <- mymap
    res <- seqlevels(x)
    orig <- seqlevels(edb)
    idx <- match(mymap$Ensembl, orig)
    expect_equal(res[idx], mymap$myway)
    seqlevelsStyle(x) <- "Ensembl"
    expect_equal(orig, seqlevels(x))
})

test_that("seqlevelsStyle,EnsDb works with custom mapping", {
    x <- edb
    mymap <- data.frame(myway = c("a", "b", "c", "e", "D"),
                        Ensembl = c(1, 5, 8, 9, "MT"))
    expect_error(seqlevelsStyle(x) <- "myway", "not known")
    expect_error(seqlevelsStyle(x) <- data.frame(a = 1:10, b = 1:10),
                 "column needs")
    expect_error(seqlevelsStyle(x) <- mymap[, "Ensembl", drop = FALSE],
                 "at least two")

    seqlevelsStyle(x) <- mymap
    expect_equal(x@.properties$genomeStyle, mymap)
    expect_equal(x@.properties$seqlevelsStyle, "myway")
})
