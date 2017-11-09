test_that(".proteinCoordsToTx works", {
    prts <- IRanges(start = c(1, 2), end = c(1, 2))
    res <- .proteinCoordsToTx(prts)
    expect_equal(start(res), c(1, 4))
    expect_equal(end(res), c(3, 6))

    res <- .proteinCoordsToTx(IRanges(start = 3, end = 5))
    expect_equal(start(res), 7)
    expect_equal(end(res), 15)
})
