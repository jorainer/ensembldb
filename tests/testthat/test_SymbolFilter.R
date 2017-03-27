
test_that("SymbolFilter works for gene", {
    sf <- SymbolFilter("SKA2")
    gnf <- GenenameFilter("SKA2")
    returnFilterColumns(edb) <- FALSE
    gns_sf <- genes(edb, filter = sf)
    gns_gnf <- genes(edb, filter = gnf)
    expect_equal(gns_sf, gns_gnf)
    returnFilterColumns(edb) <- TRUE
    gns_sf <- genes(edb, filter=sf)
    expect_equal(gns_sf$gene_name, gns_sf$symbol)
    ## Hm, what happens if we use both?
    gns <- genes(edb, filter=list(sf, gnf))
    ## All fine.
})

test_that("SymbolFilter works for tx", {
    sf <- SymbolFilter("SKA2")
    gnf <- GenenameFilter("SKA2")
    returnFilterColumns(edb) <- FALSE
    tx_sf <- transcripts(edb, filter=sf)
    tx_gnf <- transcripts(edb, filter=gnf)
    expect_equal(tx_sf, tx_gnf)
    returnFilterColumns(edb) <- TRUE
    tx_sf <- transcripts(edb, filter=sf, columns=c("gene_name"))
    expect_equal(tx_sf$gene_name, tx_sf$symbol)
})

test_that("SymbolFilter works for exons", {
    sf <- SymbolFilter("SKA2")
    gnf <- GenenameFilter("SKA2")
    returnFilterColumns(edb) <- FALSE
    ex_sf <- exons(edb, filter=sf)
    ex_gnf <- exons(edb, filter=gnf)
    expect_equal(ex_sf, ex_gnf)
    returnFilterColumns(edb) <- TRUE
    ex_sf <- exons(edb, filter=sf, columns=c("gene_name"))
    expect_equal(ex_sf$gene_name, ex_sf$symbol)
})

test_that("SymbolFilter works", {
    sf <- SymbolFilter("SKA2")
    res <- genes(edb, filter = sf, return.type = "data.frame")
    expect_equal(res$gene_id, "ENSG00000182628")
    ## We need now also a column "symbol"!
    expect_equal(res$symbol, res$gene_name)
    ## Asking explicitely for symbol
    res <- genes(edb, filter = sf, return.type = "data.frame",
                 columns = c("symbol", "gene_id"))
    expect_equal(colnames(res), c("symbol", "gene_id"))
    ## Some more stuff, also shuffling the order.
    res <- genes(edb, filter = sf, return.type = "data.frame",
                 columns = c("gene_name", "symbol", "gene_id"))
    expect_equal(colnames(res), c("gene_name", "symbol", "gene_id"))
    res <- genes(edb, filter = sf, return.type = "data.frame",
                 columns = c("gene_id", "gene_name", "symbol"))
    expect_equal(colnames(res), c("gene_id", "gene_name", "symbol"))
    ## And with GRanges as return type.
    res <- genes(edb, filter = sf, return.type = "GRanges",
                 columns = c("gene_id", "gene_name", "symbol"))
    expect_equal(colnames(mcols(res)), c("gene_id", "gene_name", "symbol"))

    ## Combine tx_name and symbol
    res <- genes(edb, filter = sf, columns = c("tx_name", "symbol"),
                 return.type = "data.frame")
    expect_equal(colnames(res), c("tx_name", "symbol", "gene_id"))
    expect_true(all(res$symbol == "SKA2"))

    ## Test for transcripts
    res <- transcripts(edb, filter=sf, return.type="data.frame")
    expect_true(all(res$symbol == "SKA2"))
    res <- transcripts(edb, filter = sf, return.type = "data.frame",
                       columns = c("symbol", "tx_id", "gene_name"))
    expect_true(all(res$symbol == "SKA2"))
    expect_equal(res$symbol, res$gene_name)
    expect_equal(colnames(res), c("symbol", "tx_id", "gene_name"))

    ## Test for exons
    res <- exons(edb, filter=sf, return.type="data.frame")
    expect_true(all(res$symbol == "SKA2"))
    res <- exons(edb, filter = c(sf, TxBiotypeFilter("nonsense_mediated_decay")),
                 return.type = "data.frame",
                 columns = c("symbol", "tx_id", "gene_name"))
    expect_true(all(res$symbol == "SKA2"))
    expect_equal(res$symbol, res$gene_name)
    expect_equal(colnames(res), c("symbol", "tx_id", "gene_name", "exon_id",
                                  "tx_biotype"))

    ## Test for exonsBy
    res <- exonsBy(edb, filter=sf)
    expect_true(all(unlist(res)$symbol == "SKA2"))
    res <- exonsBy(edb, filter = c(sf, TxBiotypeFilter("nonsense_mediated_decay")),
                   columns = c("symbol", "tx_id", "gene_name"))
    expect_true(all(unlist(res)$symbol == "SKA2"))

    expect_equal(unlist(res)$symbol, unlist(res)$gene_name)
})


