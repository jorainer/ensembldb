############################################################
## Testing the SymbolFilter.
library(EnsDb.Hsapiens.v75)
edb <- EnsDb.Hsapiens.v75

test_sf_on_genes <- function(){
    sf <- SymbolFilter("SKA2")
    gnf <- GenenameFilter("SKA2")

    returnFilterColumns(edb) <- FALSE
    gns_sf <- genes(edb, filter=sf)
    gns_gnf <- genes(edb, filter=gnf)
    checkEquals(gns_sf, gns_gnf)

    returnFilterColumns(edb) <- TRUE
    gns_sf <- genes(edb, filter=sf)
    checkEquals(gns_sf$gene_name, gns_sf$symbol)

    ## Hm, what happens if we use both?
    gns <- genes(edb, filter=list(sf, gnf))
    ## All fine.
}


test_sf_on_tx <- function(){
    sf <- SymbolFilter("SKA2")
    gnf <- GenenameFilter("SKA2")

    returnFilterColumns(edb) <- FALSE
    tx_sf <- transcripts(edb, filter=sf)
    tx_gnf <- transcripts(edb, filter=gnf)
    checkEquals(tx_sf, tx_gnf)

    returnFilterColumns(edb) <- TRUE
    tx_sf <- transcripts(edb, filter=sf, columns=c("gene_name"))
    checkEquals(tx_sf$gene_name, tx_sf$symbol)
}


test_sf_on_exons <- function(){
    sf <- SymbolFilter("SKA2")
    gnf <- GenenameFilter("SKA2")

    returnFilterColumns(edb) <- FALSE
    ex_sf <- exons(edb, filter=sf)
    ex_gnf <- exons(edb, filter=gnf)
    checkEquals(ex_sf, ex_gnf)

    returnFilterColumns(edb) <- TRUE
    ex_sf <- exons(edb, filter=sf, columns=c("gene_name"))
    checkEquals(ex_sf$gene_name, ex_sf$symbol)
}

test_SymbolFilter <- function() {
    edb <- EnsDb.Hsapiens.v75
    sf <- SymbolFilter("SKA2")

    Res <- genes(edb, filter = sf, return.type = "data.frame")
    checkEquals(Res$gene_id, "ENSG00000182628")
    ## We need now also a column "symbol"!
    checkEquals(Res$symbol, Res$gene_name)
    ## Asking explicitely for symbol
    Res <- genes(edb, filter = sf, return.type = "data.frame",
                 columns = c("symbol", "gene_id"))
    checkEquals(colnames(Res), c("symbol", "gene_id"))
    ## Some more stuff, also shuffling the order.
    Res <- genes(edb, filter = sf, return.type = "data.frame",
                 columns = c("gene_name", "symbol", "gene_id"))
    checkEquals(colnames(Res), c("gene_name", "symbol", "gene_id"))
    Res <- genes(edb, filter = sf, return.type = "data.frame",
                 columns = c("gene_id", "gene_name", "symbol"))
    checkEquals(colnames(Res), c("gene_id", "gene_name", "symbol"))
    ## And with GRanges as return type.
    Res <- genes(edb, filter = sf, return.type = "GRanges",
                 columns = c("gene_id", "gene_name", "symbol"))
    checkEquals(colnames(mcols(Res)), c("gene_id", "gene_name", "symbol"))

    ## Combine tx_name and symbol
    Res <- genes(edb, filter = sf, columns = c("tx_name", "symbol"),
                 return.type = "data.frame")
    checkEquals(colnames(Res), c("tx_name", "symbol", "gene_id"))
    checkTrue(all(Res$symbol == "SKA2"))

    ## Test for transcripts
    Res <- transcripts(edb, filter=sf, return.type="data.frame")
    checkTrue(all(Res$symbol == "SKA2"))
    Res <- transcripts(edb, filter = sf, return.type = "data.frame",
                       columns = c("symbol", "tx_id", "gene_name"))
    checkTrue(all(Res$symbol == "SKA2"))
    checkEquals(Res$symbol, Res$gene_name)
    checkEquals(colnames(Res), c("symbol", "tx_id", "gene_name"))

    ## Test for exons
    Res <- exons(edb, filter=sf, return.type="data.frame")
    checkTrue(all(Res$symbol == "SKA2"))
    Res <- exons(edb, filter = c(sf, TxBiotypeFilter("nonsense_mediated_decay")),
                 return.type = "data.frame",
                 columns = c("symbol", "tx_id", "gene_name"))
    checkTrue(all(Res$symbol == "SKA2"))
    checkEquals(Res$symbol, Res$gene_name)
    checkEquals(colnames(Res), c("symbol", "tx_id", "gene_name", "exon_id",
                                 "tx_biotype"))

    ## Test for exonsBy
    Res <- exonsBy(edb, filter=sf)
    checkTrue(all(unlist(Res)$symbol == "SKA2"))
    Res <- exonsBy(edb, filter = c(sf, TxBiotypeFilter("nonsense_mediated_decay")),
                 columns = c("symbol", "tx_id", "gene_name"))
    checkTrue(all(unlist(Res)$symbol == "SKA2"))

    checkEquals(unlist(Res)$symbol, unlist(Res)$gene_name)
}

## Here we want to test if we get always also the filter columns back.
test_multiFilterReturnCols <- function() {
    cols <- ensembldb:::addFilterColumns(edb, cols = c("exon_id"),
                                         filter = SymbolFilter("SKA2"))
    checkEquals(cols, c("exon_id", "symbol"))
    ## Two filter
    cols <- ensembldb:::addFilterColumns(edb, cols = c("exon_id"),
                                         filter = list(SymbolFilter("SKA2"),
                                                       GenenameFilter("SKA2")))
    checkEquals(cols, c("exon_id", "symbol", "gene_name"))
    cols <- ensembldb:::addFilterColumns(edb, cols = c("exon_id"),
                                         filter = list(SymbolFilter("SKA2"),
                                                       GenenameFilter("SKA2"),
                                                       GRangesFilter(
                                                           GRanges("3",
                                                                   IRanges(3, 5)
                                                                   ))))
    checkEquals(cols, c("exon_id", "symbol", "gene_name", "gene_seq_start",
                        "gene_seq_end", "seq_name", "seq_strand"))
    cols <- ensembldb:::addFilterColumns(edb, cols = c("exon_id"),
                                         filter = list(SymbolFilter("SKA2"),
                                                       GenenameFilter("SKA2"),
                                                       GRangesFilter(
                                                           GRanges("3",
                                                                   IRanges(3, 5)
                                                                   ),
                                                           feature = "exon")))
    checkEquals(cols, c("exon_id", "symbol", "gene_name", "exon_seq_start",
                        "exon_seq_end", "seq_name", "seq_strand"))
    ## SeqStartFilter and GRangesFilter
    ssf <- TxStartFilter(123)
    cols <- ensembldb:::addFilterColumns(edb, cols = c("exon_id"),
                                         filter = list(SymbolFilter("SKA2"),
                                                       GenenameFilter("SKA2"),
                                                       GRangesFilter(
                                                           GRanges("3",
                                                                   IRanges(3, 5)
                                                                   ),
                                                           feature = "exon"),
                                                       ssf))
    checkEquals(cols, c("exon_id", "symbol", "gene_name", "exon_seq_start",
                        "exon_seq_end", "seq_name", "seq_strand", "tx_seq_start"))

}


