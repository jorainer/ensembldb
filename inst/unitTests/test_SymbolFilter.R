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


############################################################
##   select method

