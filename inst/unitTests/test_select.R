####============================================================
##  test cases for AnnotationDbi methods.
##
####------------------------------------------------------------
library(EnsDb.Hsapiens.v75)
edb <- EnsDb.Hsapiens.v75

test_columns <- function(){
    cols <- columns(edb)
    ## Don't expect to see any _ there...
    checkEquals(length(grep(cols, pattern="_")), 0)
}

test_keytypes <- function(){
    keyt <- keytypes(edb)
    checkEquals(all(c("GENEID", "EXONID", "TXID") %in% keyt), TRUE)
}

test_mapper <- function(){
    Test <- ensembldb:::ensDbColumnForColumn(edb, "GENEID")
    checkEquals(unname(Test), "gene_id")

    Test <- ensembldb:::ensDbColumnForColumn(edb, c("GENEID", "TXID"))
    checkEquals(unname(Test), c("gene_id", "tx_id"))

    Test <- ensembldb:::ensDbColumnForColumn(edb, c("GENEID", "TXID", "bla"))
    checkEquals(unname(Test), c("gene_id", "tx_id"))
}

test_keys <- function(){
    ## get all gene ids
    system.time(
        ids <- keys(edb, "GENEID")
    )
    checkEquals(length(ids), length(unique(ids)))
    ## get all tx ids
    system.time(
        ids <- keys(edb, "TXID")
    )
    ## Get the TXNAME...
    nms <- keys(edb, "TXNAME")
    checkEquals(nms, ids)
    checkEquals(length(ids), length(unique(ids)))
    ## get all gene names
    system.time(
        ids <- keys(edb, "GENENAME")
    )
    checkEquals(length(ids), length(unique(ids)))
    ## get all seq names
    system.time(
        ids <- keys(edb, "SEQNAME")
    )
    checkEquals(length(ids), length(unique(ids)))
    ## get all seq strands
    system.time(
        ids <- keys(edb, "SEQSTRAND")
    )
    checkEquals(length(ids), length(unique(ids)))
    ## get all gene biotypes
    system.time(
        ids <- keys(edb, "GENEBIOTYPE")
    )
    checkEquals(ids, listGenebiotypes(edb))
}

test_select <- function(){
    ## Test:
    ## Provide GenenameFilter.
    gf <- GenenameFilter("BCL2")
    system.time(
        Test <- select(edb, keys=gf)
    )
    ## Provide list of GenenameFilter and TxbiotypeFilter.
    Test2 <- select(edb, keys=list(gf, TxbiotypeFilter("protein_coding")))
    checkEquals(Test$EXONID[Test$TXBIOTYPE == "protein_coding"], Test2$EXONID)
    ## Choose selected columns.
    Test3 <- select(edb, keys=gf, columns=c("GENEID", "GENENAME", "SEQNAME"))
    checkEquals(unique(Test[, c("GENEID", "GENENAME", "SEQNAME")]), Test3)
    ## Provide keys.
    Test4 <- select(edb, keys="BCL2", keytype="GENENAME")
    checkEquals(Test[, colnames(Test4)], Test4)
    txs <- keys(edb, "TXID")
    ## Just get stuff from the tx table; should be faster.
    system.time(
        Test <- select(edb, keys=txs, columns=c("TXID", "TXBIOTYPE", "GENEID"), keytype="TXID")
    )
    checkEquals(all(Test$TXID==txs), TRUE)
    ## Get all lincRNA genes
    Test <- select(edb, keys="lincRNA", columns=c("GENEID", "GENEBIOTYPE", "GENENAME"),
                   keytype="GENEBIOTYPE")
    Test2 <- select(edb, keys=GenebiotypeFilter("lincRNA"),
                    columns=c("GENEID", "GENEBIOTYPE", "GENENAME"))
    checkEquals(Test[, colnames(Test2)], Test2)
    ## All on chromosome 21
    Test <- select(edb, keys="21", columns=c("GENEID", "GENEBIOTYPE", "GENENAME"),
                   keytype="SEQNAME")
    Test2 <- select(edb, keys=SeqnameFilter("21"), columns=c("GENEID", "GENEBIOTYPE", "GENENAME"))
    checkEquals(Test[, colnames(Test2)], Test2)
    ## What if we can't find it?
    Test <- select(edb, keys="bla", columns=c("GENEID", "GENENAME"), keytype="GENENAME")
    ## Run the full thing.
    ## system.time(
    ##     All <- select(edb)
    ## )
    ## Test <- select(edb, keys=txs, keytype="TXID")
    ## checkEquals(Test, All)
    Test <- select(edb, keys="ENST00000000233", columns=c("GENEID", "GENENAME"), keytype="TXNAME")
    checkEquals(Test$TXNAME, "ENST00000000233")
    ## Check what happens if we just add TXNAME and also TXID.
    Test2 <- select(edb, keys=list(gf, TxbiotypeFilter("protein_coding")), columns=c("TXID", "TXNAME",
                                                                                     "GENENAME", "GENEID"))

}

test_mapIds <- function(){
    ## Simple... map gene ids to gene names
    allgenes <- keys(edb, keytype="GENEID")
    randordergenes <- allgenes[sample(1:length(allgenes), 100)]
    system.time(
        mi <- mapIds(edb, keys=allgenes, keytype="GENEID", column = "GENENAME")
    )
    checkEquals(allgenes, names(mi))
    ## What happens if the ordering is different:
    mi <- mapIds(edb, keys=randordergenes, keytype="GENEID", column = "GENENAME")
    checkEquals(randordergenes, names(mi))

    ## Now check the different options:
    ## Handle multi mappings.
    ## first
    first <- mapIds(edb, keys=randordergenes, keytype="GENEID", column="TXID")
    checkEquals(names(first), randordergenes)
    ## list
    lis <- mapIds(edb, keys=randordergenes, keytype="GENEID", column="TXID", multiVals="list")
    checkEquals(names(lis), randordergenes)
    Test <- lapply(lis, function(z){return(z[1])})
    checkEquals(first, unlist(Test))
    ## filter
    filt <- mapIds(edb, keys=randordergenes, keytype="GENEID", column="TXID", multiVals="filter")
    checkEquals(filt, unlist(lis[unlist(lapply(lis, length)) == 1]))
    ## asNA
    asNA <- mapIds(edb, keys=randordergenes, keytype="GENEID", column="TXID", multiVals="asNA")

    ## Check what happens if we provide 2 identical keys.
    Test <- mapIds(edb, keys=c("BCL2", "BCL2L11", "BCL2"), keytype="GENENAME", column="TXID")

    ## Submit Filter:
    Test <- mapIds(edb, keys=SeqnameFilter("Y"), column="GENEID", multiVals="list")
    TestS <- select(edb, keys=Test[[1]], columns="SEQNAME", keytype="GENEID")
    checkEquals(unique(TestS$SEQNAME), "Y")
    ## Submit 2 filter.
    Test <- mapIds(edb, keys=list(SeqnameFilter("Y"), SeqstrandFilter("-")), multiVals="list",
                   column="GENEID")
    TestS <- select(edb, keys=Test[[1]], keytype="GENEID", columns=c("SEQNAME", "SEQSTRAND"))
    checkTrue(all(TestS$SEQNAME == "Y"))
    checkTrue(all(TestS$SEQSTRAND == -1))
}

## Test if the results are properly sorted if we submit a single filter or just keys.
test_select_sorted <- function() {
    ks <- c("ZBTB16", "BCL2", "SKA2", "BCL2L11")
    ## gene_name
    res <- select(edb, keys = ks, keytype = "GENENAME")
    checkEquals(unique(res$GENENAME), ks)
    res <- select(edb, keys = GenenameFilter(ks))
    checkEquals(unique(res$GENENAME), ks)

    ## Using two filters;
    res <- select(edb, keys = list(GenenameFilter(ks),
                                   TxbiotypeFilter("nonsense_mediated_decay")))
    ## We don't expect same sorting here!
    checkTrue(!all(unique(res$GENENAME) == ks[ks %in% unique(res$GENENAME)]))

    ## symbol
    res <- select(edb, keys = ks, keytype = "SYMBOL",
                  columns = c("GENENAME", "SYMBOL", "SEQNAME"))

    ## tx_biotype
    ks <- c("retained_intron", "nonsense_mediated_decay")
    res <- select(edb, keys = ks, keytype = "TXBIOTYPE",
                  columns = c("GENENAME", "TXBIOTYPE"))
    checkEquals(unique(res$TXBIOTYPE), ks)
    res <- select(edb, keys = TxbiotypeFilter(ks),
                  keytype = "TXBIOTYPE", columns = c("GENENAME", "TXBIOTYPE"))
    checkEquals(unique(res$TXBIOTYPE), ks)
}

test_select_symbol <- function() {
    ## Can I use SYMBOL as keytype?
    ks <- c("ZBTB16", "BCL2", "SKA2", "BCL2L11")
    res <- select(edb, keys = ks, keytype = "GENENAME")
    res2 <- select(edb, keys = ks, keytype = "SYMBOL")
    checkEquals(res, res2)

    ## Can I use the SymbolFilter?
    res <- select(edb, keys = GenenameFilter(ks),
                  columns = c("TXNAME", "SYMBOL", "GENEID"))
    checkEquals(colnames(res), c("TXNAME", "SYMBOL", "GENEID", "GENENAME"))

    res <- select(edb, keys = SymbolFilter(ks), columns=c("GENEID"))
    checkEquals(colnames(res), c("GENEID", "SYMBOL"))
    checkEquals(res$SYMBOL, ks)

    ## Can I ask for SYMBOL?
    res <- select(edb, keys = list(SeqnameFilter("Y"),
                                   GenebiotypeFilter("lincRNA")),
                  columns = c("GENEID", "SYMBOL"))
    checkEquals(colnames(res), c("GENEID", "SYMBOL", "SEQNAME", "GENEBIOTYPE"))
}

test_select_symbol_n_txname <- function() {
    ks <- c("ZBTB16", "BCL2", "SKA2")
    ## Symbol allowed in keytype
    res <- select(edb, keys = ks, keytype = "SYMBOL", columns = "GENENAME")
    checkEquals(colnames(res), c("SYMBOL", "GENENAME"))
    checkEquals(res$SYMBOL, ks)

    ## Symbol using SymbolFilter
    res <- select(edb, keys = SymbolFilter(ks), columns = "GENENAME")
    checkEquals(colnames(res), c("GENENAME", "SYMBOL"))
    checkEquals(res$SYMBOL, ks)

    ## Symbol as a column.
    res <- select(edb, keys = ks, keytype = "GENENAME", columns = "SYMBOL")
    checkEquals(colnames(res), c("GENENAME", "SYMBOL"))

    ## TXNAME as a column
    res <- select(edb, keys = ks, keytype = "GENENAME", columns = c("TXNAME"))
    checkEquals(colnames(res), c("GENENAME", "TXNAME"))
}
