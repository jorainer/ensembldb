library("EnsDb.Hsapiens.v75")

## testing GeneidFilter
test_GeneidFilter <- function(){
    GF <- GeneidFilter("ENSG0000001")
    ## check if column matches the present database.
    checkEquals(column(GF, EnsDb.Hsapiens.v75), "gene.gene_id")
    ## check error if value is not as expected.
    checkException(GeneidFilter("ENSG000001", ">"))
    ## expect the filter to change the condition if lenght of values
    ## is > 1
    checkMultiValsIn(GeneidFilter(c("a", "b"), "="))
    checkMultiValsNotIn(GeneidFilter(c("a", "b"), "!="))
}

test_GenebiotypeFilter <- function(){
    Filt <- GenebiotypeFilter("protein_coding")
    checkEquals(column(Filt, EnsDb.Hsapiens.v75), "gene.gene_biotype")
    checkException(GenebiotypeFilter("protein_coding", ">"))
    ## expect the filter to change the condition if lenght of values
    ## is > 1
    checkMultiValsIn(GenebiotypeFilter(c("a", "b"), "="))
    checkMultiValsNotIn(GenebiotypeFilter(c("a", "b"), "!="))

}

test_GenenameFilter <- function(){
    Filt <- GenenameFilter("genename")
    checkEquals(column(Filt, EnsDb.Hsapiens.v75), "gene.gene_name")
    checkException(GenenameFilter("genename", ">"))
    ## expect the filter to change the condition if lenght of values
    ## is > 1
    checkMultiValsIn(GenenameFilter(c("a", "b"), "="))
    checkMultiValsNotIn(GenenameFilter(c("a", "b"), "!="))
    ## check if we're escaping correctly!
    Filt <- GenenameFilter("I'm a gene")
    checkEquals(where(Filt, EnsDb.Hsapiens.v75), "gene.gene_name = 'I''m a gene'")
}

test_TxidFilter <- function(){
    Filt <- TxidFilter("a")
    checkEquals(column(Filt, EnsDb.Hsapiens.v75), "tx.tx_id")
    checkException(TxidFilter("a", ">"))
    ## expect the filter to change the condition if lenght of values
    ## is > 1
    checkMultiValsIn(TxidFilter(c("a", "b"), "="))
    checkMultiValsNotIn(TxidFilter(c("a", "b"), "!="))
}

test_TxbiotypeFilter <- function(){
    Filt <- TxbiotypeFilter("a")
    checkEquals(column(Filt, EnsDb.Hsapiens.v75), "tx.tx_biotype")
    checkException(TxbiotypeFilter("a", ">"))
    ## expect the filter to change the condition if lenght of values
    ## is > 1
    checkMultiValsIn(TxbiotypeFilter(c("a", "b"), "="))
    checkMultiValsNotIn(TxbiotypeFilter(c("a", "b"), "!="))
}

test_ExonidFilter <- function(){
    Filt <- ExonidFilter("a")
    checkEquals(column(Filt, EnsDb.Hsapiens.v75), "tx2exon.exon_id")
    checkException(ExonidFilter("a", ">"))
    ## expect the filter to change the condition if lenght of values
    ## is > 1
    checkMultiValsIn(ExonidFilter(c("a", "b"), "="))
    checkMultiValsNotIn(ExonidFilter(c("a", "b"), "!="))
}

## SeqnameFilter
test_SeqnameFilter <- function(){
    Filt <- SeqnameFilter("a")
    checkEquals(column(Filt, EnsDb.Hsapiens.v75), "gene.seq_name")
    checkException(SeqnameFilter("a", ">"))
}

## SeqstrandFilter
test_SeqstrandFilter <- function(){
    checkException(SeqstrandFilter("a"))
    Filt <- SeqstrandFilter("-")
    checkEquals(column(Filt, EnsDb.Hsapiens.v75), "gene.seq_strand")
}

## SeqstartFilter, feature
test_SeqstartFilter <- function(){
    Filt <- SeqstartFilter(123, feature="gene")
    checkEquals(column(Filt, EnsDb.Hsapiens.v75), "gene.gene_seq_start")
    Filt <- SeqstartFilter(123, feature="transcript")
    checkEquals(column(Filt, EnsDb.Hsapiens.v75), "tx.tx_seq_start")
}

## SeqendFilter
test_SeqendFilter <- function(){
    Filt <- SeqendFilter(123, feature="gene")
    checkEquals(column(Filt, EnsDb.Hsapiens.v75), "gene.gene_seq_end")
    Filt <- SeqendFilter(123, feature="transcript")
    checkEquals(column(Filt, EnsDb.Hsapiens.v75), "tx.tx_seq_end")
}



## checks if "condition" of the filter is "in"
checkMultiValsIn <- function(filt){
    checkEquals(condition(filt), "in")
}
## checks if "condition" of the filter is "in"
checkMultiValsNotIn <- function(filt){
    checkEquals(condition(filt), "not in")
}

test_ExonrankFilter <- function(){
    Filt <- ExonrankFilter(123)
    checkException(ExonrankFilter("a"))

    edb <- EnsDb.Hsapiens.v75
    checkException(value(Filt) <- "b")

    checkEquals(column(Filt), "exon_idx")
    checkEquals(column(Filt, edb), "tx2exon.exon_idx")
    where(Filt, edb)
}

## SymbolFilter
test_SymbolFilter <- function(){
    edb <- EnsDb.Hsapiens.v75
    sf <- SymbolFilter("SKA2")

    ## Check the column method.
    checkEquals(column(sf), "symbol")
    ## For EnsDb we want it to link to gene_name
    checkEquals(column(sf, edb), "gene.gene_name")
    checkException(column(sf, edb, with.tables = c("tx", "exon")))

    ## Check the where method.
    checkEquals(where(sf), "symbol = 'SKA2'")
    condition(sf) <- "!="
    checkEquals(where(sf, edb), "gene.gene_name != 'SKA2'")

    ## Test if we can use it:
    condition(sf) <- "="
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
    checkEquals(colnames(Res), c("tx_name", "symbol"))
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
    Res <- exons(edb, filter = c(sf, TxbiotypeFilter("nonsense_mediated_decay")),
                 return.type = "data.frame",
                 columns = c("symbol", "tx_id", "gene_name"))
    checkTrue(all(Res$symbol == "SKA2"))
    checkEquals(Res$symbol, Res$gene_name)
    checkEquals(colnames(Res), c("symbol", "tx_id", "gene_name"))

    ## Test for exonsBy
    Res <- exonsBy(edb, filter=sf)
    checkTrue(all(unlist(Res)$symbol == "SKA2"))
    Res <- exonsBy(edb, filter = c(sf, TxbiotypeFilter("nonsense_mediated_decay")),
                 columns = c("symbol", "tx_id", "gene_name"))
    checkTrue(all(unlist(Res)$symbol == "SKA2"))

    checkEquals(unlist(Res)$symbol, unlist(Res)$gene_name)

    ## Test for transcriptsBy too
}


## Here we want to test if we get always also the filter columns back.
test_multiFilterReturnCols <- function() {
    cols <- ensembldb:::addSymlinkFilterCols(cols=c("exon_id"), filter = SymbolFilter("SKA2"))
    checkEquals(cols, c("exon_id", "symbol"))
    ## Two filter
    cols <- ensembldb:::addSymlinkFilterCols(cols=c("exon_id"),
                                             filter = list(SymbolFilter("SKA2"),
                                                           GenenameFilter("SKA2")))
    checkEquals(cols, c("exon_id", "symbol", "gene_name"))
    cols <- ensembldb:::addSymlinkFilterCols(cols=c("exon_id"),
                                             filter = list(SymbolFilter("SKA2"),
                                                           GenenameFilter("SKA2"),
                                                           GRangesFilter(GRanges("3",
                                                                                 IRanges(3, 5)
                                                                                 ))))
    checkEquals(cols, c("exon_id", "symbol", "gene_name", "gene_seq_start",
                        "gene_seq_end", "seq_name", "seq_strand"))
    cols <- ensembldb:::addSymlinkFilterCols(cols=c("exon_id"),
                                             filter = list(SymbolFilter("SKA2"),
                                                           GenenameFilter("SKA2"),
                                                           GRangesFilter(GRanges("3",
                                                                                 IRanges(3, 5)
                                                                                 ))),
                                             feature = "exon")
    checkEquals(cols, c("exon_id", "symbol", "gene_name", "exon_seq_start",
                        "exon_seq_end", "seq_name", "seq_strand"))
    ## SeqstartFilter and GRangesFilter
    ssf <- SeqstartFilter(123, feature="tx")
    cols <- ensembldb:::addSymlinkFilterCols(cols=c("exon_id"),
                                             filter = list(SymbolFilter("SKA2"),
                                                           GenenameFilter("SKA2"),
                                                           GRangesFilter(GRanges("3",
                                                                                 IRanges(3, 5)
                                                                                 )),
                                                           ssf),
                                             feature = "exon")
    checkEquals(cols, c("exon_id", "symbol", "gene_name", "exon_seq_start",
                        "exon_seq_end", "seq_name", "seq_strand", "tx_seq_start"))

}

