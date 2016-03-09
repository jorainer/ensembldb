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



