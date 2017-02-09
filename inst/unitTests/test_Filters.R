library("EnsDb.Hsapiens.v75")
edb <- EnsDb.Hsapiens.v75

## test_validateConditionFilter <- function() {
##     GF <- GeneidFilter("ENSG0000001")
##     GF@condition = "*"
##     checkTrue(is.character(ensembldb:::validateConditionFilter(GF)))
##     F <- SeqstartFilter(14)
##     F@condition = "*"
##     checkTrue(is.character(ensembldb:::validateConditionFilter(F)))
##     F@condition = ">"
##     F@value = "a"
##     checkTrue(is.character(ensembldb:::validateConditionFilter(F)))
## }

## test_dotWhere <- function() {
##     GF <- GeneidFilter("ENSG0000001")
##     checkEquals(ensembldb:::.where(GF), "= 'ENSG0000001'")
## }

## test_condition <- function() {
##     GF <- GeneidFilter("ENSG0000001")
##     checkException(condition(GF) <- "*")
##     checkEquals(condition(GF), "=")
## }

## test_value <- function() {
##     GF <- GeneidFilter("ENSG0000001")
##     value(GF) <- "a"
##     checkEquals(value(GF), "a")
##     F <- SeqstartFilter(14)
##     checkException(value(F) <- "a")
## }

test_where_list <- function() {
    filts <- list(GenenameFilter("BCL2L11"), GenenameFilter("BCL2"))
    res <- where(filts)
    checkEquals(res, " where gene_name = 'BCL2L11' and gene_name = 'BCL2'")
    res <- where(filts, edb)
    checkEquals(res, paste0(" where gene.gene_name = 'BCL2L11'",
                            " and gene.gene_name = 'BCL2'"))
    res <- where(filts, edb, with.tables = "gene")
    checkEquals(res, paste0(" where gene.gene_name = 'BCL2L11'",
                            " and gene.gene_name = 'BCL2'"))

}

## testing GeneidFilter
test_GeneidFilter <- function(){
    GF <- GeneidFilter("ENSG0000001")
    ## check if column matches the present database.
    checkEquals(column(GF, edb), "gene.gene_id")
    checkEquals(column(GF, edb, with.tables = "tx"),
                "tx.gene_id")
    ## check error if value is not as expected.
    checkException(GeneidFilter("ENSG000001", ">"))
    ## expect the filter to change the condition if lenght of values
    ## is > 1
    checkMultiValsIn(GeneidFilter(c("a", "b"), "="))
    checkMultiValsNotIn(GeneidFilter(c("a", "b"), "!="))
    ## where and column
    checkEquals(where(GF, edb), "gene.gene_id = 'ENSG0000001'")
    checkEquals(where(GF), "gene_id = 'ENSG0000001'")
    checkEquals(where(GF, edb, with.tables = "tx"),
                "tx.gene_id = 'ENSG0000001'")
}

test_GenebiotypeFilter <- function(){
    Filt <- GenebiotypeFilter("protein_coding")
    ## column
    checkEquals(column(Filt), "gene_biotype")
    checkEquals(column(Filt, edb), "gene.gene_biotype")
    checkEquals(column(Filt, edb, with.tables = "gene"),
                "gene.gene_biotype")
    checkException(GenebiotypeFilter("protein_coding", ">"))
    ## expect the filter to change the condition if lenght of values
    ## is > 1
    checkMultiValsIn(GenebiotypeFilter(c("a", "b"), "="))
    checkMultiValsNotIn(GenebiotypeFilter(c("a", "b"), "!="))

    ## where
    checkEquals(where(Filt), "gene_biotype = 'protein_coding'")
    checkEquals(where(Filt, edb), "gene.gene_biotype = 'protein_coding'")
    checkEquals(where(Filt, edb, with.tables = "gene"),
                "gene.gene_biotype = 'protein_coding'")
}

test_GenenameFilter <- function(){
    Filt <- GenenameFilter("genename")
    ## column
    checkEquals(column(Filt, edb), "gene.gene_name")
    checkEquals(column(Filt), "gene_name")
    checkEquals(column(Filt, edb, with.tables = "gene"), "gene.gene_name")
    checkException(GenenameFilter("genename", ">"))
    ## expect the filter to change the condition if lenght of values
    ## is > 1
    checkMultiValsIn(GenenameFilter(c("a", "b"), "="))
    checkMultiValsNotIn(GenenameFilter(c("a", "b"), "!="))

    ## where
    checkEquals(where(Filt, edb), "gene.gene_name = 'genename'")
    checkEquals(where(Filt), "gene_name = 'genename'")
    checkEquals(where(Filt, edb, with.tables = "gene"),
                "gene.gene_name = 'genename'")

    ## check if we're escaping correctly!
    Filt <- GenenameFilter("I'm a gene")
    checkEquals(where(Filt, EnsDb.Hsapiens.v75), "gene.gene_name = 'I''m a gene'")
}

test_EntrezidFilter <- function(){
    Filt <- EntrezidFilter("123")
    ## column
    checkEquals(column(Filt, edb), "gene.entrezid")
    checkEquals(column(Filt), "entrezid")
    checkEquals(column(Filt, edb, with.tables = "gene"), "gene.entrezid")
    ## where
    checkEquals(where(Filt, edb), "gene.entrezid = '123'")
    checkEquals(where(Filt), "entrezid = '123'")
    checkEquals(where(Filt, edb, with.tables = "gene"),
                "gene.entrezid = '123'")

}

test_TxidFilter <- function(){
    Filt <- TxidFilter("a")
    ## column
    checkEquals(column(Filt), "tx_id")
    checkEquals(column(Filt, edb), "tx.tx_id")
    checkEquals(column(Filt, edb, with.tables = "tx2exon"), "tx2exon.tx_id")
    checkException(TxidFilter("a", ">"))
    ## expect the filter to change the condition if lenght of values
    ## is > 1
    checkMultiValsIn(TxidFilter(c("a", "b"), "="))
    checkMultiValsNotIn(TxidFilter(c("a", "b"), "!="))
    ## where
    checkEquals(where(Filt), "tx_id = 'a'")
    checkEquals(where(Filt, edb), "tx.tx_id = 'a'")
    checkEquals(where(Filt, edb, with.tables = "tx2exon"),
                "tx2exon.tx_id = 'a'")

}

test_TxbiotypeFilter <- function(){
    Filt <- TxbiotypeFilter("a")
    ## column
    checkEquals(column(Filt), "tx_biotype")
    checkEquals(column(Filt, edb), "tx.tx_biotype")
    checkEquals(column(Filt, edb, with.tables = "tx"), "tx.tx_biotype")
    checkException(TxbiotypeFilter("a", ">"))
    ## expect the filter to change the condition if lenght of values
    ## is > 1
    checkMultiValsIn(TxbiotypeFilter(c("a", "b"), "="))
    checkMultiValsNotIn(TxbiotypeFilter(c("a", "b"), "!="))
    ## where
    checkEquals(where(Filt), "tx_biotype = 'a'")
    checkEquals(where(Filt, edb), "tx.tx_biotype = 'a'")
    checkEquals(where(Filt, edb, with.tables = "tx"), "tx.tx_biotype = 'a'")
}

test_ExonidFilter <- function(){
    Filt <- ExonidFilter("a")
    ## column
    checkEquals(column(Filt), "exon_id")
    checkEquals(column(Filt, edb), "tx2exon.exon_id")
    checkEquals(column(Filt, edb, with.tables = "exon"), "exon.exon_id")
    checkException(ExonidFilter("a", ">"))
    ## expect the filter to change the condition if lenght of values
    ## is > 1
    checkMultiValsIn(ExonidFilter(c("a", "b"), "="))
    checkMultiValsNotIn(ExonidFilter(c("a", "b"), "!="))
    ## where
    checkEquals(where(Filt), "exon_id = 'a'")
    checkEquals(where(Filt, edb), "tx2exon.exon_id = 'a'")
    checkEquals(where(Filt, edb, with.tables = "exon"), "exon.exon_id = 'a'")
}

test_GRangesFilter <- function() {
    gr <- GRanges(seqnames = "a",
                  ranges = IRanges(start = 1, end = 5))
    F <- GRangesFilter(value = gr)
    checkEquals(unname(column(F)), c("gene_seq_start", "gene_seq_end",
                                     "seq_name", "seq_strand"))
    checkEquals(unname(column(F, edb)),
                c("gene.gene_seq_start", "gene.gene_seq_end",
                  "gene.seq_name", "gene.seq_strand"))

    checkEquals(where(F), paste0("gene_seq_start >= 1 and gene_seq_end",
                                 " <= 5 and seq_name == 'a'"))
    checkEquals(where(F, edb), paste0("gene.gene_seq_start >= 1 and",
                                      " gene.gene_seq_end",
                                      " <= 5 and gene.seq_name == 'a'"))
    ## tx
    F@feature <- "tx"
    checkEquals(unname(column(F)), c("tx_seq_start", "tx_seq_end",
                                     "seq_name", "seq_strand"))
    checkEquals(unname(column(F, edb)),
                c("tx.tx_seq_start", "tx.tx_seq_end",
                  "gene.seq_name", "gene.seq_strand"))
    ## exon
    F@feature <- "exon"
    checkEquals(unname(column(F)), c("exon_seq_start", "exon_seq_end",
                                     "seq_name", "seq_strand"))
    checkEquals(unname(column(F, edb)),
                c("exon.exon_seq_start", "exon.exon_seq_end",
                  "gene.seq_name", "gene.seq_strand"))
    ## Check the buildWhere
    res <- ensembldb:::buildWhereForGRanges(F, columns = column(F))
    checkEquals(res,
                "exon_seq_start >= 1 and exon_seq_end <= 5 and seq_name == 'a'")
}

test_strand2something <- function() {
    checkEquals(ensembldb:::strand2num("+"), 1)
    checkEquals(ensembldb:::strand2num("-"), -1)
    checkException(ensembldb:::strand2num("a"))
    checkEquals(ensembldb:::num2strand(-1), "-")
    checkEquals(ensembldb:::num2strand(1), "+")
}

## SeqnameFilter
test_SeqnameFilter <- function(){
    Filt <- SeqnameFilter("a")
    checkEquals(column(Filt), "seq_name")
    checkEquals(column(Filt, edb), "gene.seq_name")
    checkEquals(column(Filt, edb, with.tables = "chromosome"),
                "chromosome.seq_name")
    checkException(SeqnameFilter("a", ">"))
    checkEquals(where(Filt), "seq_name = 'a'")
    checkEquals(where(Filt, edb), "gene.seq_name = 'a'")
    checkEquals(where(Filt, edb, with.tables = "chromosome"),
                "chromosome.seq_name = 'a'")
}

## SeqstrandFilter
test_SeqstrandFilter <- function(){
    checkException(SeqstrandFilter("a"))
    Filt <- SeqstrandFilter("-")
    checkEquals(column(Filt), "seq_strand")
    checkEquals(column(Filt, edb), "gene.seq_strand")
    checkEquals(column(Filt, edb, with.tables = "gene"), "gene.seq_strand")
    ## where
    checkEquals(where(Filt), "seq_strand = -1")
    checkEquals(where(Filt, edb), "gene.seq_strand = -1")
    checkEquals(where(Filt, edb, with.tables = "gene"), "gene.seq_strand = -1")

}

## SeqstartFilter, feature
test_SeqstartFilter <- function(){
    Filt <- SeqstartFilter(123, feature="gene")
    checkEquals(column(Filt), "gene_seq_start")
    checkEquals(column(Filt, edb), "gene.gene_seq_start")
    checkEquals(column(Filt, edb, with.tables = "gene"), "gene.gene_seq_start")
    ## where
    checkEquals(where(Filt), "gene_seq_start = 123")
    checkEquals(where(Filt, edb), "gene.gene_seq_start = 123")
    checkEquals(where(Filt, edb, with.tables = "gene"),
                "gene.gene_seq_start = 123")

    Filt <- SeqstartFilter(123, feature="transcript")
    ## column
    checkEquals(column(Filt), "tx_seq_start")
    checkEquals(column(Filt, edb), "tx.tx_seq_start")
    checkEquals(column(Filt, edb, with.tables = "tx"), "tx.tx_seq_start")
    ## where
    checkEquals(where(Filt), "tx_seq_start = 123")
    checkEquals(where(Filt, edb), "tx.tx_seq_start = 123")
    checkEquals(where(Filt, edb, with.tables = "tx"), "tx.tx_seq_start = 123")

    Filt <- SeqstartFilter(123, feature="exon")
    ## column
    checkEquals(column(Filt), "exon_seq_start")
    checkEquals(column(Filt, edb), "exon.exon_seq_start")
    checkEquals(column(Filt, edb, with.tables = "exon"), "exon.exon_seq_start")
    ## where
    checkEquals(where(Filt), "exon_seq_start = 123")
    checkEquals(where(Filt, edb), "exon.exon_seq_start = 123")
    checkEquals(where(Filt, edb, with.tables = "exon"),
                "exon.exon_seq_start = 123")
}

## SeqendFilter
test_SeqendFilter <- function(){
    Filt <- SeqendFilter(123, feature="gene")
    ## column
    checkEquals(column(Filt), "gene_seq_end")
    checkEquals(column(Filt, edb), "gene.gene_seq_end")
    checkEquals(column(Filt, edb, with.tables = "gene"), "gene.gene_seq_end")
    ## where
    checkEquals(where(Filt), "gene_seq_end = 123")
    checkEquals(where(Filt, edb), "gene.gene_seq_end = 123")
    checkEquals(where(Filt, edb, with.tables = "gene"),
                "gene.gene_seq_end = 123")

    Filt <- SeqendFilter(123, feature="transcript")
    ## column
    checkEquals(column(Filt), "tx_seq_end")
    checkEquals(column(Filt, edb), "tx.tx_seq_end")
    checkEquals(column(Filt, edb, with.tables = "tx"), "tx.tx_seq_end")
    ## where
    checkEquals(where(Filt), "tx_seq_end = 123")
    checkEquals(where(Filt, edb), "tx.tx_seq_end = 123")
    checkEquals(where(Filt, edb, with.tables = "tx"), "tx.tx_seq_end = 123")

    Filt <- SeqendFilter(123, feature="exon")
    ## column
    checkEquals(column(Filt), "exon_seq_end")
    checkEquals(column(Filt, edb), "exon.exon_seq_end")
    checkEquals(column(Filt, edb, with.tables = "exon"), "exon.exon_seq_end")
    ## where
    checkEquals(where(Filt), "exon_seq_end = 123")
    checkEquals(where(Filt, edb), "exon.exon_seq_end = 123")
    checkEquals(where(Filt, edb, with.tables = "exon"),
                "exon.exon_seq_end = 123")
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
    checkException(value(Filt) <- "b")

    ## column
    checkEquals(column(Filt), "exon_idx")
    checkEquals(column(Filt, edb), "tx2exon.exon_idx")
    checkEquals(column(Filt, edb, with.tables = "tx2exon"),
                "tx2exon.exon_idx")
    ## where
    checkEquals(where(Filt), "exon_idx = 123")
    checkEquals(where(Filt, edb), "tx2exon.exon_idx = 123")
    checkEquals(where(Filt, edb, with.tables = "tx2exon"),
                "tx2exon.exon_idx = 123")
}

## SymbolFilter
test_SymbolFilter <- function() {
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
    Res <- exons(edb, filter = c(sf, TxbiotypeFilter("nonsense_mediated_decay")),
                 return.type = "data.frame",
                 columns = c("symbol", "tx_id", "gene_name"))
    checkTrue(all(Res$symbol == "SKA2"))
    checkEquals(Res$symbol, Res$gene_name)
    checkEquals(colnames(Res), c("symbol", "tx_id", "gene_name", "exon_id", "tx_biotype"))

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
                                                       GRangesFilter(GRanges("3",
                                                                             IRanges(3, 5)
                                                                             ))))
    checkEquals(cols, c("exon_id", "symbol", "gene_name", "gene_seq_start",
                        "gene_seq_end", "seq_name", "seq_strand"))
    cols <- ensembldb:::addFilterColumns(edb, cols = c("exon_id"),
                                         filter = list(SymbolFilter("SKA2"),
                                                       GenenameFilter("SKA2"),
                                                       GRangesFilter(GRanges("3",
                                                                             IRanges(3, 5)
                                                                             ),
                                                                     feature = "exon")))
    checkEquals(cols, c("exon_id", "symbol", "gene_name", "exon_seq_start",
                        "exon_seq_end", "seq_name", "seq_strand"))
    ## SeqstartFilter and GRangesFilter
    ssf <- SeqstartFilter(123, feature="tx")
    cols <- ensembldb:::addFilterColumns(edb, cols = c("exon_id"),
                                         filter = list(SymbolFilter("SKA2"),
                                                       GenenameFilter("SKA2"),
                                                       GRangesFilter(GRanges("3",
                                                                             IRanges(3, 5)
                                                                             ),
                                                                     feature = "exon"),
                                                       ssf))
    checkEquals(cols, c("exon_id", "symbol", "gene_name", "exon_seq_start",
                        "exon_seq_end", "seq_name", "seq_strand", "tx_seq_start"))

}

