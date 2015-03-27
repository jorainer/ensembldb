## that's just a plain simple R-script calling the standard methods.

library( "EnsDb.Hsapiens.v75" )
DB <- EnsDb.Hsapiens.v75

## testing genes method.
test_genes <- function(){
    Gns <- genes(DB)
    ## checkEquals(class(genes(DB, return.type="DataFrame",
    ##                         filter=list(SeqnameFilter("Y")))), "DataFrame" )
}

test_transcripts <- function(){
    Tns <- transcripts(DB)
}

test_transcriptsBy <- function(){
    TnsBy <- transcriptsBy(DB, filter=list(SeqnameFilter("Y")), by="gene")
}

test_exons <- function(){
    Exns <- exons(DB, filter=list(SeqnameFilter("X")))
}

test_exonsBy <- function(){
    ExnsBy <- exonsBy(DB, filter=list(SeqnameFilter("X")), by="tx")
}

test_dbfunctionality <- function(){
    GBT <- listGenebiotypes(DB)
    TBT <- listTxbiotypes(DB)
}



