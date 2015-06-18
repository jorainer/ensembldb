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
    Exns <- exons(DB, filter=SeqnameFilter("X"))
    Exns <- exons(DB, filter=list(SeqnameFilter("X")))
}

test_exonsBy <- function(){
    ExnsBy <- exonsBy(DB, filter=list(SeqnameFilter("X")), by="tx")
}

test_dbfunctionality <- function(){
    GBT <- listGenebiotypes(DB)
    TBT <- listTxbiotypes(DB)
}

## test if we get the expected exceptions if we're not submitting
## correct filter objects
test_filterExceptions <- function(){
    checkException(genes(DB, filter="d"))
    checkException(genes(DB, filter=list(SeqnameFilter("X"),
                                 "z")))
    checkException(transcripts(DB, filter="d"))
    checkException(transcripts(DB, filter=list(SeqnameFilter("X"),
                                 "z")))
    checkException(exons(DB, filter="d"))
    checkException(exons(DB, filter=list(SeqnameFilter("X"),
                                 "z")))
    checkException(exonsBy(DB, filter="d"))
    checkException(exonsBy(DB, filter=list(SeqnameFilter("X"),
                                 "z")))
    checkException(transcriptsBy(DB, filter="d"))
    checkException(transcriptsBy(DB, filter=list(SeqnameFilter("X"),
                                 "z")))
}

test_promoters <- function(){
    promoters(EnsDb.Hsapiens.v75, filter=GeneidFilter(c("ENSG00000184895",
                                                    "ENSG00000092377")))
}


