###================================================
##  Here we check functionality to use EnsDbs with
##  UCSC chromosome names
###------------------------------------------------
library(EnsDb.Hsapiens.v75)
edb <- EnsDb.Hsapiens.v75

test_updateEnsDb <- function(){
    edb2 <- updateEnsDb(edb)
    checkEquals(edb2@tables, edb@tables)
    checkTrue(.hasSlot(edb2, ".properties"))
}

test_properties <- function(){
    checkEquals(ensembldb:::getProperty(edb, "foo"), NA)

    checkException(ensembldb:::setProperty(edb, "foo"))

    edb <- ensembldb:::setProperty(edb, foo="bar")
    checkEquals(ensembldb:::getProperty(edb, "foo"), "bar")
    checkEquals(length(ensembldb:::properties(edb)), 1)
}

test_check_SeqnameFilter <- function(){
    Orig <- getOption("ucscChromosomeNames", FALSE)
    options(ucscChromosomeNames=TRUE)
    snf <- SeqnameFilter(c("chrX", "chr3"))
    checkEquals(value(snf), c("chrX", "chr3"))
    checkEquals(value(snf, edb), c("X", "3"))

    options(ucscChromosomeNames=FALSE)
    checkEquals(value(snf, edb), c("X", "3"))

    ## No matter what, where has to return names without chr!
    checkEquals(where(snf, edb), "gene.seq_name in ('X','3')")

    ## GRangesFilter:
    grf <- GRangesFilter(GRanges("chrX", IRanges(123, 345)))
    checkEqualsNumeric(length(grep(where(grf), pattern="seq_name == 'chrX'")), 1)
    checkEqualsNumeric(length(grep(where(grf, edb), pattern="seq_name == 'X'")), 1)

    ## Check chromosome MT/chrM
    options(ucscChromosomeNames=FALSE)
    snf <- SeqnameFilter("MT")
    checkEquals(where(snf, edb), "gene.seq_name = 'MT'")
    snf <- SeqnameFilter("chrM")
    checkEquals(where(snf, edb), "gene.seq_name = 'MT'")
    options(ucscChromosomeNames=TRUE)
    snf <- SeqnameFilter("MT")
    checkEquals(where(snf, edb), "gene.seq_name = 'MT'")
    snf <- SeqnameFilter("chrM")
    checkEquals(where(snf, edb), "gene.seq_name = 'MT'")

    options(ucscChromosomeNames=Orig)
}

test_check_retrieve_data <- function(){
    Orig <- getOption("ucscChromosomeNames", FALSE)

    options(ucscChromosomeNames=FALSE)
    genes <- genes(edb, filter=SeqnameFilter(c("21", "Y", "X")))
    checkEquals(all(seqlevels(genes) %in% c("21", "X", "Y")), TRUE)
    options(ucscChromosomeNames=TRUE)
    genes <- genes(edb, filter=SeqnameFilter(c("21", "Y", "X")))
    checkEquals(all(seqlevels(genes) %in% c("chr21", "chrX", "chrY")), TRUE)

    ## Check chromosome MT
    options(ucscChromosomeNames=FALSE)
    exons <- exons(edb, filter=SeqnameFilter("MT"))
    checkEquals(seqlevels(exons), "MT")
    options(ucscChromosomeNames=TRUE)
    exons <- exons(edb, filter=SeqnameFilter("MT"))
    checkEquals(seqlevels(exons), "chrM")

    options(ucscChromosomeNames=Orig)
}

## Use the stuff from GenomeInfoDb!
notrun_test_newstuff <- function(){
    library(GenomeInfoDb)
    Map <- mapSeqlevels(seqlevels(edb), style="Ensembl")
    Map <- mapSeqlevels(seqlevels(edb), style="UCSC")
    ## just check what's out there
    genomeStyles()
}


