###================================================
##  Here we check functionality to use EnsDbs with
##  UCSC chromosome names
###------------------------------------------------
library(EnsDb.Hsapiens.v75)
edb <- EnsDb.Hsapiens.v75

## library(EnsDb.Hsapiens.v83)
## edb <- EnsDb.Hsapiens.v83
## library(EnsDb.Hsapiens.v81)

test_seqlevels <- function(){
    orig <- getOption("ensembldb.seqnameNotFound")
    options(ensembldb.seqnameNotFound=NA)
    edb <- EnsDb.Hsapiens.v75
    SL <- seqlevels(edb)
    ucscs <- paste0("chr", c(1:22, "X", "Y", "M"))
    seqlevelsStyle(edb) <- "UCSC"
    SL2 <- seqlevels(edb)
    checkEquals(sort(ucscs), sort(SL2[!is.na(SL2)]))
    ## Check if we throw an error message
    options(ensembldb.seqnameNotFound="MISSING")
    checkException(seqlevels(edb))
    ## Check if returning original names works.
    options(ensembldb.seqnameNotFound="ORIGINAL")
    SL3 <- seqlevels(edb)
    idx <- which(SL3 %in% ucscs)
    checkEquals(sort(SL[-idx]), sort(SL3[-idx]))
    options(ensembldb.seqnameNotFound=orig)
}

test_seqinfo <- function(){
    edb <- EnsDb.Hsapiens.v75
    orig <- getOption("ensembldb.seqnameNotFound")
    options(ensembldb.seqnameNotFound="MISSING")
    seqlevelsStyle(edb) <- "UCSC"
    checkException(seqinfo(edb))
    options(ensembldb.seqnameNotFound="ORIGINAL")
    si <- seqinfo(edb)
    options(ensembldb.seqnameNotFound=orig)
}

## Testing if getWhat returns what we expect.
test_getWhat_seqnames <- function(){
    orig <- getOption("ensembldb.seqnameNotFound")
    edb <- EnsDb.Hsapiens.v75
    seqlevelsStyle(edb) <- "Ensembl"
    ensRes <- ensembldb:::getWhat(edb, columns=c("seq_name", "seq_strand"))
    seqlevelsStyle(edb) <- "UCSC"
    ucscRes <- ensembldb:::getWhat(edb, columns=c("seq_name", "seq_strand"))
    seqlevelsStyle(edb) <- "NCBI"
    ncbiRes <- ensembldb:::getWhat(edb, columns=c("seq_name", "seq_strand"))
    options(ensembldb.seqnameNotFound=orig)
}

test_SeqnameFilter_seqnames <- function(){
    orig <- getOption("ensembldb.seqnameNotFound")
    options(ensembldb.seqnameNotFound="MISSING")
    edb <- EnsDb.Hsapiens.v75
    seqlevelsStyle(edb) <- "Ensembl"
    snf <- SeqnameFilter("chrX")
    snfEns <- SeqnameFilter(c("X", "Y"))
    snfNo <- SeqnameFilter(c("bla", "blu"))
    snfSomeNo <- SeqnameFilter(c("bla", "X"))

    seqlevelsStyle(edb) <- "Ensembl"
    checkEquals(value(snf), "chrX")
    ## That makes no sense for a query though.
    checkEquals(value(snf, edb), "chrX")
    checkEquals(value(snfEns, edb), c("X", "Y"))
    seqlevelsStyle(edb) <- "UCSC"
    checkEquals(value(snf, edb), "X")
    checkException(value(snfEns, edb))
    checkException(value(snfNo, edb))
    checkException(value(snfSomeNo, edb))

    ## Setting the options to "ORIGINAL"
    options(ensembldb.seqnameNotFound="ORIGINAL")
    checkEquals(value(snf, edb), "X")
    checkEquals(value(snfEns, edb), c("X", "Y"))
    checkEquals(value(snfNo, edb), c("bla", "blu"))
    checkEquals(value(snfSomeNo, edb), c("bla", "X"))
    ##
    snf <- SeqnameFilter(c("chrX", "Y"))
    checkEquals(value(snf, edb), c("X", "Y"))

    options(ensembldb.seqnameNotFound=orig)
}

test_genes_seqnames <- function(){
    orig <- getOption("ensembldb.seqnameNotFound")
    edb <- EnsDb.Hsapiens.v75
    ## Here we want to test whether the result returned by the function does really
    ## work when changing the seqnames.
    seqlevelsStyle(edb) <- "Ensembl"
    ensAll <- genes(edb)
    ens21Y <- genes(edb, filter=SeqnameFilter(c("Y", "21")))
    checkEquals(sort(as.character(unique(seqnames(ens21Y)))), c("21", "Y"))
    gr <- GRanges(seqnames="Y", ranges=IRanges(start=1, end=59373566), strand="+")
    ensY <- genes(edb, filter=GRangesFilter(gr))
    checkEquals(seqlevels(ensY), "Y")
    checkEquals(unique(as.character(strand(ensY))), "+")

    ## Check UCSC stuff
    seqlevelsStyle(edb) <- "UCSC"
    options(ensembldb.seqnameNotFound="ORIGINAL")
    ## Just visually inspect the seqinfo and seqnames for the "all" query.
    ucscAll <- genes(edb)
    as.character(unique(seqnames(ucscAll)))
    ucsc21Y <- genes(edb, filter=SeqnameFilter(c("chrY", "chr21")))
    checkEquals(sort(seqlevels(ucsc21Y)), c("chr21", "chrY"))
    checkEquals(sort(names(ens21Y)), sort(names(ucsc21Y)))
    ## GRangesFilter.
    gr <- GRanges(seqnames="chrY", ranges=IRanges(start=1, end=59373566), strand="+")
    ucscY <- genes(edb, filter=GRangesFilter(gr))
    checkEquals(seqlevels(ucscY), "chrY")
    checkEquals(unique(as.character(strand(ucscY))), "+")
    checkEquals(sort(names(ensY)), sort(names(ucscY)))
    options(ensembldb.seqnameNotFound=orig)
}

test_transcripts_seqnames <- function(){
    orig <- getOption("ensembldb.seqnameNotFound")
    edb <- EnsDb.Hsapiens.v75
    seqlevelsStyle(edb) <- "Ensembl"
    ens21Y <- transcripts(edb, filter=SeqnameFilter(c("Y", "21")))
    checkEquals(sort(seqlevels(ens21Y)), c("21", "Y"))
    gr <- GRanges(seqnames="Y", ranges=IRanges(start=1, end=59373566), strand="+")
    ensY <- transcripts(edb, filter=GRangesFilter(gr))
    checkEquals(seqlevels(ensY), "Y")
    checkEquals(unique(as.character(strand(ensY))), "+")

    ## Check UCSC stuff
    seqlevelsStyle(edb) <- "UCSC"
    options(ensembldb.seqnameNotFound="ORIGINAL")
    ucsc21Y <- transcripts(edb, filter=SeqnameFilter(c("chrY", "chr21")))
    checkEquals(sort(seqlevels(ucsc21Y)), c("chr21", "chrY"))
    checkEquals(sort(names(ens21Y)), sort(names(ucsc21Y)))
    ## GRangesFilter.
    gr <- GRanges(seqnames="chrY", ranges=IRanges(start=1, end=59373566), strand="+")
    ucscY <- transcripts(edb, filter=GRangesFilter(gr))
    checkEquals(seqlevels(ucscY), "chrY")
    checkEquals(unique(as.character(strand(ucscY))), "+")
    checkEquals(sort(names(ensY)), sort(names(ucscY)))
    options(ensembldb.seqnameNotFound=orig)
}

test_transcriptsBy_seqnames <- function(){
    orig <- getOption("ensembldb.seqnameNotFound")
    edb <- EnsDb.Hsapiens.v75
    seqlevelsStyle(edb) <- "Ensembl"
    ens21Y <- transcriptsBy(edb, filter=SeqnameFilter(c("Y", "21")))
    checkEquals(sort(seqlevels(ens21Y)), c("21", "Y"))
    gr <- GRanges(seqnames="Y", ranges=IRanges(start=1, end=59373566), strand="+")
    ensY <- transcriptsBy(edb, filter=GRangesFilter(gr))
    checkEquals(seqlevels(ensY), "Y")
    checkEquals(unique(as.character(unlist(strand(ensY)))), "+")

    ## Check UCSC stuff
    seqlevelsStyle(edb) <- "UCSC"
    options(ensembldb.seqnameNotFound="ORIGINAL")
    ucsc21Y <- transcriptsBy(edb, filter=SeqnameFilter(c("chrY", "chr21")))
    checkEquals(sort(seqlevels(ucsc21Y)), c("chr21", "chrY"))
    checkEquals(sort(names(ens21Y)), sort(names(ucsc21Y)))
    ## GRangesFilter.
    gr <- GRanges(seqnames="chrY", ranges=IRanges(start=1, end=59373566), strand="+")
    ucscY <- transcriptsBy(edb, filter=GRangesFilter(gr))
    checkEquals(seqlevels(ucscY), "chrY")
    checkEquals(unique(as.character(unlist(strand(ucscY)))), "+")
    checkEquals(sort(names(ensY)), sort(names(ucscY)))
    options(ensembldb.seqnameNotFound=orig)
}

test_exons_seqnames <- function(){
    orig <- getOption("ensembldb.seqnameNotFound")
    edb <- EnsDb.Hsapiens.v75
    seqlevelsStyle(edb) <- "Ensembl"
    ens21Y <- exons(edb, filter=SeqnameFilter(c("Y", "21")))
    checkEquals(sort(seqlevels(ens21Y)), c("21", "Y"))
    gr <- GRanges(seqnames="Y", ranges=IRanges(start=1, end=59373566), strand="+")
    ensY <- exons(edb, filter=GRangesFilter(gr))
    checkEquals(seqlevels(ensY), "Y")
    checkEquals(unique(as.character(strand(ensY))), "+")

    ## Check UCSC stuff
    seqlevelsStyle(edb) <- "UCSC"
    options(ensembldb.seqnameNotFound="ORIGINAL")
    ucsc21Y <- exons(edb, filter=SeqnameFilter(c("chrY", "chr21")))
    checkEquals(sort(seqlevels(ucsc21Y)), c("chr21", "chrY"))
    checkEquals(sort(names(ens21Y)), sort(names(ucsc21Y)))
    ## GRangesFilter.
    gr <- GRanges(seqnames="chrY", ranges=IRanges(start=1, end=59373566), strand="+")
    ucscY <- exons(edb, filter=GRangesFilter(gr))
    checkEquals(seqlevels(ucscY), "chrY")
    checkEquals(unique(as.character(strand(ucscY))), "+")
    checkEquals(sort(names(ensY)), sort(names(ucscY)))
    options(ensembldb.seqnameNotFound=orig)
}

test_exonsBy_seqnames <- function(){
    orig <- getOption("ensembldb.seqnameNotFound")
    edb <- EnsDb.Hsapiens.v75
    seqlevelsStyle(edb) <- "Ensembl"
    ens21Y <- exonsBy(edb, filter=SeqnameFilter(c("Y", "21")))
    checkEquals(sort(seqlevels(ens21Y)), c("21", "Y"))
    gr <- GRanges(seqnames="Y", ranges=IRanges(start=1, end=59373566), strand="+")
    ensY <- exonsBy(edb, filter=GRangesFilter(gr))
    checkEquals(seqlevels(ensY), "Y")
    checkEquals(unique(as.character(unlist(strand(ensY)))), "+")

    ## Check UCSC stuff
    seqlevelsStyle(edb) <- "UCSC"
    options(ensembldb.seqnameNotFound="ORIGINAL")
    ucsc21Y <- exonsBy(edb, filter=SeqnameFilter(c("chrY", "chr21")))
    checkEquals(sort(seqlevels(ucsc21Y)), c("chr21", "chrY"))
    checkEquals(sort(names(ens21Y)), sort(names(ucsc21Y)))
    ## GRangesFilter.
    gr <- GRanges(seqnames="chrY", ranges=IRanges(start=1, end=59373566), strand="+")
    ucscY <- exonsBy(edb, filter=GRangesFilter(gr))
    checkEquals(seqlevels(ucscY), "chrY")
    checkEquals(unique(as.character(unlist(strand(ucscY)))), "+")
    checkEquals(sort(names(ensY)), sort(names(ucscY)))
    options(ensembldb.seqnameNotFound=orig)
}


test_cdsBy_seqnames <- function(){
    orig <- getOption("ensembldb.seqnameNotFound")
    edb <- EnsDb.Hsapiens.v75
    seqlevelsStyle(edb) <- "Ensembl"
    ens21Y <- cdsBy(edb, filter=SeqnameFilter(c("Y", "21")))
    checkEquals(sort(seqlevels(ens21Y)), c("21", "Y"))
    gr <- GRanges(seqnames="Y", ranges=IRanges(start=1, end=59373566), strand="+")
    ensY <- cdsBy(edb, filter=GRangesFilter(gr))
    checkEquals(seqlevels(ensY), "Y")
    checkEquals(unique(as.character(unlist(strand(ensY)))), "+")

    ## Check UCSC stuff
    seqlevelsStyle(edb) <- "UCSC"
    options(ensembldb.seqnameNotFound="ORIGINAL")
    ucsc21Y <- cdsBy(edb, filter=SeqnameFilter(c("chrY", "chr21")))
    checkEquals(sort(seqlevels(ucsc21Y)), c("chr21", "chrY"))
    checkEquals(sort(names(ens21Y)), sort(names(ucsc21Y)))
    ## GRangesFilter.
    gr <- GRanges(seqnames="chrY", ranges=IRanges(start=1, end=59373566), strand="+")
    ucscY <- cdsBy(edb, filter=GRangesFilter(gr))
    checkEquals(seqlevels(ucscY), "chrY")
    checkEquals(unique(as.character(unlist(strand(ucscY)))), "+")
    checkEquals(sort(names(ensY)), sort(names(ucscY)))
    options(ensembldb.seqnameNotFound=orig)
}

test_threeUTRsByTranscript_seqnames <- function(){
    orig <- getOption("ensembldb.seqnameNotFound")
    edb <- EnsDb.Hsapiens.v75
    seqlevelsStyle(edb) <- "Ensembl"
    ens21Y <- threeUTRsByTranscript(edb, filter=SeqnameFilter(c("Y", "21")))
    checkEquals(sort(seqlevels(ens21Y)), c("21", "Y"))
    gr <- GRanges(seqnames="Y", ranges=IRanges(start=1, end=59373566), strand="+")
    ensY <- threeUTRsByTranscript(edb, filter=GRangesFilter(gr))
    checkEquals(seqlevels(ensY), "Y")
    checkEquals(unique(as.character(unlist(strand(ensY)))), "+")

    ## Check UCSC stuff
    seqlevelsStyle(edb) <- "UCSC"
    options(ensembldb.seqnameNotFound="ORIGINAL")
    ucsc21Y <- threeUTRsByTranscript(edb, filter=SeqnameFilter(c("chrY", "chr21")))
    checkEquals(sort(seqlevels(ucsc21Y)), c("chr21", "chrY"))
    checkEquals(sort(names(ens21Y)), sort(names(ucsc21Y)))
    ## GRangesFilter.
    gr <- GRanges(seqnames="chrY", ranges=IRanges(start=1, end=59373566), strand="+")
    ucscY <- threeUTRsByTranscript(edb, filter=GRangesFilter(gr))
    checkEquals(seqlevels(ucscY), "chrY")
    checkEquals(unique(as.character(unlist(strand(ucscY)))), "+")
    checkEquals(sort(names(ensY)), sort(names(ucscY)))
    options(ensembldb.seqnameNotFound=orig)
}

test_fiveUTRsByTranscript_seqnames <- function(){
    orig <- getOption("ensembldb.seqnameNotFound")
    edb <- EnsDb.Hsapiens.v75
    seqlevelsStyle(edb) <- "Ensembl"
    ens21Y <- fiveUTRsByTranscript(edb, filter=SeqnameFilter(c("Y", "21")))
    checkEquals(sort(seqlevels(ens21Y)), c("21", "Y"))
    gr <- GRanges(seqnames="Y", ranges=IRanges(start=1, end=59373566), strand="+")
    ensY <- fiveUTRsByTranscript(edb, filter=GRangesFilter(gr))
    checkEquals(seqlevels(ensY), "Y")
    checkEquals(unique(as.character(unlist(strand(ensY)))), "+")

    ## Check UCSC stuff
    seqlevelsStyle(edb) <- "UCSC"
    options(ensembldb.seqnameNotFound="ORIGINAL")
    ucsc21Y <- fiveUTRsByTranscript(edb, filter=SeqnameFilter(c("chrY", "chr21")))
    checkEquals(sort(seqlevels(ucsc21Y)), c("chr21", "chrY"))
    checkEquals(sort(names(ens21Y)), sort(names(ucsc21Y)))
    ## GRangesFilter.
    gr <- GRanges(seqnames="chrY", ranges=IRanges(start=1, end=59373566), strand="+")
    ucscY <- fiveUTRsByTranscript(edb, filter=GRangesFilter(gr))
    checkEquals(seqlevels(ucscY), "chrY")
    checkEquals(unique(as.character(unlist(strand(ucscY)))), "+")
    checkEquals(sort(names(ensY)), sort(names(ucscY)))
    options(ensembldb.seqnameNotFound=orig)
}


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
    checkEquals(length(ensembldb:::properties(edb)), 2)
}

test_set_get_seqlevelsStyle <- function(){
    edb <- EnsDb.Hsapiens.v75
    ## Testing the getter/setter for the seqlevelsStyle.
    checkEquals(seqlevelsStyle(edb), "Ensembl")
    checkEquals(NA, ensembldb:::getProperty(edb, "seqlevelsStyle"))

    seqlevelsStyle(edb) <- "Ensembl"
    checkEquals(seqlevelsStyle(edb), "Ensembl")
    checkEquals("Ensembl", ensembldb:::getProperty(edb, "seqlevelsStyle"))

    ## Try NCBI.
    seqlevelsStyle(edb) <- "NCBI"
    checkEquals(seqlevelsStyle(edb), "NCBI")

    ## Try UCSC.
    seqlevelsStyle(edb) <- "UCSC"
    checkEquals(seqlevelsStyle(edb), "UCSC")

    ## Error checking:
    checkException(seqlevelsStyle(edb) <- "bla")
}

## Just dry run this without any actual query.
test_formatSeqnamesForQuery <- function(){
    ## Testing if the formating/mapping between seqnames works as expected
    ## We want to map anything TO Ensembl.
    ## Check also the warning messages!
    ucscs <- c("chr1", "chr3", "chr1", "chr9", "chrM", "chr1", "chrX")
    enses <- c("1", "3", "1", "9", "MT", "1", "X")
    ## reset
    edb <- EnsDb.Hsapiens.v75
    ## Shouldn't do anything here.
    seqlevelsStyle(edb)
    ensembldb:::dbSeqlevelsStyle(edb)
    got <- ensembldb:::formatSeqnamesForQuery(edb, enses)
    checkEquals(got, enses)
    ## Change the seqlevels to UCSC
    seqlevelsStyle(edb) <- "UCSC"
    ## If ifNotFound is not specified we suppose to get an error.
    options(ensembldb.seqnameNotFound="MISSING")
    checkException(ensembldb:::formatSeqnamesForQuery(edb, enses))
    ## With specifying ifNotFound
    got <- ensembldb:::formatSeqnamesForQuery(edb, enses, ifNotFound=NA)
    checkEquals(all(is.na(got)), TRUE)
    ## Same by setting the option
    options(ensembldb.seqnameNotFound=NA)
    got <- ensembldb:::formatSeqnamesForQuery(edb, enses)
    checkEquals(all(is.na(got)), TRUE)

    ## Now the working example:
    got <- ensembldb:::formatSeqnamesForQuery(edb, ucscs)
    checkEquals(got, enses)
    ## What if one is not mappable:
    got <- ensembldb:::formatSeqnamesForQuery(edb, c(ucscs, "asdfd"), ifNotFound=NA)
    checkEquals(got, c(enses, NA))
}

## Just dry run this without any actual query
test_formatSeqnamesFromQuery <- function(){
    ucscs <- c("chr1", "chr3", "chr1", "chr9", "chrM", "chr1", "chrX")
    enses <- c("1", "3", "1", "9", "MT", "1", "X")
    edb <- EnsDb.Hsapiens.v75
    ## Shouldn't do anything here.
    seqlevelsStyle(edb)
    ensembldb:::dbSeqlevelsStyle(edb)
    got <- ensembldb:::formatSeqnamesFromQuery(edb, enses)
    checkEquals(got, enses)
    ## Change the seqlevels to UCSC
    seqlevelsStyle(edb) <- "UCSC"
    ## If ifNotFound is not specified we suppose to get an error.
    options(ensembldb.seqnameNotFound="MISSING")
    checkException(ensembldb:::formatSeqnamesFromQuery(edb, ucsc))
    ## With specifying ifNotFound
    got <- ensembldb:::formatSeqnamesFromQuery(edb, ucscs, ifNotFound=NA)
    checkEquals(all(is.na(got)), TRUE)
    ## Same using options
    options(ensembldb.seqnameNotFound=NA)
    got <- ensembldb:::formatSeqnamesFromQuery(edb, ucscs, ifNotFound=NA)
    checkEquals(all(is.na(got)), TRUE)
    ## Now the working example:
    got <- ensembldb:::formatSeqnamesFromQuery(edb, enses)
    checkEquals(got, ucscs)
    ## What if one is not mappable:
    got <- ensembldb:::formatSeqnamesFromQuery(edb, c(enses, "asdfd"), ifNotFound=NA)
    checkEquals(got, c(ucscs, NA))
    got <- ensembldb:::formatSeqnamesFromQuery(edb, c(enses, "asdfd"))
    checkEquals(got, c(ucscs, NA))
}

notrun_test_set_seqlevels <- function(){
    ## To test what happens if no mapping is available
    ##gff <- "/Users/jo/Projects/EnsDbs/83/gadus_morhua/Gadus_morhua.gadMor1.83.gff3.gz"
    library(AnnotationHub)
    ah <- AnnotationHub()
    ah <- ah["AH47962"]
    fromG <- ensDbFromAH(ah, outfile=tempfile())
    edb <- EnsDb(fromG)
    seqlevelsStyle(edb)
    checkException(seqlevelsStyle(edb) <- "UCSC")
}





deprecated_test_check_SeqnameFilter <- function(){
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

deprecated_test_check_retrieve_data <- function(){
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


notrun_check_get_sequence_bsgenome <- function(){
    edb <- EnsDb.Hsapiens.v75
    ## Using first the Ensembl fasta stuff.
    ensSeqs <- extractTranscriptSeqs(getGenomeFaFile(edb),
                                     exonsBy(edb, "tx", filter=SeqnameFilter("Y")))
    ## Now the same using the BSgenome stuff.
    seqlevelsStyle(edb) <- "UCSC"
    options(ensembldb.seqnameNotFound="ORIGINAL")
    exs <- exonsBy(edb, "tx", filter=SeqnameFilter("chrY"))
    library(BSgenome.Hsapiens.UCSC.hg19)
    bsg <- BSgenome.Hsapiens.UCSC.hg19
    ucscSeqs <- extractTranscriptSeqs(bsg, exs)

    checkEquals(as.character(ensSeqs), as.character(ucscSeqs))
}


## Use the stuff from GenomeInfoDb!
notrun_test_newstuff <- function(){
    library(GenomeInfoDb)
    Map <- mapSeqlevels(seqlevels(edb), style="Ensembl")
    Map <- mapSeqlevels(seqlevels(edb), style="UCSC")
    ## just check what's out there
    genomeStyles()
}


