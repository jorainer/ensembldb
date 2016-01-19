library(EnsDb.Hsapiens.v75)
edb <- EnsDb.Hsapiens.v75

notrun_test_extractTranscriptSeqs <- function(){
    ## That's how we want to test the transcript seqs.
    ZBTB <- extractTranscriptSeqs(getGenomeFaFile(edb), edb, filter=GenenameFilter("ZBTB16"))

    library(TxDb.Hsapiens.UCSC.hg19.knownGene)
    txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
    library(BSgenome.Hsapiens.UCSC.hg19)
    genome <- BSgenome.Hsapiens.UCSC.hg19
    All <- extractTranscriptSeqs(genome, txdb)
}

notrun_test_getCdsSequence <- function(){
    ## That's when we like to get the sequence from the coding region.
    genome <- getGenomeFaFile(edb)
    tx <- extractTranscriptSeqs(genome, edb, filter=SeqnameFilter("Y"))
    cdsSeq <- extractTranscriptSeqs(genome, cdsBy(edb, filter=SeqnameFilter("Y")))
    ## that's basically to get the CDS sequence.
    ## UTR sequence:
    tutr <- extractTranscriptSeqs(genome, threeUTRsByTranscript(edb, filter=SeqnameFilter("Y")))
    futr <- extractTranscriptSeqs(genome, fiveUTRsByTranscript(edb, filter=SeqnameFilter("Y")))
    theTx <- "ENST00000602770"
    fullSeq <- as.character(tx[theTx])
    ## build the one from 5', cds and 3'
    compSeq <- ""
    if(any(names(futr) == theTx))
        compSeq <- paste0(compSeq, as.character(futr[theTx]))
    if(any(names(cdsSeq) == theTx))
        compSeq <- paste0(compSeq, as.character(cdsSeq[theTx]))
    if(any(names(tutr) == theTx))
        compSeq <- paste(compSeq, as.character(tutr[theTx]))
    checkEquals(unname(fullSeq), compSeq)
}

notrun_test_cds <- function(){
    library(TxDb.Hsapiens.UCSC.hg19.knownGene)
    txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
    cds <- cds(txdb)
    cby <- cdsBy(txdb, by="tx")

    gr <- cby[[7]][1]
    seqlevels(gr) <- sub(seqlevels(gr), pattern="chr", replacement="")
    tx <- transcripts(edb, filter=GRangesFilter(gr, condition="overlapping"))
    cby[[7]]

    ## Note: so that fits! And we've to include the stop_codon feature for GTF import!
    ## Make an TxDb from GTF:
    gtf <- "/Users/jo/Projects/EnsDbs/75/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz"
    library(GenomicFeatures)
    Test <- makeTxDbFromGFF(gtf, format="gtf", organism="Homo sapiens")
    scds <- cdsBy(Test, by="tx")
    gr <- scds[[7]][1]
    tx <- transcripts(edb, filter=GRangesFilter(gr, condition="overlapping"))
    scds[[7]]
    ## Compare:
    ## TxDb form GTF has: 865692 879533
    ## EnsDb: 865692 879533

    ## Next test:
    gr <- scds[[2]][1]
    tx <- transcripts(edb, filter=GRangesFilter(gr, condition="overlapping"))
    tx
    scds[[2]]
    ## start_codon: 367659 367661, stop_codon: 368595 368597 CDS: 367659 368594.
    ## TxDb from GTF includes the stop_codon!
}

