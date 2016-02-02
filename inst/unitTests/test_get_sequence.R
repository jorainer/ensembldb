library(EnsDb.Hsapiens.v75)
edb <- EnsDb.Hsapiens.v75

## That's now using the BSGenome package...
test_extractTranscriptSeqs_with_BSGenome <- function(){
    library(BSgenome.Hsapiens.UCSC.hg19)
    bsg <- BSgenome.Hsapiens.UCSC.hg19

    ## Changing the seqlevels tyle to UCSC
    seqlevelsStyle(edb) <- "UCSC"
    ZBTB <- extractTranscriptSeqs(bsg, edb, filter=GenenameFilter("ZBTB16"))
    ## Load the sequences for one ZBTB16 transcript from FA.
    faf <- system.file("txt/ENST00000335953.fa.gz", package="ensembldb")
    Seqs <- readDNAStringSet(faf)
    tx <- "ENST00000335953"
    ## cDNA
    checkEquals(unname(as.character(ZBTB[tx])),
                unname(as.character(Seqs[grep(names(Seqs), pattern="cdna")])))
    ## CDS
    cBy <- cdsBy(edb, "tx", filter=TxidFilter(tx))
    CDS <- extractTranscriptSeqs(bsg, cBy)
    checkEquals(unname(as.character(CDS)),
                unname(as.character(Seqs[grep(names(Seqs), pattern="cds")])))
    ## 5' UTR
    fBy <- fiveUTRsByTranscript(edb, filter=TxidFilter(tx))
    UTR <- extractTranscriptSeqs(bsg, fBy)
    checkEquals(unname(as.character(UTR)),
                unname(as.character(Seqs[grep(names(Seqs), pattern="utr5")])))
    ## 3' UTR
    tBy <- threeUTRsByTranscript(edb, filter=TxidFilter(tx))
    UTR <- extractTranscriptSeqs(bsg, tBy)
    checkEquals(unname(as.character(UTR)),
                unname(as.character(Seqs[grep(names(Seqs), pattern="utr3")])))


    ## Another gene on the reverse strand:
    faf <- system.file("txt/ENST00000200135.fa.gz", package="ensembldb")
    Seqs <- readDNAStringSet(faf)
    tx <- "ENST00000200135"
    ## cDNA
    cDNA <- extractTranscriptSeqs(bsg, edb, filter=TxidFilter(tx))
    checkEquals(unname(as.character(cDNA)),
                unname(as.character(Seqs[grep(names(Seqs), pattern="cdna")])))
    ## do the same, but from other strand
    exns <- exonsBy(edb, "tx", filter=TxidFilter(tx))
    cDNA <- extractTranscriptSeqs(bsg, exns)
    checkEquals(unname(as.character(cDNA)),
                unname(as.character(Seqs[grep(names(Seqs), pattern="cdna")])))
    strand(exns) <- "+"
    cDNA <- extractTranscriptSeqs(bsg, exns)
    checkTrue(unname(as.character(cDNA)) !=
              unname(as.character(Seqs[grep(names(Seqs), pattern="cdna")])))
    ## CDS
    cBy <- cdsBy(edb, "tx", filter=TxidFilter(tx))
    CDS <- extractTranscriptSeqs(bsg, cBy)
    checkEquals(unname(as.character(CDS)),
                unname(as.character(Seqs[grep(names(Seqs), pattern="cds")])))
    ## 5' UTR
    fBy <- fiveUTRsByTranscript(edb, filter=TxidFilter(tx))
    UTR <- extractTranscriptSeqs(bsg, fBy)
    checkEquals(unname(as.character(UTR)),
                unname(as.character(Seqs[grep(names(Seqs), pattern="utr5")])))
    ## 3' UTR
    tBy <- threeUTRsByTranscript(edb, filter=TxidFilter(tx))
    UTR <- extractTranscriptSeqs(bsg, tBy)
    checkEquals(unname(as.character(UTR)),
                unname(as.character(Seqs[grep(names(Seqs), pattern="utr3")])))
}


notrun_test_extractTranscriptSeqs <- function(){
    ## Note: we can't run that by default as we can not assume everybody has
    ## AnnotationHub and the required ressource installed.
    ## That's how we want to test the transcript seqs.
    genome <- getGenomeFaFile(edb)
    ZBTB <- extractTranscriptSeqs(genome, edb, filter=GenenameFilter("ZBTB16"))
    ## Load the sequences for one ZBTB16 transcript from FA.
    faf <- system.file("txt/ENST00000335953.fa.gz", package="ensembldb")
    Seqs <- readDNAStringSet(faf)
    tx <- "ENST00000335953"
    ## cDNA
    checkEquals(unname(as.character(ZBTB[tx])),
                unname(as.character(Seqs[grep(names(Seqs), pattern="cdna")])))
    ## CDS
    cBy <- cdsBy(edb, "tx", filter=TxidFilter(tx))
    CDS <- extractTranscriptSeqs(genome, cBy)
    checkEquals(unname(as.character(CDS)),
                unname(as.character(Seqs[grep(names(Seqs), pattern="cds")])))
    ## 5' UTR
    fBy <- fiveUTRsByTranscript(edb, filter=TxidFilter(tx))
    UTR <- extractTranscriptSeqs(genome, fBy)
    checkEquals(unname(as.character(UTR)),
                unname(as.character(Seqs[grep(names(Seqs), pattern="utr5")])))
    ## 3' UTR
    tBy <- threeUTRsByTranscript(edb, filter=TxidFilter(tx))
    UTR <- extractTranscriptSeqs(genome, tBy)
    checkEquals(unname(as.character(UTR)),
                unname(as.character(Seqs[grep(names(Seqs), pattern="utr3")])))


    ## Another gene on the reverse strand:
    faf <- system.file("txt/ENST00000200135.fa.gz", package="ensembldb")
    Seqs <- readDNAStringSet(faf)
    tx <- "ENST00000200135"
    ## cDNA
    cDNA <- extractTranscriptSeqs(genome, edb, filter=TxidFilter(tx))
    checkEquals(unname(as.character(cDNA)),
                unname(as.character(Seqs[grep(names(Seqs), pattern="cdna")])))
    ## do the same, but from other strand
    exns <- exonsBy(edb, "tx", filter=TxidFilter(tx))
    cDNA <- extractTranscriptSeqs(genome, exns)
    checkEquals(unname(as.character(cDNA)),
                unname(as.character(Seqs[grep(names(Seqs), pattern="cdna")])))
    strand(exns) <- "+"
    cDNA <- extractTranscriptSeqs(genome, exns)
    checkTrue(unname(as.character(cDNA)) !=
              unname(as.character(Seqs[grep(names(Seqs), pattern="cdna")])))
    ## CDS
    cBy <- cdsBy(edb, "tx", filter=TxidFilter(tx))
    CDS <- extractTranscriptSeqs(genome, cBy)
    checkEquals(unname(as.character(CDS)),
                unname(as.character(Seqs[grep(names(Seqs), pattern="cds")])))
    ## 5' UTR
    fBy <- fiveUTRsByTranscript(edb, filter=TxidFilter(tx))
    UTR <- extractTranscriptSeqs(genome, fBy)
    checkEquals(unname(as.character(UTR)),
                unname(as.character(Seqs[grep(names(Seqs), pattern="utr5")])))
    ## 3' UTR
    tBy <- threeUTRsByTranscript(edb, filter=TxidFilter(tx))
    UTR <- extractTranscriptSeqs(genome, tBy)
    checkEquals(unname(as.character(UTR)),
                unname(as.character(Seqs[grep(names(Seqs), pattern="utr3")])))
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

