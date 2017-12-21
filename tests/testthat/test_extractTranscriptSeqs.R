
test_that("extractTranscriptSeqs works with BSGenome", {
    library(BSgenome.Hsapiens.NCBI.GRCh38)
    bsg <- BSgenome.Hsapiens.NCBI.GRCh38

    ## ## Changing the seqlevels tyle to UCSC
    ## seqlevelsStyle(edb) <- "UCSC"
    ZBTB <- extractTranscriptSeqs(bsg, edb, filter=GenenameFilter("ZBTB16"))
    ## Load the sequences for one ZBTB16 transcript from FA.
    faf <- system.file("txt/ENST00000335953.fa.gz", package="ensembldb")
    Seqs <- readDNAStringSet(faf)
    tx <- "ENST00000335953"
    ## cDNA
    expect_equal(unname(as.character(ZBTB[tx])),
                 unname(as.character(Seqs[grep(names(Seqs), pattern="cdna")])))
    ## CDS
    cBy <- cdsBy(edb, "tx", filter=TxIdFilter(tx))
    CDS <- extractTranscriptSeqs(bsg, cBy)
    expect_equal(unname(as.character(CDS)),
                 unname(as.character(Seqs[grep(names(Seqs), pattern="cds")])))
    ## 5' UTR
    fBy <- fiveUTRsByTranscript(edb, filter=TxIdFilter(tx))
    UTR <- extractTranscriptSeqs(bsg, fBy)
    expect_equal(unname(as.character(UTR)),
                 unname(as.character(Seqs[grep(names(Seqs), pattern="utr5")])))
    ## 3' UTR
    tBy <- threeUTRsByTranscript(edb, filter=TxIdFilter(tx))
    UTR <- extractTranscriptSeqs(bsg, tBy)
    expect_equal(unname(as.character(UTR)),
                 unname(as.character(Seqs[grep(names(Seqs), pattern="utr3")])))

    ## Another gene on the reverse strand:
    faf <- system.file("txt/ENST00000200135.fa.gz", package="ensembldb")
    Seqs <- readDNAStringSet(faf)
    tx <- "ENST00000200135"
    ## cDNA
    cDNA <- extractTranscriptSeqs(bsg, edb, filter=TxIdFilter(tx))
    expect_equal(unname(as.character(cDNA)),
                 unname(as.character(Seqs[grep(names(Seqs), pattern="cdna")])))
    ## do the same, but from other strand
    exns <- exonsBy(edb, "tx", filter=TxIdFilter(tx))
    cDNA <- extractTranscriptSeqs(bsg, exns)
    expect_equal(unname(as.character(cDNA)),
                 unname(as.character(Seqs[grep(names(Seqs), pattern="cdna")])))
    strand(exns) <- "+"
    cDNA <- extractTranscriptSeqs(bsg, exns)
    expect_true(unname(as.character(cDNA)) !=
                unname(as.character(Seqs[grep(names(Seqs), pattern="cdna")])))
    ## CDS
    cBy <- cdsBy(edb, "tx", filter=TxIdFilter(tx))
    CDS <- extractTranscriptSeqs(bsg, cBy)
    expect_equal(unname(as.character(CDS)),
                 unname(as.character(Seqs[grep(names(Seqs), pattern="cds")])))
    ## 5' UTR
    fBy <- fiveUTRsByTranscript(edb, filter=TxIdFilter(tx))
    UTR <- extractTranscriptSeqs(bsg, fBy)
    expect_equal(unname(as.character(UTR)),
                 unname(as.character(Seqs[grep(names(Seqs), pattern="utr5")])))
    ## 3' UTR
    tBy <- threeUTRsByTranscript(edb, filter=TxIdFilter(tx))
    UTR <- extractTranscriptSeqs(bsg, tBy)
    expect_equal(unname(as.character(UTR)),
                 unname(as.character(Seqs[grep(names(Seqs), pattern="utr3")])))
})

