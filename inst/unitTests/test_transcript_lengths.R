####============================================================
##  Tests related to transcript/feature length calculations.
##
##
####------------------------------------------------------------
## Loading data and stuff
library(EnsDb.Hsapiens.v75)
edb <- EnsDb.Hsapiens.v75

test_transcriptLengths <- function(){

    ## With filter.
    daFilt <- SeqnameFilter("Y")
    allTxY <- transcripts(edb, filter=daFilt)
    txLenY <- transcriptLengths(edb, filter=daFilt)
    checkEquals(names(allTxY), rownames(txLenY))

    ## Check if lengths are OK:
    txLenY2 <- lengthOf(edb, "tx", filter=daFilt)
    checkEquals(unname(txLenY2[rownames(txLenY)]), txLenY$tx_len)

    ## Include the cds, 3' and 5' UTR
    txLenY <- transcriptLengths(edb, with.cds_len = TRUE, with.utr5_len = TRUE,
                                with.utr3_len = TRUE,
                                filter=daFilt)
    ## sum of 5' CDS and 3' has to match tx_len:
    txLen <- rowSums(txLenY[, c("cds_len", "utr5_len", "utr3_len")])
    checkEquals(txLenY[!is.na(txLen), "tx_len"], unname(txLen[!is.na(txLen)]))
    ## just to be sure...
    checkEquals(txLenY[!is.na(txLenY$utr3_len), "tx_len"],
                unname(txLen[!is.na(txLenY$utr3_len)]))
    ## Seems to be OK.

    ## Next check the 5' UTR lengths: that also verifies the fiveUTR call.
    futr <- fiveUTRsByTranscript(edb, filter=daFilt)
    futrLen <- sum(width(futr))
    checkEquals(unname(futrLen), txLenY[names(futrLen), "utr5_len"])
    ## 3'
    tutr <- threeUTRsByTranscript(edb, filter=daFilt)
    tutrLen <- sum(width(tutr))
    checkEquals(unname(tutrLen), txLenY[names(tutrLen), "utr3_len"])
}

notrun_compare_full <- function(){
    ## That's on the full thing.
    ## Test if the result has the same ordering than the transcripts call.
    allTx <- transcripts(edb)
    txLen <- transcriptLengths(edb, with.cds_len=TRUE, with.utr5_len=TRUE,
                               with.utr3_len=TRUE)
    checkEquals(names(allTx), rownames(txLen))
    system.time(
        futr <- fiveUTRsByTranscript(edb)
    )
    ## 23 secs.
    futrLen <- sum(width(futr))  ## do I need reduce???
    checkEquals(unname(futrLen), txLen[names(futrLen), "utr5_len"])
    ## 3'
    system.time(
        tutr <- threeUTRsByTranscript(edb)
    )
    system.time(
        tutrLen <- sum(width(tutr))
    )
    checkEquals(unname(tutrLen), txLen[names(tutrLen), "utr3_len"])
}

notrun_compare_to_genfeat <- function(){
    library(TxDb.Hsapiens.UCSC.hg19.knownGene)
    txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

    system.time(
        Len <- transcriptLengths(edb)
    )
    ## Woa, 52 sec
    system.time(
        txLen <- lengthOf(edb, "tx")
    )
    ## Faster, 31 sec
    checkEquals(Len$tx_len, unname(txLen[rownames(Len)]))
    system.time(
        Len2 <- transcriptLengths(txdb)
    )
    ## :) 2.5 sec.
    ## Next.
    system.time(
        Len <- transcriptLengths(edb, with.cds_len = TRUE)
    )
    ## 56 sec
    system.time(
        Len2 <- transcriptLengths(txdb, with.cds_len=TRUE)
    )
    ## 4 sec.
}


