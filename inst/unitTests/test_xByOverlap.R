####============================================================
##  tests for exonsByOverlaps, transcriptsByOverlaps
##
####------------------------------------------------------------
library(EnsDb.Hsapiens.v75)
edb <- EnsDb.Hsapiens.v75

test_transcriptsByOverlaps <- function(){
    ir2 <- IRanges(start=c(2654890, 2709520, 28111770),
                   end=c(2654900, 2709550, 28111790))
    gr2 <- GRanges(rep("Y", length(ir2)), ir2)
    grf2 <- GRangesFilter(gr2, condition="overlapping")
    Test <- transcripts(edb, filter=grf2)

    Test2 <- transcriptsByOverlaps(edb, gr2)
    checkEquals(names(Test), names(Test2))

    ## on one strand.
    gr2 <- GRanges(rep("Y", length(ir2)), ir2, strand=rep("-", length(ir2)))
    grf2 <- GRangesFilter(gr2, condition="overlapping")
    Test <- transcripts(edb, filter=grf2)
    Test2 <- transcriptsByOverlaps(edb, gr2)
    checkEquals(names(Test), names(Test2))

    ## Combine with filter...
    gr2 <- GRanges(rep("Y", length(ir2)), ir2)
    Test3 <- transcriptsByOverlaps(edb, gr2, filter=SeqstrandFilter("-"))
    checkEquals(names(Test), names(Test3))
}

test_exonsByOverlaps <- function(){
    ir2 <- IRanges(start=c(2654890, 2709520, 28111770),
                   end=c(2654900, 2709550, 28111790))
    gr2 <- GRanges(rep("Y", length(ir2)), ir2)
    grf2 <- GRangesFilter(gr2, condition="overlapping")
    Test <- exons(edb, filter=grf2)

    Test2 <- exonsByOverlaps(edb, gr2)
    checkEquals(names(Test), names(Test2))

    ## on one strand.
    gr2 <- GRanges(rep("Y", length(ir2)), ir2, strand=rep("-", length(ir2)))
    grf2 <- GRangesFilter(gr2, condition="overlapping")
    Test <- exons(edb, filter=grf2)
    Test2 <- exonsByOverlaps(edb, gr2)
    checkEquals(names(Test), names(Test2))

    ## Combine with filter...
    gr2 <- GRanges(rep("Y", length(ir2)), ir2)
    Test3 <- exonsByOverlaps(edb, gr2, filter=SeqstrandFilter("-"))
    checkEquals(names(Test), names(Test3))
}


testing_txByOverlap <- function(){
    ## Apparently, a combination between transcripts and findoverlaps.
    grf <- GRangesFilter(GRanges(seqname="Y", IRanges(start=2655145, end=2655500)),
                         condition="overlapping")
    grf2 <- GRangesFilter(GRanges(seqname="Y", IRanges(start=28740998, end=28741998)),
                          condition="overlapping")
    transcripts(edb, filter=list(SeqnameFilter("Y"), GenebiotypeFilter("protein_coding")))
    where(grf)
    con <- dbconn(edb)
    library(RSQLite)
    q <- paste0("select * from gene where (", where(grf, edb),
                ") or (", where(grf2), ")")
    Test <- dbGetQuery(con, q)

    ## Here we go...
    ir <- IRanges(start=c(142999, 231380, 27635900),
                  end=c(143300, 231800, 27636200))
    gr <- GRanges(seqname=rep("Y", length(ir)), ir)
    grf <- GRangesFilter(gr, condition="overlapping")
    where(grf)
    where(grf, edb)
    Test <- transcripts(edb, filter=grf)
    ## ?? Nothing ??
    ir2 <- IRanges(start=c(2654890, 2709520, 28111770),
                   end=c(2654900, 2709550, 28111790))
    grf2 <- GRangesFilter(GRanges(rep("Y", length(ir2)), ir2), condition="overlapping")
    Test <- transcripts(edb, filter=grf2)
    checkEquals(names(Test), c("ENST00000383070", "ENST00000250784", "ENST00000598545"))
    ## ## TxDb...
    ## library(TxDb.Hsapiens.UCSC.hg19.knownGene)
    ## txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
    ## gr <- GRanges(seqname=c("chrY", "chrY", "chrY", "chrY"),
    ##               IRanges(start=c(2655145, 28740998, 2709990, 28111770),
    ##                       end=c(2655200, 28741998, 2709999, 28112800)))
    ## transcriptsByOverlaps(txdb, GRanges(seqname=rep("chrY", length(ir)), ir))
    ## transcriptsByOverlaps(txdb, GRanges(seqname=rep("chrY", length(ir2)), ir2))

}

notrun_txdb <- function(){
    txdb <- loadDb(system.file("extdata", "hg19_knownGene_sample.sqlite",
                               package="GenomicFeatures"))
    gr <- GRanges(seqnames = rep("chr1",2),
                  ranges = IRanges(start=c(500,10500), end=c(10000,30000)),
                  strand = strand(rep("-",2)))
    transcriptsByOverlaps(txdb, gr)
}

