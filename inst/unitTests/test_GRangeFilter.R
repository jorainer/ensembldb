###============================================================
##  Testing the GRangesFilter
###------------------------------------------------------------
library(EnsDb.Hsapiens.v75)
edb <- EnsDb.Hsapiens.v75

test_GRangesFilterValidity <- function(){
    checkException(GRangesFilter(value="bla"))
    checkException(GRangesFilter(GRanges(seqnames="X", ranges=IRanges(4, 6)),
                                 condition=">"))
    ## Testing slots
    gr <- GRanges("X", ranges=IRanges(123, 234), strand="-")
    grf <- GRangesFilter(gr, condition="within")
    ## Now check some stuff
    checkEquals(start(grf), start(gr))
    checkEquals(end(grf), end(gr))
    checkEquals(as.character(strand(gr)), strand(grf))
    checkEquals(as.character(seqnames(gr)), seqnames(grf))

    ## Test column:
    ## filter alone.
    tocomp <- c(start="gene_seq_start", end="gene_seq_end", seqname="seq_name",
                strand="seq_strand")
    checkEquals(column(grf), tocomp)
    grf@feature <- "tx"
    tocomp <- c(start="tx_seq_start", end="tx_seq_end", seqname="seq_name",
                strand="seq_strand")
    checkEquals(column(grf), tocomp)
    grf@feature <- "exon"
    tocomp <- c(start="exon_seq_start", end="exon_seq_end", seqname="seq_name",
                strand="seq_strand")
    checkEquals(column(grf), tocomp)
    ## filter and ensdb.
    tocomp <- c(start="exon.exon_seq_start", end="exon.exon_seq_end", seqname="gene.seq_name",
                strand="gene.seq_strand")
    checkEquals(column(grf, edb), tocomp)
    grf@feature <- "tx"
    tocomp <- c(start="tx.tx_seq_start", end="tx.tx_seq_end", seqname="gene.seq_name",
                strand="gene.seq_strand")
    checkEquals(column(grf, edb), tocomp)
    grf@feature <- "gene"
    tocomp <- c(start="gene.gene_seq_start", end="gene.gene_seq_end", seqname="gene.seq_name",
                strand="gene.seq_strand")
    checkEquals(column(grf, edb), tocomp)

    ## Test where:
    ## filter alone.
    tocomp <- "gene_seq_start >= 123 and gene_seq_end <= 234 and seq_name == 'X' and seq_strand = -1"
    checkEquals(where(grf), tocomp)
    ## what if we set strand to *
    grf2 <- GRangesFilter(GRanges("1", IRanges(123, 234)))
    tocomp <- "gene.gene_seq_start >= 123 and gene.gene_seq_end <= 234 and gene.seq_name == '1'"
    checkEquals(where(grf2, edb), tocomp)

    ## Now, using overlapping.
    grf@location <- "overlapping"
    grf@feature <- "transcript"
    tocomp <- "tx.tx_seq_start <= 234 and tx.tx_seq_end >= 123 and gene.seq_name = 'X' and gene.seq_strand = -1"
    checkEquals(where(grf, edb), tocomp)
}

## Here we check if we fetch what we expect from the database.
test_GRangesFilterQuery <- function(){
    do.plot <- FALSE
    zbtb <- genes(edb, filter=GenenameFilter("ZBTB16"))
    txs <- transcripts(edb, filter=GenenameFilter("ZBTB16"))

    ## Now use the GRangesFilter to fetch all tx
    txs2 <- transcripts(edb, filter=GRangesFilter(zbtb))
    checkEquals(txs$tx_id, txs2$tx_id)

    ## Exons:
    exs <- exons(edb, filter=GenenameFilter("ZBTB16"))
    exs2 <- exons(edb, filter=GRangesFilter(zbtb))
    checkEquals(exs$exon_id, exs2$exon_id)

    ## Now check the filter with "overlapping".
    intr <- GRanges("11", ranges=IRanges(114000000, 114000050), strand="+")
    gns <- genes(edb, filter=GRangesFilter(intr, condition="overlapping"))
    checkEquals(gns$gene_name, "ZBTB16")

    txs <- transcripts(edb, filter=GRangesFilter(intr, condition="overlapping"))
    if(do.plot){
        plot(3, 3, pch=NA, xlim=c(start(zbtb), end(zbtb)), ylim=c(0, length(txs2)))
        rect(xleft=start(intr), xright=end(intr), ybottom=0, ytop=length(txs2), col="red", border="red")
        for(i in 1:length(txs2)){
            current <- txs2[i]
            rect(xleft=start(current), xright=end(current), ybottom=i-0.975, ytop=i-0.125, border="grey")
            text(start(current), y=i-0.5,pos=4, cex=0.75, labels=current$tx_id)
        }
        ## OK, that' OK.
    }

    ## OK, now for a GRangesFilter with more than one GRanges.
    ir2 <- IRanges(start=c(2654890, 2709520, 28111770),
                   end=c(2654900, 2709550, 28111790))
    grf2 <- GRangesFilter(GRanges(rep("Y", length(ir2)), ir2), condition="overlapping")
    Test <- transcripts(edb, filter=grf2)
    checkEquals(names(Test), c("ENST00000383070", "ENST00000250784", "ENST00000598545"))

}

