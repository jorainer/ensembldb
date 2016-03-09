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

test_return_columns_gene <- function(){
    cols <- c("gene_name", "seq_name", "tx_id")
    Resu <- genes(DB, filter=SeqnameFilter("Y"), columns=cols, return.type="data.frame")
    checkEquals(cols, colnames(Resu))
    Resu <- genes(DB, filter=SeqnameFilter("Y"), columns=cols, return.type="DataFrame")
    checkEquals(cols, colnames(Resu))
}
test_return_columns_tx <- function(){
    cols <- c("tx_id", "exon_id", "tx_biotype")
    Resu <- transcripts(DB, filter=SeqnameFilter("Y"), columns=cols, return.type="data.frame")
    checkEquals(cols, colnames(Resu))
    Resu <- transcripts(DB, filter=SeqnameFilter("Y"), columns=cols, return.type="DataFrame")
    checkEquals(cols, colnames(Resu))
}
test_return_columns_exon <- function(){
    cols <- c("tx_id", "exon_id", "tx_biotype", "seq_name")
    Resu <- exons(DB, filter=SeqnameFilter("Y"), columns=cols, return.type="data.frame")
    checkEquals(cols, colnames(Resu))
    Resu <- exons(DB, filter=SeqnameFilter("Y"), columns=cols, return.type="DataFrame")
    checkEquals(cols, colnames(Resu))
}

test_cdsBy <- function(){
    do.plot <- FALSE
    ## By tx
    cs <- cdsBy(DB, filter=list(SeqnameFilter("Y"), SeqstrandFilter("+")))
    tx <- exonsBy(DB, filter=list(SeqnameFilter("Y"), SeqstrandFilter("+")))
    ## Check for the first if it makes sense:
    whichTx <- names(cs)[1]
    whichCs <- cs[[1]]
    tx <- transcripts(DB, filter=TxidFilter(whichTx),
                      columns=c("tx_seq_start", "tx_seq_end", "tx_cds_seq_start",
                                "tx_cds_seq_end", "exon_seq_start", "exon_seq_end",
                                "exon_idx", "exon_id", "seq_strand"),
                      return.type="data.frame")
    checkSingleTx(tx=tx, cds=whichCs, do.plot=do.plot)
    ## Next one:
    whichTx <- names(cs)[2]
    tx <- transcripts(DB, filter=TxidFilter(whichTx),
                      columns=c("tx_seq_start", "tx_seq_end", "tx_cds_seq_start",
                                "tx_cds_seq_end", "exon_seq_start", "exon_seq_end",
                                "exon_idx", "exon_id"), return.type="data.frame")
    checkSingleTx(tx=tx, cds=cs[[2]], do.plot=do.plot)

    ## Now for reverse strand:
    cs <- cdsBy(DB, filter=list(SeqnameFilter("Y"), SeqstrandFilter("-")))
    whichTx <- names(cs)[1]
    whichCs <- cs[[1]]
    tx <- transcripts(DB, filter=TxidFilter(whichTx),
                      columns=c("tx_seq_start", "tx_seq_end", "tx_cds_seq_start",
                                "tx_cds_seq_end", "exon_seq_start", "exon_seq_end",
                                "exon_idx", "exon_id"), return.type="data.frame")
    ## order the guys by seq_start
    whichCs <- whichCs[order(start(whichCs))]
    checkSingleTx(tx=tx, cds=whichCs, do.plot=do.plot)
    ## Next one:
    whichTx <- names(cs)[2]
    whichCs <- cs[[2]]
    tx <- transcripts(DB, filter=TxidFilter(whichTx),
                      columns=c("tx_seq_start", "tx_seq_end", "tx_cds_seq_start",
                                "tx_cds_seq_end", "exon_seq_start", "exon_seq_end",
                                "exon_idx", "exon_id"), return.type="data.frame")
    ## order the guys by seq_start
    whichCs <- whichCs[order(start(whichCs))]
    checkSingleTx(tx=tx, cds=whichCs, do.plot=do.plot)

    ## Check adding columns
    Test <- cdsBy(DB, filter=list(SeqnameFilter("Y")),
                  columns=c("gene_biotype", "gene_name"))
}

test_cdsByGene <- function(){
    do.plot <- FALSE
    ## By gene.
    cs <- cdsBy(DB, filter=list(SeqnameFilter("Y"), SeqstrandFilter("+")),
                by="gene", columns=NULL)
    checkSingleGene(cs[[1]], gene=names(cs)[[1]], do.plot=do.plot)
    checkSingleGene(cs[[2]], gene=names(cs)[[2]], do.plot=do.plot)
    ## - strand
    cs <- cdsBy(DB, filter=list(SeqnameFilter("Y"), SeqstrandFilter("-")),
                by="gene", columns=NULL)
    checkSingleGene(cs[[1]], gene=names(cs)[[1]], do.plot=do.plot)
    checkSingleGene(cs[[2]], gene=names(cs)[[2]], do.plot=do.plot)

    ## looks good!
    cs2 <- cdsBy(DB, filter=list(SeqnameFilter("Y"), SeqstrandFilter("+")),
                by="gene", use.names=TRUE)
}

test_UTRs <- function(){
    do.plot <- FALSE
    fUTRs <- fiveUTRsByTranscript(DB, filter=list(SeqnameFilter("Y"), SeqstrandFilter("+")))
    tUTRs <- threeUTRsByTranscript(DB, filter=list(SeqnameFilter("Y"), SeqstrandFilter("+")))
    cds <- cdsBy(DB, "tx", filter=list(SeqnameFilter("Y"), SeqstrandFilter("+")))
    ## Check a TX:
    tx <- names(fUTRs)[1]
    checkGeneUTRs(fUTRs[[tx]], tUTRs[[tx]], cds[[tx]], tx=tx, do.plot=do.plot)
    tx <- names(fUTRs)[2]
    checkGeneUTRs(fUTRs[[tx]], tUTRs[[tx]], cds[[tx]], tx=tx, do.plot=do.plot)
    tx <- names(fUTRs)[3]
    checkGeneUTRs(fUTRs[[tx]], tUTRs[[tx]], cds[[tx]], tx=tx, do.plot=do.plot)

    ## Reverse strand
    fUTRs <- fiveUTRsByTranscript(DB, filter=list(SeqnameFilter("Y"), SeqstrandFilter("-")))
    tUTRs <- threeUTRsByTranscript(DB, filter=list(SeqnameFilter("Y"), SeqstrandFilter("-")))
    cds <- cdsBy(DB, "tx", filter=list(SeqnameFilter("Y"), SeqstrandFilter("-")))
    ## Check a TX:
    tx <- names(fUTRs)[1]
    checkGeneUTRs(fUTRs[[tx]], tUTRs[[tx]], cds[[tx]], tx=tx, do.plot=do.plot)
    tx <- names(fUTRs)[2]
    checkGeneUTRs(fUTRs[[tx]], tUTRs[[tx]], cds[[tx]], tx=tx, do.plot=do.plot)
    tx <- names(fUTRs)[3]
    checkGeneUTRs(fUTRs[[tx]], tUTRs[[tx]], cds[[tx]], tx=tx, do.plot=do.plot)
}

checkGeneUTRs <- function(f, t, c, tx, do.plot=FALSE){
    if(any(strand(c) == "+")){
        ## End of five UTR has to be smaller than any start of cds
        checkTrue(max(end(f)) < min(start(c)))
        ## 3'
        checkTrue(min(start(t)) > max(end(c)))
    }else{
        ## 5'
        checkTrue(min(start(f)) > max(end(c)))
        ## 3'
        checkTrue(max(end(t)) < min(start(c)))
    }
    ## just plot...
    if(do.plot){
        tx <- transcripts(DB, filter=TxidFilter(tx), columns=c("exon_seq_start", "exon_seq_end"),
                          return.type="data.frame")
        XL <- range(c(start(f), start(c), start(t), end(f), end(c), end(t)))
        YL <- c(0, 4)
        plot(4, 4, pch=NA, xlim=XL, ylim=YL, yaxt="n", ylab="", xlab="")
        ## five UTR
        rect(xleft=start(f), xright=end(f), ybottom=0.1, ytop=0.9, col="blue")
        ## cds
        rect(xleft=start(c), xright=end(c), ybottom=1.1, ytop=1.9)
        ## three UTR
        rect(xleft=start(t), xright=end(t), ybottom=2.1, ytop=2.9, col="red")
        ## all exons
        rect(xleft=tx$exon_seq_start, xright=tx$exon_seq_end, ybottom=3.1, ytop=3.9)
    }
}

checkSingleGene <- function(whichCs, gene, do.plot=FALSE){
    tx <- transcripts(DB, filter=GeneidFilter(gene),
                      columns=c("tx_seq_start", "tx_seq_end", "tx_cds_seq_start", "tx_cds_seq_end", "tx_id",
                                "exon_id", "exon_seq_start", "exon_seq_end"), return.type="data.frame")
    XL <- range(tx[, c("tx_seq_start", "tx_seq_end")])
    tx <- split(tx, f=tx$tx_id)
    if(do.plot){
        ##XL <- range(c(start(whichCs), end(whichCs)))
        YL <- c(0, length(tx) + 1)
        plot(4, 4, pch=NA, xlim=XL, ylim=YL, yaxt="n", ylab="", xlab="")
        ## plot the txses
        for(i in 1:length(tx)){
            current <- tx[[i]]
            rect(xleft=current$exon_seq_start, xright=current$exon_seq_end,
                 ybottom=rep((i-1+0.1), nrow(current)), ytop=rep((i-0.1), nrow(current)))
            ## coding:
            rect(xleft=current$tx_cds_seq_start, xright=current$tx_cds_seq_end,
                 ybottom=rep((i-1+0.1), nrow(current)), ytop=rep((i-0.1), nrow(current)),
                 border="blue")
        }
        rect(xleft=start(whichCs), xright=end(whichCs), ybottom=rep(length(tx)+0.1, length(whichCs)),
             ytop=rep(length(tx)+0.9, length(whichCs)), border="red")
    }
}

checkSingleTx <- function(tx, cds, do.plot=FALSE){
    rownames(tx) <- tx$exon_id
    tx <- tx[cds$exon_id, ]
    ## cds start and end have to be within the correct range.
    checkTrue(all(start(cds) >= min(tx$tx_cds_seq_start)))
    checkTrue(all(end(cds) <= max(tx$tx_cds_seq_end)))
    ## For all except the first and the last we have to assume that exon_seq_start
    ## is equal to start of cds.
    checkTrue(all(start(cds)[-1] == tx$exon_seq_start[-1]))
    checkTrue(all(end(cds)[-nrow(tx)] == tx$exon_seq_end[-nrow(tx)]))
    ## just plotting the stuff...
    if(do.plot){
        XL <- range(tx[, c("exon_seq_start", "exon_seq_end")])
        YL <- c(0, 4)
        plot(3, 3, pch=NA, xlim=XL, ylim=YL, xlab="", yaxt="n", ylab="")
        ## plotting the "real" exons:
        rect(xleft=tx$exon_seq_start, xright=tx$exon_seq_end, ybottom=rep(0, nrow(tx)),
             ytop=rep(1, nrow(tx)))
        ## plotting the cds:
        rect(xleft=start(cds), xright=end(cds), ybottom=rep(1.2, nrow(tx)),
             ytop=rep(2.2, nrow(tx)), col="blue")
    }
}


##*****************************************************************
## Gviz stuff
notrun_test_genetrack_df <- function(){
    do.plot <- FALSE
    if(do.plot){
        library(Gviz)
        options(ucscChromosomeNames=FALSE)
        data(geneModels)
        geneModels$chromosome <- 7
        chr <- 7
        start <- min(geneModels$start)
        end <- max(geneModels$end)
        myGeneModels <- getGeneRegionTrackForGviz(DB, chromosome=chr, start=start,
                                                  end=end)
        ## chromosome has to be the same....
        gtrack <- GenomeAxisTrack()
        gvizTrack <- GeneRegionTrack(geneModels, name="Gviz")
        ensdbTrack <- GeneRegionTrack(myGeneModels, name="ensdb")
        plotTracks(list(gtrack, gvizTrack, ensdbTrack))
        plotTracks(list(gtrack, gvizTrack, ensdbTrack), from=26700000, to=26780000)
        ## Looks very nice...
    }
    ## Put the stuff below into the vignette:
    ## Next we get all lincRNAs on chromosome Y
    Lncs <- getGeneRegionTrackForGviz(DB,
                                      filter=list(SeqnameFilter("Y"),
                                                  GenebiotypeFilter("lincRNA")))
    Prots <- getGeneRegionTrackForGviz(DB,
                                       filter=list(SeqnameFilter("Y"),
                                                   GenebiotypeFilter("protein_coding")))
    if(do.plot){
        plotTracks(list(gtrack, GeneRegionTrack(Lncs, name="lincRNAs"),
                        GeneRegionTrack(Prots, name="proteins")))
        plotTracks(list(gtrack, GeneRegionTrack(Lncs, name="lincRNAs"),
                        GeneRegionTrack(Prots, name="proteins")),
                   from=5000000, to=7000000, transcriptAnnotation="symbol")
    }
    ## is that the same than:
    TestL <- getGeneRegionTrackForGviz(DB,
                                      filter=list(GenebiotypeFilter("lincRNA")),
                                      chromosome="Y", start=5000000, end=7000000)
    TestP <- getGeneRegionTrackForGviz(DB,
                                      filter=list(GenebiotypeFilter("protein_coding")),
                                      chromosome="Y", start=5000000, end=7000000)
    if(do.plot){
        plotTracks(list(gtrack, GeneRegionTrack(Lncs, name="lincRNAs"),
                        GeneRegionTrack(Prots, name="proteins"),
                        GeneRegionTrack(TestL, name="compareL"),
                        GeneRegionTrack(TestP, name="compareP")),
                   from=5000000, to=7000000, transcriptAnnotation="symbol")
    }
    checkTrue(all(TestL$exon %in% Lncs$exon))
    checkTrue(all(TestP$exon %in% Prots$exon))
    ## Crazy amazing stuff
    ## system.time(
    ##     All <- getGeneRegionTrackForGviz(DB)
    ## )
}

####============================================================
##  length stuff
##
####------------------------------------------------------------
test_lengthOf <- function(){
    system.time(
        lenY <- lengthOf(DB, "tx", filter=SeqnameFilter("Y"))
    )
    ## Check what would happen if we do it ourselfs...
    system.time(
        lenY2 <- sum(width(reduce(exonsBy(DB, "tx", filter=SeqnameFilter("Y")))))
    )
    checkEquals(lenY, lenY2)
    ## Same for genes.
    system.time(
        lenY <- lengthOf(DB, "gene", filter=SeqnameFilter("Y"))
    )
    ## Check what would happen if we do it ourselfs...
    system.time(
        lenY2 <- sum(width(reduce(exonsBy(DB, "gene", filter=SeqnameFilter("Y")))))
    )
    checkEquals(lenY, lenY2)
    ## Just using the transcriptLengths


}

####============================================================
##  ExonrankFilter
##
####------------------------------------------------------------
test_ExonrankFilter <- function(){
    txs <- transcripts(DB, columns=c("exon_id", "exon_idx"),
                       filter=SeqnameFilter(c("Y")))
    txs <- txs[order(names(txs))]

    txs2 <- transcripts(DB, columns=c("exon_id"),
                        filter=list(SeqnameFilter(c("Y")),
                                    ExonrankFilter(3)))
    txs2 <- txs[order(names(txs2))]
    ## hm, that's weird somehow.
    exns <- exons(DB, columns=c("tx_id", "exon_idx"),
                  filter=list(SeqnameFilter("Y"),
                              ExonrankFilter(3)))
    checkTrue(all(exns$exon_idx == 3))
    exns <- exons(DB, columns=c("tx_id", "exon_idx"),
                  filter=list(SeqnameFilter("Y"),
                              ExonrankFilter(3, condition="<")))
    checkTrue(all(exns$exon_idx < 3))
}


notrun_lengthOf <- function(){
    ## How does TxDb do that?s
    library(TxDb.Hsapiens.UCSC.hg19.knownGene)
    txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
    Test <- transcriptLengths(txdb)
    head(Test)
}




