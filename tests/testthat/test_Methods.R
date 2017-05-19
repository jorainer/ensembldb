
## testing genes method.
test_that("genes method works", {
    Gns <- genes(edb, filter = ~ genename == "BCL2")
    expect_identical(Gns$gene_name, "BCL2")
    Gns <- genes(edb, filter = SeqNameFilter("Y"), return.type = "DataFrame")
    expect_identical(sort(colnames(Gns)),
                     sort(unique(c(listColumns(edb, "gene"), "entrezid"))))
    Gns <- genes(edb, filter = ~ seq_name == "Y" & gene_id == "ENSG00000012817",
                 return.type = "DataFrame",
                 columns = c("gene_id", "tx_name"))
    expect_identical(colnames(Gns), c("gene_id", "tx_name", "seq_name"))
    expect_true(all(Gns$seq_name == "Y"))
    expect_true(all(Gns$gene_id == "ENSG00000012817"))
    Gns <- genes(edb,
                 filter = AnnotationFilterList(SeqNameFilter("Y"),
                                               GeneIdFilter("ENSG00000012817")),
                 columns = c("gene_id", "gene_name"))
    ## Here we don't need the seqnames in mcols!
    expect_identical(colnames(mcols(Gns)), c("gene_id", "gene_name"))
    expect_true(all(Gns$seq_name == "Y"))
    expect_true(all(Gns$gene_id == "ENSG00000012817"))

    Gns <- genes(edb, filter = ~ seq_name == "Y" | genename == "BCL2",
                 return.type = "DataFrame")
    expect_true(all(Gns$seq_name %in% c("18", "Y")))
    Gns <- genes(edb,
                 filter = AnnotationFilterList(SeqNameFilter("Y"),
                                               GenenameFilter("BCL2")),
                 return.type = "DataFrame")
    expect_true(nrow(Gns) == 0)
    Gns <- genes(edb,
                 filter = AnnotationFilterList(SeqNameFilter("Y"),
                                               GenenameFilter("BCL2"),
                                               logOp = "|"),
                 return.type = "DataFrame")
    expect_true(all(Gns$seq_name %in% c("18", "Y")))

    afl <- AnnotationFilterList(GenenameFilter(c("BCL2", "BCL2L11")),
                                SeqNameFilter(18), logOp = "&")
    afl2 <- AnnotationFilterList(SeqNameFilter("Y"), afl, logOp = "|")
    Gns <- genes(edb, filter = afl2, columns = "gene_name",
                 return.type = "DataFrame")
    expect_identical(colnames(Gns), c("gene_name", "gene_id", "seq_name"))
    expect_true(!any(Gns$gene_name == "BCL2L11"))
    expect_true(any(Gns$gene_name == "BCL2"))
    expect_true(all(Gns$seq_name %in% c("Y", "18")))
})

test_that("transcripts method works", {
    Tns <- transcripts(edb, filter = SeqNameFilter("Y"),
                       return.type = "DataFrame")
    expect_identical(sort(colnames(Tns)), sort(c(listColumns(edb, "tx"),
                                                 "seq_name")))
    Tns <- transcripts(edb, columns = c("tx_id", "tx_name"),
                       filter = list(SeqNameFilter("Y"),
                                     TxIdFilter("ENST00000435741")))
    expect_identical(sort(colnames(mcols(Tns))), sort(c("tx_id", "tx_name")))
    expect_true(all(Tns$seq_name == "Y"))
    expect_true(all(Tns$tx_id == "ENST00000435741"))
    ## Check the default ordering.
    Tns <- transcripts(edb, filter = list(TxBiotypeFilter("protein_coding"),
                                          SeqNameFilter("X")),
                       return.type = "data.frame",
                       columns = c("seq_name", listColumns(edb, "tx")))
    expect_identical(order(Tns$seq_name, method = "radix"), 1:nrow(Tns))
})

test_that("promoters works", {
    res <- promoters(edb, filter = ~ genename == "ZBTB16")
    res_2 <- transcripts(edb, filter = GenenameFilter("ZBTB16"))
    expect_identical(length(res), length(res_2))
    expect_true(all(width(res) == 2200))
})

test_that("transcriptsBy works", {
    ## Expect results on the forward strand to be ordered by tx_seq_start
    res <- transcriptsBy(edb, filter = ~ seq_name == "Y" & seq_strand == "-",
                         by = "gene")
    fw <- res[[3]]
    expect_identical(order(start(fw)), 1:length(fw))
    ## Expect results on the reverse strand to be ordered by -tx_seq_end
    res <- transcriptsBy(edb, filter = list(SeqNameFilter("Y"),
                                            SeqStrandFilter("-")), by = "gene")
    rv <- res[[3]]
    expect_identical(order(start(rv), decreasing = TRUE), 1:length(rv))
})

test_that("exons works", {
    Exns <- exons(edb, filter = SeqNameFilter("Y"), return.type = "DataFrame")
    expect_identical(sort(colnames(Exns)),
                     sort(c(listColumns(edb, "exon"), "seq_name")))
    ## Check correct ordering.
    Exns <- exons(edb, return.type = "data.frame", filter = SeqNameFilter(20:22))
    expect_identical(order(Exns$seq_name, method = "radix"), 1:nrow(Exns))
})

test_that("exonsBy works", {
    ##ExnsBy <- exonsBy(edb, filter=list(SeqNameFilter("X")), by="tx")
    ExnsBy <- exonsBy(edb, filter = list(SeqNameFilter("Y")), by = "tx",
                      columns = c("tx_name"))
    expect_identical(sort(colnames(mcols(ExnsBy[[1]]))),
                     sort(c("exon_id", "exon_rank", "tx_name")))
    suppressWarnings(
        ExnsBy <- exonsBy(edb, filter = list(SeqNameFilter("Y")), by = "tx",
                          columns = c("tx_name"), use.names = TRUE)
    )
    expect_identical(sort(colnames(mcols(ExnsBy[[1]]))),
                     sort(c("exon_id", "exon_rank", "tx_name")))
    
    ## Check what happens if we specify tx_id.
    ExnsBy <- exonsBy(edb, filter=list(SeqNameFilter("Y")), by="tx",
                      columns=c("tx_id"))
    expect_identical(sort(colnames(mcols(ExnsBy[[1]]))),
                     sort(c("exon_id", "exon_rank", "tx_id")))
    ExnsBy <- exonsBy(edb, filter=list(SeqNameFilter("Y"), SeqStrandFilter("+")),
                      by="gene")
    ## Check that ordering is on start on the forward strand.
    fw <- ExnsBy[[3]]
    expect_identical(order(start(fw)), 1:length(fw))
    ##
    ExnsBy <- exonsBy(edb, filter=list(SeqNameFilter("Y"), SeqStrandFilter("-")),
                      by="gene")
    ## Check that ordering is on start on the forward strand.
    rv <- ExnsBy[[3]]
    expect_identical(order(end(rv), decreasing = TRUE), 1:length(rv))
})

test_that("listGenebiotypes works", {
    GBT <- listGenebiotypes(edb)
    TBT <- listTxbiotypes(edb)
})

## test if we get the expected exceptions if we're not submitting
## correct filter objects
test_that("Filter errors work in methods", {
    expect_error(genes(edb, filter="d"))
    expect_error(genes(edb, filter=list(SeqNameFilter("X"), "z")))
    expect_error(transcripts(edb, filter="d"))
    expect_error(transcripts(edb, filter=list(SeqNameFilter("X"), "z")))
    expect_error(exons(edb, filter="d"))
    expect_error(exons(edb, filter=list(SeqNameFilter("X"), "z")))
    expect_error(exonsBy(edb, filter="d"))
    expect_error(exonsBy(edb, filter=list(SeqNameFilter("X"), "z")))
    expect_error(transcriptsBy(edb, filter="d"))
    expect_error(transcriptsBy(edb, filter=list(SeqNameFilter("X"), "z")))
    expect_error(transcripts(edb, filter = ~ other_filter == "b"))
})

test_that("genes returns correct columns", {
    cols <- c("gene_name", "tx_id")
    Resu <- genes(edb, filter=SeqNameFilter("Y"), columns=cols,
                  return.type = "data.frame")
    expect_identical(sort(c(cols, "seq_name", "gene_id")), sort(colnames(Resu)))

    Resu <- genes(edb, filter=SeqNameFilter("Y"), columns=cols,
                  return.type = "DataFrame")
    expect_identical(sort(c(cols, "seq_name", "gene_id")), sort(colnames(Resu)))

    Resu <- genes(edb, filter=SeqNameFilter("Y"), columns=cols)
    expect_identical(sort(c(cols, "gene_id")), sort(colnames(mcols(Resu))))
})

test_that("transcripts return correct columns", {
    cols <- c("tx_id", "exon_id", "tx_biotype")
    Resu <- transcripts(edb, filter=SeqNameFilter("Y"), columns=cols,
                        return.type = "data.frame")
    expect_identical(sort(c(cols, "seq_name")), sort(colnames(Resu)))
    Resu <- transcripts(edb, filter=SeqNameFilter("Y"), columns=cols,
                        return.type = "DataFrame")
    expect_identical(sort(c(cols, "seq_name")), sort(colnames(Resu)))
    Resu <- transcripts(edb, filter=SeqNameFilter("Y"), columns=cols)
    expect_identical(sort(cols), sort(colnames(mcols(Resu))))
})

test_that("exons returns correct columns", {
    cols <- c("tx_id", "exon_id", "tx_biotype")
    Resu <- exons(edb, filter=SeqNameFilter("Y"), columns=cols,
                  return.type = "data.frame")
    expect_identical(sort(c(cols, "seq_name")), sort(colnames(Resu)))
    Resu <- exons(edb, filter=SeqNameFilter("Y"), columns=cols,
                  return.type = "DataFrame")
    expect_identical(sort(c(cols, "seq_name")), sort(colnames(Resu)))
    Resu <- exons(edb, filter=SeqNameFilter("Y"), columns=cols)
    expect_identical(sort(cols), sort(colnames(mcols(Resu))))
})

test_that("cdsBy works", {
    checkSingleTx <- function(tx, cds, do.plot=FALSE){
        rownames(tx) <- tx$exon_id
        tx <- tx[cds$exon_id, ]
        ## cds start and end have to be within the correct range.
        expect_true(all(start(cds) >= min(tx$tx_cds_seq_start)))
        expect_true(all(end(cds) <= max(tx$tx_cds_seq_end)))
        ## For all except the first and the last we have to assume that
        ## exon_seq_start
        ## is equal to start of cds.
        expect_true(all(start(cds)[-1] == tx$exon_seq_start[-1]))
        expect_true(all(end(cds)[-nrow(tx)] == tx$exon_seq_end[-nrow(tx)]))
        ## just plotting the stuff...
        if(do.plot){
            XL <- range(tx[, c("exon_seq_start", "exon_seq_end")])
            YL <- c(0, 4)
            plot(3, 3, pch=NA, xlim=XL, ylim=YL, xlab="", yaxt="n", ylab="")
            ## plotting the "real" exons:
            rect(xleft=tx$exon_seq_start, xright=tx$exon_seq_end,
                 ybottom=rep(0, nrow(tx)),
                 ytop=rep(1, nrow(tx)))
            ## plotting the cds:
            rect(xleft=start(cds), xright=end(cds), ybottom=rep(1.2, nrow(tx)),
                 ytop=rep(2.2, nrow(tx)), col="blue")
        }
    }

    ## Just checking if we get also tx_name
    cs <- cdsBy(edb, filter = SeqNameFilter("Y"), column="tx_name")
    expect_true(any(colnames(mcols(cs[[1]])) == "tx_name"))

    do.plot <- FALSE
    ## By tx
    cs <- cdsBy(edb, filter=list(SeqNameFilter("Y"), SeqStrandFilter("+")))
    tx <- exonsBy(edb, filter=list(SeqNameFilter("Y"), SeqStrandFilter("+")))
    ## Check for the first if it makes sense:
    whichTx <- names(cs)[1]
    whichCs <- cs[[1]]
    tx <- transcripts(edb, filter=TxIdFilter(whichTx),
                      columns=c("tx_seq_start", "tx_seq_end",
                                "tx_cds_seq_start", "tx_cds_seq_end",
                                "exon_seq_start", "exon_seq_end",
                                "exon_idx", "exon_id", "seq_strand"),
                      return.type="data.frame")
    checkSingleTx(tx=tx, cds=whichCs, do.plot=do.plot)
    ## Next one:
    whichTx <- names(cs)[2]
    tx <- transcripts(edb, filter=TxIdFilter(whichTx),
                      columns=c("tx_seq_start", "tx_seq_end",
                                "tx_cds_seq_start", "tx_cds_seq_end",
                                "exon_seq_start", "exon_seq_end",
                                "exon_idx", "exon_id"),
                      return.type="data.frame")
    checkSingleTx(tx=tx, cds=cs[[2]], do.plot=do.plot)

    ## Now for reverse strand:
    cs <- cdsBy(edb, filter=list(SeqNameFilter("Y"), SeqStrandFilter("-")))
    whichTx <- names(cs)[1]
    whichCs <- cs[[1]]
    tx <- transcripts(edb, filter=TxIdFilter(whichTx),
                      columns=c("tx_seq_start", "tx_seq_end",
                                "tx_cds_seq_start", "tx_cds_seq_end",
                                "exon_seq_start", "exon_seq_end",
                                "exon_idx", "exon_id"),
                      return.type="data.frame")
    ## order the guys by seq_start
    whichCs <- whichCs[order(start(whichCs))]
    checkSingleTx(tx=tx, cds=whichCs, do.plot=do.plot)
    ## Next one:
    whichTx <- names(cs)[2]
    whichCs <- cs[[2]]
    tx <- transcripts(edb, filter=TxIdFilter(whichTx),
                      columns=c("tx_seq_start", "tx_seq_end",
                                "tx_cds_seq_start", "tx_cds_seq_end",
                                "exon_seq_start", "exon_seq_end",
                                "exon_idx", "exon_id"),
                      return.type="data.frame")
    ## order the guys by seq_start
    whichCs <- whichCs[order(start(whichCs))]
    checkSingleTx(tx=tx, cds=whichCs, do.plot=do.plot)
    ## Check adding columns
    Test <- cdsBy(edb, filter=list(SeqNameFilter("Y")),
                  columns=c("gene_biotype", "gene_name"))
})

test_that("cdsBy with gene works", {
    checkSingleGene <- function(whichCs, gene, do.plot=FALSE){
        tx <- transcripts(edb, filter=GeneIdFilter(gene),
                          columns=c("tx_seq_start", "tx_seq_end",
                                    "tx_cds_seq_start", "tx_cds_seq_end",
                                    "tx_id", "exon_id", "exon_seq_start",
                                    "exon_seq_end"),
                          return.type="data.frame")
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
                     ybottom=rep((i-1+0.1), nrow(current)),
                     ytop=rep((i-0.1), nrow(current)))
                ## coding:
                rect(xleft = current$tx_cds_seq_start,
                     xright = current$tx_cds_seq_end,
                     ybottom = rep((i-1+0.1), nrow(current)),
                     ytop = rep((i-0.1), nrow(current)),
                     border = "blue")
            }
            rect(xleft=start(whichCs), xright=end(whichCs),
                 ybottom=rep(length(tx)+0.1, length(whichCs)),
                 ytop=rep(length(tx)+0.9, length(whichCs)), border="red")
        }
    }
    do.plot <- FALSE
    ## By gene.
    cs <- cdsBy(edb, filter=list(SeqNameFilter("Y"), SeqStrandFilter("+")),
                by="gene", columns=NULL)
    checkSingleGene(cs[[1]], gene=names(cs)[[1]], do.plot=do.plot)
    checkSingleGene(cs[[2]], gene=names(cs)[[2]], do.plot=do.plot)
    ## - strand
    cs <- cdsBy(edb, filter=list(SeqNameFilter("Y"), SeqStrandFilter("-")),
                by="gene", columns=NULL)
    checkSingleGene(cs[[1]], gene=names(cs)[[1]], do.plot=do.plot)
    checkSingleGene(cs[[2]], gene=names(cs)[[2]], do.plot=do.plot)
    ## looks good!
    cs2 <- cdsBy(edb, filter=list(SeqNameFilter("Y"), SeqStrandFilter("+")),
                 by="gene", use.names=TRUE)
})

test_that("UTRs work", {
    checkGeneUTRs <- function(f, t, c, tx, do.plot=FALSE){
        if(any(strand(c) == "+")){
            ## End of five UTR has to be smaller than any start of cds
            expect_true(max(end(f)) < min(start(c)))
            ## 3'
            expect_true(min(start(t)) > max(end(c)))
        }else{
            ## 5'
            expect_true(min(start(f)) > max(end(c)))
            ## 3'
            expect_true(max(end(t)) < min(start(c)))
        }
        ## just plot...
        if(do.plot){
            tx <- transcripts(edb, filter=TxIdFilter(tx),
                              columns=c("exon_seq_start", "exon_seq_end"),
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
            rect(xleft=tx$exon_seq_start, xright=tx$exon_seq_end,
                 ybottom=3.1, ytop=3.9)
        }
    }
    ## check presence of tx_name
    fUTRs <- fiveUTRsByTranscript(edb,
                                  filter = TxIdFilter("ENST00000155093"),
                                  column = "tx_name")
    expect_true(any(colnames(mcols(fUTRs[[1]])) == "tx_name"))

    do.plot <- FALSE
    fUTRs <- fiveUTRsByTranscript(edb, filter = list(SeqNameFilter("Y"),
                                                     SeqStrandFilter("+")))
    tUTRs <- threeUTRsByTranscript(edb, filter = list(SeqNameFilter("Y"),
                                                      SeqStrandFilter("+")))
    cds <- cdsBy(edb, "tx", filter = list(SeqNameFilter("Y"),
                                          SeqStrandFilter("+")))
    ## Check a TX:
    tx <- names(fUTRs)[1]
    checkGeneUTRs(fUTRs[[tx]], tUTRs[[tx]], cds[[tx]], tx = tx,
                  do.plot = do.plot)
    tx <- names(fUTRs)[2]
    checkGeneUTRs(fUTRs[[tx]], tUTRs[[tx]], cds[[tx]], tx = tx,
                  do.plot = do.plot)
    tx <- names(fUTRs)[3]
    checkGeneUTRs(fUTRs[[tx]], tUTRs[[tx]], cds[[tx]], tx = tx,
                  do.plot = do.plot)

    ## Reverse strand
    fUTRs <- fiveUTRsByTranscript(edb, filter = list(SeqNameFilter("Y"),
                                                     SeqStrandFilter("-")))
    tUTRs <- threeUTRsByTranscript(edb, filter = list(SeqNameFilter("Y"),
                                                      SeqStrandFilter("-")))
    cds <- cdsBy(edb, "tx", filter = list(SeqNameFilter("Y"),
                                          SeqStrandFilter("-")))
    ## Check a TX:
    tx <- names(fUTRs)[1]
    checkGeneUTRs(fUTRs[[tx]], tUTRs[[tx]], cds[[tx]], tx = tx,
                  do.plot = do.plot)
    tx <- names(fUTRs)[2]
    checkGeneUTRs(fUTRs[[tx]], tUTRs[[tx]], cds[[tx]], tx = tx,
                  do.plot = do.plot)
    tx <- names(fUTRs)[3]
    checkGeneUTRs(fUTRs[[tx]], tUTRs[[tx]], cds[[tx]], tx = tx,
                  do.plot = do.plot)

    res_1 <- ensembldb:::getUTRsByTranscript(edb, what = "five",
                                             filter = TxIdFilter("ENST00000335953"))
    res_2 <- fiveUTRsByTranscript(edb, filter = TxIdFilter("ENST00000335953"))
    expect_identical(res_1, res_2)
    res_1 <- ensembldb:::getUTRsByTranscript(edb, what = "three",
                                             filter = TxIdFilter("ENST00000335953"))
    res_2 <- threeUTRsByTranscript(edb, filter = TxIdFilter("ENST00000335953"))
    expect_identical(res_1, res_2)
})

test_that("lengthOf works", {
    system.time(
        lenY <- lengthOf(edb, "tx", filter=SeqNameFilter("Y"))
    )
    ## Check what would happen if we do it ourselfs...
    system.time(
        lenY2 <- sum(width(reduce(exonsBy(edb, "tx",
                                          filter=SeqNameFilter("Y")))))
    )
    expect_identical(lenY, lenY2)
    ## Same for genes.
    system.time(
        lenY <- lengthOf(edb, "gene", filter= ~ seq_name == "Y")
    )
    ## Check what would happen if we do it ourselfs...
    system.time(
        lenY2 <- sum(width(reduce(exonsBy(edb, "gene",
                                          filter=SeqNameFilter("Y")))))
    )
    expect_identical(lenY, lenY2)
    ## Just using the transcriptLengths

    res <- ensembldb:::.transcriptLengths(edb, filter = GenenameFilter("ZBTB16"))
    res_2 <- lengthOf(edb, "tx", filter = GenenameFilter("ZBTB16"))
    expect_identical(sort(res$tx_len), unname(sort(res_2)))
    ## also cds lengths etc.
    res <- ensembldb:::.transcriptLengths(edb, filter = GenenameFilter("ZBTB16"),
                                          with.cds_len = TRUE,
                                          with.utr5_len = TRUE,
                                          with.utr3_len = TRUE)
    expect_identical(colnames(res), c("tx_id", "gene_id", "nexon", "tx_len",
                                      "cds_len", "utr5_len", "utr3_len"))
    tx <- transcripts(edb, filter = list(GenenameFilter("ZBTB16"),
                                         TxBiotypeFilter("protein_coding")))
    expect_true(all(!is.na(res[names(tx), "cds_len"])))
    expect_equal(unname(res[names(tx), "tx_len"]),
                 unname(rowSums(res[names(tx),
                                    c("utr5_len", "cds_len", "utr3_len")])))
})


####============================================================
##  ExonRankFilter
##
####------------------------------------------------------------
test_that("ExonRankFilter works with methods", {
    txs <- transcripts(edb, columns=c("exon_id", "exon_idx"),
                       filter=SeqNameFilter(c("Y")))
    txs <- txs[order(names(txs))]

    txs2 <- transcripts(edb, columns=c("exon_id"),
                        filter=list(SeqNameFilter(c("Y")),
                                    ExonRankFilter(3)))
    txs2 <- txs[order(names(txs2))]
    ## hm, that's weird somehow.
    exns <- exons(edb, columns=c("tx_id", "exon_idx"),
                  filter=list(SeqNameFilter("Y"),
                              ExonRankFilter(3)))
    expect_true(all(exns$exon_idx == 3))
    exns <- exons(edb, columns=c("tx_id", "exon_idx"),
                  filter=list(SeqNameFilter("Y"),
                              ExonRankFilter(3, condition="<")))
    expect_true(all(exns$exon_idx < 3))
})

test_that("buildQuery and getWhat works", {
    library(RSQLite)
    Q <- buildQuery(edb, columns = c("gene_name", "gene_id"))
    expect_identical(Q, "select distinct gene.gene_name,gene.gene_id from gene")

    gf <- GeneIdFilter("ENSG00000000005")
    Q <- buildQuery(edb, columns = c("gene_name", "exon_idx"),
                    filter = AnnotationFilterList(gf))
    res <- dbGetQuery(dbconn(edb), Q)
    Q_2 <- paste0("select * from gene join tx on (gene.gene_id=tx.gene_id)",
                  " join tx2exon on (tx.tx_id=tx2exon.tx_id) where",
                  " gene.gene_id = 'ENSG00000000005'")
    res_2 <- dbGetQuery(dbconn(edb), Q_2)
    expect_identical(res, unique(res_2[, colnames(res)]))
    res_3 <- ensembldb:::getWhat(edb, columns = c("gene_name", "exon_idx"),
                                 filter = AnnotationFilterList(gf))
    expect_identical(res_3, unique(res_2[, colnames(res_3)]))
})

test_that("toSaf works", {
    txs <- transcriptsBy(edb, filter = GenenameFilter("ZBTB16"))
    saf <- ensembldb:::.toSaf(txs)
    expect_identical(nrow(saf), sum(lengths(txs)))
    saf2 <- toSAF(txs)
    expect_identical(saf2, saf)
})

test_that("disjointExons works", {
    dje <- disjointExons(edb, filter = GenenameFilter("ZBTB16"))
    exns <- exons(edb, filter = GenenameFilter("ZBTB16"))
    ## Expect that dje is shorter than exns, since overlapping exon parts have
    ## been fused.
    expect_true(length(dje) < length(exns))
    dje <- disjointExons(edb, filter = GenenameFilter("ZBTB16"),
                         aggregateGenes = TRUE)
    expect_true(length(dje) < length(exns))
})

test_that("getGeneRegionTrackForGviz works", {
    res <- getGeneRegionTrackForGviz(edb, filter = GenenameFilter("ZBTB16"))
    expect_true(all(res$feature %in% c("protein_coding", "utr5", "utr3")))
    ## Do the same without a filter:
    ## LLLLL
    res2 <- getGeneRegionTrackForGviz(edb, chromosome = "11", start = 113930000,
                                      end = 113935000)
    expect_true(all(res2$symbol == "ZBTB16"))
})

test_that("filter columns are correctly added in methods", {
    filtList <- AnnotationFilterList(GenenameFilter("a"),
                                     ExonStartFilter(123),
                                     SymbolFilter("b"), TxIdFilter("c"))
    res <- ensembldb:::addFilterColumns(cols = c("a"), filter = filtList,
                                        edb = edb)
    expect_identical(res, c("a", "gene_name", "exon_seq_start", "symbol", "tx_id"))
    res <- ensembldb:::addFilterColumns(filter = filtList,
                                        edb = edb)
    expect_identical(res, c("gene_name", "exon_seq_start", "symbol", "tx_id"))
    ## New filts
    filtList <- AnnotationFilterList(GenenameFilter("a"), ExonStartFilter(123),
                                     SymbolFilter("b"), TxIdFilter("c"))
    res <- ensembldb:::addFilterColumns(cols = c("a"), filter = filtList,
                                        edb = edb)
    expect_identical(res, c("a", "gene_name", "exon_seq_start", "symbol", "tx_id"))
    res <- ensembldb:::addFilterColumns(filter = filtList,
                                        edb = edb)
    expect_identical(res, c("gene_name", "exon_seq_start", "symbol", "tx_id"))
})

test_that("supportedFilters works", {
    res <- ensembldb:::.supportedFilters(edb)
    if (!hasProteinData(edb))
        expect_equal(length(res), 19)
    else 
        expect_equal(length(res), 24)
    res <- supportedFilters(edb)
    if (!hasProteinData(edb))
        expect_equal(length(res), 19)
    else 
        expect_equal(length(res), 24)
})

## Here we check if we fetch what we expect from the database.
test_that("GRangesFilter works in queries", {
    do.plot <- FALSE
    zbtb <- genes(edb, filter = GenenameFilter("ZBTB16"))
    txs <- transcripts(edb, filter = GenenameFilter("ZBTB16"))
    ## Now use the GRangesFilter to fetch all tx
    txs2 <- transcripts(edb, filter = GRangesFilter(zbtb))
    expect_equal(txs$tx_id, txs2$tx_id)
    ## Exons:
    exs <- exons(edb, filter = GenenameFilter("ZBTB16"))
    exs2 <- exons(edb, filter = GRangesFilter(zbtb))
    expect_equal(exs$exon_id, exs2$exon_id)
    ## Now check the filter with "overlapping".
    intr <- GRanges("11", ranges = IRanges(114000000, 114000050), strand = "+")
    gns <- genes(edb, filter = GRangesFilter(intr, type = "any"))
    expect_equal(gns$gene_name, "ZBTB16")
    ##
    txs <- transcripts(edb, filter = GRangesFilter(intr, type = "any"))
    expect_equal(sort(txs$tx_id), sort(c("ENST00000335953", "ENST00000541602",
                                         "ENST00000392996", "ENST00000539918")))
    if(do.plot){
        plot(3, 3, pch=NA, xlim=c(start(zbtb), end(zbtb)),
             ylim=c(0, length(txs2)))
        rect(xleft=start(intr), xright=end(intr), ybottom=0, ytop=length(txs2),
             col="red", border="red")
        for(i in 1:length(txs2)){
            current <- txs2[i]
            rect(xleft=start(current), xright=end(current), ybottom=i-0.975,
                 ytop=i-0.125, border="grey")
            text(start(current), y=i-0.5,pos=4, cex=0.75, labels=current$tx_id)
        }
        ## OK, that' OK.
    }

    ## OK, now for a GRangesFilter with more than one GRanges.
    ir2 <- IRanges(start=c(2654890, 2709520, 28111770),
                   end=c(2654900, 2709550, 28111790))
    grf2 <- GRangesFilter(GRanges(rep("Y", length(ir2)), ir2),
                          type = "any")
    Test <- transcripts(edb, filter = grf2)
    expect_equal(names(Test), c("ENST00000383070", "ENST00000250784",
                                "ENST00000598545"))
})

test_that("show works", {
    res <- capture.output(show(edb))
    expect_equal(res[1], "EnsDb for Ensembl:")
    expect_equal(res[9], "|ensembl_version: 75")
})

test_that("organism method works", {
    res <- organism(edb)
    expect_equal(res, "Homo sapiens")
})

test_that("metadata method works", {
    res <- metadata(edb)
    expect_equal(res, dbGetQuery(dbconn(edb), "select * from metadata"))
})

test_that("ensemblVersion works", {
    expect_equal(ensemblVersion(edb), "75")
})

test_that("getMetadataValue works", {
    expect_error(ensembldb:::getMetadataValue(edb))
})

test_that("seqinfo and seqlevels work", {
    si <- seqinfo(edb)
    expect_true(is(si, "Seqinfo"))
    sl <- seqlevels(edb)
    library(RSQLite)
    chrs <- dbGetQuery(dbconn(edb), "select seq_name from chromosome")[, 1]
    expect_true(all(sl %in% chrs))
    expect_true(all(seqlevels(si) %in% chrs))
})

test_that("ensVersionFromSourceUrl works", {
    res <- ensembldb:::.ensVersionFromSourceUrl(
                           "ftp://ftp.ensembl.org/release-85/gtf")
    expect_equal(res, 85)
})

test_that("listBiotypes works", {
    res <- listTxbiotypes(edb)
    library(RSQLite)
    res_2 <- dbGetQuery(dbconn(edb), "select distinct tx_biotype from tx")[, 1]
    expect_true(all(res %in% res_2))
    res <- listGenebiotypes(edb)
    res_2 <- dbGetQuery(dbconn(edb), "select distinct gene_biotype from gene")[, 1]
    expect_true(all(res %in% res_2))
})

test_that("listTables works", {
    res <- listTables(edb)
    schema_version <- ensembldb:::dbSchemaVersion(edb)
    if (!hasProteinData(edb)) {
        expect_equal(names(res),
                     names(ensembldb:::.ensdb_tables(schema_version)))
    } else {
        expect_equal(
            sort(names(res)),
            sort(unique(c(names(ensembldb:::.ensdb_tables(schema_version)),
                          names(ensembldb:::.ensdb_protein_tables(
                                                schema_version))))))
    }
    ## Repeat with deleting the cached tables
    edb@tables <- list()
    res <- listTables(edb)
    if (!hasProteinData(edb)) {
        expect_equal(names(res),
                     names(ensembldb:::.ensdb_tables(schema_version)))
    } else {
        expect_equal(
            sort(names(res)),
            sort(unique(c(names(ensembldb:::.ensdb_tables(schema_version)),
                          names(ensembldb:::.ensdb_protein_tables(
                                                schema_version))))))
    }
})

test_that("listColumns works", {
    res <- listColumns(edb, table = "gene")
    expect_equal(res, c(ensembldb:::.ensdb_tables(ensembldb:::dbSchemaVersion(dbconn(edb)))$gene, "symbol"))
    res <- listColumns(edb, table = "tx")
    expect_equal(res, c(ensembldb:::.ensdb_tables(ensembldb:::dbSchemaVersion(dbconn(edb)))$tx, "tx_name"))
    res <- listColumns(edb, table = "exon")
    expect_equal(res, c(ensembldb:::.ensdb_tables(ensembldb:::dbSchemaVersion(dbconn(edb)))$exon))
    res <- listColumns(edb, table = "chromosome")
    expect_equal(res, c(ensembldb:::.ensdb_tables(ensembldb:::dbSchemaVersion(dbconn(edb)))$chromosome))
    res <- listColumns(edb, table = "tx2exon")
    expect_equal(res, c(ensembldb:::.ensdb_tables(ensembldb:::dbSchemaVersion(dbconn(edb)))$tx2exon))
    if (hasProteinData(edb)) {
        res <- listColumns(edb, table = "protein")
        expect_equal(res, ensembldb:::.ensdb_protein_tables(ensembldb:::dbSchemaVersion(dbconn(edb)))$protein)
        res <- listColumns(edb, table = "uniprot")
        expect_equal(res, ensembldb:::.ensdb_protein_tables(ensembldb:::dbSchemaVersion(dbconn(edb)))$uniprot)
        res <- listColumns(edb, table = "protein_domain")
        expect_equal(res, ensembldb:::.ensdb_protein_tables(ensembldb:::dbSchemaVersion(dbconn(edb)))$protein_domain)
    }
    ## Repeat with deleting the cached tables
    edb@tables <- list()
    res <- listColumns(edb, table = "gene")
    expect_equal(res, c(ensembldb:::.ensdb_tables(ensembldb:::dbSchemaVersion(dbconn(edb)))$gene, "symbol"))
    res <- listColumns(edb, table = "tx")
    expect_equal(res, c(ensembldb:::.ensdb_tables(ensembldb:::dbSchemaVersion(dbconn(edb)))$tx, "tx_name"))
    res <- listColumns(edb, table = "exon")
    expect_equal(res, c(ensembldb:::.ensdb_tables(ensembldb:::dbSchemaVersion(dbconn(edb)))$exon))
    res <- listColumns(edb, table = "chromosome")
    expect_equal(res, c(ensembldb:::.ensdb_tables(ensembldb:::dbSchemaVersion(dbconn(edb)))$chromosome))
    res <- listColumns(edb, table = "tx2exon")
    expect_equal(res, c(ensembldb:::.ensdb_tables(ensembldb:::dbSchemaVersion(dbconn(edb)))$tx2exon))
    if (hasProteinData(edb)) {
        res <- listColumns(edb, table = "protein")
        expect_equal(res, ensembldb:::.ensdb_protein_tables(ensembldb:::dbSchemaVersion(dbconn(edb)))$protein)
        res <- listColumns(edb, table = "uniprot")
        expect_equal(res, ensembldb:::.ensdb_protein_tables(ensembldb:::dbSchemaVersion(dbconn(edb)))$uniprot)
        res <- listColumns(edb, table = "protein_domain")
        expect_equal(res, ensembldb:::.ensdb_protein_tables(ensembldb:::dbSchemaVersion(dbconn(edb)))$protein_domain)
    }
})

test_that("cleanColumns works", {
    cols <- c("gene_id", "tx_id", "tx_name")
    res <- ensembldb:::cleanColumns(edb, cols)
    expect_equal(cols, res)
    cols <- c(cols, "not there")
    suppressWarnings(
        res <- ensembldb:::cleanColumns(edb, cols)
    )
    expect_equal(cols[1:3], res)
    cols <- c("gene_id", "protein_id", "tx_id", "protein_sequence")
    suppressWarnings(
        res <- ensembldb:::cleanColumns(edb, cols)
    )
    if (hasProteinData(edb)) {
        expect_equal(res, cols)
    } else {
        expect_equal(res, cols[c(1, 3)])
    }
    ## with full names:
    cols <- c("gene.gene_id", "protein.protein_id", "tx.tx_id",
              "protein.protein_sequence")
    suppressWarnings(
        res <- ensembldb:::cleanColumns(edb, cols)
    )
    if (hasProteinData(edb)) {
        expect_equal(res, cols)
    } else {
        expect_equal(res, cols[c(1, 3)])
    }
})

test_that("tablesForColumns works", {
    expect_error(ensembldb:::tablesForColumns(edb))
    res <- ensembldb:::tablesForColumns(edb, columns = "tx_id")
    if (hasProteinData(edb))
        expect_equal(res, c("tx", "tx2exon", "protein"))
    else
        expect_equal(res, c("tx", "tx2exon"))
    res <- ensembldb:::tablesForColumns(edb, columns = "seq_name")
    expect_equal(res, c("gene", "chromosome"))
    if (hasProteinData(edb)) {
        res <- ensembldb:::tablesForColumns(edb, columns = "protein_id")
        expect_equal(res, c("protein", "uniprot", "protein_domain"))
    }
})

test_that("tablesByDegree works", {
    res <- ensembldb:::tablesByDegree(edb,
                                      tab = c("chromosome", "gene", "tx"))
    expect_equal(res, c("gene", "tx", "chromosome"))
})

test_that("updateEnsDb works", {
    edb2 <- updateEnsDb(edb)
    expect_equal(edb2@tables, edb@tables)
    expect_true(.hasSlot(edb2, ".properties"))
})

test_that("properties work", {
    origProps <- ensembldb:::properties(edb)
    expect_equal(ensembldb:::getProperty(edb, "foo"), NA)
    expect_error(ensembldb:::setProperty(edb, "foo"))
    edb <- ensembldb:::setProperty(edb, foo="bar")
    expect_equal(ensembldb:::getProperty(edb, "foo"), "bar")
    expect_equal(length(ensembldb:::properties(edb)),
                 length(origProps) + 1)
    expect_true(any(names(ensembldb:::properties(edb)) == "foo"))
    edb <- ensembldb:::dropProperty(edb, "foo")
    expect_true(all(names(ensembldb:::properties(edb)) != "foo"))
})

## Compare the results for genes call with and without ordering in R
test_that("ordering works in genes calls", {
    orig <- ensembldb:::orderResultsInR(edb)
    ensembldb:::orderResultsInR(edb) <- FALSE
    res_sql <- genes(edb, return.type = "data.frame")
    ensembldb:::orderResultsInR(edb) <- TRUE
    res_r <- genes(edb, return.type = "data.frame")
    rownames(res_sql) <- NULL
    rownames(res_r) <- NULL
    expect_equal(res_sql, res_r)
    ## Join tx table
    ensembldb:::orderResultsInR(edb) <- FALSE
    res_sql <- genes(edb, columns = c("gene_id", "tx_id"),
                     return.type = "data.frame")
    ensembldb:::orderResultsInR(edb) <- TRUE
    res_r <- genes(edb, columns = c("gene_id", "tx_id"),
                   return.type = "data.frame")
    rownames(res_sql) <- NULL
    rownames(res_r) <- NULL
    expect_equal(res_sql, res_r)
    ## Join tx table and use an SeqNameFilter
    ensembldb:::orderResultsInR(edb) <- FALSE
    res_sql <- genes(edb, columns = c("gene_id", "tx_id"),
                     filter = SeqNameFilter("Y"))
    ensembldb:::orderResultsInR(edb) <- TRUE
    res_r <- genes(edb, columns = c("gene_id", "tx_id"),
                   filter = SeqNameFilter("Y"))
    expect_equal(res_sql, res_r)

    ensembldb:::orderResultsInR(edb) <- orig
})

test_that("transcriptLengths works",{
    ## With filter.
    daFilt <- SeqNameFilter("Y")
    allTxY <- transcripts(edb, filter = daFilt)
    txLenY <- transcriptLengths(edb, filter = daFilt)
    expect_equal(names(allTxY), txLenY$tx_id)
    rownames(txLenY) <- txLenY$tx_id

    ## Check if lengths are OK:
    txLenY2 <- lengthOf(edb, "tx", filter = daFilt)
    expect_equal(unname(txLenY2[txLenY$tx_id]), txLenY$tx_len)

    ## Include the cds, 3' and 5' UTR
    txLenY <- transcriptLengths(edb, with.cds_len = TRUE, with.utr5_len = TRUE,
                                with.utr3_len = TRUE,
                                filter=daFilt)
    ## sum of 5' CDS and 3' has to match tx_len:
    txLen <- rowSums(txLenY[, c("cds_len", "utr5_len", "utr3_len")])
    expect_equal(txLenY[txLenY$cds_len > 0, "tx_len"],
                 unname(txLen[txLenY$cds_len > 0]))
    ## just to be sure...
    expect_equal(txLenY[txLenY$utr3_len > 0, "tx_len"],
                unname(txLen[txLenY$utr3_len > 0]))
    ## Seems to be OK.

    ## Next check the 5' UTR lengths: that also verifies the fiveUTR call.
    futr <- fiveUTRsByTranscript(edb, filter = daFilt)
    futrLen <- sum(width(futr))
    rownames(txLenY) <- txLenY$tx_id
    expect_equal(unname(futrLen), txLenY[names(futrLen), "utr5_len"])
    ## 3'
    tutr <- threeUTRsByTranscript(edb, filter=daFilt)
    tutrLen <- sum(width(tutr))
    expect_equal(unname(tutrLen), txLenY[names(tutrLen), "utr3_len"])
})

test_that("transcriptsByOverlaps works", {
    ir2 <- IRanges(start = c(2654890, 2709520, 28111770),
                   end = c(2654900, 2709550, 28111790))
    gr2 <- GRanges(rep("Y", length(ir2)), ir2)
    grf2 <- GRangesFilter(gr2, type = "any")
    Test <- transcripts(edb, filter = grf2)
    Test2 <- transcriptsByOverlaps(edb, gr2)
    expect_equal(names(Test), names(Test2))
    ## on one strand.
    gr2 <- GRanges(rep("Y", length(ir2)), ir2, strand = rep("-", length(ir2)))
    grf2 <- GRangesFilter(gr2, type = "any")
    Test <- transcripts(edb, filter = grf2)
    Test2 <- transcriptsByOverlaps(edb, gr2)
    expect_equal(names(Test), names(Test2))

    ## Combine with filter...
    gr2 <- GRanges(rep("Y", length(ir2)), ir2)
    Test3 <- transcriptsByOverlaps(edb, gr2, filter = SeqStrandFilter("-"))
    expect_equal(names(Test), names(Test3))
})

test_that("exonsByOverlaps works", {
    ir2 <- IRanges(start=c(2654890, 2709520, 28111770),
                   end=c(2654900, 2709550, 28111790))
    gr2 <- GRanges(rep("Y", length(ir2)), ir2)
    grf2 <- GRangesFilter(gr2, type = "any")
    Test <- exons(edb, filter = grf2)
    Test2 <- exonsByOverlaps(edb, gr2)
    expect_equal(names(Test), names(Test2))
    ## on one strand.
    gr2 <- GRanges(rep("Y", length(ir2)), ir2, strand=rep("-", length(ir2)))
    grf2 <- GRangesFilter(gr2, type = "any")
    Test <- exons(edb, filter = grf2)
    Test2 <- exonsByOverlaps(edb, gr2)
    expect_equal(names(Test), names(Test2))
    ## Combine with filter...
    gr2 <- GRanges(rep("Y", length(ir2)), ir2)
    Test3 <- exonsByOverlaps(edb, gr2, filter=SeqStrandFilter("-"))
    expect_equal(names(Test), names(Test3))
})
