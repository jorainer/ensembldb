
test_that("columns works", {
    cols <- columns(edb)
    ## Don't expect to see any _ there...
    expect_equal(length(grep(cols, pattern="_")), 0)
})

test_that("keytypes works", {
    keyt <- keytypes(edb)
    expect_equal(all(c("GENEID", "EXONID", "TXID") %in% keyt), TRUE)
})

test_that("ensDbColumnForColumn works", {
    Test <- ensembldb:::ensDbColumnForColumn(edb, "GENEID")
    expect_equal(unname(Test), "gene_id")
    Test <- ensembldb:::ensDbColumnForColumn(edb, c("GENEID", "TXID"))
    expect_equal(unname(Test), c("gene_id", "tx_id"))
    suppressWarnings(
        Test <- ensembldb:::ensDbColumnForColumn(edb, c("GENEID", "TXID", "bla"))
    )
    expect_equal(unname(Test), c("gene_id", "tx_id"))
})

test_that("keys works", {
    ## get all gene ids
    ids <- keys(edb, "GENEID")
    expect_true(length(ids) > 0)
    expect_equal(length(ids), length(unique(ids)))
    ## get all tx ids
    ids <- keys(edb, "TXID")
    expect_true(length(ids) > 0)
    ## Get the TXNAME...
    nms <- keys(edb, "TXNAME")
    expect_equal(nms, ids)
    expect_equal(length(ids), length(unique(ids)))
    ## get all gene names
    ids <- keys(edb, "GENENAME")
    expect_true(length(ids) > 0)
    expect_equal(length(ids), length(unique(ids)))
    ## get all seq names
    ids <- keys(edb, "SEQNAME")
    expect_true(length(ids) > 0)
    expect_equal(length(ids), length(unique(ids)))
    ## get all seq strands
    ids <- keys(edb, "SEQSTRAND")
    expect_true(length(ids) > 0)
    expect_equal(length(ids), length(unique(ids)))
    ## get all gene biotypes
    ids <- keys(edb, "GENEBIOTYPE")
    expect_true(length(ids) > 0)
    expect_equal(ids, listGenebiotypes(edb))
    ## Now with protein data.
    if (hasProteinData(edb)) {
        library(RSQLite)
        ls <- keys(edb, "PROTEINID")
        ls_2 <- dbGetQuery(dbconn(edb),
                           "select distinct protein_id from protein")$protein_id
        expect_equal(sort(ls), sort(ls_2))
        ##
        ks <- keys(edb, "UNIPROTID")
        ks_2 <- dbGetQuery(dbconn(edb),
                           "select distinct uniprot_id from uniprot")$uniprot_id
        expect_equal(sort(ks), sort(ks_2))
        ##
        ks <- keys(edb, "PROTEINDOMAINID")
        ks_2 <- dbGetQuery(dbconn(edb),
                           paste0("select distinct protein_domain_id from",
                                  " protein_domain"))$protein_domain_id
        expect_equal(sort(ks), sort(ks_2))
    }
    ## keys with filter:
    res <- keys(edb, "GENENAME", filter = ~ genename == "BCL2")
    expect_equal(res, "BCL2")
})

test_that("select method works", {
    .comprehensiveCheckForGene <- function(x) {
        ##   Check if we've got all of the transcripts.
        idx_rem <- grep("LRG", x$TXID)
        if (length(idx_rem))
            x <- x[-idx_rem, ]
        txs <- dbGetQuery(
            dbconn(edb),
            paste0("select tx_id from tx where gene_id = '",
                   x$GENEID[1], "';"))
        expect_equal(sort(txs$tx_id), sort(unique(x$TXID)))
        ##   Check if we've got all exons.
        exs <- dbGetQuery(
            dbconn(edb),
            paste0("select exon_id from tx2exon where tx_id in (",
                   paste0("'", txs$tx_id, "'", collapse = ", "),")"))
        a <- sort(unique(exs$exon_id))
        b <- sort(unique(x$EXONID))
        a <- a[!is.na(a)]
        b <- b[!is.na(b)]
        expect_equal(a, b)
        if (hasProteinData(edb)) {
            ##  Check if we've got all proteins
            prt <- dbGetQuery(
                dbconn(edb),
                paste0("select protein_id from protein where tx_id in (",
                       paste0("'", txs$tx_id, "'", collapse = ", "), ")"))
            a <- sort(prt$protein_id)
            b <- sort(unique(x$PROTEINID))
            a <- a[!is.na(a)]
            b <- b[!is.na(b)]
            expect_equal(a, b)
            ##  Check if we've got all uniprots.
            res <- dbGetQuery(
                dbconn(edb),
                paste0("select uniprot_id from uniprot where ",
                       "protein_id in (", paste0("'", prt$protein_id,
                                                 "'", collapse = ", ") ,")"))
            a <- sort(unique(res$uniprot_id))
            b <- sort(unique(x$UNIPROTID))
            a <- a[!is.na(a)]
            b <- b[!is.na(b)]
            expect_equal(a, b)
            ##  Check if we've got all protein_domains.
            res <- dbGetQuery(
                dbconn(edb),
                paste0("select protein_domain_id from protein_domain ",
                       "where protein_id in (",
                       paste0("'", prt$protein_id,
                              "'", collapse = ", "), ")"))
            a <- sort(unique(res$protein_domain_id))
            b <- sort(unique(x$PROTEINDOMAINID))
            a <- a[!is.na(a)]
            b <- b[!is.na(b)]
            expect_equal(a, b)
        }
    }

    library(RSQLite)
    edb <- filter(edb, filter = ~ seq_name %in% c(18, 11, 2, 5, 7, 21))
    ## 1) Test:
    ##   Provide GenenameFilter.
    gf <- GenenameFilter("BCL2")
    Test <- select(edb, keys = gf)
    expect_true(all(Test$GENENAME == "BCL2"))
    .comprehensiveCheckForGene(Test)
    Test2 <- select(edb, keys = ~ symbol == "BCL2")
    expect_equal(Test, Test2)
    ## ZBTB16
    tmp <- select(edb, keys = GenenameFilter("ZBTB16"))
    .comprehensiveCheckForGene(tmp)
    ## BCL2L11
    tmp <- select(edb, keys = GenenameFilter("BCL2L11"))
    .comprehensiveCheckForGene(tmp)
    ## NR3C1
    tmp <- select(edb, keys = GenenameFilter("NR3C1"))
    .comprehensiveCheckForGene(tmp)
    ## Combine GenenameFilter and TxBiotypeFilter.
    Test2 <- select(edb, keys = ~ symbol == "BCL2" &
                             tx_biotype == "protein_coding")
    expect_equal(Test$EXONID[Test$TXBIOTYPE == "protein_coding"], Test2$EXONID)
    ## Choose selected columns.
    Test3 <- select(edb, keys = gf, columns = c("GENEID", "GENENAME", "SEQNAME"))
    expect_equal(unique(Test[, c("GENEID", "GENENAME", "SEQNAME")]), Test3)
    ## Provide keys.
    Test4 <- select(edb, keys = "BCL2", keytype = "GENENAME")
    expect_equal(Test[, colnames(Test4)], Test4)
    gns <- keys(edb, "GENEID")
    ## Just get stuff from the tx table; should be faster.
    Test <- select(edb, keys = gns, columns = c("GENEID", "SEQNAME"),
                   keytype = "GENEID")
    expect_equal(all(Test$GENEID == gns), TRUE)
    ## Get all lincRNA genes
    Test <- select(edb, keys = "lincRNA", columns = c("GENEID", "GENEBIOTYPE",
                                                      "GENENAME"),
                   keytype = "GENEBIOTYPE")
    Test2 <- select(edb, keys = GeneBiotypeFilter("lincRNA"),
                    columns = c("GENEID", "GENEBIOTYPE", "GENENAME"))
    expect_equal(Test[, colnames(Test2)], Test2)
    ## All on chromosome 21
    Test <- select(edb, keys = "21", columns = c("GENEID", "GENEBIOTYPE",
                                                 "GENENAME"),
                   keytype = "SEQNAME")
    Test2 <- select(edb, keys = ~ seq_name == "21",
                    columns = c("GENEID", "GENEBIOTYPE", "GENENAME"))
    expect_equal(Test[, colnames(Test2)], Test2)
    ## What if we can't find it?
    Test <- select(edb, keys = "bla", columns = c("GENEID", "GENENAME"),
                   keytype = "GENENAME")
    expect_equal(colnames(Test), c("GENEID", "GENENAME"))
    expect_true(nrow(Test) == 0)
    ## TXNAME
    Test <- select(edb, keys = "ENST00000000233",
                   columns = c("GENEID", "GENENAME"), keytype = "TXNAME")
    expect_equal(Test$TXNAME, "ENST00000000233")
    ## Check what happens if we just add TXNAME and also TXID.
    Test2 <- select(edb, keys = list(gf, TxBiotypeFilter("protein_coding")),
                    columns = c("TXID", "TXNAME", "GENENAME", "GENEID"))
    expect_equal(colnames(Test2), c("TXID", "TXNAME", "GENENAME", "GENEID",
                                    "TXBIOTYPE"))
    ## Protein stuff.
    if (hasProteinData(edb)) {
        ## Test:
        ## o if we're fetching with PROTEINID keys we're just getting protein
        ##   coding tx, i.e. those with a tx_cds_seq_start not NULL AND we get
        ##   also those with a uniprot ID null.
        pids <- c("ENSP00000338157", "ENSP00000437716", "ENSP00000443013",
                  "ENSP00000376721", "ENSP00000445047")
        res <- select(edb, keys = pids, keytype = "PROTEINID",
                      columns = c("TXID", "TXCDSSEQSTART", "TXBIOTYPE",
                                  "PROTEINID", "UNIPROTID", "PROTEINDOMAINID"))
        expect_equal(sort(pids), sort(unique(res$PROTEINID)))
        res_2 <- select(edb, keys = ProteinIdFilter(pids),
                        columns = c("TXID", "TXCDSSEQSTART", "TXBIOTYPE",
                                    "PROTEINID", "UNIPROTID", "PROTEINDOMAINID"))
        expect_equal(sort(pids), sort(unique(res_2$PROTEINID)))
        expect_equal(res, res_2)
        expect_true(all(!is.na(res$TXCDSSEQSTART)))
        ## Do we have all of the uniprot ids?
        tmp <- dbGetQuery(dbconn(edb),
                          paste0("select uniprot_id from uniprot where ",
                                 "protein_id in (",
                                 paste0("'", pids,"'", collapse = ", "),")"))
        a <- sort(unique(res$UNIPROTID))
        b <- sort(unique(tmp$uniprot_id))
        a <- a[!is.na(a)]
        b <- b[!is.na(b)]
        expect_equal(a, b)
        ## Do we have all protein domain ids?
        tmp <- dbGetQuery(dbconn(edb),
                          paste0("select protein_domain_id from protein_domain ",
                                 "where protein_id in (",
                                 paste0("'", pids,"'", collapse = ", "),")"))
        a <- sort(unique(res$PROTEINDOMAINID))
        b <- sort(unique(tmp$protein_domain_id))
        a <- a[!is.na(a)]
        b <- b[!is.na(b)]
        expect_equal(a, b)

        ## o if we're fetching with uniprot and protein id filter we get all
        ##   even if they don't have a protein domain.
        upids <- c("ZBT16_HUMAN", "Q71UL7_HUMAN", "Q71UL6_HUMAN", "Q71UL5_HUMAN")
        res <- select(edb, keys = upids, keytype = "UNIPROTID",
                      columns = c("PROTEINID", "UNIPROTID", "PROTEINDOMAINID"))
    }
    edb <- dropFilter(edb)
})

test_that("mapIds works", {
    ## Simple... map gene ids to gene names
    allgenes <- keys(edb, keytype = "GENEID")
    randordergenes <- allgenes[sample(1:length(allgenes), 100)]
    mi <- mapIds(edb, keys = allgenes, keytype = "GENEID", column = "GENENAME")
    expect_equal(allgenes, names(mi))
    ## Ordering should always match the ordering of the input:
    mi <- mapIds(edb, keys = randordergenes, keytype = "GENEID",
                 column = "GENENAME")
    expect_equal(randordergenes, names(mi))
    ## Handle multi mappings.
    ## o first
    first <- mapIds(edb, keys = randordergenes, keytype = "GENEID",
                    column = "TXID")
    expect_equal(names(first), randordergenes)
    ## o list
    lis <- mapIds(edb, keys = randordergenes, keytype = "GENEID",
                  column = "TXID", multiVals = "list")
    expect_equal(names(lis), randordergenes)
    Test <- lapply(lis, function(z){return(z[1])})
    expect_equal(first, unlist(Test))
    ## o filter
    filt <- mapIds(edb, keys = randordergenes, keytype = "GENEID",
                   column = "TXID", multiVals = "filter")
    expect_equal(filt, unlist(lis[unlist(lapply(lis, length)) == 1]))
    ## o asNA
    asNA <- mapIds(edb, keys = randordergenes, keytype = "GENEID",
                   column = "TXID", multiVals = "asNA")
    ## Check what happens if we provide 2 identical keys.
    Test <- mapIds(edb, keys = c("BCL2", "BCL2L11", "BCL2"),
                   keytype = "GENENAME", column = "TXID")
    expect_equal(names(Test), c("BCL2", "BCL2L11", "BCL2"))
    expect_true(length(unique(Test)) == 2)
    ## Submit Filter:
    Test <- mapIds(edb, keys = SeqNameFilter("Y"), column = "GENEID",
                   multiVals = "list")
    TestS <- select(edb, keys = Test[[1]], columns = "SEQNAME",
                    keytype = "GENEID")
    expect_equal(unique(TestS$SEQNAME), "Y")
    ## Submit 2 filter.LLLLL
    Test <- mapIds(edb, keys = ~ seq_name == "Y" & seq_strand == "-",
                   multiVals = "list", column = "GENEID")
    TestS <- select(edb, keys = Test[[1]], keytype = "GENEID",
                    columns = c("SEQNAME", "SEQSTRAND"))
    expect_true(all(TestS$SEQNAME == "Y"))
    expect_true(all(TestS$SEQSTRAND == -1))

    ## Now using protein annotations:
    if (hasProteinData(edb)) {
        library(RSQLite)
        txids <- keys(edb, keytype = "TXID", filter = GenenameFilter("ZBTB16"))
        mapd <- mapIds(edb, keys = txids, keytype = "TXID", column = "GENENAME")
        expect_equal(names(mapd), txids)
        expect_true(all(mapd == "ZBTB16"))
        ## Map to protein ids.
        mapd <- mapIds(edb, keys = txids, keytype = "TXID", column = "PROTEINID")
        res <- dbGetQuery(dbconn(edb),
                          paste0("select protein_id from protein where tx_id in",
                                 " (", paste0("'", txids,"'", collapse = ", "),
                                 ")"))
        pids <- mapd[!is.na(mapd)]
        expect_true(all(pids %in% res$protein_id))
        ## multi-mapping:
        ## proteins and uniprot.
        mapd <- mapIds(edb, keys = pids, keytype = "PROTEINID",
                       column = "UNIPROTID", multiVals = "list")
        mapd <- mapd[!is.na(mapd)]
        res <- dbGetQuery(dbconn(edb),
                          paste0("select protein_id, uniprot_id from uniprot ",
                                 "where protein_id in (",
                                 paste0("'", pids, "'", collapse = ", "), ")"))
        res <- split(res$uniprot_id, res$protein_id)
        expect_equal(mapd, res[names(mapd)])
        ## Just to ensure:
        tmp <- proteins(edb, filter = ProteinIdFilter(pids),
                        columns = c("uniprot_id", "protein_id"))
        upids <- tmp$uniprot_id[!is.na(tmp$uniprot_id)]
        expect_true(all(res$uniprot_id %in% upids))
        ## map protein ids to gene name
        mapd <- mapIds(edb, keys = pids, keytype = "PROTEINID",
                       column = "GENENAME")
        expect_true(all(mapd == "ZBTB16"))
    }
})

## Test if the results are properly sorted if we submit a single filter or just keys.
test_that("select results are properly sorted", {
    ks <- c("ZBTB16", "BCL2", "SKA2", "BCL2L11")
    ## gene_id
    res <- select(edb, keys = ks, keytype = "GENENAME")
    expect_equal(unique(res$GENENAME), ks)
    res <- select(edb, keys = GenenameFilter(ks))
    expect_equal(unique(res$GENENAME), ks)
    ## Using two filters;
    res <- select(edb, keys = list(GenenameFilter(ks),
                                   TxBiotypeFilter("nonsense_mediated_decay")))
    ## We don't expect same sorting here!
    expect_true(!all(unique(res$GENENAME) == ks[ks %in% unique(res$GENENAME)]))
    res2 <- select(edb, keys = ~ genename == ks &
                            tx_biotype == "nonsense_mediated_decay")
    expect_equal(res, res2)
    ## symbol
    res <- select(edb, keys = ks, keytype = "SYMBOL",
                  columns = c("GENENAME", "SYMBOL", "SEQNAME"))
    expect_equal(res$SYMBOL, ks)
    expect_equal(res$GENENAME, ks)
    ## tx_biotype
    ks <- c("retained_intron", "nonsense_mediated_decay")
    res <- select(edb, keys = ks, keytype = "TXBIOTYPE",
                  columns = c("GENENAME", "TXBIOTYPE"))
    expect_equal(unique(res$TXBIOTYPE), ks)
    res <- select(edb, keys = TxBiotypeFilter(ks),
                  keytype = "TXBIOTYPE", columns = c("GENENAME", "TXBIOTYPE"))
    expect_equal(unique(res$TXBIOTYPE), ks)
})

test_that("select works with symbol as keytype", {
    ks <- c("ZBTB16", "BCL2", "SKA2", "BCL2L11")
    res <- select(edb, keys = ks, keytype = "GENENAME")
    res2 <- select(edb, keys = ks, keytype = "SYMBOL")
    expect_equal(res, res2)
    res <- select(edb, keys = GenenameFilter(ks),
                  columns = c("TXNAME", "SYMBOL", "GENEID"))
    expect_equal(colnames(res), c("TXNAME", "SYMBOL", "GENEID", "GENENAME"))
    res <- select(edb, keys = ~ symbol == ks, columns=c("GENEID"))
    expect_equal(colnames(res), c("GENEID", "SYMBOL"))
    expect_equal(unique(res$SYMBOL), ks)
    res <- select(edb, keys = list(SeqNameFilter("Y"),
                                   GeneBiotypeFilter("lincRNA")),
                  columns = c("GENEID", "SYMBOL"))
    expect_equal(colnames(res), c("GENEID", "SYMBOL", "SEQNAME", "GENEBIOTYPE"))
    expect_true(all(res$SEQNAME == "Y"))
})

test_that("select works with txname", {
    ## TXNAME as a column
    ks <- c("ZBTB16", "BCL2", "SKA2")
    res <- select(edb, keys = ks, keytype = "GENENAME", columns = c("TXNAME"))
    expect_equal(colnames(res), c("GENENAME", "TXNAME"))
})

test_that(".keytype2FilterMapping works", {
    ## Check whether or not we're getting protein columns.
    res <- ensembldb:::.keytype2FilterMapping()
    expect_equal(names(res),
                 c("ENTREZID", "GENEID", "GENEBIOTYPE", "GENENAME", "TXID",
                   "TXBIOTYPE", "EXONID", "SEQNAME", "SEQSTRAND", "TXNAME",
                   "SYMBOL"))
    res <- ensembldb:::.keytype2FilterMapping(TRUE)
    expect_true(all(c("PROTEINID", "UNIPROTID", "PROTEINDOMAINID") %in%
                    names(res)))
})

test_that("keytypes works", {
    keyt <- c("ENTREZID", "GENEID", "GENEBIOTYPE", "GENENAME", "TXID",
              "TXBIOTYPE", "EXONID", "SEQNAME", "SEQSTRAND", "TXNAME",
              "SYMBOL")
    res <- keytypes(edb)
    if (hasProteinData(edb)) {
        expect_equal(res, sort(c(keyt, "PROTEINID", "UNIPROTID",
                                 "PROTEINDOMAINID")))
    } else {
        expect_equal(res, sort(keyt))
    }
})

test_that("filterForKeytype works", {
    res <- ensembldb:::filterForKeytype("SYMBOL")
    expect_true(is(res, "SymbolFilter"))
    if (hasProteinData(edb)) {
        res <- ensembldb:::filterForKeytype("PROTEINDOMAINID", edb)
        expect_true(is(res, "ProtDomIdFilter"))
    }
    res <- ensembldb:::filterForKeytype("TXID")
    expect_true(is(res, "TxIdFilter"))
})
