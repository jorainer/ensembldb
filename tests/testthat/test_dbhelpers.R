
test_that("prefixColumns works", {
    res <- ensembldb:::prefixColumns(edb, columns = "a")
    expect_true(is.null(res))
    expect_error(ensembldb:::prefixColumns(edb, columns = "a", clean = FALSE))
    res <- ensembldb:::prefixColumns(edb, columns = c("gene_id", "a"),
                                     clean = FALSE)
    expect_equal(names(res), "gene")
    expect_equal(res$gene, "gene.gene_id")
    ## The "new" prefixColumns function ALWAYS returns the first table in which
    ## a column was found; tables are ordered as in listTables
    res <- ensembldb:::prefixColumns(edb, columns = c("tx_id", "gene_id",
                                                      "tx_biotype"))
    want <- list(gene = "gene.gene_id",
                 tx = c("tx.tx_id", "tx.tx_biotype"))
    expect_equal(res, want)
    ##
    res <- ensembldb:::prefixColumns(edb, columns = c("exon_idx", "seq_name",
                                                      "gene_id"))
    want <- list(gene = c("gene.gene_id", "gene.seq_name"),
                 tx2exon = "tx2exon.exon_idx")
    expect_equal(res, want)
    ##
    res <- ensembldb:::prefixColumns(edb, columns = c("exon_idx", "seq_name",
                                                      "gene_id", "exon_id"))
    want <- list(gene = c("gene.gene_id", "gene.seq_name"),
                 tx2exon = c("tx2exon.exon_id", "tx2exon.exon_idx"))
    expect_equal(res, want)

    if (hasProteinData(edb)) {
        res <- ensembldb:::prefixColumns(edb,
                                         columns = c("tx_id", "protein_id"))
        want <- list(tx = "tx.tx_id", protein = "protein.protein_id")
        expect_equal(res, want)
        ##
        res <- ensembldb:::prefixColumns(edb,
                                         columns = c("uniprot_id",
                                                     "protein_domain_id"))
        want <- list(uniprot = "uniprot.uniprot_id",
                     protein_domain = "protein_domain.protein_domain_id")
        expect_equal(res, want)
        ##
        res <- ensembldb:::prefixColumns(edb,
                                         columns = c("uniprot_id",
                                                     "protein_domain_id",
                                                     "protein_id", "tx_id"))
        want = list(tx = "tx.tx_id", protein = "protein.protein_id",
                    uniprot = "uniprot.uniprot_id",
                    protein_domain = "protein_domain.protein_domain_id")
        expect_equal(res, want)
    }
})

############################################################
## Test the new join engine.
## o use the startWith argument.
## o change the join argument.
test_that("joinTwoTables works", {
    ## Check errors:
    expect_error(ensembldb:::joinTwoTables(a = "gene", b = "dont exist"))
    expect_error(ensembldb:::joinTwoTables(a = c("a", "b"), b = "gene"))
    ## Working example:
    res <- ensembldb:::joinTwoTables(a = c("a", "gene"), b = "tx")
    expect_equal(sort(res[1:2]), c("gene", "tx"))
    expect_equal(res[3], "on (gene.gene_id=tx.gene_id)")
    ## Error
    expect_error(ensembldb:::joinTwoTables(a = "tx", b = "exon"))
    ## Working example:
    res <- ensembldb:::joinTwoTables(a = c("tx"), b = c("exon", "tx2exon"))
    expect_equal(sort(res[1:2]), c("tx", "tx2exon"))
    expect_equal(res[3], "on (tx.tx_id=tx2exon.tx_id)")
    res <- ensembldb:::joinTwoTables(a = c("chromosome", "gene", "tx"),
                                     b = c("exon", "protein", "tx2exon"))
    expect_equal(sort(res[1:2]), c("tx", "tx2exon"))
    expect_equal(res[3], "on (tx.tx_id=tx2exon.tx_id)")
})

test_that("joinQueryOnTables2 and joinQueryOnColumns2 work", {
    ## exceptions
    expect_error(ensembldb:::joinQueryOnTables2(edb, tab = c("a", "exon")))
    res <- ensembldb:::joinQueryOnTables2(edb, tab = c("gene", "exon"))
    want <- paste0("gene join tx on (gene.gene_id=tx.gene_id) join",
                   " tx2exon on (tx.tx_id=tx2exon.tx_id) join",
                   " exon on (tx2exon.exon_id=exon.exon_id)")
    ## The "default" order is gene->tx->tx2exon->exon
    expect_equal(res, want)
    res <- ensembldb:::joinQueryOnColumns2(edb, columns = c("exon_seq_start",
                                                            "gene_name"))
    expect_equal(res, want)
    ## Same but in the order: exon->tx2exon->tx->gene
    res <- ensembldb:::joinQueryOnTables2(edb, tab = c("gene", "exon"),
                                          startWith = "exon")
    want <- paste0("exon join tx2exon on (tx2exon.exon_id=exon.exon_id)",
                   " join tx on (tx.tx_id=tx2exon.tx_id) join",
                   " gene on (gene.gene_id=tx.gene_id)")
    expect_equal(res, want)
    res <- ensembldb:::joinQueryOnColumns2(edb, columns = c("exon_seq_start",
                                                            "gene_name"),
                                           startWith = "exon")
    expect_equal(res, want)
    ## That would be less expensive, but with "startWith" we force it to start
    ## from table exon, instead of just using tx2exon and tx.
    res <- ensembldb:::joinQueryOnColumns2(edb, columns = c("exon_id",
                                                            "gene_id"),
                                           startWith = "exon")
    expect_equal(res, want)
    ## Check proteins too.
    if (hasProteinData(edb)) {
        res <- ensembldb:::joinQueryOnTables2(edb, tab = c("protein", "gene",
                                                           "exon"))
        ## That should be: gene->tx->tx2exon->exon->protein
        want <- paste0("gene join tx on (gene.gene_id=tx.gene_id) join",
                       " tx2exon on (tx.tx_id=tx2exon.tx_id) join",
                       " exon on (tx2exon.exon_id=exon.exon_id) left outer join",
                       " protein on (tx.tx_id=protein.tx_id)")
        expect_equal(res, want)
        res <- ensembldb:::joinQueryOnColumns2(edb,
                                               columns = c("protein_id",
                                                           "gene_name",
                                                           "exon_seq_start"))
        expect_equal(res, want)
        res <- ensembldb:::joinQueryOnTables2(edb, tab = c("protein", "gene"),
                                              startWith = "protein")
        want <- paste0("protein left outer join tx on (tx.tx_id=protein.tx_id)",
                       " join gene on (gene.gene_id=tx.gene_id)")
        expect_equal(res, want)
        res <- ensembldb:::joinQueryOnColumns2(edb, columns = c("protein_id",
                                                                "gene_name"),
                                               startWith = "protein")
        expect_equal(res, want)
    }
})

test_that("addRequiredTables works", {
    have <- c("exon", "gene")
    need <- c("exon", "gene", "tx2exon", "tx")
    expect_equal(sort(need), sort(ensembldb:::addRequiredTables(edb, have)))

    have <- c("exon", "chromosome")
    need <- c("exon", "tx2exon", "tx", "gene", "chromosome")
    expect_equal(sort(need), sort(ensembldb:::addRequiredTables(edb, have)))

    have <- c("chromosome", "tx")
    need <- c("chromosome", "tx", "gene")
    expect_equal(sort(need), sort(ensembldb:::addRequiredTables(edb, have)))

    if (hasProteinData(edb)) {
        have <- c("uniprot", "exon")
        need <- c("uniprot", "exon", "protein", "tx", "tx2exon")
        expect_equal(sort(need),
                     sort(ensembldb:::addRequiredTables(edb, have)))

        have <- c("uniprot", "chromosome")
        need <- c("uniprot", "chromosome", "protein", "tx", "gene")
        expect_equal(sort(need),
                     sort(ensembldb:::addRequiredTables(edb, have)))

        have <- c("protein_domain", "gene")
        need <- c("protein_domain", "gene", "protein", "tx")
        expect_equal(sort(need),
                     sort(ensembldb:::addRequiredTables(edb, have)))

        have <- c("protein", "exon")
        need <- c("protein", "exon", "tx", "tx2exon")
        expect_equal(sort(need),
                     sort(ensembldb:::addRequiredTables(edb, have)))
    }
})

test_that(".buildQuery with filter works", {
    columns <- c("gene_id", "gene_name", "exon_id")
    gnf <- GenenameFilter("BCL2")
    Q <- ensembldb:::.buildQuery(edb, columns = columns,
                                 filter = AnnotationFilterList(gnf))
    want <- paste0("select distinct gene.gene_id,gene.gene_name,",
                   "tx2exon.exon_id from gene join tx on (gene.gene_id",
                   "=tx.gene_id) join tx2exon on (tx.tx_id=tx2exon.tx_id)",
                   " where (gene.gene_name = 'BCL2')")
    expect_equal(Q, want)
    library(RSQLite)
    res <- dbGetQuery(dbconn(edb), Q)
    expect_equal(unique(res$gene_name), "BCL2")
    ## Two GeneNameFilters combined with or
    gnf2 <- GenenameFilter("BCL2L11")
    columns <- c("gene_id", "gene_name", "exon_id")
    Q <- ensembldb:::.buildQuery(edb, columns = columns,
                                 filter = AnnotationFilterList(gnf, gnf2,
                                                               logOp = "|"))
    want <- paste0("select distinct gene.gene_id,gene.gene_name,",
                   "tx2exon.exon_id from gene join tx on (gene.gene_id",
                   "=tx.gene_id) join tx2exon on (tx.tx_id=tx2exon.tx_id)",
                   " where (gene.gene_name = 'BCL2' or gene.gene_name = ",
                   "'BCL2L11')")
    expect_equal(Q, want)
    res <- dbGetQuery(dbconn(edb), Q)
    expect_true(all(res$gene_name %in% c("BCL2", "BCL2L11")))
    ## Combine with a SeqnameFilter.
    snf <- SeqNameFilter(2)
    flt <- AnnotationFilterList(gnf, gnf2, snf, logOp = c("|", "&"))
    Q <- ensembldb:::.buildQuery(edb, columns = columns, filter = flt)
    want <- paste0("select distinct gene.gene_id,gene.gene_name,",
                   "tx2exon.exon_id,gene.seq_name from gene join tx on (",
                   "gene.gene_id=tx.gene_id) join tx2exon on (tx.tx_id=",
                   "tx2exon.tx_id) where (gene.gene_name = 'BCL2' or ",
                   "gene.gene_name = 'BCL2L11' and gene.seq_name = '2')")
    expect_equal(Q, want)
    res <- dbGetQuery(dbconn(edb), Q)
    expect_true(all(res$gene_name %in% c("BCL2", "BCL2L11")))
    ## now with a nested AnnotationFilterList:
    flt <- AnnotationFilterList(AnnotationFilterList(gnf, gnf2, logOp = "|"),
                                snf, logOp = "&")
    Q <- ensembldb:::.buildQuery(edb, columns = columns, filter = flt)
    want <- paste0("select distinct gene.gene_id,gene.gene_name,",
                   "tx2exon.exon_id,gene.seq_name from gene join tx on (",
                   "gene.gene_id=tx.gene_id) join tx2exon on (tx.tx_id=",
                   "tx2exon.tx_id) where ((gene.gene_name = 'BCL2' or ",
                   "gene.gene_name = 'BCL2L11') and gene.seq_name = '2')")
    expect_equal(Q, want)
    res <- dbGetQuery(dbconn(edb), Q)
    expect_true(all(res$gene_name %in% c("BCL2L11")))
    ## If we only want to get BCL2L11 back:
    flt <- AnnotationFilterList(GenenameFilter(c("BCL2", "BCL2L11")), snf,
                                logOp = "&")
    Q <- ensembldb:::.buildQuery(edb, columns = columns, filter = flt)
    want <- paste0("select distinct gene.gene_id,gene.gene_name,",
                   "tx2exon.exon_id,gene.seq_name from gene join tx on (",
                   "gene.gene_id=tx.gene_id) join tx2exon on (tx.tx_id=",
                   "tx2exon.tx_id) where (gene.gene_name in ('BCL2','BCL2L11'",
                   ") and gene.seq_name = '2')")
    expect_equal(Q, want)
    res <- dbGetQuery(dbconn(edb), Q)
    expect_true(all(res$gene_name == "BCL2L11"))
    
    ## Check with a GRangesFilter.
    grf <- GRangesFilter(GRanges(seqnames = 18, IRanges(60790600, 60790700)))
    flt <- AnnotationFilterList(grf)
    Q <- ensembldb:::.buildQuery(edb, columns = columns, filter = flt)
    want <- paste0("select distinct gene.gene_id,gene.gene_name,tx2exon.",
                   "exon_id,gene.gene_seq_start,gene.gene_seq_end,gene.seq_name",
                   ",gene.seq_strand from gene join tx on (gene.gene_id",
                   "=tx.gene_id) join tx2exon on (tx.tx_id=tx2exon.tx_id) ",
                   "where ((gene.gene_seq_start<=60790700 and gene.gene_seq",
                   "_end>=60790600 and gene.seq_name='18'))")
    expect_equal(Q, want)
    res <- dbGetQuery(dbconn(edb), Q)
    expect_true(all(res$gene_name == "BCL2"))
})

test_that("buildQuery with startWith works", {
    columns <- c("gene_id", "gene_name", "exon_id")
    Q <- ensembldb:::.buildQuery(edb, columns = columns)
    want <- paste0("select distinct gene.gene_id,gene.gene_name,",
                   "tx2exon.exon_id from gene join tx on (gene.gene_id",
                   "=tx.gene_id) join tx2exon on (tx.tx_id=tx2exon.tx_id)")
    expect_equal(Q, want)
    ## Different if we use startWith = exon
    Q <- ensembldb:::.buildQuery(edb, columns = columns, startWith = "exon")
    want <- paste0("select distinct gene.gene_id,gene.gene_name,",
                   "tx2exon.exon_id from exon join tx2exon on (tx2exon.exon_id",
                   "=exon.exon_id) join tx on (tx.tx_id=tx2exon.tx_id)",
                   " join gene on (gene.gene_id=tx.gene_id)")
    expect_equal(Q, want)
    Q <- ensembldb:::.buildQuery(edb, columns = c("gene_id", "tx_biotype"))
    want <- paste0("select distinct gene.gene_id,tx.tx_biotype from gene ",
                   "join tx on (gene.gene_id=tx.gene_id)")
    expect_equal(Q, want)
    Q <- ensembldb:::.buildQuery(edb, columns = c("gene_id", "tx_biotype"),
                                 startWith = "exon")
    want <- paste0("select distinct gene.gene_id,tx.tx_biotype from exon ",
                   "join tx2exon on (tx2exon.exon_id=exon.exon_id) join ",
                   "tx on (tx.tx_id=tx2exon.tx_id) join ",
                   "gene on (gene.gene_id=tx.gene_id)")
    expect_equal(Q, want)
    if (hasProteinData(edb)) {
        ## Protein columns.
        Q <- ensembldb:::.buildQuery(edb,
                                     columns = c("protein_id", "uniprot_id",
                                                 "protein_domain_id"))
        want <- paste0("select distinct protein.protein_id,uniprot.uniprot_id,",
                       "protein_domain.protein_domain_id from protein left ",
                       "outer join protein_domain on (protein.protein_id=",
                       "protein_domain.protein_id) left outer join ",
                       "uniprot on (protein.protein_id=uniprot.protein_id)")
        expect_equal(Q, want)
        ## start at protein
        Q <- ensembldb:::.buildQuery(edb,
                                     columns = c("protein_id", "uniprot_id",
                                                 "protein_domain_id"),
                                     startWith = "protein")
        want <- paste0("select distinct protein.protein_id,uniprot.uniprot_id,",
                       "protein_domain.protein_domain_id from protein left ",
                       "outer join protein_domain on (protein.protein_id=",
                       "protein_domain.protein_id) left outer join ",
                       "uniprot on (protein.protein_id=uniprot.protein_id)")
        expect_equal(Q, want)
        ## start at uniprot.
        Q <- ensembldb:::.buildQuery(edb,
                                     columns = c("protein_id", "uniprot_id",
                                                 "protein_domain_id"),
                                     startWith = "uniprot")
        want <- paste0("select distinct protein.protein_id,uniprot.uniprot_id,",
                       "protein_domain.protein_domain_id from uniprot left ",
                       "outer join protein on (protein.protein_id=",
                       "uniprot.protein_id) left outer join",
                       " protein_domain on (protein.protein_id=",
                       "protein_domain.protein_id)")
        expect_equal(Q, want)
        ## join with tx.
        Q <- ensembldb:::.buildQuery(edb, columns = c("tx_id", "protein_id",
                                                      "uniprot_id", "gene_id"))
        want <- paste0("select distinct tx.tx_id,protein.protein_id,",
                       "uniprot.uniprot_id,gene.gene_id from gene join ",
                       "tx on (gene.gene_id=tx.gene_id) left outer join protein",
                       " on (tx.tx_id=protein.tx_id) left outer join uniprot on",
                       " (protein.protein_id=uniprot.protein_id)")
        expect_equal(Q, want)
        ## if we started from protein:
        Q <- ensembldb:::.buildQuery(edb, columns = c("tx_id", "protein_id",
                                                      "uniprot_id", "gene_id"),
                                     startWith = "protein")
        want <- paste0("select distinct tx.tx_id,protein.protein_id,",
                       "uniprot.uniprot_id,gene.gene_id from protein left outer",
                       " join tx on (tx.tx_id=protein.tx_id) join gene on",
                       " (gene.gene_id=tx.gene_id) left outer join uniprot on",
                       " (protein.protein_id=uniprot.protein_id)")
        expect_equal(Q, want)
    }
})



## This test is an important one as it checks that we don't miss any entries
## from the database, e.g. if we query gene and join with protein that we don't
## miss any non-coding transcripts, or if we join protein with uniprot or
## protein_domain that we don't miss any values.
test_that("query is valid", {
    ## Check RNA/DNA tables; shouldn't be a problem there, though.
    Ygns <- genes(edb, filter = SeqNameFilter("Y"), return.type = "data.frame")
    Ytxs <- transcripts(edb, filter = SeqNameFilter("Y"),
                        return.type = "data.frame",
                        columns = c("gene_id", "tx_id", "tx_biotype"))
    Yexns <- exons(edb, filter = SeqNameFilter("Y"), return.type = "data.frame",
                   columns = c("exon_id", "gene_id"))
    expect_true(all(unique(Ygns$gene_id) %in% unique(Yexns$gene_id)))
    expect_true(all(unique(Ygns$gene_id) %in% unique(Ytxs$gene_id)))
    ## Check gene with protein
    if (hasProteinData(edb)) {
        library(RSQLite)
        ## Simulate what a simple join would do:
        gns_f <- dbGetQuery(dbconn(edb),
                            paste0("select gene.gene_id, tx.tx_id, tx_biotype, ",
                                   "protein_id from gene join tx on ",
                                   "(gene.gene_id=tx.gene_id) join protein on ",
                                   "(tx.tx_id=protein.tx_id) ",
                                   "where seq_name = 'Y'"))
        ## We expect that gns_f is smaller, but that all protein_coding tx are
        ## there.
        expect_true(length(unique(gns_f$gene_id)) < length(unique(Ygns$gene_id)))
        expect_true(all(unique(Ytxs[Ytxs$tx_biotype == "protein_coding", "tx_id"])
                        %in% unique(gns_f$tx_id)))
        ## Now test the "real" query:
        Ygns_2 <- genes(edb, filter = SeqNameFilter("Y"),
                        return.type = "data.frame",
                        columns = c("gene_id", "tx_id", "tx_biotype",
                                    "protein_id"))
        ## We expect that ALL genes are present and ALL tx:
        expect_true(all(unique(Ygns$gene_id) %in% unique(Ygns_2$gene_id)))
        expect_true(all(unique(Ygns$tx_id) %in% unique(Ygns_2$tx_id)))

        ## Get all the tx with protein_id
        txs <- transcripts(edb, columns = c("tx_id", "protein_id"),
                           return.type = "data.frame")
        txids <- dbGetQuery(dbconn(edb), "select tx_id from tx;")[, "tx_id"]
        protids <- dbGetQuery(dbconn(edb),
                              "select protein_id from protein;")[, "protein_id"]
        expect_true(all(txids %in% txs$tx_id))
        expect_true(all(protids %in% txs$protein_id))

        ## Check protein with uniprot
        uniprotids <- dbGetQuery(dbconn(edb),
                                 "select uniprot_id from uniprot")$uniprot_id

        ## Check protein with protein domain
        ## Check protein_domain with uniprot
    }
})

test_that(".getWhat works", {
    library(RSQLite)
    Q_2 <- paste0("select * from gene join tx on (gene.gene_id=tx.gene_id)",
                  " join tx2exon on (tx.tx_id=tx2exon.tx_id) where",
                  " gene.gene_id = 'ENSG00000000005'")
    res_2 <- dbGetQuery(dbconn(edb), Q_2)
    gf <- GeneIdFilter("ENSG00000000005")
    res_3 <- ensembldb:::.getWhat(edb, columns = c("gene_name", "exon_idx"),
                                 filter = AnnotationFilterList(gf))
    expect_identical(res_3, unique(res_2[, colnames(res_3)]))
})

test_that(".logOp2SQL works", {
    expect_equal(ensembldb:::.logOp2SQL("|"), "or")
    expect_equal(ensembldb:::.logOp2SQL("&"), "and")
    expect_equal(ensembldb:::.logOp2SQL("dfdf"), NULL)
})

