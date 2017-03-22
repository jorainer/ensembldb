#p###########################################################
## Tests for the dbhelpers.R functions
library(EnsDb.Hsapiens.v75)
edb <- EnsDb.Hsapiens.v75

test_addRequiredTables <- function() {
    have <- c("exon", "gene")
    need <- c("exon", "gene", "tx2exon", "tx")
    checkEquals(sort(need), sort(ensembldb:::addRequiredTables(edb, have)))

    have <- c("exon", "chromosome")
    need <- c("exon", "tx2exon", "tx", "gene", "chromosome")
    checkEquals(sort(need), sort(ensembldb:::addRequiredTables(edb, have)))

    have <- c("chromosome", "tx")
    need <- c("chromosome", "tx", "gene")
    checkEquals(sort(need), sort(ensembldb:::addRequiredTables(edb, have)))

    if (hasProteinData(edb)) {
        have <- c("uniprot", "exon")
        need <- c("uniprot", "exon", "protein", "tx", "tx2exon")
        checkEquals(sort(need),
                    sort(ensembldb:::addRequiredTables(edb, have)))

        have <- c("uniprot", "chromosome")
        need <- c("uniprot", "chromosome", "protein", "tx", "gene")
        checkEquals(sort(need),
                    sort(ensembldb:::addRequiredTables(edb, have)))

        have <- c("protein_domain", "gene")
        need <- c("protein_domain", "gene", "protein", "tx")
        checkEquals(sort(need),
                    sort(ensembldb:::addRequiredTables(edb, have)))

        have <- c("protein", "exon")
        need <- c("protein", "exon", "tx", "tx2exon")
        checkEquals(sort(need),
                    sort(ensembldb:::addRequiredTables(edb, have)))
    }
}

test_prefixColumns <- function() {
    res <- ensembldb:::prefixColumns(edb, columns = "a")
    checkTrue(is.null(res))
    checkException(ensembldb:::prefixColumns(edb, columns = "a", clean = FALSE))
    res <- ensembldb:::prefixColumns(edb, columns = c("gene_id", "a"),
                                      clean = FALSE)
    checkEquals(names(res), "gene")
    checkEquals(res$gene, "gene.gene_id")
    ## The "new" prefixColumns function ALWAYS returns the first table in which
    ## a column was found; tables are ordered as in listTables
    res <- ensembldb:::prefixColumns(edb, columns = c("tx_id", "gene_id",
                                                       "tx_biotype"))
    want <- list(gene = "gene.gene_id",
                 tx = c("tx.tx_id", "tx.tx_biotype"))
    checkEquals(res, want)
    ##
    res <- ensembldb:::prefixColumns(edb, columns = c("exon_idx", "seq_name",
                                                       "gene_id"))
    want <- list(gene = c("gene.gene_id", "gene.seq_name"),
                 tx2exon = "tx2exon.exon_idx")
    checkEquals(res, want)
    ##
    res <- ensembldb:::prefixColumns(edb, columns = c("exon_idx", "seq_name",
                                                       "gene_id", "exon_id"))
    want <- list(gene = c("gene.gene_id", "gene.seq_name"),
                 tx2exon = c("tx2exon.exon_id", "tx2exon.exon_idx"))
    checkEquals(res, want)

    if (hasProteinData(edb)) {
        res <- ensembldb:::prefixColumns(edb,
                                          columns = c("tx_id", "protein_id"))
        want <- list(tx = "tx.tx_id", protein = "protein.protein_id")
        checkEquals(res, want)
        ##
        res <- ensembldb:::prefixColumns(edb,
                                          columns = c("uniprot_id",
                                                      "protein_domain_id"))
        want <- list(uniprot = "uniprot.uniprot_id",
                     protein_domain = "protein_domain.protein_domain_id")
        checkEquals(res, want)
        ##
        res <- ensembldb:::prefixColumns(edb,
                                          columns = c("uniprot_id",
                                                      "protein_domain_id",
                                                      "protein_id", "tx_id"))
        want = list(tx = "tx.tx_id", protein = "protein.protein_id",
                    uniprot = "uniprot.uniprot_id",
                    protein_domain = "protein_domain.protein_domain_id")
        checkEquals(res, want)
    }
}

test_buildQuery_filter <- function() {
    columns <- c("gene_id", "gene_name", "exon_id")
    gnf <- GenenameFilter("BCL2")
    Q <- ensembldb:::.buildQuery(edb, columns = columns,
                                 filter = AnnotationFilterList(gnf))
    want <- paste0("select distinct gene.gene_id,gene.gene_name,",
                   "tx2exon.exon_id from gene join tx on (gene.gene_id",
                   "=tx.gene_id) join tx2exon on (tx.tx_id=tx2exon.tx_id)",
                   " where (gene.gene_name = 'BCL2')")
    checkEquals(Q, want)
    library(RSQLite)
    res <- dbGetQuery(dbconn(edb), Q)
    checkEquals(unique(res$gene_name), "BCL2")
    ## Two GeneNameFilters combined with or
    gnf2 <- GenenameFilter("BCL2L11")
    columns <- c("gene_id", "gene_name", "exon_id")
    Q <- ensembldb:::.buildQuery(edb, columns = columns,
                                 filter = AnnotationFilterList(gnf, gnf2,
                                                               logOp = "|"))
    want <- paste0("select distinct gene.gene_id,gene.gene_name,",
                   "tx2exon.exon_id from gene join tx on (gene.gene_id",
                   "=tx.gene_id) join tx2exon on (tx.tx_id=tx2exon.tx_id)",
                   " where (gene.gene_name = 'BCL2') or (gene.gene_name = ",
                   "'BCL2L11')")
    checkEquals(Q, want)
    res <- dbGetQuery(dbconn(edb), Q)
    checkTrue(all(res$gene_name %in% c("BCL2", "BCL2L11")))
    ## Combine with a SeqnameFilter.
    snf <- SeqNameFilter(2)
    flt <- AnnotationFilterList(gnf, gnf2, snf, logOp = c("|", "&"))
    Q <- ensembldb:::.buildQuery(edb, columns = columns, filter = flt)
    want <- paste0("select distinct gene.gene_id,gene.gene_name,",
                   "tx2exon.exon_id,gene.seq_name from gene join tx on (",
                   "gene.gene_id=tx.gene_id) join tx2exon on (tx.tx_id=",
                   "tx2exon.tx_id) where (gene.gene_name = 'BCL2') or (",
                   "gene.gene_name = 'BCL2L11') and (gene.seq_name = '2')")
    checkEquals(Q, want)
    res <- dbGetQuery(dbconn(edb), Q)
    checkTrue(all(res$gene_name %in% c("BCL2", "BCL2L11")))
    ## If we only want to get BCL2L11 back:
    flt <- AnnotationFilterList(GenenameFilter(c("BCL2", "BCL2L11")), snf,
                                logOp = "&")
    Q <- ensembldb:::.buildQuery(edb, columns = columns, filter = flt)
    want <- paste0("select distinct gene.gene_id,gene.gene_name,",
                   "tx2exon.exon_id,gene.seq_name from gene join tx on (",
                   "gene.gene_id=tx.gene_id) join tx2exon on (tx.tx_id=",
                   "tx2exon.tx_id) where (gene.gene_name in ('BCL2','BCL2L11'",
                   ")) and (gene.seq_name = '2')")
    checkEquals(Q, want)
    res <- dbGetQuery(dbconn(edb), Q)
    checkTrue(all(res$gene_name == "BCL2L11"))
    
    ## Check with a GRangesFilter.
    grf <- GRangesFilter(GRanges(seqnames = 18, IRanges(60790600, 60790700)))
    flt <- AnnotationFilterList(grf)
    Q <- ensembldb:::.buildQuery(edb, columns = columns, filter = flt)
    want <- paste0("select distinct gene.gene_id,gene.gene_name,tx2exon.",
                   "exon_id,gene.seq_name,gene.gene_seq_start,gene.gene_seq",
                   "_end,gene.seq_strand from gene join tx on (gene.gene_id",
                   "=tx.gene_id) join tx2exon on (tx.tx_id=tx2exon.tx_id) ",
                   "where (gene.gene_seq_start<=60790700 and gene.gene_seq",
                   "_end>=60790600 and gene.seq_name='18')")
    checkEquals(Q, want)
    res <- dbGetQuery(dbconn(edb), Q)
    checkTrue(all(res$gene_name == "BCL2"))
}

test_buildQuery_startWith <- function() {
    columns <- c("gene_id", "gene_name", "exon_id")
    Q <- ensembldb:::.buildQuery(edb, columns = columns)
    want <- paste0("select distinct gene.gene_id,gene.gene_name,",
                   "tx2exon.exon_id from gene join tx on (gene.gene_id",
                   "=tx.gene_id) join tx2exon on (tx.tx_id=tx2exon.tx_id)")
    checkEquals(Q, want)
    ## Different if we use startWith = exon
    Q <- ensembldb:::.buildQuery(edb, columns = columns, startWith = "exon")
    want <- paste0("select distinct gene.gene_id,gene.gene_name,",
                   "tx2exon.exon_id from exon join tx2exon on (tx2exon.exon_id",
                   "=exon.exon_id) join tx on (tx.tx_id=tx2exon.tx_id)",
                   " join gene on (gene.gene_id=tx.gene_id)")
    checkEquals(Q, want)
    Q <- ensembldb:::.buildQuery(edb, columns = c("gene_id", "tx_biotype"))
    want <- paste0("select distinct gene.gene_id,tx.tx_biotype from gene ",
                   "join tx on (gene.gene_id=tx.gene_id)")
    checkEquals(Q, want)
    Q <- ensembldb:::.buildQuery(edb, columns = c("gene_id", "tx_biotype"),
                                 startWith = "exon")
    want <- paste0("select distinct gene.gene_id,tx.tx_biotype from exon ",
                   "join tx2exon on (tx2exon.exon_id=exon.exon_id) join ",
                   "tx on (tx.tx_id=tx2exon.tx_id) join ",
                   "gene on (gene.gene_id=tx.gene_id)")
    checkEquals(Q, want)
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
        checkEquals(Q, want)
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
        checkEquals(Q, want)
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
        checkEquals(Q, want)
        ## join with tx.
        Q <- ensembldb:::.buildQuery(edb, columns = c("tx_id", "protein_id",
                                                      "uniprot_id", "gene_id"))
        want <- paste0("select distinct tx.tx_id,protein.protein_id,",
                       "uniprot.uniprot_id,gene.gene_id from gene join ",
                       "tx on (gene.gene_id=tx.gene_id) left outer join protein",
                       " on (tx.tx_id=protein.tx_id) left outer join uniprot on",
                       " (protein.protein_id=uniprot.protein_id)")
        checkEquals(Q, want)
        ## if we started from protein:
        Q <- ensembldb:::.buildQuery(edb, columns = c("tx_id", "protein_id",
                                                      "uniprot_id", "gene_id"),
                                     startWith = "protein")
        want <- paste0("select distinct tx.tx_id,protein.protein_id,",
                       "uniprot.uniprot_id,gene.gene_id from protein left outer",
                       " join tx on (tx.tx_id=protein.tx_id) join gene on",
                       " (gene.gene_id=tx.gene_id) left outer join uniprot on",
                       " (protein.protein_id=uniprot.protein_id)")
        checkEquals(Q, want)
    }
}

############################################################
## Test the new join engine.
## o use the startWith argument.
## o change the join argument.
test_joinTwoTables <- function() {
    ## Check errors:
    checkException(ensembldb:::joinTwoTables(a = "gene", b = "dont exist"))
    checkException(ensembldb:::joinTwoTables(a = c("a", "b"), b = "gene"))
    ## Working example:
    res <- ensembldb:::joinTwoTables(a = c("a", "gene"), b = "tx")
    checkEquals(sort(res[1:2]), c("gene", "tx"))
    checkEquals(res[3], "on (gene.gene_id=tx.gene_id)")
    ## Error
    checkException(ensembldb:::joinTwoTables(a = "tx", b = "exon"))
    ## Working example:
    res <- ensembldb:::joinTwoTables(a = c("tx"), b = c("exon", "tx2exon"))
    checkEquals(sort(res[1:2]), c("tx", "tx2exon"))
    checkEquals(res[3], "on (tx.tx_id=tx2exon.tx_id)")
    res <- ensembldb:::joinTwoTables(a = c("chromosome", "gene", "tx"),
                                     b = c("exon", "protein", "tx2exon"))
    checkEquals(sort(res[1:2]), c("tx", "tx2exon"))
    checkEquals(res[3], "on (tx.tx_id=tx2exon.tx_id)")
}

test_joinQueryOnTables2_joinQueryOnColumns2 <- function() {
    ## exceptions
    checkException(ensembldb:::joinQueryOnTables2(edb, tab = c("a", "exon")))
    res <- ensembldb:::joinQueryOnTables2(edb, tab = c("gene", "exon"))
    want <- paste0("gene join tx on (gene.gene_id=tx.gene_id) join",
                   " tx2exon on (tx.tx_id=tx2exon.tx_id) join",
                   " exon on (tx2exon.exon_id=exon.exon_id)")
    ## The "default" order is gene->tx->tx2exon->exon
    checkEquals(res, want)
    res <- ensembldb:::joinQueryOnColumns2(edb, columns = c("exon_seq_start",
                                                            "gene_name"))
    checkEquals(res, want)
    ## Same but in the order: exon->tx2exon->tx->gene
    res <- ensembldb:::joinQueryOnTables2(edb, tab = c("gene", "exon"),
                                          startWith = "exon")
    want <- paste0("exon join tx2exon on (tx2exon.exon_id=exon.exon_id)",
                   " join tx on (tx.tx_id=tx2exon.tx_id) join",
                   " gene on (gene.gene_id=tx.gene_id)")
    checkEquals(res, want)
    res <- ensembldb:::joinQueryOnColumns2(edb, columns = c("exon_seq_start",
                                                            "gene_name"),
                                           startWith = "exon")
    checkEquals(res, want)
    ## That would be less expensive, but with "startWith" we force it to start
    ## from table exon, instead of just using tx2exon and tx.
    res <- ensembldb:::joinQueryOnColumns2(edb, columns = c("exon_id",
                                                            "gene_id"),
                                           startWith = "exon")
    checkEquals(res, want)
    ## Check proteins too.
    if (hasProteinData(edb)) {
        res <- ensembldb:::joinQueryOnTables2(edb, tab = c("protein", "gene",
                                                           "exon"))
        ## That should be: gene->tx->tx2exon->exon->protein
        want <- paste0("gene join tx on (gene.gene_id=tx.gene_id) join",
                       " tx2exon on (tx.tx_id=tx2exon.tx_id) join",
                       " exon on (tx2exon.exon_id=exon.exon_id) left outer join",
                       " protein on (tx.tx_id=protein.tx_id)")
        checkEquals(res, want)
        res <- ensembldb:::joinQueryOnColumns2(edb,
                                               columns = c("protein_id",
                                                           "gene_name",
                                                           "exon_seq_start"))
        checkEquals(res, want)
        res <- ensembldb:::joinQueryOnTables2(edb, tab = c("protein", "gene"),
                                              startWith = "protein")
        want <- paste0("protein left outer join tx on (tx.tx_id=protein.tx_id)",
                       " join gene on (gene.gene_id=tx.gene_id)")
        checkEquals(res, want)
        res <- ensembldb:::joinQueryOnColumns2(edb, columns = c("protein_id",
                                                                "gene_name"),
                                               startWith = "protein")
        checkEquals(res, want)
    }
}

## This test is an important one as it checks that we don't miss any entries
## from the database, e.g. if we query gene and join with protein that we don't
## miss any non-coding transcripts, or if we join protein with uniprot or
## protein_domain that we don't miss any values.
test_query_validity <- function() {
    ## Check RNA/DNA tables; shouldn't be a problem there, though.
    Ygns <- genes(edb, filter = SeqNameFilter("Y"), return.type = "data.frame")
    Ytxs <- transcripts(edb, filter = SeqNameFilter("Y"),
                        return.type = "data.frame",
                        columns = c("gene_id", "tx_id", "tx_biotype"))
    Yexns <- exons(edb, filter = SeqNameFilter("Y"), return.type = "data.frame",
                   columns = c("exon_id", "gene_id"))
    checkTrue(all(unique(Ygns$gene_id) %in% unique(Yexns$gene_id)))
    checkTrue(all(unique(Ygns$gene_id) %in% unique(Ytxs$gene_id)))
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
        checkTrue(length(unique(gns_f$gene_id)) < length(unique(Ygns$gene_id)))
        checkTrue(all(unique(Ytxs[Ytxs$tx_biotype == "protein_coding", "tx_id"])
                      %in% unique(gns_f$tx_id)))
        ## Now test the "real" query:
        Ygns_2 <- genes(edb, filter = SeqNameFilter("Y"),
                        return.type = "data.frame",
                        columns = c("gene_id", "tx_id", "tx_biotype",
                                    "protein_id"))
        ## We expect that ALL genes are present and ALL tx:
        checkTrue(all(unique(Ygns$gene_id) %in% unique(Ygns_2$gene_id)))
        checkTrue(all(unique(Ygns$tx_id) %in% unique(Ygns_2$tx_id)))

        ## Get all the tx with protein_id
        txs <- transcripts(edb, columns = c("tx_id", "protein_id"),
                           return.type = "data.frame")
        txids <- dbGetQuery(dbconn(edb), "select tx_id from tx;")[, "tx_id"]
        protids <- dbGetQuery(dbconn(edb),
                              "select protein_id from protein;")[, "protein_id"]
        checkTrue(all(txids %in% txs$tx_id))
        checkTrue(all(protids %in% txs$protein_id))

        ## Check protein with uniprot
        uniprotids <- dbGetQuery(dbconn(edb),
                                 "select uniprot_id from uniprot")$uniprot_id

        ## Check protein with protein domain
        ## Check protein_domain with uniprot
    }
}


## Compare the performance of inner join with left outer join.
dontrun_test_outer_join_performance <- function() {
    Q_1 <- ensembldb:::joinQueryOnTables2(edb, tab = c("gene", "exon"))
    Q_2 <- ensembldb:::joinQueryOnTables2(edb, tab = c("gene", "exon"),
                                          startWith = "exon")
    Q_3 <- ensembldb:::joinQueryOnTables2(edb, tab = c("gene", "exon"),
                                          startWith = "exon",
                                          join = "left outer join")
    library(microbenchmark)
    library(RSQLite)
    microbenchmark(dbGetQuery(dbconn(edb), paste0("select * from ", Q_1)),
                   dbGetQuery(dbconn(edb), paste0("select * from ", Q_2)),
                   dbGetQuery(dbconn(edb), paste0("select * from ", Q_3)),
                   times = 10)
    ## Result: Q_1 is a second faster (13 instead of 14).
    ## Check performance joining tx and genes.
    Q_1 <- ensembldb:::joinQueryOnTables2(edb, tab = c("tx", "gene"))
    Q_2 <- ensembldb:::joinQueryOnTables2(edb, tab = c("tx", "gene"),
                                          startWith = "tx")
    Q_3 <- ensembldb:::joinQueryOnTables2(edb, tab = c("tx", "gene"),
                                          startWith = "tx",
                                          join = "left outer join")
    microbenchmark(dbGetQuery(dbconn(edb), paste0("select * from ", Q_1)),
                   dbGetQuery(dbconn(edb), paste0("select * from ", Q_2)),
                   dbGetQuery(dbconn(edb), paste0("select * from ", Q_3)),
                   times = 10)
    ## No difference.
}
