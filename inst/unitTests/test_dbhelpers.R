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

## test_prefixColumns <- function() {
##     ## gene
##     res <- ensembldb:::prefixColumns(edb, "gene_id")
##     checkEquals(unname(unlist(res)), "gene.gene_id")
##     res <- ensembldb:::prefixColumns(edb, "gene_id",
##                                      with.tables = c("tx", "exon"))
##     checkEquals(unname(unlist(res)), "tx.gene_id")
##     ## tx
##     res <- ensembldb:::prefixColumns(edb, "tx_id")
##     checkEquals(unname(unlist(res)), "tx.tx_id")
##     res <- ensembldb:::prefixColumns(edb, "tx_id",
##                                      with.tables = "tx2exon")
##     checkEquals(unname(unlist(res)), "tx2exon.tx_id")
##     ## exon
##     res <- ensembldb:::prefixColumns(edb, "exon_id")
##     checkEquals(unname(unlist(res)), "tx2exon.exon_id")
##     res <- ensembldb:::prefixColumns(edb, "exon_id",
##                                      with.tables = c("exon"))
##     checkEquals(unname(unlist(res)), "exon.exon_id")

##     if (hasProteinData(edb)) {
##         res <- ensembldb:::prefixColumns(edb, "protein_id")
##         checkEquals(unname(unlist(res)), "protein.protein_id")
##         res <- ensembldb:::prefixColumns(edb, "protein_id",
##                                          with.tables = c("uniprot",
##                                                          "protein_domain"))
##         checkEquals(unname(unlist(res)), "uniprot.protein_id")
##     }
## }

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

## test_joinQueryOnTables <- function() {
##     res <- ensembldb:::joinQueryOnTables(edb, tab = c("exon", "gene"))
##     checkEquals(res, paste0("gene join tx on (gene.gene_id=tx.gene_id) ",
##                             "join tx2exon on (tx.tx_id=tx2exon.tx_id) ",
##                             "join exon on (tx2exon.exon_id=exon.exon_id)"))
##     res <- ensembldb:::joinQueryOnTables(edb, tab = c("gene", "chromosome"))
##     checkEquals(res, "gene join chromosome on (gene.seq_name=chromosome.seq_name)")
##     res <- ensembldb:::joinQueryOnTables(edb, tab = c("exon", "tx"))
##     checkEquals(res, paste0("tx join tx2exon on (tx.tx_id=tx2exon.tx_id) ",
##                             "join exon on (tx2exon.exon_id=exon.exon_id)"))
##     res <- ensembldb:::joinQueryOnTables(edb, tab = c("chromosome", "tx"))
##     checkEquals(res, paste0("gene join tx on (gene.gene_id=tx.gene_id) ",
##                             "join chromosome on (gene.seq_name=chromosome.seq_name)"))
##     if (hasProteinData(edb)) {
##         res <- ensembldb:::joinQueryOnTables(edb, tab = c("tx", "uniprot"))
##         checkEquals(res, paste0("tx join protein on (tx.tx_id=protein.tx_id) ",
##                                 "join uniprot on (protein.protein_id=uniprot.protein_id)"))
##         res <- ensembldb:::joinQueryOnTables(edb, tab = c("protein", "gene"))
##         checkEquals(res, paste0("gene join tx on (gene.gene_id=tx.gene_id) ",
##                                 "join protein on (tx.tx_id=protein.tx_id)"))
##         res <- ensembldb:::joinQueryOnTables(edb,
##                                              tab = c("uniprot", "protein_domain"))
##         checkEquals(res, "uniprot join protein_domain on (uniprot.protein_id=protein_domain.protein_id)")
##         res <- ensembldb:::joinQueryOnTables(edb, tab = c("uniprot", "exon"))
##         checkEquals(res, paste0("tx join tx2exon on (tx.tx_id=tx2exon.tx_id) ",
##                                 "join exon on (tx2exon.exon_id=exon.exon_id) ",
##                                 "join protein on (tx.tx_id=protein.tx_id) ",
##                                 "join uniprot on (protein.protein_id=uniprot.protein_id)"))
##     }
## }

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
                       "protein_domain.protein_domain_id from protein join ",
                       "protein_domain on (protein.protein_id=",
                       "protein_domain.protein_id) join ",
                       "uniprot on (protein.protein_id=uniprot.protein_id)")
        checkEquals(Q, want)
        ## start at protein
        Q <- ensembldb:::.buildQuery(edb,
                                     columns = c("protein_id", "uniprot_id",
                                                 "protein_domain_id"),
                                     startWith = "protein")
        want <- paste0("select distinct protein.protein_id,uniprot.uniprot_id,",
                       "protein_domain.protein_domain_id from protein join ",
                       "protein_domain on (protein.protein_id=",
                       "protein_domain.protein_id) join ",
                       "uniprot on (protein.protein_id=uniprot.protein_id)")
        checkEquals(Q, want)
        ## start at uniprot.
        Q <- ensembldb:::.buildQuery(edb,
                                     columns = c("protein_id", "uniprot_id",
                                                 "protein_domain_id"),
                                     startWith = "uniprot")
        want <- paste0("select distinct protein.protein_id,uniprot.uniprot_id,",
                       "protein_domain.protein_domain_id from uniprot join ",
                       "protein on (protein.protein_id=uniprot.protein_id) join",
                       " protein_domain on (protein.protein_id=",
                       "protein_domain.protein_id)")
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
                       " exon on (tx2exon.exon_id=exon.exon_id) join",
                       " protein on (tx.tx_id=protein.tx_id)")
        checkEquals(res, want)
        res <- ensembldb:::joinQueryOnColumns2(edb,
                                               columns = c("protein_id",
                                                           "gene_name",
                                                           "exon_seq_start"))
        checkEquals(res, want)
        res <- ensembldb:::joinQueryOnTables2(edb, tab = c("protein", "gene"),
                                              startWith = "protein",
                                              join = "left outer join")
        want <- paste0("protein left outer join tx on (tx.tx_id=protein.tx_id)",
                       " left outer join gene on (gene.gene_id=tx.gene_id)")
        checkEquals(res, want)
        res <- ensembldb:::joinQueryOnColumns2(edb, columns = c("protein_id",
                                                                "gene_name"),
                                               startWith = "protein",
                                               join = "left outer join")
        checkEquals(res, want)
    }
}

## This test is an important one as it checks that we don't miss any entries
## from the database, e.g. if we query gene and join with protein that we don't
## miss any non-coding transcripts, or if we join protein with uniprot or
## protein_domain that we don't miss any values.
test_outer_join_validity <- function() {
    ## Check gene with protein
    if (hasProteinData(edb)) {
        library(RSQLite)
        gns <- genes(edb, filter = SeqnameFilter("Y"),
                     return.type = "data.frame",
                     columns = c("gene_id", "tx_id", "tx_biotype"))
        gns_f <- dbGetQuery(dbconn(edb),
                            paste0("select gene.gene_id, tx.tx_id, tx_biotype, ",
                                   "protein_id from gene join tx on ",
                                   "(gene.gene_id=tx.gene_id) join protein on ",
                                   "(tx.tx_id=protein.tx_id) ",
                                   "where seq_name = 'Y'"))
        gns_2 <- genes(edb, filter = SeqnameFilter("Y"),
                       return.type = "data.frame",
                       columns = c("gene_id", "tx_id", "tx_biotype",
                                   "protein_id"))
        checkTrue(nrow(gns_f) < nrow(gns_2))
        checkTrue(all(unique(gns$gene_id) %in% unique(gns_2$gene_id)))

        ## Check transcript with protein
        txs <- transcripts(edb, filter = SeqnameFilter("Y"),
                           return.type = "data.frame",
                           columns = c("tx_id", "tx_biotype"))
        txs_f <- dbGetQuery(dbconn(edb),
                            paste0("select tx.tx_id, tx_biotype, ",
                                   "protein_id from ",
                                   "tx join gene on (tx.gene_id=gene.gene_id) ",
                                   "join protein on (tx.tx_id=protein.tx_id)",
                                   " where seq_name = 'Y'"))
        txs_2 <- transcripts(edb, filter = SeqnameFilter("Y"),
                             return.type = "data.frame",
                             columns = c("tx_id", "tx_biotype", "protein_id"))
        checkTrue(nrow(txs_f) < nrow(txs_2))
        checkTrue(all(unique(txs$tx_id) %in% unique(txs_2$tx_id)))
        ## Just getting tx_id and prote_id;
        txs_2 <- transcripts(edb, filter = SeqnameFilter("Y"),
                             return.type = "data.frame",
                             columns = c("tx_id", "protein_id"))
        checkTrue(all(unique(txs$tx_id) %in% unique(txs_2$tx_id)))

        ## Check protein with uniprot
        prts <- dbGetQuery(dbconn(edb), "select protein_id from protein")
        prts_2 <- dbGetQuery(dbconn(edb),
                             paste0("select protein.protein_id,uniprot_id from",
                                    " protein left outer join",
                                    " uniprot on (protein.protein_id=",
                                    "uniprot.protein_id)"))
        prts_f <- dbGetQuery(dbconn(edb),
                             paste0("select * from protein join",
                                    " uniprot on (protein.protein_id=",
                                    "uniprot.protein_id)"))
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
