############################################################
## Tests related to protein data.
library(EnsDb.Hsapiens.v75)
edb <- EnsDb.Hsapiens.v75

test_listProteinColumns <- function() {
    if (hasProteinData(edb)) {
        res <- listProteinColumns(edb)
        checkTrue(any(res == "protein_id"))
        checkTrue(any(res == "uniprot_id"))
        checkTrue(any(res == "protein_domain_id"))
        ## That's new columns fetched for Uniprot:
        checkTrue(any(res == "uniprot_db"))
        checkTrue(any(res == "uniprot_mapping_type"))
    } else {
        checkException(listProteinColumns(edb))
    }
}

############################################################
## Getting protein data in other methods.
test_genes_with_proteins <- function() {
    if (hasProteinData(edb)) {
        res <- genes(edb, columns = c("gene_name", "gene_id", "protein_id",
                                      "uniprot_id", "tx_id", "tx_biotype"),
                     filter = GenenameFilter("ZBTB16"),
                     return.type = "data.frame")
        ## have a 1:n mapping of protein_id to uniprot id:
        checkTrue(length(unique(res$protein_id)) <
                  nrow(res))
        checkEquals(colnames(res), c("gene_name", "gene_id", "protein_id",
                                     "uniprot_id", "tx_id", "tx_biotype"))
        ## All protein_coding have an uniprot_id
        checkTrue(all(!is.na(res[res$tx_biotype == "protein_coding",
                                 "uniprot_id"])))
        ## combine with cdsBy:
        cds <- cdsBy(edb, columns = c("tx_biotype", "protein_id"),
                     filter = GenenameFilter("ZBTB16"))
        codingTx <- unique(res[!is.na(res$protein_id), "tx_id"])
        checkEquals(sort(names(cds)), sort(codingTx))
        ## Next one fetching also protein domain data.
        res <- genes(edb, columns = c("gene_name", "tx_id", "protein_id",
                                      "protein_domain_id"),
                     filter = GenenameFilter("ZBTB16"),
                     return.type = "data.frame")
        checkEquals(colnames(res), c("gene_name", "tx_id", "protein_id",
                                     "protein_domain_id", "gene_id"))
        checkTrue(nrow(res) > length(unique(res$protein_id)))
        checkTrue(nrow(res) > length(unique(res$tx_id)))
    }
}

## transcripts
test_transcripts_with_proteins <- function() {
    if (hasProteinData(edb)) {
        res <- transcripts(edb, columns = c("tx_biotype", "protein_id",
                                            "uniprot_id"),
                           filter = TxidFilter("ENST00000335953"),
                           return.type = "data.frame")
        ## 1:1 mapping for tx_id <-> protein_id
        checkTrue(nrow(unique(res[, c("tx_id", "protein_id")])) == 1)
        ## Mapping tx_id -> uniprot_id is (0,1):n
        checkTrue(nrow(res) > length(unique(res$tx_id)))
        ## Add protein domains.
        res <- transcripts(edb, columns = c("tx_biotype", "protein_id",
                                            "uniprot_id",
                                            "protein_domain_id"),
                           filter = TxidFilter("ENST00000335953"),
                           return.type = "data.frame")
        resL <- split(res, f = res$uniprot_id)
        ## All have the same protein domains:
        resM <- do.call(rbind, lapply(resL, function(z) z$protein_domain_id))
        checkEquals(nrow(unique(resM)), 1)
    }
}

## exons
test_exons_with_proteins <- function() {
    if (hasProteinData(edb)) {
        ## Check if a call that includes a protein_id returns same data than one
        ## without.
        exns <- exons(edb, filter = GenenameFilter("BCL2L11"),
                      return.type = "data.frame")
        exns_2 <- exons(edb, filter = GenenameFilter("BCL2L11"),
                        return.type = "data.frame",
                        columns = c("exon_id", "protein_id"))
        checkEquals(sort(unique(exns$exon_id)),
                    sort(unique(exns_2$exon_id)))
        ## ZBTB16
        exns <- exons(edb, filter = GenenameFilter("ZBTB16"),
                      return.type = "data.frame")
        exns_2 <- exons(edb, filter = GenenameFilter("ZBTB16"),
                        return.type = "data.frame",
                        columns = c("exon_id", "protein_id"))
        checkEquals(sort(unique(exns$exon_id)),
                    sort(unique(exns_2$exon_id)))
    }
}

## exonsBy
test_exonsBy_with_proteins <- function() {
    if (hasProteinData(edb)) {
        exns <- exonsBy(edb, filter = GenenameFilter("ZBTB16"),
                        columns = "tx_biotype")
        exns_2 <- exonsBy(edb, filter = GenenameFilter("ZBTB16"),
                          columns = c("protein_id", "tx_biotype"))
        checkEquals(names(exns), names(exns_2))
        exns <- unlist(exns)
        exns_2 <- unlist(exns_2)
        checkTrue(any(is.na(exns_2$protein_id)))
        checkEquals(exns$exon_id, exns_2$exon_id)
    }
}

## transcriptsBy
test_transcriptsBy_with_proteins <- function() {
    if (hasProteinData(edb)) {
        txs <- transcriptsBy(edb, filter = GenenameFilter("ZBTB16"),
                             columns = "gene_biotype")
        txs_2 <- transcriptsBy(edb, filter = GenenameFilter("ZBTB16"),
                               columns = c("protein_id", "gene_biotype"))
        checkEquals(names(txs), names(txs_2))
        txs <- unlist(txs)
        txs_2 <- unlist(txs_2)
        checkTrue(any(is.na(txs_2$protein_id)))
        checkEquals(start(txs), start(txs_2))
    }
}

## cdsBy
test_cdsBy_with_proteins <- function() {
    if (hasProteinData(edb)) {
        cds <- cdsBy(edb, filter = GenenameFilter("ZBTB16"),
                     columns = "gene_biotype")
        cds_2 <- cdsBy(edb, filter = GenenameFilter("ZBTB16"),
                       columns = c("protein_id", "gene_biotype"))
        checkEquals(names(cds), names(cds_2))
        cds <- unlist(cds)
        cds_2 <- unlist(cds_2)
        checkTrue(all(!is.na(cds_2$protein_id)))
        checkEquals(start(cds), start(cds_2))
    }
}

## fiveUTRsByTranscript
test_fiveUTRsByTranscript_with_proteins <- function() {
    if (hasProteinData(edb)) {
        utrs <- fiveUTRsByTranscript(edb, filter = GenenameFilter("ZBTB16"),
                                     columns = "tx_biotype")
        utrs_2 <- fiveUTRsByTranscript(edb, filter = GenenameFilter("ZBTB16"),
                                       columns = c("protein_id", "gene_biotype"))
        checkEquals(names(utrs), names(utrs_2))
        utrs <- unlist(utrs)
        utrs_2 <- unlist(utrs_2)
        checkTrue(all(!is.na(utrs_2$protein_id)))
        checkEquals(start(utrs), start(utrs_2))
    }
}

############################################################
## Tests using protein filters
test_genes_with_protein_filters <- function() {
    ## o ProteinidFilter
    pif <- ProteinidFilter("ENSP00000376721")
    if (hasProteinData(edb)) {
        gns <- genes(edb, filter = pif, return.type = "data.frame")
        checkEquals(gns$gene_name, "ZBTB16")
    }
    ## o UniprotidFilter
    uif <- UniprotidFilter("Q71UL7_HUMAN")
    if (hasProteinData(edb)) {
        gns <- genes(edb, filter = uif, return.type = "data.frame",
                     columns = c("protein_id", "gene_name", "tx_id"))
        checkTrue("ENSP00000376721" %in% gns$protein_id)
        checkTrue(nrow(gns) == 2)
    }
    ## o ProtdomidFilter
    pdif <- ProtdomidFilter("PF00096")
    if (hasProteinData(edb)) {
        gns <- genes(edb, filter = list(pdif, GenenameFilter("ZBTB%", "like")),
                     return.type = "data.frame",
                     column = c("gene_name", "gene_biotype"))
        checkTrue(all(gns$gene_biotype == "protein_coding"))
    }
}

############################################################
## Using protein data based filters.
test_ProteinidFilter <- function() {
    pf <- ProteinidFilter("ABC")
    checkEquals(value(pf), "ABC")
    checkEquals(column(pf), "protein_id")
    checkEquals(where(pf), "protein_id = 'ABC'")
    if (hasProteinData(edb)) {
        checkEquals(column(pf, edb), "protein.protein_id")
        checkEquals(column(pf, edb, with.tables = "protein_domain"),
                    "protein_domain.protein_id")
        checkEquals(column(pf, edb, with.tables = "uniprot"),
                    "uniprot.protein_id")
        checkEquals(where(pf, edb), "protein.protein_id = 'ABC'")
        checkEquals(where(pf, edb, with.tables = "uniprot"),
                    "uniprot.protein_id = 'ABC'")
    } else {
        checkException(column(pf, edb))
        checkException(where(pf, edb))
        checkException(column(pf, edb, with.tables = "uniprot"))
        checkException(where(pf, edb, with.tables = "uniprot"))
    }
    pf <- ProteinidFilter(c("A", "B"))
    checkEquals(where(pf), "protein_id in ('A','B')")
    checkException(ProteinidFilter("B", condition = ">"))
}

test_UniprotidFilter <- function() {
    pf <- UniprotidFilter("ABC")
    checkEquals(value(pf), "ABC")
    checkEquals(column(pf), "uniprot_id")
    checkEquals(where(pf), "uniprot_id = 'ABC'")
    if (hasProteinData(edb)) {
        checkEquals(column(pf, edb), "uniprot.uniprot_id")
        checkEquals(column(pf, edb, with.tables = "uniprot"),
                    "uniprot.uniprot_id")
        checkEquals(where(pf, edb), "uniprot.uniprot_id = 'ABC'")
        checkEquals(where(pf, edb, with.tables = "uniprot"),
                    "uniprot.uniprot_id = 'ABC'")
    } else {
        checkException(column(pf, edb))
        checkException(where(pf, edb))
        checkException(column(pf, edb, with.tables = "uniprot"))
        checkException(where(pf, edb, with.tables = "uniprot"))
    }
    pf <- UniprotidFilter(c("A", "B"))
    checkEquals(where(pf), "uniprot_id in ('A','B')")
    checkException(UniprotidFilter("B", condition = ">"))
}

test_ProtdomidFilter <- function() {
    pf <- ProtdomidFilter("ABC")
    checkEquals(value(pf), "ABC")
    checkEquals(column(pf), "protein_domain_id")
    checkEquals(where(pf), "protein_domain_id = 'ABC'")
    if (hasProteinData(edb)) {
        checkEquals(column(pf, edb), "protein_domain.protein_domain_id")
        checkEquals(column(pf, edb, with.tables = "protein_domain"),
                    "protein_domain.protein_domain_id")
        checkEquals(where(pf, edb), "protein_domain.protein_domain_id = 'ABC'")
        checkEquals(where(pf, edb, with.tables = "protein_domain"),
                    "protein_domain.protein_domain_id = 'ABC'")
    } else {
        checkException(column(pf, edb))
        checkException(where(pf, edb))
        checkException(column(pf, edb, with.tables = "protein_domain"))
        checkException(where(pf, edb, with.tables = "protein_domain"))
    }
    pf <- ProtdomidFilter(c("A", "B"))
    checkEquals(where(pf), "protein_domain_id in ('A','B')")
    checkException(ProtdomidFilter("B", condition = ">"))
}

test_UniprotdbFilter <- function() {
    pf <- UniprotdbFilter("ABC")
    checkEquals(value(pf), "ABC")
    checkEquals(column(pf), "uniprot_db")
    checkEquals(where(pf), "uniprot_db = 'ABC'")
    if (hasProteinData(edb)) {
        checkEquals(column(pf, edb), "uniprot.uniprot_db")
        checkEquals(column(pf, edb, with.tables = "uniprot"),
                    "uniprot.uniprot_db")
        checkEquals(where(pf, edb), "uniprot.uniprot_db = 'ABC'")
        checkEquals(where(pf, edb, with.tables = "uniprot"),
                    "uniprot.uniprot_db = 'ABC'")
    } else {
        checkException(column(pf, edb))
        checkException(where(pf, edb))
        checkException(column(pf, edb, with.tables = "uniprot"))
        checkException(where(pf, edb, with.tables = "uniprot"))
    }
    pf <- UniprotdbFilter(c("A", "B"))
    checkEquals(where(pf), "uniprot_db in ('A','B')")
    checkException(UniprotdbFilter("B", condition = ">"))
}

test_UniprotmappingtypeFilter <- function() {
    pf <- UniprotmappingtypeFilter("ABC")
    checkEquals(value(pf), "ABC")
    checkEquals(column(pf), "uniprot_mapping_type")
    checkEquals(where(pf), "uniprot_mapping_type = 'ABC'")
    if (hasProteinData(edb)) {
        checkEquals(column(pf, edb), "uniprot.uniprot_mapping_type")
        checkEquals(column(pf, edb, with.tables = "uniprot"),
                    "uniprot.uniprot_mapping_type")
        checkEquals(where(pf, edb), "uniprot.uniprot_mapping_type = 'ABC'")
        checkEquals(where(pf, edb, with.tables = "uniprot"),
                    "uniprot.uniprot_mapping_type = 'ABC'")
    } else {
        checkException(column(pf, edb))
        checkException(where(pf, edb))
        checkException(column(pf, edb, with.tables = "uniprot"))
        checkException(where(pf, edb, with.tables = "uniprot"))
    }
    pf <- UniprotmappingtypeFilter(c("A", "B"))
    checkEquals(where(pf), "uniprot_mapping_type in ('A','B')")
    checkException(UniprotmappingtypeFilter("B", condition = ">"))
}



############################################################
## The dedicated methods to fetch protein data.
test_proteins <- function() {
    if (hasProteinData(edb)) {
        ## Check return type.
        prts_DF <- proteins(edb, filter = GenenameFilter("ZBTB16"))
        checkTrue(is(prts_DF, "DataFrame"))
        prts_df <- proteins(edb, filter = GenenameFilter("ZBTB16"),
                            return.type = "data.frame")
        checkTrue(is(prts_df, "data.frame"))
        prts_aa <- proteins(edb, filter = GenenameFilter("ZBTB16"),
                            return.type = "AAStringSet")
        checkTrue(is(prts_aa, "AAStringSet"))
        ## Check content.
        library(RSQLite)
        res_q <- dbGetQuery(dbconn(edb),
                            paste0("select tx.tx_id, protein_id, gene_name from ",
                                   "protein left outer join tx on (protein.tx_id=",
                                   "tx.tx_id) join gene on (gene.gene_id=",
                                   "tx.gene_id) where gene_name = 'ZBTB16'"))
        checkEquals(res_q$tx_id, prts_df$tx_id)
        checkEquals(res_q$protein_id, prts_df$protein_id)
        checkEquals(prts_df$protein_id, names(prts_aa))

        ## Add protein domain information to the proteins.
        prts_df <- proteins(edb, filter = ProteinidFilter(c("ENSP00000338157",
                                                            "ENSP00000443013")),
                            columns = c("protein_id", "protein_domain_id",
                                        "uniprot_id"),
                            return.type = "data.frame")
        ## Check if we have all data that we expect:
        uniprots <- dbGetQuery(dbconn(edb),
                               paste0("select uniprot_id from uniprot where",
                                      " protein_id in ('ENSP00000338157',",
                                      "'ENSP00000443013')"))$uniprot_id
        checkTrue(all(uniprots %in% prts_df$uniprot_id))
        protdoms <- dbGetQuery(dbconn(edb),
                               paste0("select protein_domain_id from",
                                      " protein_domain where protein_id",
                                      " in ('ENSP00000338157',",
                                      "'ENSP00000443013')"))$protein_domain_id
        checkTrue(all(protdoms %in% prts_df$protein_domain_id))
    }
}

## Testing protein to uniprot mappings.
test_proteins_uniprot <- function() {
    ## ZBTB16 and the mapping of 1 protein to two Uniprot IDs, one with DIRECT
    ## mapping type.
}

test_isProteinFilter <- function() {
    ## TRUE
    checkTrue(ensembldb:::isProteinFilter(ProteinidFilter("a")))
    checkTrue(ensembldb:::isProteinFilter(UniprotidFilter("a")))
    checkTrue(ensembldb:::isProteinFilter(ProtdomidFilter("a")))
    checkTrue(ensembldb:::isProteinFilter(UniprotdbFilter("a")))
    checkTrue(ensembldb:::isProteinFilter(UniprotmappingtypeFilter("a")))
    ## FALSE
    checkTrue(!ensembldb:::isProteinFilter(GeneidFilter("a")))
    checkTrue(!ensembldb:::isProteinFilter(SymbolFilter("a")))
    checkTrue(!ensembldb:::isProteinFilter(3))
    checkTrue(!ensembldb:::isProteinFilter("dfdf"))
}

notrun_test_protein_domains <- function() {
    res <- ensembldb:::getWhat(edb, columns = c("protein_id", "tx_id", "gene_id",
                                                "gene_name"),
                               filter = list(ProtdomidFilter("PF00096")))
}

## test_ProteinsFromDataframe <- function() {
##     if (hasProteinData(edb)) {
##         prt_df <- proteins(edb, filter = GenenameFilter("ZBTB16"),
##                            return.type = "data.frame")
##         prt <- ensembldb:::.ProteinsFromDataframe(prt_df)
##     }
## }



## What would be nice for Proteins: be a little more like a GRanges.
## o Access mcols of the aa directly.
## o Access columns from the metadata with $
## o inherit directly from AAStringSet?
## o names() returns result from seqnames?
## o show method be a little more like a GRanges.
## o do 'names' really have to be unique? what with 1:n mapping between
##   protein_id and Uniprot ID? would be nice if the binding would NOT be by
##   name but rather by index!

## How should a result object look like:
## ProteinResult:
## o extend an AAStringSet (that's fine!)
## o transcript_id is required.
## o domain data as IRangesList
## Functionality:
## o getGenomeMapping get a GRanges with the cds for this (cdsBy tx).

