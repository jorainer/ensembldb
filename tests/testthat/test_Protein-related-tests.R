
############################################################
## Getting protein data in other methods.
test_that("genes works with proteins", {
    if (hasProteinData(edb)) {
        res <- genes(edb, columns = c("gene_name", "gene_id", "protein_id",
                                      "uniprot_id", "tx_id", "tx_biotype"),
                     filter = GenenameFilter("ZBTB16"),
                     return.type = "data.frame")
        ## have a 1:n mapping of protein_id to uniprot id:
        expect_true(length(unique(res$protein_id)) <
                  nrow(res))
        expect_equal(colnames(res), c("gene_name", "gene_id", "protein_id",
                                     "uniprot_id", "tx_id", "tx_biotype"))
        ## All protein_coding have an uniprot_id
        expect_true(all(!is.na(res[res$tx_biotype == "protein_coding",
                                 "uniprot_id"])))
        ## combine with cdsBy:
        cds <- cdsBy(edb, columns = c("tx_biotype", "protein_id"),
                     filter = GenenameFilter("ZBTB16"))
        codingTx <- unique(res[!is.na(res$protein_id), "tx_id"])
        expect_equal(sort(names(cds)), sort(codingTx))
        ## Next one fetching also protein domain data.
        res <- genes(edb, columns = c("gene_name", "tx_id", "protein_id",
                                      "protein_domain_id"),
                     filter = GenenameFilter("ZBTB16"),
                     return.type = "data.frame")
        expect_equal(colnames(res), c("gene_name", "tx_id", "protein_id",
                                     "protein_domain_id", "gene_id"))
        expect_true(nrow(res) > length(unique(res$protein_id)))
        expect_true(nrow(res) > length(unique(res$tx_id)))
    }
})

test_that("transcripts works with proteins", {
    if (hasProteinData(edb)) {
        res <- transcripts(edb, columns = c("tx_biotype", "protein_id",
                                            "uniprot_id"),
                           filter = TxIdFilter("ENST00000335953"),
                           return.type = "data.frame")
        ## 1:1 mapping for tx_id <-> protein_id
        expect_true(nrow(unique(res[, c("tx_id", "protein_id")])) == 1)
        ## Mapping tx_id -> uniprot_id is (0,1):n
        expect_true(nrow(res) > length(unique(res$tx_id)))
        ## Add protein domains.
        res <- transcripts(edb, columns = c("tx_biotype", "protein_id",
                                            "uniprot_id",
                                            "protein_domain_id"),
                           filter = TxIdFilter("ENST00000335953"),
                           return.type = "data.frame")
        resL <- split(res, f = res$uniprot_id)
        ## All have the same protein domains:
        resM <- do.call(rbind, lapply(resL, function(z) z$protein_domain_id))
        expect_equal(nrow(unique(resM)), 1)
    }
})

## exons
test_that("exons works with proteins",  {
    if (hasProteinData(edb)) {
        ## Check if a call that includes a protein_id returns same data than one
        ## without.
        exns <- exons(edb, filter = GenenameFilter("BCL2L11"),
                      return.type = "data.frame")
        exns_2 <- exons(edb, filter = GenenameFilter("BCL2L11"),
                        return.type = "data.frame",
                        columns = c("exon_id", "protein_id"))
        expect_equal(sort(unique(exns$exon_id)),
                    sort(unique(exns_2$exon_id)))
        ## ZBTB16
        exns <- exons(edb, filter = GenenameFilter("ZBTB16"),
                      return.type = "data.frame")
        exns_2 <- exons(edb, filter = GenenameFilter("ZBTB16"),
                        return.type = "data.frame",
                        columns = c("exon_id", "protein_id"))
        expect_equal(sort(unique(exns$exon_id)),
                    sort(unique(exns_2$exon_id)))
    }
})

## exonsBy
test_that("exonsBy works with proteins", {
    if (hasProteinData(edb)) {
        exns <- exonsBy(edb, filter = GenenameFilter("ZBTB16"),
                        columns = "tx_biotype")
        exns_2 <- exonsBy(edb, filter = GenenameFilter("ZBTB16"),
                          columns = c("protein_id", "tx_biotype"))
        expect_equal(names(exns), names(exns_2))
        exns <- unlist(exns)
        exns_2 <- unlist(exns_2)
        expect_true(any(is.na(exns_2$protein_id)))
        expect_equal(exns$exon_id, exns_2$exon_id)
    }
})

## transcriptsBy
test_that("transcriptsBy works with proteins", {
    if (hasProteinData(edb)) {
        txs <- transcriptsBy(edb, filter = GenenameFilter("ZBTB16"),
                             columns = "gene_biotype")
        txs_2 <- transcriptsBy(edb, filter = GenenameFilter("ZBTB16"),
                               columns = c("protein_id", "gene_biotype"))
        expect_equal(names(txs), names(txs_2))
        txs <- unlist(txs)
        txs_2 <- unlist(txs_2)
        expect_true(any(is.na(txs_2$protein_id)))
        expect_equal(start(txs), start(txs_2))
    }
})

## cdsBy
test_that("cdsBy works with proteins", {
    if (hasProteinData(edb)) {
        cds <- cdsBy(edb, filter = GenenameFilter("ZBTB16"),
                     columns = "gene_biotype")
        cds_2 <- cdsBy(edb, filter = GenenameFilter("ZBTB16"),
                       columns = c("protein_id", "gene_biotype"))
        expect_equal(names(cds), names(cds_2))
        cds <- unlist(cds)
        cds_2 <- unlist(cds_2)
        expect_true(all(!is.na(cds_2$protein_id)))
        expect_equal(start(cds), start(cds_2))
    }
})

## fiveUTRsByTranscript
test_that("fiveUTRsByTranscript works with proteins", {
    if (hasProteinData(edb)) {
        utrs <- fiveUTRsByTranscript(edb, filter = GenenameFilter("ZBTB16"),
                                     columns = "tx_biotype")
        utrs_2 <- fiveUTRsByTranscript(edb, filter = GenenameFilter("ZBTB16"),
                                       columns = c("protein_id", "gene_biotype"))
        expect_equal(names(utrs), names(utrs_2))
        utrs <- unlist(utrs)
        utrs_2 <- unlist(utrs_2)
        expect_true(all(!is.na(utrs_2$protein_id)))
        expect_equal(start(utrs), start(utrs_2))
    }
})

test_that("genes works with protein filters", {
    ## o ProteinIdFilter
    pif <- ProteinIdFilter("ENSP00000376721")
    if (hasProteinData(edb)) {
        gns <- genes(edb, filter = pif, return.type = "data.frame")
        expect_equal(gns$gene_name, "ZBTB16")
    }
    ## o UniprotFilter
    uif <- UniprotFilter("Q05516")
    if (hasProteinData(edb)) {
        gns <- genes(edb, filter = uif, return.type = "data.frame",
                     columns = c("protein_id", "gene_name", "tx_id"))
        expect_true("ENSP00000376721" %in% gns$protein_id)
        expect_true(nrow(gns) == 2)
    }
    ## o ProtDomIdFilter
    pdif <- ProtDomIdFilter("PF00096")
    if (hasProteinData(edb)) {
        gns <- genes(edb, filter = list(pdif,
                                        GenenameFilter("ZBTB%", "startsWith")),
                     return.type = "data.frame",
                     column = c("gene_name", "gene_biotype"))
        expect_true(all(gns$gene_biotype == "protein_coding"))
    }
})

test_that("proteins works", {
    if (hasProteinData(edb)) {
        ## Check return type.
        prts_DF <- proteins(edb, filter = GenenameFilter("ZBTB16"))
        expect_true(is(prts_DF, "DataFrame"))
        prts_df <- proteins(edb, filter = GenenameFilter("ZBTB16"),
                            return.type = "data.frame")
        expect_true(is(prts_df, "data.frame"))
        prts_aa <- proteins(edb, filter = GenenameFilter("ZBTB16"),
                            return.type = "AAStringSet")
        expect_true(is(prts_aa, "AAStringSet"))
        ## Check content.
        library(RSQLite)
        res_q <- dbGetQuery(
            dbconn(edb),
            paste0("select tx.tx_id, protein_id, gene_name from ",
                   "protein left outer join tx on (protein.tx_id=",
                   "tx.tx_id) join gene on (gene.gene_id=",
                   "tx.gene_id) where gene_name = 'ZBTB16'"))
        expect_equal(res_q$tx_id, prts_df$tx_id)
        expect_equal(res_q$protein_id, prts_df$protein_id)
        expect_equal(prts_df$protein_id, names(prts_aa))
        ## Add protein domain information to the proteins.
        prts_df <- proteins(edb, filter = ProteinIdFilter(c("ENSP00000338157",
                                                            "ENSP00000443013")),
                            columns = c("protein_id", "protein_domain_id",
                                        "uniprot_id"),
                            return.type = "data.frame")
        ## Check if we have all data that we expect:
        uniprots <- dbGetQuery(dbconn(edb),
                               paste0("select uniprot_id from uniprot where",
                                      " protein_id in ('ENSP00000338157',",
                                      "'ENSP00000443013')"))$uniprot_id
        expect_true(all(uniprots %in% prts_df$uniprot_id))
        protdoms <- dbGetQuery(dbconn(edb),
                               paste0("select protein_domain_id from",
                                      " protein_domain where protein_id",
                                      " in ('ENSP00000338157',",
                                      "'ENSP00000443013')"))$protein_domain_id
        expect_true(all(protdoms %in% prts_df$protein_domain_id))
    }
})

test_that("proteins works with uniprot mapping", {
    ## ZBTB16 and the mapping of 1 protein to two Uniprot IDs, one with DIRECT
    ## mapping type.
    if (hasProteinData(edb)) {
        prts <- proteins(edb, filter = GenenameFilter("ZBTB16"),
                         columns = c("uniprot_id", "uniprot_db",
                                     "uniprot_mapping_type"),
                         return.type = "DataFrame")
        ## NOTE: this is true for Ensembl 86, but might not be the case for 75!
        ## Here we have the n:m mapping:
        ## Q05516 is assigned to ENSP00000338157 and ENSP00000376721,
        ## Each of the two proteins is however also annotated to a second
        ## Uniprot ID: A0A024R3C6
        ## If we use the UniprotMappingTypeFilter with only DIRECT mapping
        ## we expect to reduce it to the 1:n mapping between Uniprot and Ensembl
        prts <- proteins(edb, filter = list(GenenameFilter("ZBTB16"),
                                            UniprotMappingTypeFilter("DIRECT")),
                         columns = c("uniprot_id", "uniprot_db",
                                     "uniprot_mapping_type"),
                         return.type = "DataFrame")
        expect_true(all(prts$uniprot_mapping_type == "DIRECT"))
        ## Check the UniprotDbFilter
        prts <- proteins(edb, filter = list(GenenameFilter("ZBTB16"),
                                            UniprotDbFilter("SPTREMBL")),
                         columns = c("uniprot_id", "uniprot_db", "protein_id"),
                         return.type = "DataFrame")
        expect_true(all(prts$uniprot_db == "SPTREMBL"))
    }
})

test_that("isProteinFilter works", {
    ## TRUE
    expect_true(ensembldb:::isProteinFilter(ProteinIdFilter("a")))
    expect_true(ensembldb:::isProteinFilter(UniprotFilter("a")))
    expect_true(ensembldb:::isProteinFilter(ProtDomIdFilter("a")))
    expect_true(ensembldb:::isProteinFilter(UniprotDbFilter("a")))
    expect_true(ensembldb:::isProteinFilter(UniprotMappingTypeFilter("a")))
    ## FALSE
    expect_true(!ensembldb:::isProteinFilter(GeneIdFilter("a")))
    expect_true(!ensembldb:::isProteinFilter(SymbolFilter("a")))
    expect_true(!ensembldb:::isProteinFilter(3))
    expect_true(!ensembldb:::isProteinFilter("dfdf"))
})

test_that("ProteinDomainId and ProtDomId work", {
    pdid <- "PS50063"
    res <- proteins(edb, filter = ~ protein_domain_id == pdid)
    res_2 <- proteins(edb, filter = ~ prot_dom_id == pdid)
    expect_equal(res, res_2)
    res <- proteins(edb, filter = ~ genename == "BCL2" &
                             protein_domain_source == "pfam")
    expect_equal(nrow(res), 3)
})
