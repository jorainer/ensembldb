############################################################
## Tests related to protein data.
library(EnsDb.Hsapiens.v75)
edb <- EnsDb.Hsapiens.v75

############################################################
## Getting protein data in other methods.
test_genes_with_proteins <- function() {
    suppressWarnings(
        res <- genes(edb, columns = c("gene_name", "gene_id", "protein_id",
                                      "uniprot_id"),
                     filter = GenenameFilter("ZBTB16"),
                     return.type = "data.frame")
    )
    if (hasProteinData(edb)) {
        ## have a 1:n mapping of protein_id to uniprot id:
        checkTrue(length(unique(res$protein_id)) <
                  nrow(res))
        checkEquals(colnames(res), c("gene_name", "gene_id", "protein_id",
                                     "uniprot_id"))
    } else {
        checkEquals(colnames(res), c("gene_name", "gene_id"))
        checkEquals(nrow(res), length(unique(res$gene_id)))
    }
    ## Next one fetching also protein domain data.
    suppressWarnings(
        res <- genes(edb, columns = c("gene_name", "tx_id", "protein_id",
                                      "protein_domain_id"),
                     filter = GenenameFilter("ZBTB16"),
                     return.type = "data.frame")
    )
    if (hasProteinData(edb)) {
        checkEquals(colnames(res), c("gene_name", "tx_id", "protein_id",
                                     "protein_domain_id", "gene_id"))
        checkTrue(nrow(res) > length(unique(res$protein_id)))
        checkTrue(nrow(res) > length(unique(res$tx_id)))
    } else {
        checkEquals(colnames(res), c("gene_name", "tx_id", "gene_id"))
        checkTrue(nrow(res) == length(unique(res$tx_id)))
    }
}

############################################################
## Using protein data based filters.

############################################################
## The dedicated methods to fetch protein data.
