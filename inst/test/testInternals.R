detachem <- function(x){
    NS <- loadedNamespaces()
    if(any(NS==x)){
        pkgn <- paste0("package:", x)
        detach(pkgn, unload=TRUE, character.only=TRUE)
    }
}
Pkgs <- c("EnsDb.Hsapiens.v75", "ensembldb")
tmp <- sapply(Pkgs, detachem)
tmp <- sapply(Pkgs, library, character.only=TRUE)
DB <- EnsDb.Hsapiens.v75


#######################################################
##
##  add required tables if needed.
##
## check if we get what we want...
Expect <- c("exon", "tx2exon", "tx")
Get <- ensembldb:::addRequiredTables(EnsDb.Hsapiens.v75, c("exon", "tx"))
Get
if(sum(Get %in% Expect)!=length(Expect))
    stop("Didn't get what I expected!")


Expect <- c("exon", "tx2exon", "tx", "gene")
Get <- ensembldb:::addRequiredTables(EnsDb.Hsapiens.v75, c("exon", "gene"))
Get
if(sum(Get %in% Expect)!=length(Expect))
    stop("Didn't get what I expected!")



Expect <- c("exon", "tx2exon", "tx", "gene")
Get <- ensembldb:::addRequiredTables(EnsDb.Hsapiens.v75, c("exon", "gene", "tx"))
Get
if(sum(Get %in% Expect)!=length(Expect))
    stop("Didn't get what I expected!")


#######################################################
##
##  join queries
##
ensembldb:::joinQueryOnTables(EnsDb.Hsapiens.v75, c("exon", "t2exon", "tx"))


ensembldb:::joinQueryOnTables(EnsDb.Hsapiens.v75, c("exon"))


ensembldb:::joinQueryOnTables(EnsDb.Hsapiens.v75, c("exon", "t2exon", "tx", "gene"))


ensembldb:::joinQueryOnTables(EnsDb.Hsapiens.v75, c("tx", "gene"))


ensembldb:::joinQueryOnTables(EnsDb.Hsapiens.v75, c("chromosome", "gene"))




#######################################################
##
##  join queries on column names
##
## for that query we don't need the exon table
ensembldb:::cleanColumns(EnsDb.Hsapiens.v75, c("gene_id","tx_id", "bla", "value"))

## don't require the exon table here, exon_id is also in tx2exon.
ensembldb:::joinQueryOnColumns(EnsDb.Hsapiens.v75, c("gene_id", "tx_id", "gene_name", "exon_id"))

##
ensembldb:::joinQueryOnColumns(EnsDb.Hsapiens.v75, c("gene_id", "tx_id", "gene_name", "exon_idx"))


ensembldb:::joinQueryOnColumns(EnsDb.Hsapiens.v75, c("gene_id", "tx_id", "gene_name", "exon_id", "exon_seq_start"))



#######################################################
##
##  clean columns
##
ensembldb:::cleanColumns(EnsDb.Hsapiens.v75, c("gene_id" ,"bma", "gene.gene_biotype"))

ensembldb:::cleanColumns(EnsDb.Hsapiens.v75, c("gene_id" ,"gene.gene_name", "gene.gene_biotype"))



#######################################################
##
##  check built queries
##
ensembldb:::.buildQuery(EnsDb.Hsapiens.v75, columns=c("gene_id", "gene_name", "tx_id", "exon_id"), filter=list(SeqnameFilter("Y"), SeqstrandFilter("-")))


## throws a warning
ensembldb:::.buildQuery(EnsDb.Hsapiens.v75, columns=c("gene_id", "gene_name", "tx_id", "exon_id"), filter=list(SeqnameFilter("Y"), SeqstrandFilter("-")), order.by="exon_seq_end", order.type="desc")


## works
ensembldb:::.buildQuery(EnsDb.Hsapiens.v75, columns=c("gene_id", "gene_name", "tx_id", "exon_id", "exon_seq_end"), filter=list(SeqnameFilter("Y"), SeqstrandFilter("-")), order.by="exon_seq_end", order.type="desc")


ensembldb:::.buildQuery(EnsDb.Hsapiens.v75, columns=c("tx_id", "exon_id", "exon_seq_end"), filter=list(SeqnameFilter("Y"), SeqstrandFilter("-")), order.by="exon_seq_end", order.type="desc")


ensembldb:::.buildQuery(EnsDb.Hsapiens.v75, columns=c("tx_id", "gene_id"))


## check the new filter thingy.
GF <- GeneidFilter("a")
where(GF)
column(GF)

## with db
column(GF, DB)
where(GF, DB)

## with db and with.tables
column(GF, DB, with.tables="tx")
where(GF, DB, with.tables="tx")


column(GF, DB, with.tables=c("gene", "tx"))
where(GF, DB, with.tables=c("gene", "tx"))

## does throw an error!
##column(GF, DB, with.tables="exon")

## silently drops the submitted ones.
column(GF, DB, with.tables="blu")

##
ensembldb:::.buildQuery(DB, columns=c("tx_id", "gene_id"))
## with filter
ensembldb:::.buildQuery(DB, columns=c("tx_id", "gene_id"),
                        filter=list(GeneidFilter("a")))
ensembldb:::.buildQuery(DB, columns=c("tx_id", "gene_id"),
                        filter=list(GeneidFilter("a"),
                                    SeqnameFilter(1)))

ensembldb:::.buildQuery(DB, columns=c("tx_id", "gene_id", "exon_idx"),
                        filter=list(GeneidFilter("a"),
                                    SeqnameFilter(1)))

