## check namespace.
detachem <- function( x ){
    NS <- loadedNamespaces()
    if( any( NS==x ) ){
        pkgn <- paste0( "package:", x )
        detach( pkgn, unload=TRUE, character.only=TRUE )
    }
}
Pkgs <- c( "EnsDb.Hsapiens.v75", "ensembldb" )
tmp <- sapply( Pkgs, detachem )
tmp <- sapply( Pkgs, library, character.only=TRUE )

###

## just get all genes.
cat( "getting all genes..." )
Gns <- genes( EnsDb.Hsapiens.v75 )
Gns
cat("done\n")

cat( "getting all transcripts..." )
Gns <- transcripts( EnsDb.Hsapiens.v75 )
Gns
cat("done\n")

cat( "getting all exons..." )
Gns <- exons( EnsDb.Hsapiens.v75 )
Gns
cat("done\n")

## get exons, sort by exon_seq_start
Gns <- exons( EnsDb.Hsapiens.v75, columns=c( "exon_id", "tx_id" ), filter=list( TxidFilter( "a" ) ) )
ensembldb:::.buildQuery( EnsDb.Hsapiens.v75, columns=c( "exon_id", "tx_id" ), filter=list( TxidFilter( "a" ) ))

cat( "all transcripts by..." )
tmp <- transcriptsBy( EnsDb.Hsapiens.v75 )
tmp
cat("done\n")

cat( "all exons by..." )
tmp <- exonsBy( EnsDb.Hsapiens.v75 )
tmp
cat("done\n")


###########
## getWhat... generic query interface to the database.
Test <- ensembldb:::getWhat( EnsDb.Hsapiens.v75, columns=c( "gene_id", "gene_biotype", "gene_name", "seq_name" ), filter=list( SeqnameFilter( "Y" ) ) )
head(Test)
dim(Test)

## now let's joind exon...
Test <- ensembldb:::getWhat( EnsDb.Hsapiens.v75, columns=c( "gene_id", "gene_biotype", "gene_name", "seq_name", "exon_id", "exon_seq_start", "exon_seq_end" ), order.by="exon_seq_end", order.type="desc", filter=list( SeqnameFilter( "Y" ) ) )
head(Test)
dim(Test)


## throws a warning since exon_chrom_end is not valid.
Test <- ensembldb:::getWhat( EnsDb.Hsapiens.v75, columns=c( "gene_id", "gene_biotype", "gene_name", "seq_name", "exon_id", "exon_seq_start", "exon_seq_end" ), order.by="exon_chrom_end", order.type="desc", filter=list( SeqnameFilter( "Y" ) ) )
head(Test)
dim(Test)


## add a Txid Filter.
Test <- ensembldb:::getWhat( EnsDb.Hsapiens.v75, columns=c( "gene_id", "gene_biotype", "gene_name", "seq_name", "exon_id", "exon_seq_start", "exon_seq_end", "tx_id" ), order.by="exon_seq_end", order.type="desc", filter=list( TxidFilter( "ENST00000028008" ) ) )
Test

Test <- ensembldb:::getWhat( EnsDb.Hsapiens.v75, columns=c( "gene_id", "gene_biotype", "gene_name", "seq_name" ), filter=list( TxidFilter( "ENST00000028008" ) ) )
Test



######
## exonsBy
## get all Exons by gene for genes encoded on chromosomes 1, 2, 4
Test <- exonsBy( EnsDb.Hsapiens.v75, by="gene", columns=c( "gene_id", "gene_name", "gene_biotype" ), filter=list( SeqnameFilter( c( 1, 2,4 ) ), SeqstrandFilter( "-" ) ) )
Test

## tx_biotype and tx_id have been removed.
Test <- exonsBy( EnsDb.Hsapiens.v75, by="gene", columns=c( "gene_id", "gene_name", "gene_biotype", "tx_biotype", "tx_id" ), filter=list( SeqnameFilter( c( 1, 2,4 ) ), SeqstrandFilter( "-" ) ) )
Test

Test <- exonsBy( EnsDb.Hsapiens.v75, by="tx", columns=c( "gene_id", "tx_id", "tx_biotype" ), filter=list( SeqnameFilter( c( 1, 2,4 ) ) ) )
Test

## exons for a specific transcript
Test <- exonsBy( EnsDb.Hsapiens.v75, by="tx", columns=c( "gene_id", "tx_id", "tx_biotype" ), filter=list( TxidFilter( "ENST00000028008" ) ) )
Test

## that also works, albeit throwing an warning.
Test <- exonsBy( EnsDb.Hsapiens.v75, by="gene", columns=c( "gene_id", "tx_id", "tx_biotype" ), filter=list( TxidFilter( "ENST00000028008" ) ) )
Test



########
## transcriptsBy
Test <- transcriptsBy( EnsDb.Hsapiens.v75, by="gene", filter=list( SeqstrandFilter( "+" ), SeqnameFilter( "X" ) ) )
Test

## that should throw a warning
Test <- transcriptsBy( EnsDb.Hsapiens.v75, by="gene", filter=list( SeqstrandFilter( "+" ), SeqnameFilter( "X" ) ), columns=c( "exon_id", "exon_seq_start" ) )
Test


Test <- transcriptsBy( EnsDb.Hsapiens.v75, by="exon", filter=list( SeqstrandFilter( "+" ), SeqnameFilter( "X" ) ), columns="tx_biotype" )
Test

## that should throw a warning
Test <- transcriptsBy( EnsDb.Hsapiens.v75, by="exon", filter=list( SeqstrandFilter( "+" ), SeqnameFilter( "X" ) ), columns=c( "exon_id", "exon_seq_start", "tx_biotype" ) )
Test


######
## genes
Test <- genes( EnsDb.Hsapiens.v75, filter=list( GenebiotypeFilter( "lincRNA" ) ) )
head( Test )
length( Test )

## adding tx properties along with gene columns; this will return a data.frame with the
## additional information; gene columns can however no longer be unique in the data.frame
Test <- genes( EnsDb.Hsapiens.v75, filter=list( GenebiotypeFilter( "lincRNA" ) ), columns=c( listColumns( EnsDb.Hsapiens.v75, "gene"), "tx_id", "tx_biotype" ) )
head( Test )
length( Test )

######
## transcripts
## get all transcripts that are target to nonsense mediated decay
Test <- transcripts( EnsDb.Hsapiens.v75, filter=list( TxbiotypeFilter( "nonsense_mediated_decay" ) ) )
head( Test )
length( Test )

## order the transcripts by seq_name; this does not work.
Test <- transcripts( EnsDb.Hsapiens.v75, filter=list( TxbiotypeFilter( "nonsense_mediated_decay" ) ), order.by="seq_name" )
head( Test )
nrow( Test )

## order the transcripts by seq_name; have to explicitely add seq_name to the columns.
Test <- transcripts( EnsDb.Hsapiens.v75, filter=list( TxbiotypeFilter( "nonsense_mediated_decay" ) ), order.by="seq_name", columns=c( listColumns( EnsDb.Hsapiens.v75, "tx" ), "seq_name" ) )
head( Test )
nrow( Test )

## get in addition the gene_name and gene_id
Test <- transcripts( EnsDb.Hsapiens.v75, filter=list( TxbiotypeFilter( "nonsense_mediated_decay" ) ), columns=c( listColumns( EnsDb.Hsapiens.v75, "tx" ), "gene_id", "gene_name" ) )
head( Test )
nrow( Test )

## get in addition the gene_name and gene_id and also exon_id and exon_idx
Test <- transcripts( EnsDb.Hsapiens.v75, filter=list( TxbiotypeFilter( "nonsense_mediated_decay" ) ), columns=c( listColumns( EnsDb.Hsapiens.v75, "tx" ), "gene_id", "gene_name", "exon_id", "exon_idx" ) )
head( Test )
nrow( Test )


#####
## exons
##
Test <- exons( EnsDb.Hsapiens.v75, filter=list( TxidFilter( "ENST00000028008" ) ), columns=c( "gene_id","gene_name", "gene_biotype" ) )
Test




##################
## examples from EnsDb-class:

## display some information:
EnsDb.Hsapiens.v75

organism( EnsDb.Hsapiens.v75 )

seqinfo( EnsDb.Hsapiens.v75 )

## show the tables
listTables( EnsDb.Hsapiens.v75 )


######    buildQuery
##
## join tables gene and transcript and return gene_id and tx_id
buildQuery( EnsDb.Hsapiens.v75, columns=c( "gene_id", "tx_id" ) )


## get all exon_ids and transcript ids of genes encoded on chromosome Y.
buildQuery( EnsDb.Hsapiens.v75, columns=c( "exon_id", "tx_id" ), filter=list( SeqnameFilter(  "Y") ) )


######   genes
##
## get all genes coded on chromosome Y
AllY <- genes( EnsDb.Hsapiens.v75, filter=list( SeqnameFilter( "Y" ) ) )
head( AllY )

## return result as GRanges.
AllY.granges <- genes( EnsDb.Hsapiens.v75, filter=list( SeqnameFilter(
  "Y" ) ), return.type="GRanges" )
AllY.granges

## include all transcripts of the gene and their chromosomal
## coordinates, sort by chrom start of transcripts and return as
## GRanges.
AllY.granges.tx <- genes( EnsDb.Hsapiens.v75, filter=list(
  SeqnameFilter( "Y" ) ), return.type="GRanges", columns=c(
  "gene_id", "seq_name", "seq_strand", "tx_id", "tx_biotype",
  "tx_seq_start", "tx_seq_end" ), order.by="tx_seq_start" )
AllY.granges.tx



######   transcripts
##
## get all transcripts of a gene
Tx <- transcripts( EnsDb.Hsapiens.v75, filter=list( GeneidFilter(
  "ENSG00000184895" ) ), order.by="tx_seq_start" )
Tx

## get all transcripts of two genes along with some information on the
## gene and transcript
Tx.granges <- transcripts( EnsDb.Hsapiens.v75, filter=list(
  GeneidFilter( c( "ENSG00000184895", "ENSG00000092377" ),
  condition="in" )), return.type="GRanges", order.by="tx_seq_start",
  columns=c( "gene_id", "gene_seq_start", "gene_seq_end",
  "gene_biotype", "tx_biotype" ) )
Tx.granges



######   exons
##
## get all exons of the provided genes
Exon.granges <- exons( EnsDb.Hsapiens.v75, filter=list( GeneidFilter( c(
  "ENSG00000184895", "ENSG00000092377" ) )),
  return.type="GRanges", order.by="exon_seq_start", columns=c(
  "gene_id", "gene_seq_start", "gene_seq_end", "gene_biotype" ) )
Exon.granges



#####    exonsBy
##
## get all exons for transcripts encoded on chromosomes 1 to 22, X and Y.
ETx <- exonsBy( EnsDb.Hsapiens.v75, by="tx", filter=list( SeqnameFilter(
  c( 1:22, "X", "Y" ) ) ) )
ETx
## get all exons for genes encoded on chromosome 1 to 22, X and Y and
## include additional annotation columns in the result
EGenes <- exonsBy( EnsDb.Hsapiens.v75, by="gene", filter=list(
  SeqnameFilter( c( 1:22, "X", "Y" ) ) ), columns=c( "gene_biotype",
  "gene_name" ) )
EGenes

## Note that this might also contain "LRG" genes.
sum( grep( names( EGenes ), pattern="LRG" ) )
## fetch just Ensembl genes:
EGenes <- exonsBy( EnsDb.Hsapiens.v75, by="gene", filter=list(
  SeqnameFilter( c( 1:22, "X", "Y" ) ), GeneidFilter( "ENS%", "like" ) ), columns=c( "gene_biotype",
  "gene_name" ) )

sum( grep( names( EGenes ), pattern="LRG" ) )



#####    transcriptsBy
##
TGenes <- transcriptsBy( EnsDb.Hsapiens.v75, by="gene", filter=list(
  SeqnameFilter( c( 1:22, "X", "Y" ) ) ) )
TGenes



#####    lengthOf
##
## length of a specific gene.
lengthOf( EnsDb.Hsapiens.v75, filter=list( GeneidFilter(
  "ENSG00000000003" ) ) )

## length of a transcript
lengthOf( EnsDb.Hsapiens.v75, of="tx", filter=list( TxidFilter(
  "ENST00000494424" ) ) )

## average length of all protein coding genes
mean( lengthOf( EnsDb.Hsapiens.v75, of="gene", filter=list(
  GenebiotypeFilter( "protein_coding" ),
  SeqnameFilter( c( 1:22, "X", "Y" ) ) ) ) )

## average length of all snoRNAs
mean( lengthOf( EnsDb.Hsapiens.v75, of="gene", filter=list(
  GenebiotypeFilter( "snoRNA" ),
  SeqnameFilter( c( 1:22, "X", "Y" ) ) ) ) )

listGenebiotypes(EnsDb.Hsapiens.v75)

listTxbiotypes(EnsDb.Hsapiens.v75)

