##***********************************************************************
##
##     Generic methods
##
##***********************************************************************
if(!isGeneric("column"))
    setGeneric("column", function(object, db, with.tables, ...)
        standardGeneric("column"))
if(!isGeneric("buildQuery"))
    setGeneric("buildQuery", function(x, ...)
        standardGeneric("buildQuery"))
if(!isGeneric("cleanColumns"))
    setGeneric("cleanColumns", function(x, columns, ...)
        starndardGeneric("cleanColumns"))
if(!isGeneric("condition"))
    setGeneric("condition", function(x, ...)
        standardGeneric("condition"))
if(!isGeneric("genes"))
    setGeneric("genes", function(x, ...)
        standardGeneric("genes"))
if(!isGeneric("getWhat"))
    setGeneric("getWhat", function(x, ...)
        standardGeneric("getWhat"))
if(!isGeneric("ensemblVersion"))
    setGeneric("ensemblVersion", function(x)
        standardGeneric("ensemblVersion"))
if(!isGeneric("exons"))
    setGeneric("exons", function(x, ...)
        standardGeneric("exons"))
if(!isGeneric("exonsBy"))
    setGeneric("exonsBy", function(x, ...)
        standardGeneric("exonsBy"))
if(!isGeneric("getGenomeFaFile"))
    setGeneric("getGenomeFaFile", function(x, ...)
        standardGeneric("getGenomeFaFile"))
if(!isGeneric("getMetadataValue"))
    setGeneric("getMetadataValue", function(x, name)
        standardGeneric("getMetadataValue"))
if(!isGeneric("listColumns")){
    setGeneric("listColumns", function(x, ...)
        standardGeneric("listColumns"))
}
if(!isGeneric("listGenebiotypes")){
    setGeneric("listGenebiotypes", function(x, ...)
        standardGeneric("listGenebiotypes"))
}
if(!isGeneric("listTxbiotypes")){
    setGeneric("listTxbiotypes", function(x, ...)
        standardGeneric("listTxbiotypes"))
}
if(!isGeneric("lengthOf"))
    setGeneric("lengthOf", function(x, ...)
        standardGeneric("lengthOf"))
if(!isGeneric("print"))
    setGeneric("print", function(x, ...)
        standardGeneric("print"))
if(!isGeneric("requireTable"))
    setGeneric("requireTable", function(x, db, ...)
        standardGeneric("requireTable"))
if(!isGeneric("seqinfo"))
    setGeneric("seqinfo", function(x)
        standardGeneric("seqinfo"))
if(!isGeneric("show"))
    setGeneric("show", function(object, ...)
        standardGeneric("show"))
if(!isGeneric("toSAF"))
    setGeneric("toSAF", function(x, ...)
        standardGeneric("toSAF"))
if(!isGeneric("listTables")){
    setGeneric("listTables", function(x, ...)
        standardGeneric("listTables"))
}
if(!isGeneric("tablesByDegree")){
    setGeneric("tablesByDegree", function(x, ...)
        standardGeneric("tablesByDegree"))
}
if(!isGeneric("tablesForColumns"))
    setGeneric("tablesForColumns", function(x, attributes, ...)
        standardGeneric("tablesForColumns"))
if(!isGeneric("transcripts"))
    setGeneric("transcripts", function(x, ...)
        standardGeneric("transcripts"))
if(!isGeneric("transcriptsBy"))
    setGeneric("transcriptsBy", function(x, ...)
        standardGeneric("transcriptsBy"))
if(!isGeneric("where"))
    setGeneric("where", function(object, db, with.tables, ...)
        standardGeneric("where"))


##***********************************************************************
##
##     Methods for EnsDb classes
##
##***********************************************************************
setMethod("show", "EnsDb", function(object){
    if(is.null(object@ensdb)){
        cat("Dash it! Got an empty thing!\n")
    }else{
        info <- dbGetQuery(object@ensdb, "select * from metadata")
        cat("EnsDb for Ensembl:\n")
        for(i in 1:nrow(info)){
            cat(paste0("|", info[ i, "name" ], ": ",
                       info[ i, "value" ], "\n"))
        }
        ## gene and transcript info.
        cat(paste0("| No. of genes: ", dbGetQuery(object@ensdb, "select count(distinct gene_id) from gene")[ 1,1 ], ".\n"))
        cat(paste0("| No. of transcripts: ", dbGetQuery(object@ensdb, "select count(distinct tx_id) from tx")[ 1,1 ], ".\n"))
    }
})

setMethod("organism", "EnsDb", function(object){
    Species <- .getMetaDataValue(object@ensdb, "Organism")
    ## reformat the e.g. homo_sapiens string into Homo sapiens
                                        #
    Species <- gsub(Species, pattern="_", replacement=" ", fixed=TRUE)
    Species <- .organismName(Species)
    return(Species)
})

setMethod("metadata", "EnsDb", function(x, ...){
    Res <- dbGetQuery(dbconn(x), "select * from metadata")
    return(Res)
})
#####
## Validation
##
validateEnsDb <- function(object){
    ## check if the database contains all required tables...
    if(!is.null(object@ensdb)){
        have <- dbListTables(object@ensdb)
        need <- c("chromosome", "exon", "gene", "metadata", "tx",
                  "tx2exon")
        notthere <- need[ !(need %in% have) ]
        if(length(notthere) > 0){
            return(paste("Required tables", notthere,
                         "not found in the database!"))
        }
    }
    return(TRUE)
}
setValidity("EnsDb", validateEnsDb)
setMethod("initialize", "EnsDb", function(.Object,...){
    OK <- validateEnsDb(.Object)
    if(class(OK)=="character"){
        stop(OK)
    }
    callNextMethod(.Object, ...)
})

### connection:
## returns the connection object to the SQL database
setMethod("dbconn", "EnsDb", function(x){
    return(x@ensdb)
})

### ensemblVersion
## returns the ensembl version of the package.
setMethod("ensemblVersion", "EnsDb", function(x){
    eVersion <- getMetadataValue(x, "ensembl_version")
    return(eVersion)
})
### getMetadataValue
## returns the metadata value for the specified name/key
setMethod("getMetadataValue", "EnsDb", function(x, name){
    if(missing(name))
        stop("Argument name has to be specified!")
    return(metadata(x)[metadata(x)$name==name, "value"])
})

### seqinfo
## returns the sequence/chromosome information from the database.
setMethod("seqinfo", "EnsDb", function(x){
    Chrs <- dbGetQuery(dbconn(x), "select * from chromosome")
    Chr.build <- .getMetaDataValue(dbconn(x), "genome_build")
    SI <- Seqinfo(seqnames=Chrs$seq_name,
                  seqlengths=Chrs$seq_length,
                  isCircular=Chrs$is_circular==1, genome=Chr.build)
    return(SI)
})

### getGenomeFaFile
## queries the dna.toplevel.fa file from AnnotationHub matching the current
## Ensembl version
setMethod("getGenomeFaFile", "EnsDb", function(x, pattern="dna.toplevel.fa"){
    ah <- AnnotationHub()
    queryCol <- "title"  ## the column we want to query
    eData <- query(ah, c(organism(x), paste0("release-", ensemblVersion(x))))
    idx <- grep(mcols(eData)[, queryCol], pattern=pattern)
    if(length(idx) == 0){
        stop(paste0("No genomic fasta file found in AnnotationHub for organism ",
                    organism(x), " and Ensembl version ", ensemblVersion(x), "!\n",
                    "Available resources are:",
                    paste0(mcols(eData)[, queryCol], collapse="\n"), "\n"))
    }
    if(length(idx) > 1){
        warning(paste0("Found more than one matching file: ",
                       paste0(mcols(eData)[idx, queryCol], collapse="\n"),
                       "\nUsing only the first." ))
        idx <- idx[1]
    }
    Dna <- ah[[names(eData)[idx]]]
    ## generate an index if none is available
    if(is.na(index(Dna))){
        indexFa(Dna)
        Dna <- FaFile(path(Dna))
    }
    return(Dna)
})

### listTables
## returns a named list with database table columns
setMethod("listTables", "EnsDb", function(x, ...){
    if(length(x@tables)==0){
        tables <- dbListTables(dbconn(x))
        ## read the columns for these tables.
        Tables <- vector(length=length(tables), "list")
        for(i in 1:length(Tables)){
            Tables[[ i ]] <- colnames(dbGetQuery(dbconn(x),
                                                 paste0("select * from ",
                                                        tables[ i ],
                                                        " limit 1")))
        }
        names(Tables) <- tables
        x@tables <- Tables
    }
    Tab <- x@tables
    Tab <- Tab[ tablesByDegree(x, tab=names(Tab)) ]
    return(Tab)
})

### listColumns
## lists all columns.
setMethod("listColumns", "EnsDb", function(x,
                                           table,
                                           skip.keys=TRUE, ...){
    if(length(x@tables)==0){
        tables <- dbListTables(dbconn(x))
        ## read the columns for these tables.
        Tables <- vector(length=length(tables), "list")
        for(i in 1:length(Tables)){
            Tables[[ i ]] <- colnames(dbGetQuery(dbconn(x),
                                                 paste0("select * from ",
                                                        tables[ i ],
                                                        " limit 1")))
        }
        names(Tables) <- tables
        x@tables <- Tables
    }
    if(!missing(table)){
        columns <- x@tables[[ table ]]
    }else{
        columns <- unlist(x@tables, use.names=FALSE)
    }
    if(skip.keys){
        ## remove everything that has a _pk or _fk...
        idx <- grep(columns, pattern="_fk$")
        if(length(idx) > 0)
            columns <- columns[ -idx ]
        idx <- grep(columns, pattern="_pk$")
        if(length(idx) > 0)
            columns <- columns[ -idx ]
    }
    return(columns)
})

setMethod("listGenebiotypes", "EnsDb", function(x, ...){
    return(dbGetQuery(dbconn(x), "select distinct gene_biotype from gene")[,1])
})
setMethod("listTxbiotypes", "EnsDb", function(x, ...){
    return(dbGetQuery(dbconn(x), "select distinct tx_biotype from tx")[,1])
})

### cleanColumns
## checks columns and removes all that are not present in database tables
## the method checks internally whether the columns are in the full form,
## i.e. gene.gene_id (<table name>.<column name>)
setMethod("cleanColumns", "EnsDb", function(x,
                                            columns, ...){
    if(missing(columns))
        stop("No columns submitted!")
    ## vote of the majority
    full.name <- length(grep(columns, pattern=".", fixed=TRUE)) >
        floor(length(columns) /2)
    if(full.name){
        suppressWarnings(
            full.columns <- unlist(prefixColumns(x,
                                                 unlist(listTables(x)),
                                                 clean=FALSE),
                                   use.names=TRUE)
          )
        bm <- columns %in% full.columns
        removed <- columns[ !bm ]
    }else{
        bm <- columns %in% unlist(listTables(x)[ c("gene", "tx", "exon",
                                                   "tx2exon", "chromosome") ])
        removed <- columns[ !bm ]
    }
    if(length(removed) > 0){
        warning("Columns ", paste(removed, collapse=", "),
                " are not valid and have been removed")
    }
    return(columns[ bm ])
})

### tablesForColumns
## returns the tables for the specified columns.
setMethod("tablesForColumns", "EnsDb", function(x, columns, ...){
    if(missing(columns))
        stop("No columns submitted!")
    bm <- unlist(lapply(listTables(x), function(z){
        return(any(z %in% columns))
    }))
    if(!any(bm))
        return(NULL)
    Tables <- names(bm)[ bm ]
    Tables <- Tables[ !(Tables %in% c("metadata")) ]
    return(Tables)
})

## returns the table names ordered by degree, i.e. edges to other tables
setMethod("tablesByDegree", "EnsDb", function(x,
                                              tab=names(listTables(x)),
                                              ...){
    ## ## to do this with a graph:
    ## DBgraph <- graphNEL(nodes=c("gene", "tx", "tx2exon", "exon", "chromosome", "information"),
    ##                  edgeL=list(gene=c("tx", "chromosome"),
    ##                      tx=c("gene", "tx2exon"),
    ##                      tx2exon=c("tx", "exon"),
    ##                      exon="tx2exon",
    ##                      chromosome="gene"
    ##                          ))
    ## Tab <- names(sort(degree(DBgraph), decreasing=TRUE))
    Table.order <- c(gene=1, tx=2, tx2exon=3, exon=4, chromosome=5, metadata=6)
    ##Table.order <- c(gene=2, tx=1, tx2exon=3, exon=4, chromosome=5, metadata=6)
    Tab <- tab[ order(Table.order[ tab ]) ]
    return(Tab)
})




### genes:
## get genes from the database.
setMethod("genes", "EnsDb", function(x,
                                     columns=listColumns(x, "gene"),
                                     filter, order.by="",
                                     order.type="asc",
                                     return.type="GRanges"){
    return.type <- match.arg(return.type, c("data.frame", "GRanges", "DataFrame"))
    columns <- unique(c("gene_id", columns))
    ## if return.type is GRanges we require columns: seq_name, gene_seq_start
    ## and gene_seq_end and seq_strand
    if(return.type=="GRanges"){
        columns <- unique(c(columns, c("gene_seq_start", "gene_seq_end",
                                       "seq_name", "seq_strand")))
    }
    if(missing(filter)){
        filter=list()
    }else{
        filter <- checkFilter(filter)
    }
    Res <- getWhat(x, columns=columns, filter=filter,
                   order.by=order.by, order.type=order.type)
    if(return.type=="data.frame"){
        return(Res)
    }
    if(return.type=="GRanges"){
        metacols <- columns[ !(columns %in% c("seq_name",
                                              "seq_strand",
                                              "gene_seq_start",
                                              "gene_seq_end")) ]
        SI <- seqinfo(x)
        GR <- GRanges(seqnames=Rle(Res$seq_name),
                      ranges=IRanges(start=Res$gene_seq_start, end=Res$gene_seq_end),
                      strand=Rle(Res$seq_strand),
                      seqinfo=SI[ unique(Res$seq_name) ],
                      Res[ , metacols, drop=FALSE ]
                    )
        names(GR) <- Res$gene_id
        return(GR)
    }
    if(return.type=="DataFrame"){
        return(DataFrame(Res))
    }
})

### transcripts:
## get transcripts from the database.
setMethod("transcripts", "EnsDb", function(x, columns=listColumns(x, "tx"),
                                           filter, order.by="", order.type="asc",
                                           return.type="GRanges"){
    return.type <- match.arg(return.type, c("data.frame", "GRanges", "DataFrame"))
    columns <- unique(c("tx_id", columns))
    ## if return.type is GRanges we require columns: seq_name, gene_seq_start
    ## and gene_seq_end and seq_strand
    if(return.type=="GRanges"){
        columns <- unique(c(columns, c("tx_seq_start",
                                       "tx_seq_end",
                                       "seq_name",
                                       "seq_strand")))
    }
    if(missing(filter)){
        filter=list()
    }else{
        filter <- checkFilter(filter)
    }
    Res <- getWhat(x, columns=columns, filter=filter,
                   order.by=order.by, order.type=order.type)
    if(return.type=="data.frame"){
        return(Res)
    }
    if(return.type=="GRanges"){
        metacols <- columns[ !(columns %in% c("seq_name",
                                              "seq_strand",
                                              "tx_seq_start",
                                              "tx_seq_end")) ]
        SI <- seqinfo(x)
        GR <- GRanges(seqnames=Rle(Res$seq_name),
                      ranges=IRanges(start=Res$tx_seq_start, end=Res$tx_seq_end),
                      strand=Rle(Res$seq_strand),
                      seqinfo=SI[ unique(Res$seq_name) ],
                      Res[ , metacols, drop=FALSE ]
                    )
        names(GR) <- Res$tx_id
        return(GR)
    }
    if(return.type=="DataFrame"){
        return(DataFrame(Res))
    }
})

### promoters:
## get promoter regions from the database.
setMethod("promoters", "EnsDb",
          function(x, upstream=2000, downstream=200, ...)
          {
              gr <- transcripts(x, ...)
              trim(suppressWarnings(promoters(gr,
                                              upstream=upstream,
                                              downstream=downstream)))
          }
)

### exons:
## get exons from the database.
setMethod("exons", "EnsDb", function(x, columns=listColumns(x, "exon"), filter,
                                     order.by="", order.type="asc",
                                     return.type="GRanges"){
    return.type <- match.arg(return.type, c("data.frame", "GRanges", "DataFrame"))
    if(!any(columns %in% c(listColumns(x, "exon"), "exon_idx"))){
        ## have to have at least one column from the gene table...
        columns <- c("exon_id", columns)
    }
    columns <- unique(c("exon_id", columns))
    ## if return.type is GRanges we require columns: seq_name, gene_seq_start
    ## and gene_seq_end and seq_strand
    if(return.type=="GRanges"){
        columns <- unique(c(columns, c("exon_seq_start",
                                       "exon_seq_end",
                                       "seq_name",
                                       "seq_strand")))
    }
    if(missing(filter)){
        filter=list()
    }else{
        filter <- checkFilter(filter)
    }
    Res <- getWhat(x, columns=columns, filter=filter,
                   order.by=order.by, order.type=order.type)
    if(return.type=="data.frame"){
        return(Res)
    }
    if(return.type=="GRanges"){
        metacols <- columns[ !(columns %in% c("seq_name",
                                              "seq_strand",
                                              "exon_seq_start",
                                              "exon_seq_end")) ]
        SI <- seqinfo(x)
        GR <- GRanges(seqnames=Rle(Res$seq_name),
                      ranges=IRanges(start=Res$exon_seq_start, end=Res$exon_seq_end),
                      strand=Rle(Res$seq_strand),
                      seqinfo=SI[ unique(Res$seq_name) ],
                      Res[ , metacols, drop=FALSE ]
                    )
        names(GR) <- Res$exon_id
        return(GR)
    }
    if(return.type=="DataFrame"){
        return(DataFrame(Res))
    }
})


## should return a GRangesList
## still considerably slower than the corresponding call in the GenomicFeatures package.
setMethod("exonsBy", "EnsDb", function(x, by=c("tx", "gene"),
                                       columns=listColumns(x, "exon"), filter){
    by <- match.arg(by, c("tx", "gene"))
    ## note: if it's by gene we don't want any columns from the transcript table AND
    ## we don't want the exon_rank! rather we would like to sort by exon_chrom_start * seq_strand
    min.columns <- c(paste0(by, "_id"), "seq_name",
                     "exon_seq_start", "exon_seq_end", "exon_id", "seq_strand")
    by.id.full <- unlist(prefixColumns(x, columns=paste0(by, "_id"),
                                       clean=FALSE),
                         use.names=FALSE)
    if(by=="gene"){
        ## tx columns have to be removed, since the same exon can be part of more
        ## than one tx
        txcolumns <- listColumns(x, "tx")
        txcolumns <- txcolumns[ txcolumns != "gene_id" ]
        torem <- columns %in% txcolumns
        if(any(torem))
            warning(paste0("Columns ",
                           paste(columns[ torem ], collapse=","),
                           " have been removed as they are not allowed if exons are fetched by gene."))
        columns <- columns[ !torem ]
        ##order.by <- paste0(by.id.full, ", case when seq_strand=1 then exon_seq_start when seq_strand=-1 then (exon_seq_end * -1) end")
        order.by <- paste0("case when seq_strand=1 then exon_seq_start when seq_strand=-1 then (exon_seq_end * -1) end")
    }else{
        min.columns <- c(min.columns, "exon_idx")
        ##order.by <- paste0(by.id.full, ",exon_idx")
        order.by <- paste0("exon_idx")
    }
    ## define the minimal columns that we need...
    columns <- unique(c(columns, min.columns))
    ## get the seqinfo:
    SI <- seqinfo(x)
    ##Res <- getWhat(dbconn(x), columns=columns, filter=filter, order.by=paste0("seq_name,gene_seq_start,",by ,"_id,exon_idx"))
    if(missing(filter)){
        filter=list()
    }else{
        filter <- checkFilter(filter)
    }
    Res <- getWhat(x, columns=columns, filter=filter, order.by=order.by, skip.order.check=TRUE)
    SI <- SI[ unique(Res$seq_name) ]
    ## replace exon_idx with exon_rank
    colnames(Res)[ colnames(Res)=="exon_idx" ] <- "exon_rank"
    columns[ columns=="exon_idx" ] <- "exon_rank"
    columns.metadata <- columns[ !(columns %in% c("seq_name", "seq_strand", "exon_seq_start", "exon_seq_end", paste0(by, "_id"))) ]
    columns.metadata <- match(columns.metadata, colnames(Res))   ## presumably faster...
    GR <- GRanges(seqnames=Rle(Res$seq_name),
                  strand=Rle(Res$seq_strand),
                  ranges=IRanges(start=Res$exon_seq_start, end=Res$exon_seq_end),
                  seqinfo=SI,
                  Res[ , columns.metadata, drop=FALSE ]
                )
    ## now that GR is ordered as we wanted; once we split it it will be ordered by
    ## the value which we used for splitting.
    return(split(GR, Res[ , paste0(by, "_id") ]))
})


## should return a GRangesList
## still considerably slower than the corresponding call in the GenomicFeatures package.
setMethod("transcriptsBy", "EnsDb", function(x, by=c("gene", "exon"),
                                             columns=listColumns(x, "tx"), filter){
    by <- match.arg(by, c("gene", "exon"))
    if(by=="cds")
        stop("fetching transcripts by cds is not (yet) implemented.")
    min.columns <- c(paste0(by, "_id"), "seq_name", "tx_seq_start",
                     "tx_seq_end", "tx_id", "seq_strand")
    ## can not have exon columns!
    torem <- columns %in% c(listColumns(x, "exon"), "exon_idx")
    if(any(torem))
        warning(paste0("Columns ",
                       paste(columns[ torem ], collapse=","),
                       " have been removed as they are not allowed if transcripts are fetched."))
    columns <- columns[ !torem ]
    by.id.full <- unlist(prefixColumns(x, columns=paste0(by, "_id"),
                                       clean=FALSE), use.names=FALSE)
    order.by <- paste0(by.id.full , ", case when seq_strand=1 then tx_seq_start when seq_strand=-1 then (tx_seq_end * -1) end")
    ## define the minimal columns that we need...
    columns <- unique(c(columns, min.columns))
    ## get the seqinfo:
    SI <- seqinfo(x)
    if(missing(filter)){
        filter=list()
    }else{
        filter <- checkFilter(filter)
    }
    Res <- getWhat(x, columns=columns,
                   filter=filter,
                   order.by=order.by,
                   skip.order.check=TRUE)
    SI <- SI[ unique(Res$seq_name) ]
    columns.metadata <- columns[ !(columns %in% c("seq_name", "seq_strand", "tx_seq_start", "tx_seq_end", paste0(by, "_id"))) ]
    columns.metadata <- match(columns.metadata, colnames(Res))   ## presumably faster...
    GR <- GRanges(seqnames=Rle(Res$seq_name),
                  strand=Rle(Res$seq_strand),
                  ranges=IRanges(start=Res$tx_seq_start, end=Res$tx_seq_end),
                  seqinfo=SI,
                  Res[ , columns.metadata, drop=FALSE ]
                )
    ## now that GR is ordered as we wanted; once we split it it will be ordered by
    ## the value which we used for splitting.
    return(split(GR, Res[ , paste0(by, "_id") ]))
})


## for GRangesList...
setMethod("lengthOf", "GRangesList", function(x, ...){
    return(unlist(lapply(width(reduce(x)), sum)))
})

## return the length of genes
setMethod("lengthOf", "EnsDb", function(x, of="gene", filter=list()){
    of <- match.arg(of, c("gene", "tx"))
    ## get the exons by gene or transcript from the database...
    suppressWarnings(
        GRL <- exonsBy(x, by=of, filter=filter)
      )
    return(lengthOf(GRL))
})


## toSAF... function to transform a GRangesList into a data.frame
## corresponding to the SAF format.
## assuming the names of the GRangesList to be the GeneID and the
## element (GRanges) the start/end coordinates
## of an exon, transcript or the gene itself.
.toSaf <- function(x){
    DF <- as.data.frame(x)
    colnames(DF)[ colnames(DF)=="group_name" ] <- "GeneID"
    colnames(DF)[ colnames(DF)=="seqnames" ] <- "Chr"
    colnames(DF)[ colnames(DF)=="start" ] <- "Start"
    colnames(DF)[ colnames(DF)=="end" ] <- "End"
    colnames(DF)[ colnames(DF)=="strand" ] <- "Strand"
    return(DF[ , c("GeneID", "Chr", "Start", "End", "Strand")])
}

## for GRangesList...
setMethod("toSAF", "GRangesList", function(x, ...){
    return(.toSaf(x))
})

.requireTable <- function(db, attr){
    return(names(prefixColumns(db, columns=attr)))
}
## these function determine which tables we need for the submitted filters.
setMethod("requireTable", signature(x="GeneidFilter", db="EnsDb"),
          function(x, db, ...){
              return(.requireTable(db=db, attr="gene_id"))
          })
setMethod("requireTable", signature(x="EntrezidFilter", db="EnsDb"),
          function(x, db, ...){
              return(.requireTable(db=db, attr="entrezid"))
          })
setMethod("requireTable", signature(x="GenebiotypeFilter", db="EnsDb"),
          function(x, db, ...){
              return(.requireTable(db=db, attr="gene_biotype"))
          })
setMethod("requireTable", signature(x="GenenameFilter", db="EnsDb"),
          function(x, db, ...){
              return(.requireTable(db=db, attr="gene_name"))
          })
setMethod("requireTable", signature(x="TxidFilter", db="EnsDb"),
          function(x, db, ...){
              return(.requireTable(db=db, attr="tx_id"))
          })
setMethod("requireTable", signature(x="TxbiotypeFilter", db="EnsDb"),
          function(x, db, ...){
              return(.requireTable(db=db, attr="tx_biotype"))
          })
setMethod("requireTable", signature(x="ExonidFilter", db="EnsDb"),
          function(x, db, ...){
              return(.requireTable(db=db, attr="exon_id"))
          })
setMethod("requireTable", signature(x="SeqnameFilter", db="EnsDb"),
          function(x, db, ...){
              return(.requireTable(db=db, attr="seq_name"))
          })
setMethod("requireTable", signature(x="SeqstrandFilter", db="EnsDb"),
          function(x, db, ...){
              return(.requireTable(db=db, attr="seq_name"))
          })
setMethod("requireTable", signature(x="SeqstartFilter", db="EnsDb"),
          function(x, db, ...){
              if(x@feature=="gene")
                  return(.requireTable(db=db, attr="gene_seq_start"))
              if(x@feature=="transcript" | x@feature=="tx")
                  return(.requireTable(db=db, attr="tx_seq_start"))
              if(x@feature=="exon")
                  return(.requireTable(db=db, attr="exon_seq_start"))
              return(NA)
          })
setMethod("requireTable", signature(x="SeqendFilter", db="EnsDb"),
          function(x, db, ...){
              if(x@feature=="gene")
                  return(.requireTable(db=db, attr="gene_seq_end"))
              if(x@feature=="transcript" | x@feature=="tx")
                  return(.requireTable(db=db, attr="tx_seq_end"))
              if(x@feature=="exon")
                  return(.requireTable(db=db, attr="exon_seq_end"))
              return(NA)
          })
setMethod("buildQuery", "EnsDb",
          function(x, columns=c("gene_id", "gene_biotype", "gene_name"),
                   filter=list(), order.by="",
                   order.type="asc",
                   skip.order.check=FALSE){
              return(.buildQuery(x=x,
                                 columns=columns,
                                 filter=filter,
                                 order.by=order.by,
                                 order.type=order.type,
                                 skip.order.check=skip.order.check))
          })
setMethod("getWhat", "EnsDb",
          function(x, columns=c("gene_id", "gene_biotype", "gene_name"),
                   filter=list(), order.by="", order.type="asc",
                   group.by=NULL, skip.order.check=FALSE){
              return(.getWhat(x=x,
                              columns=columns,
                              filter=filter,
                              order.by=order.by,
                              order.type=order.type,
                              group.by=group.by,
                              skip.order.check=skip.order.check))
          })

## that's basically a copy of the code from the GenomicFeatures package.
setMethod("disjointExons", "EnsDb",
          function(x, aggregateGenes=FALSE, includeTranscripts=TRUE, filter, ...){
              if(missing(filter)){
                  filter <- list()
              }else{
                  filter <- checkFilter(filter)
              }

              exonsByGene <- exonsBy(x, by="gene", filter=filter)
              exonicParts <- disjoin(unlist(exonsByGene, use.names=FALSE))

              if (aggregateGenes) {
                  foGG <- findOverlaps(exonsByGene, exonsByGene)
                  aggregateNames <- GenomicFeatures:::.listNames(names(exonsByGene), as.list(foGG))
                  foEG <- findOverlaps(exonicParts, exonsByGene, select="first")
                  gene_id <- aggregateNames[foEG]
                  pasteNames <- GenomicFeatures:::.pasteNames(names(exonsByGene), as.list(foGG))[foEG]
                  orderByGeneName <- order(pasteNames)
                  exonic_rle <- runLength(Rle(pasteNames[orderByGeneName]))
              } else {
                  ## drop exonic parts that overlap > 1 gene
                  foEG <- findOverlaps(exonicParts, exonsByGene)
                  idxList <- as.list(foEG)
                  if (any(keep <- countQueryHits(foEG) == 1)) {
                      idxList <- idxList[keep]
                      exonicParts <- exonicParts[keep]
                  }
                  gene_id <- GenomicFeatures:::.listNames(names(exonsByGene), idxList)
                  orderByGeneName <- order(unlist(gene_id, use.names=FALSE))
                  exonic_rle <- runLength(Rle(unlist(gene_id[orderByGeneName],
                                                     use.names=FALSE)))
              }
              values <- DataFrame(gene_id)

              if (includeTranscripts) {
                  exonsByTx <- exonsBy(x, by="tx", filter=filter)
                  foET <- findOverlaps(exonicParts, exonsByTx)
                  values$tx_name <- GenomicFeatures:::.listNames(names(exonsByTx), as.list(foET))
              }
              mcols(exonicParts) <- values
              exonicParts <- exonicParts[orderByGeneName]
              exonic_part <- unlist(lapply(exonic_rle, seq_len), use.names=FALSE)
              exonicParts$exonic_part <- exonic_part
              return(exonicParts)
          }
         )


### utility functions
## checkFilter:
## checks the filter argument and ensures that a list of Filter object is returned
checkFilter <- function(x){
    if(class(x)=="list"){
        if(length(x)==0)
            return(x)
        ## check if all elements are Filter classes.
        IsAFilter <- unlist(lapply(x, function(z){
                                        return(inherits(z, what="BasicFilter"))
                                    }))
        if(any(!IsAFilter))
            stop("One of more elements in filter are not filter objects!")
    }else{
        if(inherits(x, what="BasicFilter")){
            x <- list(x)
        }else{
            stop("filter has to be a filter object or a list of filter objects!")
        }
    }
    return(x)
}

