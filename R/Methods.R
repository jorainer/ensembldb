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
    Chrs$seq_name <- formatSeqnamesFromQuery(x, Chrs$seq_name)
    SI <- Seqinfo(seqnames=Chrs$seq_name,
                  seqlengths=Chrs$seq_length,
                  isCircular=Chrs$is_circular==1, genome=Chr.build)
    return(SI)
})

### seqlevels
setMethod("seqlevels", "EnsDb", function(x){
    Chrs <- dbGetQuery(dbconn(x), "select distinct seq_name from chromosome")
    Chrs <- formatSeqnamesFromQuery(x, Chrs$seq_name)
    return(Chrs)
})

### getGenomeFaFile
## queries the dna.toplevel.fa file from AnnotationHub matching the current
## Ensembl version
## Update: if we can't find a FaFile matching the Ensembl version we suggest ones
## that might match.
setMethod("getGenomeFaFile", "EnsDb", function(x, pattern="dna.toplevel.fa"){
    ah <- AnnotationHub()
    ## Reduce the AnnotationHub to species, provider and genome version.
    ah <- .reduceAH(ah, organism=organism(x), dataprovider="Ensembl",
                    genome=unique(genome(x)))
    if(length(ah) == 0)
        stop("Can not find any ressources in AnnotationHub for organism: ",
             organism(x), ", data provider: Ensembl and genome version: ",
             unique(genome(x)), "!")
    ## Reduce to all Fasta files with toplevel or primary_assembly.
    ah <- ah[ah$rdataclass == "FaFile", ]
    if(length(ah) == 0)
        stop("No FaFiles available in AnnotationHub for organism: ",
             organism(x), ", data provider: Ensembl and genome version: ",
             unique(genome(x)), "! You might also try to use the",
             " 'getGenomeTwoBitFile' method instead.")
    ## Reduce to dna.toplevel or dna.primary_assembly.
    idx <- c(grep(ah$title, pattern="dna.toplevel"),
             grep(ah$title, pattern="dna.primary_assembly"))
    if(length(idx) == 0)
        stop("No genome assembly fasta file available for organism: ",
             organism(x), ", data provider: Ensembl and genome version: ",
             unique(genome(x)), "!")
    ah <- ah[idx, ]
    ## Get the Ensembl version from the source url.
    ensVers <- .ensVersionFromSourceUrl(ah$sourceurl)
    if(any(ensVers == ensemblVersion(x))){
        ## Got it.
        itIs <- which(ensVers == ensemblVersion(x))
    }else{
        ## Get the "closest" one.
        diffs <- abs(ensVers - as.numeric(ensemblVersion(x)))
        itIs <- which(diffs == min(diffs))[1]
        message("Returning the Fasta file for Ensembl version ", ensVers[itIs],
                " since no file for Ensembl version ", ensemblVersion(x),
                " is available.")
    }
    ## Getting the ressource.
    Dna <- ah[[names(ah)[itIs]]]
    ## generate an index if none is available
    if(is.na(index(Dna))){
        indexFa(Dna)
        Dna <- FaFile(path(Dna))
    }
    return(Dna)
})
## Just restricting the Annotation Hub to entries matching the species and the
## genome; not yet the Ensembl version.
.reduceAH <- function(ah, organism=NULL, dataprovider="Ensembl",
                      genome=NULL){
    if(!is.null(dataprovider))
        ah <- ah[ah$dataprovider == dataprovider, ]
    if(!is.null(organism))
        ah <- ah[ah$species == organism, ]
    if(!is.null(genome))
        ah <- ah[ah$genome == genome, ]
    return(ah)
}
.ensVersionFromSourceUrl <- function(url){
    url <- strsplit(url, split="/", fixed=TRUE)
    ensVers <- unlist(lapply(url, function(z){
        idx <- grep(z, pattern="^release")
        if(length(idx) == 0)
            return(-1)
        return(as.numeric(unlist(strsplit(z[idx], split="-"))[2]))
    }))
    return(ensVers)
}

####============================================================
##  getGenomeTwoBitFile
##
##  Search and retrieve a genomic DNA resource through a TwoBitFile
##  from AnnotationHub.
####------------------------------------------------------------
setMethod("getGenomeTwoBitFile", "EnsDb", function(x){
    ah <- AnnotationHub()
    ## Reduce the AnnotationHub to species, provider and genome version.
    ah <- .reduceAH(ah, organism=organism(x), dataprovider="Ensembl",
                    genome=unique(genome(x)))
    if(length(ah) == 0)
        stop("Can not find any ressources in AnnotationHub for organism: ",
             organism(x), ", data provider: Ensembl and genome version: ",
             unique(genome(x)), "!")
    ## Reduce to all Fasta files with toplevel or primary_assembly.
    ah <- ah[ah$rdataclass == "TwoBitFile", ]
    if(length(ah) == 0)
        stop("No TwoBitFile available in AnnotationHub for organism: ",
             organism(x), ", data provider: Ensembl and genome version: ",
             unique(genome(x)), "!")
    ## Reduce to dna.toplevel or dna.primary_assembly.
    idx <- c(grep(ah$title, pattern="dna.toplevel"),
             grep(ah$title, pattern="dna.primary_assembly"))
    if(length(idx) == 0)
        stop("No genome assembly fasta file available for organism: ",
             organism(x), ", data provider: Ensembl and genome version: ",
             unique(genome(x)), "!")
    ah <- ah[idx, ]
    ## Get the Ensembl version from the source url.
    ensVers <- .ensVersionFromSourceUrl(ah$sourceurl)
    if(any(ensVers == ensemblVersion(x))){
        ## Got it.
        itIs <- which(ensVers == ensemblVersion(x))
    }else{
        ## Get the "closest" one.
        diffs <- abs(ensVers - as.numeric(ensemblVersion(x)))
        itIs <- which(diffs == min(diffs))[1]
        message("Returning the TwoBit file for Ensembl version ", ensVers[itIs],
                " since no file for Ensembl version ", ensemblVersion(x),
                " is available.")
    }
    ## Getting the ressource.
    Dna <- ah[[names(ah)[itIs]]]
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
        warning("Columns ", paste(sQuote(removed), collapse=", "),
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
    retColumns <- columns
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
    filter <- setFeatureInGRangesFilter(filter, "gene")
    ## If we don't have an order.by define one.
    if(order.by == ""){
        order.by <- NULL
        if(any(columns == "gene_seq_start"))
            order.by <- "gene_seq_start"
        if(any(columns == "seq_name"))
            order.by <- paste(c("seq_name", order.by), collapse=", ")
        if(is.null(order.by))
            order.by <- ""
    }
    Res <- getWhat(x, columns=columns, filter=filter,
                   order.by=order.by, order.type=order.type)
    if(return.type=="data.frame" | return.type=="DataFrame"){
        notThere <- !(retColumns %in% colnames(Res))
        if(any(notThere))
            warning(paste0("Columns ", paste(retColumns[notThere], collapse=", "),
                           " not present in the result data.frame!"))
        retColumns <- retColumns[!notThere]
        Res <- Res[, retColumns]
        if(return.type=="DataFrame")
            Res <- DataFrame(Res)
        return(Res)
    }
    if(return.type=="GRanges"){
        metacols <- columns[ !(columns %in% c("seq_name",
                                              "seq_strand",
                                              "gene_seq_start",
                                              "gene_seq_end")) ]
        suppressWarnings(
            SI <- seqinfo(x)
        )
        SI <- SI[as.character(unique(Res$seq_name))]
        GR <- GRanges(seqnames=Rle(Res$seq_name),
                      ranges=IRanges(start=Res$gene_seq_start, end=Res$gene_seq_end),
                      strand=Rle(Res$seq_strand),
                      seqinfo=SI[as.character(unique(Res$seq_name))],
                      Res[ , metacols, drop=FALSE ]
                    )
        names(GR) <- Res$gene_id
        return(GR)
    }
})

### transcripts:
## get transcripts from the database.
setMethod("transcripts", "EnsDb", function(x, columns=listColumns(x, "tx"),
                                           filter, order.by="", order.type="asc",
                                           return.type="GRanges"){
    return.type <- match.arg(return.type, c("data.frame", "GRanges", "DataFrame"))
    retColumns <- columns
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
    filter <- setFeatureInGRangesFilter(filter, "tx")
    ## If we don't have an order.by define one.
    if(order.by == ""){
        order.by <- NULL
        if(any(columns == "tx_seq_start"))
            order.by <- "tx_seq_start"
        if(any(columns == "seq_name"))
            order.by <- paste(c("seq_name", order.by), collapse=", ")
        if(is.null(order.by))
            order.by <- ""
    }
    Res <- getWhat(x, columns=columns, filter=filter,
                   order.by=order.by, order.type=order.type)
    if(return.type=="data.frame" | return.type=="DataFrame"){
        notThere <- !(retColumns %in% colnames(Res))
        if(any(notThere))
            warning(paste0("Columns ", paste(retColumns[notThere], collapse=", "),
                           " not present in the result data.frame!"))
        retColumns <- retColumns[!notThere]
        Res <- Res[, retColumns]
        if(return.type=="DataFrame")
            Res <- DataFrame(Res)
        return(Res)
    }
    if(return.type=="GRanges"){
        notThere <- !(columns %in% colnames(Res))
        if(any(notThere))
            warning(paste0("Columns ", paste(columns[notThere], collapse=", "),
                           " not present in the result data.frame!"))
        columns <- columns[!notThere]
        metacols <- columns[ !(columns %in% c("seq_name",
                                              "seq_strand",
                                              "tx_seq_start",
                                              "tx_seq_end")) ]
        suppressWarnings(
            SI <- seqinfo(x)
        )
        SI <- SI[as.character(unique(Res$seq_name))]
        GR <- GRanges(seqnames=Rle(Res$seq_name),
                      ranges=IRanges(start=Res$tx_seq_start, end=Res$tx_seq_end),
                      strand=Rle(Res$seq_strand),
                      seqinfo=SI[as.character(unique(Res$seq_name))],
                      Res[ , metacols, drop=FALSE ]
                    )
        names(GR) <- Res$tx_id
        return(GR)
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
    retColumns <- columns
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
    ## If we don't have an order.by define one.
    if(order.by == ""){
        order.by <- NULL
        if(any(columns == "exon_seq_start"))
            order.by <- "exon_seq_start"
        if(any(columns == "seq_name"))
            order.by <- paste(c("seq_name", order.by), collapse=", ")
        if(is.null(order.by))
            order.by <- ""
    }
    filter <- setFeatureInGRangesFilter(filter, "exon")
    Res <- getWhat(x, columns=columns, filter=filter,
                   order.by=order.by, order.type=order.type)
    if(return.type=="data.frame" | return.type=="DataFrame"){
        notThere <- !(retColumns %in% colnames(Res))
        if(any(notThere))
            warning(paste0("Columns ", paste(retColumns[notThere], collapse=", "),
                           " not present in the result data.frame!"))
        retColumns <- retColumns[!notThere]
        Res <- Res[, retColumns]
        if(return.type=="DataFrame")
            Res <- DataFrame(Res)
        return(Res)
    }
    if(return.type=="GRanges"){
        notThere <- !(columns %in% colnames(Res))
        if(any(notThere))
            warning(paste0("Columns ", paste(columns[notThere], collapse=", "),
                           " not present in the result data.frame!"))
        columns <- columns[!notThere]
        metacols <- columns[ !(columns %in% c("seq_name",
                                              "seq_strand",
                                              "exon_seq_start",
                                              "exon_seq_end")) ]
        suppressWarnings(
            SI <- seqinfo(x)
        )
        SI <- SI[as.character(unique(Res$seq_name))]
        GR <- GRanges(seqnames=Rle(Res$seq_name),
                      ranges=IRanges(start=Res$exon_seq_start, end=Res$exon_seq_end),
                      strand=Rle(Res$seq_strand),
                      seqinfo=SI[as.character(unique(Res$seq_name))],
                      Res[ , metacols, drop=FALSE ]
                    )
        names(GR) <- Res$exon_id
        return(GR)
    }
})


## should return a GRangesList
## still considerably slower than the corresponding call in the GenomicFeatures package.
setMethod("exonsBy", "EnsDb", function(x, by=c("tx", "gene"),
                                       columns=listColumns(x, "exon"), filter, use.names=FALSE){
    by <- match.arg(by, c("tx", "gene"))
    bySuff <- "_id"
    if(use.names){
        if(by == "tx"){
            use.names <- FALSE
            warning("Argument use.names ignored as no transcript names are available.")
        }else{
            columns <- unique(c(columns, "gene_name"))
            bySuff <- "_name"
        }
    }
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
    suppressWarnings(
        SI <- seqinfo(x)
    )
    ##Res <- getWhat(dbconn(x), columns=columns, filter=filter, order.by=paste0("seq_name,gene_seq_start,",by ,"_id,exon_idx"))
    if(missing(filter)){
        filter=list()
    }else{
        filter <- checkFilter(filter)
    }
    ## We're applying eventual GRangesFilter to either gene or tx.
    filter <- setFeatureInGRangesFilter(filter, by)
    Res <- getWhat(x, columns=columns, filter=filter, order.by=order.by, skip.order.check=TRUE)
    SI <- SI[as.character(unique(Res$seq_name))]
    ## replace exon_idx with exon_rank
    colnames(Res)[ colnames(Res)=="exon_idx" ] <- "exon_rank"
    columns[ columns=="exon_idx" ] <- "exon_rank"
    notThere <- !(columns %in% colnames(Res))
    if(any(notThere))
        warning(paste0("Columns ", paste(columns[notThere], collapse=", "),
                       " not present in the result data.frame!"))
    columns <- columns[!notThere]
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
    return(split(GR, Res[ , paste0(by, bySuff) ]))
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
    suppressWarnings(
        SI <- seqinfo(x)
    )
    if(missing(filter)){
        filter=list()
    }else{
        filter <- checkFilter(filter)
    }
    filter <- setFeatureInGRangesFilter(filter, by)
    Res <- getWhat(x, columns=columns,
                   filter=filter,
                   order.by=order.by,
                   skip.order.check=TRUE)
    SI <- SI[as.character(unique(Res$seq_name))]
    notThere <- !(columns %in% colnames(Res))
    if(any(notThere))
        warning(paste0("Columns ", paste(columns[notThere], collapse=", "),
                       " not present in the result data.frame!"))
    columns <- columns[!notThere]
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
    return(sum(width(reduce(x))))
##    return(unlist(lapply(width(reduce(x)), sum)))
})

## return the length of genes or transcripts
setMethod("lengthOf", "EnsDb", function(x, of="gene", filter=list()){
    of <- match.arg(of, c("gene", "tx"))
    ## get the exons by gene or transcript from the database...
    suppressWarnings(
        GRL <- exonsBy(x, by=of, filter=filter)
      )
    return(lengthOf(GRL))
})

####============================================================
##  transcriptLengths
##
##  For TxDb: calls just the function (not method!) from the GenomicFeatures
##            package.
##  For EnsDb: calls the .transcriptLengths function.
####------------------------------------------------------------
## setMethod("transcriptLengths", "TxDb", function(x, with.cds_len=FALSE, with.utr5_len=FALSE,
##                                                with.utr3_len=FALSE){
##     return(GenomicFeatures::transcriptLengths(x, with.cds_len=with.cds_len,
##                                               with.utr5_len=with.utr5_len,
##                                               with.utr3_len=with.utr3_len))
## })
## setMethod("transcriptLengths", "EnsDb", function(x, with.cds_len=FALSE, with.utr5_len=FALSE,
##                                                 with.utr3_len=FALSE, filter=list()){
##     return(.transcriptLengths(x, with.cds_len=with.cds_len, with.utr5_len=with.utr3_len,
##                               with.utr3_len=with.utr3_len, filter=filter))
## })
## implement the method from the GenomicFeatures package
.transcriptLengths <- function(x, with.cds_len=FALSE, with.utr5_len=FALSE,
                               with.utr3_len=FALSE, filter=list()){
    ## First we're going to fetch the exonsBy.
    ## Or use getWhat???
    ## Dash, have to make two queries!
    allTxs <- transcripts(x, filter=filter)
    exns <- exonsBy(x, filter=filter)
    ## Match ordering
    exns <- exns[match(allTxs$tx_id, names(exns))]
    ## Calculate length of transcripts.
    txLengths <- sum(width(reduce(exns)))
    ## Calculate no. of exons.
    ## build result data frame:
    Res <- data.frame(tx_id=allTxs$tx_id, gene_id=allTxs$gene_id,
                      nexon=lengths(exns), tx_len=txLengths,
                      stringsAsFactors=FALSE)
    if(!any(c(with.cds_len, with.utr5_len, with.utr3_len))){
        ## Return what we've got thus far.
        return(Res)
    }
    if(with.cds_len)
        Res <- cbind(Res, cds_len=rep(NA, nrow(Res)))
    if(with.utr5_len)
        Res <- cbind(Res, utr5_len=rep(NA, nrow(Res)))
    if(with.utr3_len)
        Res <- cbind(Res, utr3_len=rep(NA, nrow(Res)))
    ## Otherwise do the remaining stuff...
    txs <- allTxs[!is.na(allTxs$tx_cds_seq_start)]
    if(length(txs) > 0){
        cExns <- exns[txs$tx_id]
        cReg <- GRanges(seqnames=seqnames(txs),
                             ranges=IRanges(txs$tx_cds_seq_start,
                                            txs$tx_cds_seq_end),
                             strand=strand(txs),
                             tx_id=txs$tx_id)
        cReg <- split(cReg, f=cReg$tx_id)
        ## Match order.
        cReg <- cReg[match(txs$tx_id, names(cReg))]
        cdsExns <- intersect(cReg, cExns)
        ## cExns: all exons of coding transcripts (includes untranslated
        ##        and translated region)
        ## cReg: just the start-end position of the coding region of the tx.
        ## cdsExns: the coding part of all exons of the tx.
        if(with.cds_len){
            ## Calculate CDS length
            cdsLengths <- sum(width(reduce(cdsExns)))
            Res[names(cdsLengths), "cds_len"] <- cdsLengths
        }
        if(with.utr3_len | with.utr5_len){
            ## ! UTR is the difference between the exons and the cds-exons
            ## Note: order of parameters is important!
            utrReg <- setdiff(cExns, cdsExns)
            leftOfCds <- utrReg[end(utrReg) < start(cReg)]
            rightOfCds <- utrReg[start(utrReg) > end(cReg)]
            ## Calculate lengths.
            leftOfLengths <- sum(width(reduce(leftOfCds)))
            rightOfLengths <- sum(width(reduce(rightOfCds)))
            minusTx <- which(as.character(strand(txs)) == "-" )
            if(with.utr3_len){
                ## Ordering of txs and all other stuff matches.
                tmp <- rightOfLengths
                tmp[minusTx] <- leftOfLengths[minusTx]
                Res[names(tmp), "utr3_len"] <- tmp
            }
            if(with.utr5_len){
                tmp <- leftOfLengths
                tmp[minusTx] <- rightOfLengths[minusTx]
                Res[names(tmp), "utr5_len"] <- tmp
            }
        }
    }
    return(Res)
}

## cdsBy... return coding region ranges by tx or by gene.
setMethod("cdsBy", "EnsDb", function(x, by=c("tx", "gene"),
                                     columns=NULL, filter,
                                     use.names=FALSE){
    by <- match.arg(by, c("tx", "gene"))
    if(missing(filter)){
        filter=list()
    }else{
        filter <- checkFilter(filter)
    }
    filter <- setFeatureInGRangesFilter(filter, by)
    bySuff <- "_id"
    if(by == "tx"){
        ## adding exon_id, exon_idx to the columns.
        columns <- unique(c(columns, "exon_id", "exon_idx"))
        if(use.names)
            warning("Not considering use.names as no transcript names are available.")
    }else{
        if(!is.null(columns))
            warning(paste0("Discarding argument columns as this is not supported for by='gene'."))
        columns <- NULL
        if(use.names){
            bySuff <- "_name"
            columns <- "gene_name"
        }
    }
    byId <- paste0(by, bySuff)
    byIdFull <- unlist(prefixColumns(x, columns=byId,
                                     clean=FALSE), use.names=FALSE)
    order.by <- paste0(byIdFull , ", case when seq_strand=1 then tx_seq_start when seq_strand=-1 then (tx_seq_end * -1) end")
    ## Query the data
    ## what do we need: we need columns tx_cds_seq_start and tx_cds_seq_end and exon_idx
    fetchCols <- unique(c(byId, columns, "tx_cds_seq_start", "tx_cds_seq_end",
                          "seq_name", "seq_strand", "exon_idx", "exon_id", "exon_seq_start",
                          "exon_seq_end"))
    Res <- getWhat(x, columns=fetchCols,
                   filter=filter,
                   order.by=order.by,
                   skip.order.check=TRUE)
    ## Remove rows with NA in tx_cds_seq_start
    Res <- Res[!is.na(Res$tx_cds_seq_start), ]
    ## Remove exons that are not within the cds.
    Res <- Res[Res$exon_seq_end >= Res$tx_cds_seq_start & Res$exon_seq_start <= Res$tx_cds_seq_end,
             , drop=FALSE]
    if(nrow(Res)==0){
        warning("No cds found!")
        return(NULL)
    }
    cdsStarts <- pmax.int(Res$exon_seq_start, Res$tx_cds_seq_start)
    cdsEnds <- pmin.int(Res$exon_seq_end, Res$tx_cds_seq_end)
    ## get the seqinfo:
    suppressWarnings(
        SI <- seqinfo(x)
    )
    SI <- SI[as.character(unique(Res$seq_name))]
    ## Rename columns exon_idx to exon_rank, if present
    if(any(colnames(Res) == "exon_idx")){
        colnames(Res)[colnames(Res) == "exon_idx"] <- "exon_rank"
        columns[columns == "exon_idx"] <- "exon_rank"
    }
    ## Building the result.
    if(length(columns) > 0){
        notThere <- !(columns %in% colnames(Res))
        if(any(notThere))
            warning(paste0("Columns ", paste(columns[notThere], collapse=", "),
                           " not present in the result data.frame!"))
        columns <- columns[!notThere]
        GR <- GRanges(seqnames=Rle(Res$seq_name),
                      strand=Rle(Res$seq_strand),
                      ranges=IRanges(start=cdsStarts, end=cdsEnds),
                      seqinfo=SI,
                      Res[, columns, drop=FALSE])
    }else{
        GR <- GRanges(seqnames=Rle(Res$seq_name),
                      strand=Rle(Res$seq_strand),
                      ranges=IRanges(start=cdsStarts, end=cdsEnds),
                      seqinfo=SI)
    }
    GR <- split(GR, Res[, paste0(by, bySuff)])
    ## For "by gene" we reduce the redundant ranges; that way we loose however all additional
    ## columns!
    if(by == "gene")
        GR <- reduce(GR)
    return(GR)
})


## getUTRsByTranscript
getUTRsByTranscript <- function(x, what, columns=NULL, filter){
    if(missing(filter)){
        filter=list()
    }else{
        filter <- checkFilter(filter)
    }
    filter <- setFeatureInGRangesFilter(filter, "tx")
    columns <- unique(c(columns, "exon_id", "exon_idx"))
    ## what do we need: we need columns tx_cds_seq_start and tx_cds_seq_end and exon_idx
    fetchCols <- unique(c("tx_id", columns, "tx_cds_seq_start", "tx_cds_seq_end",
                          "seq_name", "seq_strand", "exon_idx", "exon_id", "exon_seq_start",
                          "exon_seq_end"))
    by.id.full <- unlist(prefixColumns(x, columns=paste0("tx", "_id"),
                                       clean=FALSE), use.names=FALSE)
    order.by <- paste0(by.id.full ,
                       ", case when seq_strand=1 then tx_seq_start when seq_strand=-1 then (tx_seq_end * -1) end")
    ## get the seqinfo:
    suppressWarnings(
        SI <- seqinfo(x)
    )
    ## Note: doing that with a single query and some coordinate juggling is faster than
    ## calling exonsBy and GRangesList setdiff etc.
    Res <- getWhat(x, columns=fetchCols,
                   filter=filter,
                   order.by=order.by,
                   skip.order.check=TRUE)
    ## Remove rows with NA in tx_cds_seq_start
    Res <- Res[!is.na(Res$tx_cds_seq_start), ]
    ## Remove exons that are within the cds.
    Res <- Res[Res$exon_seq_start < Res$tx_cds_seq_start | Res$exon_seq_end > Res$tx_cds_seq_end,
             , drop=FALSE]
    if(nrow(Res)==0){
        warning(paste0("No ", what, "UTR found!"))
        return(NULL)
    }
    ## Rename columns exon_idx to exon_rank, if present
    if(any(colnames(Res) == "exon_idx")){
        colnames(Res)[colnames(Res) == "exon_idx"] <- "exon_rank"
        columns[columns == "exon_idx"] <- "exon_rank"
    }
    if(what == "five"){
        ## All those on the forward strand for which the exon start is smaller
        ## than the cds start and those on the reverse strand with an exon end
        ## larger than the cds end.
        Res <- Res[(Res$seq_strand > 0 & Res$exon_seq_start < Res$tx_cds_seq_start)
                   | (Res$seq_strand < 0 & Res$exon_seq_end > Res$tx_cds_seq_end), , drop=FALSE]
    }
    if(what == "three"){
        ## Other way round.
        Res <- Res[(Res$seq_strand > 0 & Res$exon_seq_end > Res$tx_cds_seq_end)
                   | (Res$seq_strand < 0 & Res$exon_seq_start < Res$tx_cds_seq_start), , drop=FALSE]
    }
    if(nrow(Res)==0){
        warning(paste0("No ", what, "UTR found!"))
        return(NULL)
    }
    ## Increase the cds end by 1 and decrease the start by 1, thus, avoiding that the UTR
    ## overlaps the cds
    Res$tx_cds_seq_end <- Res$tx_cds_seq_end+1L
    Res$tx_cds_seq_start <- Res$tx_cds_seq_start-1L
    utrStarts <- rep(0, nrow(Res))
    utrEnds <- utrStarts
    ## Distinguish between stuff which is left of and right of the CDS:
    ## Left of the CDS: can be either 5' for + strand or 3' for - strand.
    bm <- which(Res$exon_seq_start <= Res$tx_cds_seq_start)
    if(length(bm) > 0){
        if(what == "five"){
            ## 5' and left of CDS means we're having 5' CDSs
            bm <- bm[Res$seq_strand[bm] > 0]
            if(length(bm) > 0){
                utrStarts[bm] <- Res$exon_seq_start[bm]
                utrEnds[bm] <- pmin.int(Res$exon_seq_end[bm], Res$tx_cds_seq_start[bm])
            }
        }else{
            bm <- bm[Res$seq_strand[bm] < 0]
            if(length(bm) > 0){
                utrStarts[bm] <- Res$exon_seq_start[bm]
                utrEnds[bm] <- pmin.int(Res$exon_seq_end[bm], Res$tx_cds_seq_start[bm])
            }
        }
    }
    ## Right of the CDS: can be either 5' for - strand of 3' for + strand.
    bm <- which(Res$exon_seq_end >= Res$tx_cds_seq_end)
    if(length(bm) > 0){
        if(what == "five"){
            ## Right of CDS is 5' for - strand.
            bm <- bm[Res$seq_strand[bm] < 0]
            if(length(bm) > 0){
                utrStarts[bm] <- pmax.int(Res$exon_seq_start[bm], Res$tx_cds_seq_end[bm])
                utrEnds[bm] <- Res$exon_seq_end[bm]
            }
        }else{
            ## Right of CDS is 3' for + strand
            bm <- bm[Res$seq_strand[bm] > 0]
            if(length(bm) > 0){
                utrStarts[bm] <- pmax.int(Res$exon_seq_start[bm], Res$tx_cds_seq_end[bm])
                utrEnds[bm] <- Res$exon_seq_end[bm]
            }
        }
    }
    notThere <- !(columns %in% colnames(Res))
    if(any(notThere))
        warning(paste0("Columns ", paste(columns[notThere], collapse=", "),
                       " not present in the result data.frame!"))
    columns <- columns[!notThere]
    SI <- SI[as.character(unique(Res$seq_name))]
    GR <- GRanges(seqnames=Rle(Res$seq_name),
                  strand=Rle(Res$seq_strand),
                  ranges=IRanges(start=utrStarts, end=utrEnds),
                  seqinfo=SI,
                  Res[, columns, drop=FALSE])
    GR <- split(GR, Res[, "tx_id"])
    return(GR)
}

## threeUTRsByTranscript
setMethod("threeUTRsByTranscript", "EnsDb", function(x, columns=NULL, filter){
    if(missing(filter)){
        filter=list()
    }else{
        filter <- checkFilter(filter)
    }
    return(getUTRsByTranscript(x=x, what="three", columns=columns, filter=filter))
})

## fiveUTRsByTranscript
setMethod("fiveUTRsByTranscript", "EnsDb", function(x, columns=NULL, filter){
    if(missing(filter)){
        filter=list()
    }else{
        filter <- checkFilter(filter)
    }
    return(getUTRsByTranscript(x=x, what="five", columns=columns, filter=filter))
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
####
## Method that wraps the internal .getWhat function to retrieve data from the
## database. In addition, if present, we're renaming chromosome names depending
## on the ucscChromosomeNames option.
setMethod("getWhat", "EnsDb",
          function(x, columns=c("gene_id", "gene_biotype", "gene_name"),
                   filter=list(), order.by="", order.type="asc",
                   group.by=NULL, skip.order.check=FALSE){
              Res <- .getWhat(x=x,
                              columns=columns,
                              filter=filter,
                              order.by=order.by,
                              order.type=order.type,
                              group.by=group.by,
                              skip.order.check=skip.order.check)
              ## Eventually renaming seqnames according to the specified style.
              if(any(colnames(Res) == "seq_name"))
                  Res$seq_name <- formatSeqnamesFromQuery(x, Res$seq_name)
              ## ## Eventually renaming chromosome names depending on the
              ## ## value of ucscChromosomeNames.
              ## if(any(colnames(Res) == "seq_name")){
              ##     Res$seq_name <- prefixChromName(as.character(Res$seq_name))
              ## }
              return(Res)
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
    if(is(x, "list")){
        if(length(x)==0)
            return(x)
        ## check if all elements are Filter classes.
        IsAFilter <- unlist(lapply(x, function(z){
                                        return(is(z, "BasicFilter"))
                                    }))
        if(any(!IsAFilter))
            stop("One of more elements in filter are not filter objects!")
    }else{
        if(is(x, "BasicFilter")){
            x <- list(x)
        }else{
            stop("filter has to be a filter object or a list of filter objects!")
        }
    }
    return(x)
}

## Fetch data to add as a GeneTrack.
## filter ...                 Used to filter the result.
## chromosome, start, end ... Either all or none has to be specified. If specified, the function
##                            first retrieves all transcripts that have an exon in the specified
##                            range and adds them as a TranscriptidFilter to the filters. The
##                            query to fetch the "real" data is performed after.
## featureIs ...              Wheter gene_biotype or tx_biotype should be mapped to the column
##                            feature.
setMethod("getGeneRegionTrackForGviz", "EnsDb", function(x, filter=list(),
                                                         chromosome=NULL,
                                                         start=NULL,
                                                         end=NULL,
                                                         featureIs="gene_biotype"){
    featureIs <- match.arg(featureIs, c("gene_biotype", "tx_biotype"))
    filter <- checkFilter(filter)
    if(missing(chromosome))
        chromosome <- NULL
    if(missing(start))
        start <- NULL
    if(missing(end))
        end <- NULL
    ## if only chromosome is specified, create a SeqnameFilter and add it to the filter
    if(is.null(start) & is.null(end) & !is.null(chromosome)){
        filter <- c(filter, list(SeqnameFilter(chromosome)))
        chromosome <- NULL
    }
    if(any(c(!is.null(chromosome), !is.null(start), !is.null(end)))){
        ## Require however that all are defined!!!
        if(all(c(!is.null(chromosome), !is.null(start), !is.null(end)))){
            ## Fix eventually provided UCSC chromosome names:
            chromosome <- ucscToEns(chromosome)
            ## Fetch all transcripts in that region:
            tids <- dbGetQuery(dbconn(x),
                               paste0("select distinct tx.tx_id from tx join gene on",
                                      " (tx.gene_id=gene.gene_id)",
                                      " where seq_name='", chromosome, "' and (",
                                      "(tx_seq_start >=",start," and tx_seq_start <=",end,") or ",
                                      "(tx_seq_end >=",start," and tx_seq_end <=",end,") or ",
                                      "(tx_seq_start <=",start," and tx_seq_end >=",end,")",
                                      ")"))[, "tx_id"]
            if(length(tids) == 0)
                stop(paste0("Did not find any transcript on chromosome ", chromosome,
                            " from ", start, " to ", end, "!"))
            filter <- c(filter, TxidFilter(tids))
        }else{
            stop(paste0("Either all or none of arguments 'chromosome', 'start' and 'end' ",
                        " have to be specified!"))
        }
    }
    ## Return a data.frame with columns: chromosome, start, end, width, strand, feature,
    ## gene, exon, transcript and symbol.
    ## 1) Query the data as we usually would.
    ## 2) Perform an additional query to get cds and utr, remove all entries from the
    ##    first result for the same transcripts and rbind the data.frames.
    needCols <- c("seq_name", "exon_seq_start", "exon_seq_end", "seq_strand",
                  featureIs, "gene_id", "exon_id",
                  "exon_idx", "tx_id", "gene_name")
    ## That's the names to which we map the original columns from the EnsDb.
    names(needCols) <- c("chromosome", "start", "end", "strand",
                         "feature", "gene", "exon", "exon_rank", "transcript",
                         "symbol")
    txs <- transcripts(x, filter=filter,
                       columns=needCols, return.type="data.frame")
    ## Rename columns
    idx <- match(needCols, colnames(txs))
    notThere <- is.na(idx)
    idx <- idx[!notThere]
    colnames(txs)[idx] <- names(needCols)[!notThere]
    ## now processing the 5utr
    fUtr <- fiveUTRsByTranscript(x, filter=filter, columns=needCols)
    if(length(fUtr) > 0){
        fUtr <- as(unlist(fUtr, use.names=FALSE), "data.frame")
        fUtr <- fUtr[, !(colnames(fUtr) %in% c("width", "seq_name", "exon_seq_start",
                                               "exon_seq_end", "strand"))]
        colnames(fUtr)[1] <- "chromosome"
        idx <- match(needCols, colnames(fUtr))
        notThere <- is.na(idx)
        idx <- idx[!notThere]
        colnames(fUtr)[idx] <- names(needCols)[!notThere]
        ## Force being in the correct ordering:
        fUtr <- fUtr[, names(needCols)]
        fUtr$feature <- "utr5"
        ## Remove transcripts from the txs data.frame
        txs <- txs[!(txs$transcript %in% fUtr$transcript), , drop=FALSE]
    }
    tUtr <- threeUTRsByTranscript(x, filter=filter, columns=needCols)
    if(length(tUtr) > 0){
        tUtr <- as(unlist(tUtr, use.names=FALSE), "data.frame")
        tUtr <- tUtr[, !(colnames(tUtr) %in% c("width", "seq_name", "exon_seq_start",
                                               "exon_seq_end", "strand"))]
        colnames(tUtr)[1] <- "chromosome"
        idx <- match(needCols, colnames(tUtr))
        notThere <- is.na(idx)
        idx <- idx[!notThere]
        colnames(tUtr)[idx] <- names(needCols)[!notThere]
        ## Force being in the correct ordering:
        tUtr <- tUtr[, names(needCols)]
        tUtr$feature <- "utr3"
        ## Remove transcripts from the txs data.frame
        if(nrow(txs) > 0){
            txs <- txs[!(txs$transcript %in% tUtr$transcript), , drop=FALSE]
        }
    }
    cds <- cdsBy(x, filter=filter, columns=needCols)
    if(length(cds) > 0){
        cds <- as(unlist(cds, use.names=FALSE), "data.frame")
        cds <- cds[, !(colnames(cds) %in% c("width", "seq_name", "exon_seq_start",
                                            "exon_seq_end", "strand"))]
        colnames(cds)[1] <- "chromosome"
        idx <- match(needCols, colnames(cds))
        notThere <- is.na(idx)
        idx <- idx[!notThere]
        colnames(cds)[idx] <- names(needCols)[!notThere]
        ## Force being in the correct ordering:
        cds <- cds[, names(needCols)]
        ## Remove transcripts from the txs data.frame
        if(nrow(txs) > 0){
            txs <- txs[!(txs$transcript %in% cds$transcript), , drop=FALSE]
        }
    }
    if(length(fUtr) > 0){
        txs <- rbind(txs, fUtr)
    }
    if(length(tUtr) > 0){
        txs <- rbind(txs, tUtr)
    }
    if(length(cds) > 0){
        txs <- rbind(txs, cds)
    }
    ## Convert into GRanges.
    suppressWarnings(
        SI <- seqinfo(x)
    )
    SI <- SI[as.character(unique(txs$chromosome))]
    GR <- GRanges(seqnames=Rle(txs$chromosome),
                  strand=Rle(txs$strand),
                  ranges=IRanges(start=txs$start, end=txs$end),
                  seqinfo=SI,
                  txs[, c("feature", "gene", "exon", "exon_rank",
                          "transcript", "symbol"), drop=FALSE])
    return(GR)
})


## Simple helper function to set the @feature in GRangesFilter depending on the calling method.
setFeatureInGRangesFilter <- function(x, feature){
    for(i in seq(along.with=x)){
        if(is(x[[i]], "GRangesFilter")){
            x[[i]]@feature <- feature
        }
    }
    return(x)
}

####============================================================
##  properties
##
##  Get access to the "hidden" .properties slot and return it.
##  This ensures that we're not generating an error for objects that
##  do not have yet that slot.
####------------------------------------------------------------
setMethod("properties", "EnsDb", function(x, ...){
    if(.hasSlot(x, ".properties")){
        return(x@.properties)
    }else{
        warning("The present EnsDb instance has no .properties slot! ",
                "Please use 'updateEnsDb' to update the object!")
        return(list())
    }
})

####============================================================
##  getProperty
##
##  Return the value for the property with the specified name or
##  NA if not present.
####------------------------------------------------------------
setMethod("getProperty", "EnsDb", function(x, name){
    props <- properties(x)
    if(any(names(props) == name)){
        return(props[[name]])
    }else{
        return(NA)
    }
})

####============================================================
##  setProperty
##
##  Sets a property in the object. The value has to be a named vector.
####------------------------------------------------------------
setMethod("setProperty", "EnsDb", function(x, ...){
    dotL <- list(...)
    if(length(dotL) == 0){
        stop("No property specified! The property has to be submitted ",
                "in the format name=value!")
        return(x)
    }
    if(length(dotL) > 1){
        warning("'setProperty' does only support setting of a single property!",
                " Using the first submitted one.")
        dotL <- dotL[1]
    }
    if(is.null(names(dotL)) | names(dotL) == "")
        stop("A name is required! Use name=value!")
    if(.hasSlot(x, ".properties")){
        x@.properties[names(dotL)] <- dotL[[1]]
    }else{
        warning("The present EnsDb instance has no .properties slot! ",
                "Please use 'updateEnsDb' to update the object!")
    }
    return(x)
})

####============================================================
##  updateEnsDb
##
##  Update any "old" EnsDb instance to the most recent implementation.
####------------------------------------------------------------
setMethod("updateEnsDb", "EnsDb", function(x, ...){
    newE <- new("EnsDb", ensdb=x@ensdb, tables=x@tables)
    if(.hasSlot(x, ".properties"))
        newE@.properties <- x@.properties
    return(newE)
})


####============================================================
##  transcriptsByOverlaps
##
##  Just "re-implementing" the transcriptsByOverlaps methods from the
##  GenomicFeature package, finetuning and adapting it for EnsDbs
####------------------------------------------------------------
setMethod("transcriptsByOverlaps", "EnsDb", function(x, ranges, maxgap=0L, minoverlap=1L,
                                                     type=c("any", "start", "end"),
                                                     columns=listColumns(x, "tx"),
                                                     filter){
    if(missing(ranges))
        stop("Parameter 'ranges' is missing!")
    if(missing(filter)){
        filter <- list()
    }else{
        filter <- checkFilter(filter)
    }
    SLs <- unique(as.character(seqnames(ranges)))
    filter <- c(filter, SeqnameFilter(SLs))
    return(subsetByOverlaps(transcripts(x, columns=columns, filter=filter),
           ranges, maxgap=maxgap, minoverlap=minoverlap, type=match.arg(type)))
})

####============================================================
##  exonsByOverlaps
##
####------------------------------------------------------------
setMethod("exonsByOverlaps", "EnsDb", function(x, ranges, maxgap=0L, minoverlap=1L,
                                                type=c("any", "start", "end"),
                                               columns=listColumns(x, "exon"),
                                               filter){
    if(missing(ranges))
        stop("Parameter 'ranges' is missing!")
    if(missing(filter)){
        filter <- list()
    }else{
        filter <- checkFilter(filter)
    }
    SLs <- unique(as.character(seqnames(ranges)))
    filter <- c(filter, SeqnameFilter(SLs))
    return(subsetByOverlaps(exons(x, columns=columns, filter=filter),
           ranges, maxgap=maxgap, minoverlap=minoverlap, type=match.arg(type)))
})



