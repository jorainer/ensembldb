##***********************************************************************
##
##     Methods for EnsDb classes
##
##***********************************************************************
setMethod("show", "EnsDb", function(object) {
    if (is.null(object@ensdb)) {
        cat("Dash it! Got an empty thing!\n")
    } else {
        info <- dbGetQuery(object@ensdb, "select * from metadata")
        cat("EnsDb for Ensembl:\n")
        if (inherits(object@ensdb, "SQLiteConnection"))
            cat(paste0("|Backend: SQLite\n"))
        if (inherits(object@ensdb, "MySQLConnection"))
            cat(paste0("|Backend: MySQL\n"))
        for (i in 1:nrow(info)) {
            cat(paste0("|", info[ i, "name" ], ": ",
                       info[ i, "value" ], "\n"))
        }
        ## gene and transcript info.
        cat(paste0("| No. of genes: ",
                   dbGetQuery(object@ensdb,
                              "select count(distinct gene_id) from gene")[1, 1],
                   ".\n"))
        cat(paste0("| No. of transcripts: ",
                   dbGetQuery(object@ensdb,
                              "select count(distinct tx_id) from tx")[1, 1],
                   ".\n"))
        if (hasProteinData(object))
            cat("|Protein data available.\n")
    }
})

############################################################
## organism
setMethod("organism", "EnsDb", function(object){
    Species <- .getMetaDataValue(object@ensdb, "Organism")
    ## reformat the e.g. homo_sapiens string into Homo sapiens
                                        #
    Species <- gsub(Species, pattern="_", replacement=" ", fixed=TRUE)
    Species <- .organismName(Species)
    return(Species)
})

############################################################
## metadata
setMethod("metadata", "EnsDb", function(x, ...){
    Res <- dbGetQuery(dbconn(x), "select * from metadata")
    return(Res)
})

############################################################
## Validation
##
validateEnsDb <- function(object){
    ## check if the database contains all required tables...
    if(!is.null(object@ensdb)){
        msg <- validMsg(NULL, NULL)
        OK <- dbHasRequiredTables(object@ensdb)
        if (is.character(OK))
            msg <- validMsg(msg, OK)
        OK <- dbHasValidTables(object@ensdb)
        if (is.character(OK))
            msg <- validMsg(msg, OK)
        if (hasProteinData(object)) {
            OK <- dbHasRequiredTables(
                object@ensdb,
                tables = .ensdb_protein_tables(dbSchemaVersion(dbconn(object))))
            if (is.character(OK))
                msg <- validMsg(msg, OK)
            OK <- dbHasValidTables(
                object@ensdb,
                tables = .ensdb_protein_tables(dbSchemaVersion(dbconn(object))))
            if (is.character(OK))
                msg <- validMsg(msg, OK)
            cdsTx <- dbGetQuery(dbconn(object),
                                "select tx_id, tx_cds_seq_start from tx");
            if (is.character(cdsTx$tx_cds_seq_start)) {
                suppressWarnings(
                    cdsTx[, "tx_cds_seq_start"] <- as.numeric(cdsTx$tx_cds_seq_start)
                )
            }
            cdsTx <- cdsTx[!is.na(cdsTx$tx_cds_seq_start), "tx_id"]
            protTx <- dbGetQuery(dbconn(object),
                                 "select distinct tx_id from protein")$tx_id
            if (!all(cdsTx %in% protTx))
                msg <- validMsg(msg, paste0("Not all transcripts with a CDS ",
                                            "are assigned to a protein ID!"))
            if (!all(protTx %in% cdsTx))
                msg <- validMsg(msg, paste0("Not all proteins are assigned to ",
                                            "a transcript with a CDS!"))

        }
        if (is.null(msg)) TRUE
        else msg
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

############################################################
## dbconn
setMethod("dbconn", "EnsDb", function(x){
    return(x@ensdb)
})

############################################################
## ensemblVersion
##
## returns the ensembl version of the package.
setMethod("ensemblVersion", "EnsDb", function(x){
    eVersion <- getMetadataValue(x, "ensembl_version")
    return(eVersion)
})

############################################################
## getMetadataValue
##
## returns the metadata value for the specified name/key
setMethod("getMetadataValue", "EnsDb", function(x, name){
    if(missing(name))
        stop("Argument name has to be specified!")
    return(metadata(x)[metadata(x)$name==name, "value"])
})

############################################################
## seqinfo
setMethod("seqinfo", "EnsDb", function(x){
    Chrs <- dbGetQuery(dbconn(x), "select * from chromosome")
    Chr.build <- .getMetaDataValue(dbconn(x), "genome_build")
    Chrs$seq_name <- formatSeqnamesFromQuery(x, Chrs$seq_name)
    SI <- Seqinfo(seqnames=Chrs$seq_name,
                  seqlengths=Chrs$seq_length,
                  isCircular=Chrs$is_circular==1, genome=Chr.build)
    return(SI)
})

############################################################
## seqlevels
setMethod("seqlevels", "EnsDb", function(x){
    Chrs <- dbGetQuery(dbconn(x), "select distinct seq_name from chromosome")
    Chrs <- formatSeqnamesFromQuery(x, Chrs$seq_name)
    return(Chrs)
})

############################################################
## getGenomeFaFile
##
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

############################################################
##  getGenomeTwoBitFile
##
##  Search and retrieve a genomic DNA resource through a TwoBitFile
##  from AnnotationHub.
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

############################################################
## listTables
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
    Tab <- Tab[tablesByDegree(x, tab=names(Tab))]
    ## Manually add tx_name as a "virtual" column; getWhat will insert the tx_id into that.
    Tab$tx <- unique(c(Tab$tx, "tx_name"))
    ## Manually add the symbol as a "virtual" column.
    Tab$gene <- unique(c(Tab$gene, "symbol"))
    return(Tab)
})

############################################################
## listColumns
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
    Tab <- x@tables
    ## Manually add tx_name as a "virtual" column; getWhat will insert
    ## the tx_id into that.
    Tab$tx <- unique(c(Tab$tx, "tx_name"))
    ## Manually add the symbol as a "virtual" column.
    Tab$gene <- unique(c(Tab$gene, "symbol"))
    if(!missing(table)){
        columns <- unlist(Tab[names(Tab) %in% table], use.names = FALSE)
    }else{
        columns <- unlist(Tab, use.names=FALSE)
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
    return(unique(columns))
})

############################################################
## listGenebiotypes
setMethod("listGenebiotypes", "EnsDb", function(x, ...){
    return(dbGetQuery(dbconn(x), "select distinct gene_biotype from gene")[,1])
})

############################################################
## listTxbiotypes
setMethod("listTxbiotypes", "EnsDb", function(x, ...){
    return(dbGetQuery(dbconn(x), "select distinct tx_biotype from tx")[,1])
})

############################################################
## cleanColumns
##
## checks columns and removes all that are not present in database tables
## the method checks internally whether the columns are in the full form,
## i.e. gene.gene_id (<table name>.<column name>)
setMethod("cleanColumns", "EnsDb", function(x, columns, ...){
    if(missing(columns))
        stop("No columns submitted!")
    ## vote of the majority
    full.name <- length(grep(columns, pattern=".", fixed=TRUE)) >
        floor(length(columns) / 2)
    if (full.name) {
        suppressWarnings(
            full.columns <- unlist(prefixColumns(x,
                                                 unlist(listTables(x)),
                                                 clean = FALSE),
                                   use.names=TRUE)
        )
        bm <- columns %in% full.columns
        removed <- columns[ !bm ]
    } else {
        dbtabs <- names(listTables(x))
        dbtabs <- dbtabs[dbtabs != "metadata"]
        bm <- columns %in% unlist(listTables(x)[dbtabs])
        removed <- columns[!bm]
    }
    if(length(removed) > 0){
        if (length(removed) == 1)
            warning("Column ", paste(sQuote(removed), collapse=", "),
                    " is not present in the database and has been removed")
        else
            warning("Columns ", paste(sQuote(removed), collapse=", "),
                    " are not present in the database and have been removed")
    }
    return(columns[bm])
})

############################################################
## tablesForColumns
##
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

############################################################
## tablesByDegree
##
## returns the table names ordered by degree, i.e. edges to other tables
setMethod("tablesByDegree", "EnsDb", function(x,
                                              tab=names(listTables(x)),
                                              ...){
    Table.order <- c(gene = 1, tx = 2, tx2exon = 3, exon = 4, chromosome = 5,
                     protein = 6, uniprot = 7, protein_domain = 8,
                     entrezgene = 9,
                     metadata = 99)
    Tab <- tab[ order(Table.order[ tab ]) ]
    return(Tab)
})

############################################################
## hasProteinData
##
## Simply check if the database has required tables protein, uniprot
## and protein_domain.
#' @title Determine whether protein data is available in the database
#' 
#' @aliases hasProteinData
#' 
#' @description Determines whether the \code{\linkS4class{EnsDb}}
#'     provides protein annotation data.
#' 
#' @param x The \code{\linkS4class{EnsDb}} object.
#' 
#' @return A logical of length one, \code{TRUE} if protein annotations are
#'     available and \code{FALSE} otherwise.
#' 
#' @author Johannes Rainer
#' 
#' @seealso \code{\link{listTables}}
#' 
#' @examples
#' library(EnsDb.Hsapiens.v75)
#' ## Does this database/package have protein annotations?
#' hasProteinData(EnsDb.Hsapiens.v75)
setMethod("hasProteinData", "EnsDb", function(x) {
    tabs <- listTables(x)
    return(all(c("protein", "uniprot", "protein_domain") %in%
               names(tabs)))
})

############################################################
## genes
##
## get genes from the database.
setMethod("genes", "EnsDb", function(x,
                                     columns = c(listColumns(x, "gene"),
                                                 "entrezid"),
                                     filter = AnnotationFilterList(),
                                     order.by = "",
                                     order.type = "asc",
                                     return.type = "GRanges"){
    return.type <- match.arg(return.type, c("data.frame", "GRanges", "DataFrame"))
    columns <- cleanColumns(x, unique(c(columns, "gene_id")))
    ## if return.type is GRanges we require columns: seq_name, gene_seq_start
    ## and gene_seq_end and seq_strand
    if(return.type=="GRanges"){
        columns <- unique(c(columns, c("gene_seq_start", "gene_seq_end",
                                       "seq_name", "seq_strand")))
    }
    filter <- .processFilterParam(filter, x)
    filter <- setFeatureInGRangesFilter(filter, "gene")
    ## Eventually add columns for the filters:
    columns <- addFilterColumns(columns, filter, x)
    retColumns <- columns
    ## If we don't have an order.by define one.
    if(all(order.by == "")){
        order.by <- NULL
        if (any(columns == "seq_name"))
            order.by <- c(order.by, "seq_name")
        if( any(columns == "gene_seq_start"))
            order.by <- c(order.by, "gene_seq_start")
        if(is.null(order.by))
            order.by <- ""
    }
    Res <- getWhat(x, columns=columns, filter=filter,
                   order.by=order.by, order.type=order.type,
                   startWith = "gene", join = "suggested")
    ## issue #48: collapse entrezid column if dbschema 2.0 is used.
    if (as.numeric(dbSchemaVersion(x)) > 1 & any(columns == "entrezid"))
        Res <- .collapseEntrezidInTable(Res, by = "gene_id")
    if (return.type=="data.frame" | return.type=="DataFrame") {
        notThere <- !(retColumns %in% colnames(Res))
        if(any(notThere))
            warning("Columns ",
                           paste0("'", retColumns[notThere], "'", collapse=", "),
                           " not found in the database!")
        retColumns <- retColumns[!notThere]
        Res <- Res[, retColumns, drop = FALSE]
        if(return.type=="DataFrame")
            Res <- DataFrame(Res)
        return(Res)
    }
    if (return.type=="GRanges") {
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

############################################################
## transcripts:
##
## get transcripts from the database.
setMethod("transcripts", "EnsDb", function(x, columns = listColumns(x, "tx"),
                                           filter = AnnotationFilterList(),
                                           order.by = "", order.type = "asc",
                                           return.type = "GRanges"){
    return.type <- match.arg(return.type, c("data.frame", "GRanges", "DataFrame"))
    columns <- cleanColumns(x, unique(c(columns, "tx_id")))
    ## if return.type is GRanges we require columns: seq_name, gene_seq_start
    ## and gene_seq_end and seq_strand
    if(return.type=="GRanges"){
        columns <- unique(c(columns, c("tx_seq_start",
                                       "tx_seq_end",
                                       "seq_name",
                                       "seq_strand")))
    }
    filter <- .processFilterParam(filter, x)
    filter <- setFeatureInGRangesFilter(filter, "tx")
    ## Eventually add columns for the filters:
    columns <- addFilterColumns(columns, filter, x)
    retColumns <- columns
    ## If we don't have an order.by define one.
    if(all(order.by == "")){
        order.by <- NULL
        if(any(columns == "seq_name"))
            order.by <- c(order.by, "seq_name")
        if(any(columns == "tx_seq_start"))
            order.by <- c(order.by, "tx_seq_start")
        if(is.null(order.by))
            order.by <- ""
    }
    Res <- getWhat(x, columns=columns, filter = filter,
                   order.by=order.by, order.type=order.type,
                   startWith = "tx", join = "suggested")
    ## issue #48: collapse entrezid column if dbschema 2.0 is used.
    if (as.numeric(dbSchemaVersion(x)) > 1 & any(columns == "entrezid"))
        Res <- .collapseEntrezidInTable(Res, by = "tx_id")
    if(return.type=="data.frame" | return.type=="DataFrame"){
        notThere <- !(retColumns %in% colnames(Res))
        if(any(notThere))
            warning("Columns ", paste0("'", retColumns[notThere], "'",
                                       collapse=", "),
                           " not found in the database!")
        retColumns <- retColumns[!notThere]
        Res <- Res[, retColumns, drop = FALSE]
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

############################################################
## promoters:
##
setMethod("promoters", "EnsDb",
          function(x, upstream=2000, downstream=200, ...)
          {
              gr <- transcripts(x, ...)
              trim(suppressWarnings(promoters(gr,
                                              upstream=upstream,
                                              downstream=downstream)))
          }
)

############################################################
## exons
##
## get exons from the database.
setMethod("exons", "EnsDb", function(x, columns = listColumns(x, "exon"),
                                     filter = AnnotationFilterList(),
                                     order.by = "", order.type = "asc",
                                     return.type = "GRanges"){
    return.type <- match.arg(return.type, c("data.frame", "GRanges", "DataFrame"))
    if(!any(columns %in% c(listColumns(x, "exon"), "exon_idx"))){
        ## have to have at least one column from the gene table...
        columns <- c(columns, "exon_id")
    }
    columns <- cleanColumns(x, unique(c(columns, "exon_id")))
    ## if return.type is GRanges we require columns: seq_name, gene_seq_start
    ## and gene_seq_end and seq_strand
    if(return.type=="GRanges"){
        columns <- unique(c(columns, c("exon_seq_start",
                                       "exon_seq_end",
                                       "seq_name",
                                       "seq_strand")))
    }
    filter <- .processFilterParam(filter, x)
    filter <- setFeatureInGRangesFilter(filter, "exon")
    ## Eventually add columns for the filters:
    columns <- addFilterColumns(columns, filter, x)
    retColumns <- columns
    ## If we don't have an order.by define one.
    if (order.by == "") {
        order.by <- NULL
        if (any(columns == "seq_name"))
            order.by <- c(order.by, "seq_name")
        if (any(columns == "exon_seq_start"))
            order.by <- c(order.by, "exon_seq_start")
        if(is.null(order.by))
            order.by <- ""
    }
    Res <- getWhat(x, columns=columns, filter=filter,
                   order.by=order.by, order.type=order.type,
                   startWith = "exon", join = "suggested")
    ## issue #48: collapse entrezid column if dbschema 2.0 is used.
    if (as.numeric(dbSchemaVersion(x)) > 1 & any(columns == "entrezid"))
        Res <- .collapseEntrezidInTable(Res, by = "exon_id")
    if(return.type=="data.frame" | return.type=="DataFrame"){
        notThere <- !(retColumns %in% colnames(Res))
        if(any(notThere))
            warning("Columns ", paste0("'", retColumns[notThere], "'",
                                       collapse=", "),
                           " not found in the database!")
        retColumns <- retColumns[!notThere]
        Res <- Res[, retColumns, drop = FALSE]
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

############################################################
## exonsBy
##
## should return a GRangesList
setMethod("exonsBy", "EnsDb", function(x, by = c("tx", "gene"),
                                       columns = listColumns(x, "exon"),
                                       filter = AnnotationFilterList(),
                                       use.names = FALSE) {
    by <- match.arg(by, c("tx", "gene"))
    bySuff <- "_id"
    if (use.names) {
        if (by == "tx") {
            use.names <- FALSE
            warning("Argument use.names ignored as no transcript names are available.")
        } else {
            columns <- unique(c(columns, "gene_name"))
            bySuff <- "_name"
        }
    }
    filter <- .processFilterParam(filter, x)
    ## We're applying eventual GRangesFilter to either gene or tx.
    filter <- setFeatureInGRangesFilter(filter, by)
    ## Eventually add columns for the filters:
    columns <- cleanColumns(x, unique(c(columns, "exon_id")))
    columns <- addFilterColumns(columns, filter, x)
    ## Quick fix; rename any exon_rank to exon_idx.
    columns[columns == "exon_rank"] <- "exon_idx"

    ## The minimum columns we need, in addition to "columns"
    min.columns <- c(paste0(by, "_id"), "seq_name","exon_seq_start",
                     "exon_seq_end", "exon_id", "seq_strand")
    by.id.full <- unlist(prefixColumns(x, columns = paste0(by, "_id"),
                                        clean = FALSE),
                         use.names = FALSE)
    if (by == "gene") {
        ## tx columns have to be removed, since the same exon can be part of
        ## more than one tx
        txcolumns <- c(listColumns(x, "tx"), "exon_idx")
        txcolumns <- txcolumns[txcolumns != "gene_id"]
        torem <- columns %in% txcolumns
        if (any(torem))
            warning("Columns ",
                    paste(columns[ torem ], collapse = ","),
                    " have been removed as they are not allowed if exons",
                    " are fetched by gene.")
        columns <- columns[!torem]
    } else {
        min.columns <- unique(c(min.columns, "exon_idx"))
        columns <- c(columns, "exon_idx")
    }
    ## define the minimal columns that we need...
    ret_cols <- unique(columns)  ## before adding the "min.columns"
    columns <- unique(c(columns, min.columns))
    ## get the seqinfo:
    suppressWarnings(
        SI <- seqinfo(x)
    )
    ## Resolve ordering problems.
    orderR <- orderResultsInR(x)
    if (orderR) {
        order.by <- ""
    } else {
        if (by == "gene") {
            order.by <- paste0("gene.gene_id, ",
                               "case when seq_strand = 1 then exon_seq_start",
                               " when seq_strand = -1 then (exon_seq_end * -1)",
                               " end")
        } else {
            ## Funny thing is the query takes longer if I use tx2exon.tx_id!
            order.by <- "tx.tx_id, tx2exon.exon_idx"
        }
    }
    Res <- getWhat(x, columns = columns, filter = filter,
                   order.by = order.by, skip.order.check = TRUE,
                   startWith = by, join = "suggested")
    ## issue #48: collapse entrezid column if dbschema 2.0 is used.
    if (as.numeric(dbSchemaVersion(x)) > 1 & any(columns == "entrezid"))
        Res <- .collapseEntrezidInTable(Res, by = "exon_id")
    ## Now, order in R, if not already done in SQL.
    if (orderR) {
        if (by == "gene") {
            startend <- (Res$seq_strand == 1) * Res$exon_seq_start +
                (Res$seq_strand == -1) * (Res$exon_seq_end * -1)
            Res <- Res[order(Res$gene_id, startend,
                             method = "radix"), ]
        } else {
            Res <- Res[order(Res$tx_id, Res$exon_idx,
                             method = "radix"), ]
        }
    }
    SI <- SI[as.character(unique(Res$seq_name))]
    ## replace exon_idx with exon_rank
    colnames(Res)[colnames(Res) == "exon_idx"] <- "exon_rank"
    columns[columns == "exon_idx"] <- "exon_rank"
    ret_cols[ret_cols == "exon_idx"] <- "exon_rank"
    notThere <- !(ret_cols %in% colnames(Res))
    if (any(notThere))
        warning("Columns ", paste0("'", ret_cols[notThere], "'",
                                   collapse = ", "),
                " not found in the database!")
    ret_cols <- ret_cols[!notThere]
    columns.metadata <- ret_cols[!(ret_cols %in% c("seq_name", "seq_strand",
                                                   "exon_seq_start",
                                                   "exon_seq_end"))]
    columns.metadata <- match(columns.metadata, colnames(Res))
    GR <- GRanges(seqnames = Rle(Res$seq_name),
                  strand = Rle(Res$seq_strand),
                  ranges = IRanges(start = Res$exon_seq_start,
                                   end = Res$exon_seq_end),
                  seqinfo = SI,
                  Res[, columns.metadata, drop=FALSE]
                )
    return(split(GR, Res[, paste0(by, bySuff)]))
})

############################################################
## transcriptsBy
##
setMethod("transcriptsBy", "EnsDb", function(x, by = c("gene", "exon"),
                                             columns = listColumns(x, "tx"),
                                             filter = AnnotationFilterList()) {
    if (any(by == "cds"))
        stop("fetching transcripts by cds is not (yet) implemented.")
    by <- match.arg(by, c("gene", "exon"))
    byId <- paste0(by, "_id")
    min.columns <- c(paste0(by, "_id"), "seq_name", "tx_seq_start",
                     "tx_seq_end", "tx_id", "seq_strand")
    ## can not have exon columns!
    ex_cols <- c(listColumns(x, "exon"), "exon_idx")
    ex_cols <- ex_cols[ex_cols != "tx_id"]
    torem <- columns %in% ex_cols
    if (any(torem))
        warning("Columns ",
                paste(columns[ torem ], collapse=","),
                " have been removed as they are not allowed if",
                " transcripts are fetched.")
    columns <- columns[!torem]
    ## Process filters
    filter <- .processFilterParam(filter, x)
    ## GRanges filter should be based on either gene or exon coors.
    filter <- setFeatureInGRangesFilter(filter, by)
    ## Eventually add columns for the filters:
    columns <- addFilterColumns(columns, filter, x)
    ret_cols <- unique(columns)
    ## define the minimal columns that we need...
    columns <- cleanColumns(x, unique(c(columns, min.columns)))
    ## get the seqinfo:
    suppressWarnings(
        SI <- seqinfo(x)
    )
    byIdFull <- unlist(prefixColumns(x, columns = byId, clean = FALSE),
                         use.names = FALSE)
    orderR <- orderResultsInR(x)
    if (orderR) {
        order.by <- ""
    } else {
        order.by <- paste0(byIdFull ,
                           ", case when seq_strand = 1 then tx_seq_start",
                           " when seq_strand = -1 then (tx_seq_end * -1) end")
    }
    Res <- getWhat(x, columns=columns, filter=filter,
                   order.by=order.by, skip.order.check=TRUE,
                   startWith = by, join = "suggested")
    ## issue #48: collapse entrezid column if dbschema 2.0 is used.
    if (as.numeric(dbSchemaVersion(x)) > 1 & any(columns == "entrezid"))
        Res <- .collapseEntrezidInTable(Res, by = "tx_id")
    if (orderR) {
        startEnd <- (Res$seq_strand == 1) * Res$tx_seq_start +
            (Res$seq_strand == -1) * (Res$tx_seq_end * -1)
        Res <- Res[order(Res[, byId], startEnd, method = "radix"), ]
    }
    SI <- SI[as.character(unique(Res$seq_name))]
    ## Replace exon_idx with exon_rank
    colnames(Res) <- gsub(colnames(Res), pattern = "exon_idx",
                                   replacement = "exon_rank", fixed = TRUE)
    ret_cols[ret_cols == "exon_idx"] <- "exon_rank"
    notThere <- !(ret_cols %in% colnames(Res))
    if(any(notThere))
        warning("Columns ", paste0("'", ret_cols[notThere], "'", collapse=", "),
                " not found in the database!")
    ret_cols <- ret_cols[!notThere]
    columns.metadata <- ret_cols[!(ret_cols %in% c("seq_name", "seq_strand",
                                                   "tx_seq_start",
                                                   "tx_seq_end"))]
    columns.metadata <- match(columns.metadata, colnames(Res))
    GR <- GRanges(seqnames=Rle(Res$seq_name),
                  strand=Rle(Res$seq_strand),
                  ranges=IRanges(start=Res$tx_seq_start, end=Res$tx_seq_end),
                  seqinfo=SI,
                  Res[ , columns.metadata, drop=FALSE ]
                )
    return(split(GR, Res[ , byId]))
})

############################################################
## lengthOf
## for GRangesList...
setMethod("lengthOf", "GRangesList", function(x, ...){
    return(sum(width(reduce(x))))
##    return(unlist(lapply(width(reduce(x)), sum)))
})
## return the length of genes or transcripts
setMethod("lengthOf", "EnsDb", function(x, of="gene",
                                        filter=AnnotationFilterList()){
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
                               with.utr3_len=FALSE,
                               filter = AnnotationFilterList()){
    ## First we're going to fetch the exonsBy.
    ## Or use getWhat???
    ## Dash, have to make two queries!
    filter <- .processFilterParam(filter, x)
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

############################################################
## cdsBy
##
## Return coding region ranges by tx or by gene.
setMethod("cdsBy", "EnsDb", function(x, by = c("tx", "gene"),
                                     columns = NULL,
                                     filter = AnnotationFilterList(),
                                     use.names = FALSE){
    by <- match.arg(by, c("tx", "gene"))
    filter <- .processFilterParam(filter, x)
    filter <- setFeatureInGRangesFilter(filter, by)
    columns <- cleanColumns(x, columns)
    ## Eventually add columns for the filters:
    columns <- addFilterColumns(columns, filter, x)
    ## Add a filter ensuring that only coding transcripts are queried.
    filter <- AnnotationFilterList(OnlyCodingTxFilter() ,filter)
    bySuff <- "_id"
    if (by == "tx") {
        ## adding exon_id, exon_idx to the columns.
        columns <- unique(c(columns, "exon_id", "exon_idx"))
        if (use.names)
            warning("Not considering use.names as no transcript names are",
                    " available.")
    } else {
        columns <- unique(c("gene_id", columns))
        if( use.names) {
            bySuff <- "_name"
            columns <- c(columns, "gene_name")
        }
    }
    byId <- paste0(by, bySuff)
    ## Query the data
    fetchCols <- unique(c(byId, columns, "tx_cds_seq_start", "tx_cds_seq_end",
                          "seq_name", "seq_strand", "exon_idx", "exon_id",
                          "exon_seq_start", "exon_seq_end"))
    ## Ordering of the results:
    ## Force ordering in R by default here to fix issue #11
    ##orderR <- orderResultsInR(x)
    orderR <- TRUE
    if (orderR) {
        order.by <- ""
    } else {
        if (by == "tx") {
            ## Here we want to sort the exons by exon_idx
            order.by <- "tx.tx_id, tx2exon.exon_idx"
        } else {
            ## Here we want to sort the transcripts by tx start.
            order.by <- paste0("gene.gene_id, case when seq_strand = 1 then",
                               " tx_cds_seq_start when seq_strand = -1 then",
                               "(tx_cds_seq_end * -1) end")
        }
    }
    Res <- getWhat(x, columns = fetchCols,
                   filter = filter,
                   order.by = order.by,
                   skip.order.check = TRUE,
                   startWith = by, join = "suggested")
    ## issue #48: collapse entrezid column if dbschema 2.0 is used.
    if (as.numeric(dbSchemaVersion(x)) > 1 & any(columns == "entrezid"))
        Res <- .collapseEntrezidInTable(Res, by = "exon_id")
    ## Remove rows with NA in tx_cds_seq_start; that's the case for "old"
    ## databases.
    nas <- is.na(Res$tx_cds_seq_start)
    if (any(nas))
        Res <- Res[!nas, ]
    ## Remove exons that are not within the cds.
    Res <- Res[Res$exon_seq_end >= Res$tx_cds_seq_start &
               Res$exon_seq_start <= Res$tx_cds_seq_end,
             , drop = FALSE]
    if (orderR) {
        ## And finally ordering them.
        if (by == "tx") {
            Res <- Res[order(Res$tx_id, Res$exon_idx, method = "radix"), ]
        } else {
            startend <- (Res$seq_strand == 1) * Res$tx_cds_seq_start +
                (Res$seq_strand == -1) * (Res$tx_cds_seq_end * -1)
            Res <- Res[order(Res$gene_id, startend, method = "radix"), ]
        }
    }
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
    ## For "by gene" we reduce the redundant ranges;
    ## that way we loose however all additional columns!
    if(by == "gene")
        GR <- reduce(GR)
    return(GR)
})


############################################################
## getUTRsByTranscript
##
getUTRsByTranscript <- function(x, what, columns = NULL,
                                filter = AnnotationFilterList()) {
    filter <- .processFilterParam(filter, x)
    columns <- cleanColumns(x, columns)
    filter <- setFeatureInGRangesFilter(filter, "tx")
    ## Eventually add columns for the filters:
    columns <- addFilterColumns(columns, filter, x)
    columns <- unique(c(columns, "exon_id", "exon_idx"))
    ## Add the filter for coding tx only.
    filter <- AnnotationFilterList(OnlyCodingTxFilter(), filter)
    ## what do we need: tx_cds_seq_start, tx_cds_seq_end and exon_idx
    fetchCols <- unique(c("tx_id", columns, "tx_cds_seq_start",
                          "tx_cds_seq_end", "seq_name", "seq_strand",
                          "exon_seq_start", "exon_seq_end"))
    order.by <- "tx.tx_id"
    ## get the seqinfo:
    suppressWarnings(
        SI <- seqinfo(x)
    )
    ## Note: doing that with a single query and some coordinate juggling
    ## is faster than calling exonsBy and GRangesList setdiff etc.
    Res <- getWhat(x, columns=fetchCols,
                   filter=filter,
                   order.by=order.by,
                   skip.order.check=TRUE,
                   startWith = "tx", join = "suggested")
    ## issue #48: collapse entrezid column if dbschema 2.0 is used.
    if (as.numeric(dbSchemaVersion(x)) > 1 & any(columns == "entrezid"))
        Res <- .collapseEntrezidInTable(Res, by = "exon_id")
    nas <- is.na(Res$tx_cds_seq_start)
    if (any(nas))
        Res <- Res[!nas, ]
    ## Remove exons that are within the cds.
    Res <- Res[Res$exon_seq_start < Res$tx_cds_seq_start |
               Res$exon_seq_end > Res$tx_cds_seq_end, , drop=FALSE]
    if (nrow(Res) == 0) {
        warning(paste0("No ", what, "UTR found!"))
        return(NULL)
    }
    ## Rename columns exon_idx to exon_rank, if present
    if (any(colnames(Res) == "exon_idx")) {
        colnames(Res) <- sub(colnames(Res), pattern = "exon_idx",
                             replacement = "exon_rank", fixed = TRUE)
        columns[columns == "exon_idx"] <- "exon_rank"
    }
    if (what == "five") {
        ## All those on the forward strand for which the exon start is smaller
        ## than the cds start and those on the reverse strand with an exon end
        ## larger than the cds end.
        Res <- Res[(Res$seq_strand > 0 & Res$exon_seq_start < Res$tx_cds_seq_start)
                   | (Res$seq_strand < 0 & Res$exon_seq_end > Res$tx_cds_seq_end),
                 , drop=FALSE]
    } else {
        ## Other way round.
        Res <- Res[(Res$seq_strand > 0 & Res$exon_seq_end > Res$tx_cds_seq_end) |
                   (Res$seq_strand < 0 & Res$exon_seq_start < Res$tx_cds_seq_start),
                 , drop=FALSE]
    }
    if (nrow(Res) == 0) {
        warning(paste0("No ", what, "UTR found!"))
        return(NULL)
    }
    ## Increase the cds end by 1 and decrease the start by 1, thus,
    ## avoiding that the UTR overlaps the cds
    Res$tx_cds_seq_end <- Res$tx_cds_seq_end + 1L
    Res$tx_cds_seq_start <- Res$tx_cds_seq_start - 1L
    utrStarts <- rep(0, nrow(Res))
    utrEnds <- utrStarts
    ## Distinguish between stuff which is left of and right of the CDS:
    ## Left of the CDS: can be either 5' for + strand or 3' for - strand.
    bm <- which(Res$exon_seq_start <= Res$tx_cds_seq_start)
    if (length(bm) > 0) {
        if (what == "five") {
            ## 5' and left of CDS means we're having 5' CDSs
            bm <- bm[Res$seq_strand[bm] > 0]
            if(length(bm) > 0){
                utrStarts[bm] <- Res$exon_seq_start[bm]
                utrEnds[bm] <- pmin.int(Res$exon_seq_end[bm],
                                        Res$tx_cds_seq_start[bm])
            }
        } else {
            bm <- bm[Res$seq_strand[bm] < 0]
            if (length(bm) > 0) {
                utrStarts[bm] <- Res$exon_seq_start[bm]
                utrEnds[bm] <- pmin.int(Res$exon_seq_end[bm],
                                        Res$tx_cds_seq_start[bm])
            }
        }
    }
    ## Right of the CDS: can be either 5' for - strand of 3' for + strand.
    bm <- which(Res$exon_seq_end >= Res$tx_cds_seq_end)
    if (length(bm) > 0) {
        if (what == "five") {
            ## Right of CDS is 5' for - strand.
            bm <- bm[Res$seq_strand[bm] < 0]
            if (length(bm) > 0) {
                utrStarts[bm] <- pmax.int(Res$exon_seq_start[bm],
                                          Res$tx_cds_seq_end[bm])
                utrEnds[bm] <- Res$exon_seq_end[bm]
            }
        } else {
            ## Right of CDS is 3' for + strand
            bm <- bm[Res$seq_strand[bm] > 0]
            if (length(bm) > 0) {
                utrStarts[bm] <- pmax.int(Res$exon_seq_start[bm],
                                          Res$tx_cds_seq_end[bm])
                utrEnds[bm] <- Res$exon_seq_end[bm]
            }
        }
    }
    notThere <- !(columns %in% colnames(Res))
    if (any(notThere))
        warning(paste0("Columns ", paste(columns[notThere], collapse=", "),
                       " not present in the result data.frame!"))
    columns <- columns[!notThere]
    SI <- SI[as.character(unique(Res$seq_name))]
    GR <- GRanges(seqnames = Rle(Res$seq_name),
                  strand = Rle(Res$seq_strand),
                  ranges = IRanges(start=utrStarts, end=utrEnds),
                  seqinfo = SI,
                  Res[, columns, drop = FALSE])
    GR <- split(GR, Res[, "tx_id"])
    return(GR)
}

############################################################
## threeUTRsByTranscript
##
setMethod("threeUTRsByTranscript", "EnsDb",
          function(x, columns = NULL, filter = AnnotationFilterList()) {
              filter <- .processFilterParam(filter, x)
              getUTRsByTranscript(x = x, what = "three", columns = columns,
                                  filter = filter)
})

############################################################
## fiveUTRsByTranscript
##
setMethod("fiveUTRsByTranscript", "EnsDb",
          function(x, columns = NULL, filter = AnnotationFilterList()) {
    filter <- .processFilterParam(filter, x)
    getUTRsByTranscript(x = x, what = "five", columns = columns,
                        filter = filter)
})

############################################################
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

############################################################
## buildQuery
setMethod("buildQuery", "EnsDb",
          function(x, columns=c("gene_id", "gene_biotype", "gene_name"),
                   filter = AnnotationFilterList(), order.by="",
                   order.type="asc",
                   skip.order.check=FALSE){
              return(.buildQuery(x=x,
                                 columns=columns,
                                 filter=filter,
                                 order.by=order.by,
                                 order.type=order.type,
                                 skip.order.check=skip.order.check))
          })

############################################################
## getWhat
##
## Method that wraps the internal .getWhat function to retrieve data from the
## database. In addition, if present, we're renaming chromosome names depending
## on the ucscChromosomeNames option.
## Additional parameters:
## o startWith: the name of the database table from which the join should start
##   or NULL for the default behaviour (i.e. genes-> tx etc).
## o join: the type of join that should be used; one of "join",
##   "left outer join" or "suggested".
setMethod("getWhat", "EnsDb",
          function(x, columns = c("gene_id", "gene_biotype", "gene_name"),
                   filter = AnnotationFilterList(), order.by = "",
                   order.type = "asc", group.by = NULL,
                   skip.order.check = FALSE, startWith = NULL,
                   join = "suggested") {
              Res <- .getWhat(x = x,
                              columns = columns,
                              filter = filter,
                              order.by = order.by,
                              order.type = order.type,
                              group.by = group.by,
                              skip.order.check = skip.order.check,
                              startWith = startWith,
                              join = join)
              ## Eventually renaming seqnames according to the specified style.
              if(any(colnames(Res) == "seq_name"))
                  Res$seq_name <- formatSeqnamesFromQuery(x, Res$seq_name)
              return(Res)
          })

############################################################
## disjointExons
##
## that's similar to the code from the GenomicFeatures package.
setMethod("disjointExons", "EnsDb",
          function(x, aggregateGenes = FALSE, includeTranscripts = TRUE,
                   filter = AnnotationFilterList(), ...){
              filter <- .processFilterParam(filter, x)

              exonsByGene <- exonsBy(x, by = "gene", filter = filter)
              exonicParts <- disjoin(unlist(exonsByGene, use.names = FALSE))

              if (aggregateGenes) {
                  foGG <- findOverlaps(exonsByGene, exonsByGene)
                  aggregateNames <- GenomicFeatures:::.listNames(names(exonsByGene),
                                                                 as.list(foGG))
                  foEG <- findOverlaps(exonicParts, exonsByGene, select="first")
                  gene_id <- aggregateNames[foEG]
                  pasteNames <- GenomicFeatures:::.pasteNames(names(exonsByGene),
                                                              as.list(foGG))[foEG]
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
                  gene_id <- GenomicFeatures:::.listNames(names(exonsByGene),
                                                          idxList)
                  orderByGeneName <- order(unlist(gene_id, use.names=FALSE))
                  exonic_rle <- runLength(Rle(unlist(gene_id[orderByGeneName],
                                                     use.names=FALSE)))
              }
              values <- DataFrame(gene_id)

              if (includeTranscripts) {
                  exonsByTx <- exonsBy(x, by="tx", filter=filter)
                  foET <- findOverlaps(exonicParts, exonsByTx)
                  values$tx_name <- GenomicFeatures:::.listNames(names(exonsByTx),
                                                                 as.list(foET))
              }
              mcols(exonicParts) <- values
              exonicParts <- exonicParts[orderByGeneName]
              exonic_part <- unlist(lapply(exonic_rle, seq_len), use.names=FALSE)
              exonicParts$exonic_part <- exonic_part
              return(exonicParts)
          }
         )

############################################################
## getGeneRegionTrackForGviz
## Fetch data to add as a GeneTrack.
## filter ...                 Used to filter the result.
## chromosome, start, end ... Either all or none has to be specified. If
##                            specified, the function first retrieves all
##                            transcripts that have an exon in the specified
##                            range and adds them as a TranscriptidFilter to
##                            the filters. The query to fetch the "real" data
##                            is performed afterwards.
## featureIs ...              Wheter gene_biotype or tx_biotype should be
##                            mapped to the column feature.
setMethod(
    "getGeneRegionTrackForGviz",
    "EnsDb",
    function(x, filter = AnnotationFilterList(), chromosome = NULL,
             start = NULL, end = NULL, featureIs = "gene_biotype")
    {
        featureIs <- match.arg(featureIs, c("gene_biotype", "tx_biotype"))
        filter <- .processFilterParam(filter, x)
        if(missing(chromosome))
            chromosome <- NULL
        if(missing(start))
            start <- NULL
        if(missing(end))
            end <- NULL
        ## If only chromosome is specified, create a SeqNameFilter and
        ## add it to the filter
        if(is.null(start) & is.null(end) & !is.null(chromosome)){
            filter <- AnnotationFilterList(filter, SeqNameFilter(chromosome))
            chromosome <- NULL
        }
        if(any(c(!is.null(chromosome), !is.null(start), !is.null(end)))){
            ## Require however that all are defined!!!
            if(all(c(!is.null(chromosome), !is.null(start), !is.null(end)))){
                ## Fix eventually provided UCSC chromosome names:
                chromosome <- ucscToEns(chromosome)
                ## Define a GRangesFilter to include all features that overlap
                ## that region.
                grg <- GRangesFilter(GRanges(seqnames = chromosome,
                                             ranges = IRanges(start, end)),
                                     feature = "tx", type = "any")
                tids <- transcripts(x, filter = grg, columns = "tx_id")$tx_id
                filter <- AnnotationFilterList(filter, TxIdFilter(tids))
            }else{
                stop("Either all or none of arguments 'chromosome', 'start' and",
                     " 'end' have to be specified!")
            }
        }
        ## Return a data.frame with columns: chromosome, start, end, width,
        ## strand, feature,
        ## gene, exon, transcript and symbol.
        ## 1) Query the data as we usually would.
        ## 2) Perform an additional query to get cds and utr, remove all entries
        ##    from the first result for the same transcripts and rbind the
        ##    data.frames.
        needCols <- c("seq_name", "exon_seq_start", "exon_seq_end", "seq_strand",
                      featureIs, "gene_id", "exon_id",
                      "exon_idx", "tx_id", "gene_name")
        ## That's the names to which we map the original columns from the EnsDb.
        names(needCols) <- c("chromosome", "start", "end", "strand",
                             "feature", "gene", "exon", "exon_rank", "transcript",
                             "symbol")
        txs <- transcripts(x, filter = filter,
                           columns = needCols, return.type="data.frame")
        ## Rename columns
        idx <- match(needCols, colnames(txs))
        notThere <- is.na(idx)
        idx <- idx[!notThere]
        colnames(txs)[idx] <- names(needCols)[!notThere]
        ## now processing the 5utr
        fUtr <- fiveUTRsByTranscript(x, filter = filter, columns=needCols)
        if(length(fUtr) > 0){
            fUtr <- as(unlist(fUtr, use.names=FALSE), "data.frame")
            fUtr <- fUtr[, !(colnames(fUtr) %in% c("width", "seq_name",
                                                   "exon_seq_start",
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
        tUtr <- threeUTRsByTranscript(x, filter = filter, columns=needCols)
        if(length(tUtr) > 0){
            tUtr <- as(unlist(tUtr, use.names=FALSE), "data.frame")
            tUtr <- tUtr[, !(colnames(tUtr) %in% c("width", "seq_name",
                                                   "exon_seq_start",
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
        cds <- cdsBy(x, filter = filter, columns = needCols)
        if(length(cds) > 0){
            cds <- as(unlist(cds, use.names=FALSE), "data.frame")
            cds <- cds[, !(colnames(cds) %in% c("width", "seq_name",
                                                "exon_seq_start",
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
setMethod("getProperty", "EnsDb", function(x, name, default = NA){
    props <- properties(x)
    if(any(names(props) == name)){
        return(props[[name]])
    }else{
        return(default)
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

#' remove the property with the specified name.
#' @noRd
dropProperty <- function(x, name) {
    if (missing(name))
        return(x)
    prps <- x@.properties
    if (any(names(prps) == name))
        prps <- prps[names(prps) != name]
    x@.properties <- prps
    x
}

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
setMethod("transcriptsByOverlaps", "EnsDb",
          function(x, ranges, maxgap = 0L, minoverlap = 1L,
                   type = c("any", "start", "end"),
                   columns = listColumns(x, "tx"),
                   filter = AnnotationFilterList()) {
    if (missing(ranges))
        stop("Parameter 'ranges' is missing!")
    filter <- .processFilterParam(filter, x)
    SLs <- unique(as.character(seqnames(ranges)))
    filter <- AnnotationFilterList(filter, SeqNameFilter(SLs))
    columns <- cleanColumns(x, columns)
    subsetByOverlaps(transcripts(x, columns = columns, filter = filter),
                     ranges, maxgap = maxgap, minoverlap = minoverlap,
                     type = match.arg(type))
})

####============================================================
##  exonsByOverlaps
##
####------------------------------------------------------------
setMethod("exonsByOverlaps", "EnsDb",
          function(x, ranges, maxgap = 0L, minoverlap = 1L,
                   type = c("any", "start", "end"),
                   columns = listColumns(x, "exon"),
                   filter = AnnotationFilterList()) {
    if(missing(ranges))
        stop("Parameter 'ranges' is missing!")
    filter <- .processFilterParam(filter, x)
    SLs <- unique(as.character(seqnames(ranges)))
    filter <- AnnotationFilterList(filter, SeqNameFilter(SLs))
    columns <- cleanColumns(x, columns)
    subsetByOverlaps(exons(x, columns = columns, filter = filter),
                     ranges, maxgap = maxgap, minoverlap = minoverlap,
                     type = match.arg(type))
})

############################################################
## returnFilterColumns
##
## Method to set the option whether or not the filter columns should be
## returned too.
setMethod("returnFilterColumns", "EnsDb", function(x) {
    return(getProperty(x, "returnFilterColumns"))
})
setReplaceMethod("returnFilterColumns", "EnsDb", function(x, value) {
    if(!is.logical(value))
        stop("'value' has to be a logical!")
    if(length(value) > 1)
        stop("'value' has to be a logical of length 1!")
    x <- setProperty(x, returnFilterColumns=value)
    return(x)
})

############################################################
## orderResultsInR
##
## Whether the results should be ordered in R instead of in the
## SQL call
setMethod("orderResultsInR", "EnsDb", function(x) {
    return(getProperty(x, "orderResultsInR", default = FALSE))
})
setReplaceMethod("orderResultsInR", "EnsDb", function(x, value) {
    if(!is.logical(value))
        stop("'value' has to be a logical!")
    if(length(value) > 1)
        stop("'value' has to be a logical of length 1!")
    x <- setProperty(x, orderResultsInR = value)
    return(x)
})

############################################################
## useMySQL
##
## Switch from RSQlite backend to a MySQL backend.
#' @title Use a MySQL backend
#' 
#' @aliases useMySQL
#'
#' @description Change the SQL backend from \emph{SQLite} to \emph{MySQL}.
#'     When first called on an \code{\linkS4class{EnsDb}} object, the function
#'     tries to create and save all of the data into a MySQL database. All
#'     subsequent calls will connect to the already existing MySQL database.
#'
#' @details This functionality requires that the \code{RMySQL} package is
#'     installed and that the user has (write) access to a running MySQL server.
#'     If the corresponding database does already exist users without write
#'     access can use this functionality.
#'
#' @note At present the function does not evaluate whether the versions
#'     between the SQLite and MySQL database differ.
#'
#' @param x The \code{\linkS4class{EnsDb}} object.
#' 
#' @param host Character vector specifying the host on which the MySQL
#'     server runs.
#' 
#' @param port The port on which the MySQL server can be accessed.
#'
#' @param user The user name for the MySQL server.
#'
#' @param pass The password for the MySQL server.
#'
#' @return A \code{\linkS4class{EnsDb}} object providing access to the
#'      data stored in the MySQL backend.
#'
#' @author Johannes Rainer
#'
#' @examples
#' ## Load the EnsDb database (SQLite backend).
#' library(EnsDb.Hsapiens.v75)
#' edb <- EnsDb.Hsapiens.v75
#' ## Now change the backend to MySQL; my_user and my_pass should
#' ## be the user name and password to access the MySQL server.
#' \dontrun{
#' edb_mysql <- useMySQL(edb, host = "localhost", user = my_user, pass = my_pass)
#' }
setMethod("useMySQL", "EnsDb", function(x, host = "localhost",
                                        port = 3306, user, pass) {
    if (missing(user))
        stop("'user' has to be specified.")
    if (missing(pass))
        stop("'pass' has to be specified.")
    ## Check if RMySQL package is available.
    if(requireNamespace("RMySQL", quietly = TRUE)) {
        ## Check if we can connect to MySQL.
        driva <- dbDriver("MySQL")
        con <- dbConnect(driva, host = host, user = user, pass = pass,
                         port = port)
        ## Check if database is available.
        dbs <- dbGetQuery(con, "show databases;")
        ## sqliteName should be in the format EnsDb.Hsapiens.v75!
        sqliteName <- .makePackageName(dbconn(x))
        ## sqliteName <- sub(basename(dbfile(dbconn(x))),
        ##                   pattern = ".sqlite", replacement = "",
        ##                   fixed = TRUE)
        mysqlName <- SQLiteName2MySQL(sqliteName)
        if (nrow(dbs) == 0 | !any(dbs$Database == mysqlName)) {
            message("Database not available, trying to create it...",
                    appendLF = FALSE)
            dbGetQuery(con, paste0("create database ", mysqlName))
            message("OK")
        }
        dbDisconnect(con)
        ## Connect to the database and check if we've got all tables.
        con <- dbConnect(driva, host = host, user = user, pass = pass,
                         dbname = mysqlName)
        ## If we've got no tables we try to feed the SQLite database
        if (length(dbListTables(con)) == 0)
            feedEnsDb2MySQL(x, mysql = con)
        ## Check if we've got all required tables.
        OK <- dbHasRequiredTables(con)
        if (is.character(OK))
            stop(OK)
        OK <- dbHasValidTables(con)
        if (is.character(OK))
            stop(OK)
        ## Check if the versions/creation date differ.
        metadata_pkg <- metadata(x)
        ## Now store the connection into the @ensdb slot
        ## dbDisconnect(x@ensdb)
        ## x@ensdb <- NULL
        x@ensdb <- con
        metadata_db <- metadata(x)
        cre_pkg <- metadata_pkg[metadata_pkg$name == "Creation time", "value"]
        cre_db <- metadata_db[metadata_db$name == "Creation time", "value"]
        if (cre_pkg != cre_db) {
            message("Creation date between the package and the information in",
                    " the database differ:\n o package: ", cre_pkg,
                    "\n o database: ", cre_db, ".\nYou might consider to delete",
                    " the database and re-install it calling this function.")
        }
        return(x)
    } else {
        stop("Package 'RMySQL' not available.")
    }
})

############################################################
## proteins
##
## If return type is GRanges, make a seqlevel and seqinfo for each protein, i.e.
## put each protein on its own sequence.
#' @title Protein related functionality
#' 
#' @aliases proteins
#'
#' @description This help page provides information about most of the
#'     functionality related to protein annotations in \code{ensembldb}.
#'
#'     The \code{proteins} method retrieves protein related annotations from
#'     an \code{\linkS4class{EnsDb}} database.
#'
#' @details The \code{proteins} method performs the query starting from the
#'     \code{protein} tables and can hence return all annotations from the
#'     database that are related to proteins and transcripts encoding these
#'     proteins from the database. Since \code{proteins} does thus only query
#'     annotations for protein coding transcripts, the \code{\link{genes}} or
#'     \code{\link{transcripts}} methods have to be used to retrieve annotations
#'     for non-coding transcripts.
#' 
#' @param object The \code{\linkS4class{EnsDb}} object.
#'
#' @param columns For \code{proteins}: character vector defining the columns to
#'     be extracted from the database. Can be any column(s) listed by the
#'     \code{\link{listColumns}} method.
#'
#' @param filter For \code{proteins}: A filter object extending
#'     \code{AnnotationFilter} or a list of such objects to select
#'     specific entries from the database. See \code{\link{Filter-classes}} for
#'     a documentation of available filters and use
#'     \code{\link{supportedFilters}} to get the full list of supported filters.
#'
#' @param order.by For \code{proteins}: a character vector specifying the
#'     column(s) by which the result should be ordered.
#'
#' @param order.type For \code{proteins}: if the results should be ordered
#'     ascending (\code{order.type = "asc"}) or descending
#'     (\code{order.type = "desc"})
#'
#' @param return.type For \code{proteins}: character of lenght one specifying
#'     the type of the returned object. Can be either \code{"DataFrame"},
#'     \code{"data.frame"} or \code{"AAStringSet"}.
#'
#' @return The \code{proteins} method returns protein related annotations from
#'     an \code{\linkS4class{EnsDb}} object with its \code{return.type} argument
#'     allowing to define the type of the returned object. Note that if
#'     \code{return.type = "AAStringSet"} additional annotation columns are
#'     stored in a \code{DataFrame} that can be accessed with the \code{mcols}
#'     method on the returned object.
#'
#' @rdname ProteinFunctionality
#' 
#' @author Johannes Rainer
#'
#' @examples
#' library(ensembldb)
#' library(EnsDb.Hsapiens.v75)
#' edb <- EnsDb.Hsapiens.v75
#' ## Get all proteins from tha database for the gene ZBTB16, if protein
#' ## annotations are available
#' if (hasProteinData(edb))
#'     proteins(edb, filter = GenenameFilter("ZBTB16"))
setMethod("proteins", "EnsDb", function(object,
                                        columns = listColumns(object, "protein"),
                                        filter = AnnotationFilterList(),
                                        order.by = "",
                                        order.type = "asc",
                                        return.type = "DataFrame") {
    if (!hasProteinData(object))
        stop("The used EnsDb does not provide protein annotations!",
             " Thus, 'proteins' can not be used.")
    return.type <- match.arg(return.type, c("DataFrame", "AAStringSet",
                                            "data.frame"))
    columns <- cleanColumns(object, unique(c(columns, "protein_id")))
    filter <- .processFilterParam(filter, object)
    filter <- setFeatureInGRangesFilter(filter, "tx")
    ## Eventually add columns for the filters:
    columns <- addFilterColumns(columns, filter, object)
    ## Check that we don't have any exon columns here.
    ex_cols <- unique(listColumns(object, c("exon", "tx2exon")))
    ex_cols <- ex_cols[ex_cols != "tx_id"]
    if (any(columns %in% ex_cols)) {
        warning("Exon specific columns are not allowed for proteins. Columns ",
                paste0("'", columns[columns %in% ex_cols], "'", collapse = ", "),
                " have been removed.")
        columns <- columns[!(columns %in% ex_cols)]
    }
    retColumns <- columns
    ## Process order.by:
    ## If not specified we might want to order them by seq_name or tx_seq_start
    ## if present in parameter columns
    if (all(order.by == "")) {
        order.by <- NULL
        if (any(columns == "seq_name"))
            order.by <- "seq_name"
        seq_col_idx <- grep(columns, pattern = "_seq_")
        if (length(seq_col_idx) > 0)
            order.by <- c(order.by, columns[seq_col_idx[1]])
        if (is.null(order.by))
            order.by <- ""
    }
    ## If we're going to return a GRanges we need to know the length of the
    ## peptide sequence.
    if (return.type == "AAStringSet") {
        columns <- unique(c(columns, "protein_sequence"))
    }
    ## protein_id is *always* required
    columns <- unique(c(columns), "protein_id")
    ## Get the data
    Res <- getWhat(object, columns = columns, filter = filter,
                   order.by = order.by, order.type = order.type,
                   startWith = "protein", join = "suggested")
    ## issue #48: collapse entrezid column if dbschema 2.0 is used.
    if (as.numeric(dbSchemaVersion(object)) > 1 & any(columns == "entrezid"))
        Res <- .collapseEntrezidInTable(Res, by = "protein_id")
    ## Now process the result.
    cols_not_found <- !(retColumns %in% colnames(Res))
    retColumns <- retColumns[!cols_not_found]
    if (any(cols_not_found))
        warning("Columns ", paste0("'", retColumns[cols_not_found], "'",
                                   collapse = ", "),
                " not found in the database!")
    if (return.type == "AAStringSet") {
        aass <- AAStringSet(Res$protein_sequence)
        names(aass) <- Res$protein_id
        ## Add the mcols:
        retColumns <- retColumns[retColumns != "protein_sequence"]
        if (length(retColumns) > 0)
            mcols(aass) <- DataFrame(Res[, retColumns, drop = FALSE])
        return(aass)
    } else {
        Res <- Res[, retColumns, drop = FALSE]
        if (return.type == "DataFrame")
            Res <- DataFrame(Res)
        return(Res)
    }
    return(NULL)
})

############################################################
## listUniprotDbs
#' @aliases listUniprotDbs
#' 
#' @description The \code{listUniprotDbs} method lists all Uniprot database
#'     names in the \code{EnsDb}.
#' 
#' @examples
#'
#' ## List the names of all Uniprot databases from which Uniprot IDs are
#' ## available in the EnsDb
#' if (hasProteinData(edb))
#'     listUniprotDbs(edb)
#'
#' @rdname ProteinFunctionality
setMethod("listUniprotDbs", "EnsDb", function(object) {
    if (!hasProteinData(object))
        stop("The provided EnsDb database does not provide protein annotations!")
    res <- dbGetQuery(dbconn(object), "select distinct uniprot_db from uniprot")
    return(res$uniprot_db)
})

############################################################
## listUniprotMappingTypes
#' @aliases listUniprotMappingTypes
#' 
#' @description The \code{listUniprotMappingTypes} method lists all methods
#'     that were used for the mapping of Uniprot IDs to Ensembl protein IDs.
#'
#' @examples
#'
#' ## List the type of all methods that were used to map Uniprot IDs to Ensembl
#' ## protein IDs
#' if (hasProteinData(edb))
#'     listUniprotMappingTypes(edb)
#'
#' @rdname ProteinFunctionality
setMethod("listUniprotMappingTypes", "EnsDb", function(object) {
    if (!hasProteinData(object))
        stop("The provided EnsDb database does not provide protein annotations!")
    res <- dbGetQuery(dbconn(object),
                      "select distinct uniprot_mapping_type from uniprot")
    return(res$uniprot_mapping_type)
})

#' @description \code{supportedFilters} returns the names of all supported
#'     filters for the \code{EnsDb} object.
#'
#' @param object For \code{supportedFilters}: an \code{EnsDb} object.
#'
#' @param ... For \code{supportedFilters}: currently not used.
#'
#' @return For \code{supportedFilters}: the names of the supported filter
#'     classes.
#' 
#' @rdname Filter-classes
setMethod("supportedFilters", "EnsDb", function(object, ...) {
    .supportedFilters(object)
})
