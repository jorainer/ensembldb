####
## function to create a EnsDb object (or rather the SQLite database) from
## a Ensembl GTF file.
## Limitation:
## + There is no way to get the Entrezgene ID from this file.
## + Assuming that the element 2 in a row for a transcript represents its biotype, since
##   there is no explicit key transcript_biotype in element 9.
## + The CDS features in the GTF are somewhat problematic, while we're used to get just the
##   coding start and end for a transcript from the Ensembl perl API, here we get the coding
##   start and end for each exon.
ensDbFromGtf <- function(gtf, path=".", verbose=FALSE){
    options(useFancyQuotes=FALSE)
    if(verbose)
        cat("importing gtf file...")
    wanted.features <- c("gene", "transcript", "exon", "CDS")
    GTF <- import(con=gtf, format="gtf", feature.type=wanted.features, asRangedData=FALSE)
    if(verbose)
        cat("done\n")
    ## check what we've got...
    ## all wanted features?
    if(!any(wanted.features %in% levels(GTF$type)))
        stop(paste0("One or more required types are not in the gtf file. Need ",
                    paste(wanted.features, collapse=","), " but got only ",
                    paste(wanted.features[wanted.features %in% levels(GTF$type)], collapse=","),
                    "."))
    ## transcript biotype?
    if(any(colnames(mcols(GTF))=="transcript_biotype")){
        txBiotypeCol <- "transcript_biotype"
    }else{
        ## that's a little weird, but it seems that certain gtf files from Ensembl
        ## provide the transcript biotype in the element "source"
        txBiotypeCol <- "source"
    }
    ## processing the metadata:
    ## first read the header...
    tmp <- readLines(gtf, n=10)
    tmp <- tmp[grep(tmp, pattern="^#")]
    tmp <- gsub(tmp, pattern="^#", replacement="")
    tmp <- gsub(tmp, pattern="^!", replacement="")
    Header <- do.call(rbind, strsplit(tmp, split=" ", fixed=TRUE))
    colnames(Header) <- c("name", "value")
    ## getting the species name and the ensembl version from the GTF file name
    organism <- .organismName(organismFromGtfFileName(gtf))
    ensemblVersion <- ensemblVersionFromGtfFileName(gtf)
    genomeVersion <- genomeVersionFromGtfFileName(gtf)
    if(is.na(as.numeric(ensemblVersion)))
        stop("The GTF file name is not as expected: <Organism>.<genome version>.<Ensembl version>.gtf!")
    if(genomeVersion!=Header[Header[, "name"] == "genome-version", "value"]){
        stop(paste0("The GTF file name is not as expected: <Organism>.<genome version>.<Ensembl version>.gtf!",
                    " I've got genome version ", genomeVersion, " but in the header of the GTF file ",
                    Header[Header[, "name"] == "genome-version", "value"], " is specified!"))
    }
    MetaData <- data.frame(matrix(ncol=2, nrow=11))
    colnames(MetaData) <- c("name", "value")
    MetaData[1, ] <- c("Db type", "EnsDb")
    MetaData[2, ] <- c("Type of Gene ID", "Ensembl Gene ID")
    MetaData[3, ] <- c("Supporting package", "ensembldb")
    MetaData[4, ] <- c("Db created by", "ensembldb package from Bioconductor")
    MetaData[5, ] <- c("script_version", "0.0.1")
    MetaData[6, ] <- c("Creation time", date())
    MetaData[7, ] <- c("ensembl_version", ensemblVersion)
    tmp <- unlist(strsplit(gtf, split=.Platform$file.sep))
    MetaData[8, ] <- c("ensembl_host", tmp[length(tmp)])
    MetaData[9, ] <- c("Organism", organism )
    MetaData[10, ] <- c("genome_build", genomeVersion)
    MetaData[11, ] <- c("DBSCHEMAVERSION", "1.0")

    ## create database etc
    dbname <- paste0(path, .Platform$file.sep,
                     "EnsDb.",substring(organism, 1, 1),
                     unlist(strsplit(organism, split="_"))[ 2 ], ".v",
                     ensemblVersion, ".sqlite")
    con <- dbConnect(dbDriver("SQLite"), dbname=dbname)
    ## ----------------------------
    ## metadata table:
    if(verbose){
        cat("processing metadata...")
    }
    dbWriteTable(con, name="metadata", MetaData, overwrite=TRUE, row.names=FALSE)
    if(verbose)
        cat("OK\n")
    ## ----------------------------
    ##
    ## process genes
    ## we're lacking NCBI Entrezids and also the coord system, but these are not
    ## required columns anyway...
    if(verbose){
        cat("processing genes...")
    }
    ## want to have: gene_id, gene_name, entrezid, gene_biotype, gene_seq_start,
    ##               gene_seq_end, seq_name, seq_strand, seq_coord_system.
    reqCols <- c("gene_id", "gene_name", "gene_biotype")
    if(!any(reqCols %in% colnames(mcols(GTF))))
        stop(paste0("One or more required fields are not defined in the GTF! Need ",
                    paste(reqCols, collapse=","), " but got only ",
                    paste(reqCols[reqCols %in% colnames(mcols(GTF))], collapse=","),
                    "."))
    genes <- as.data.frame(GTF[GTF$type == "gene", reqCols])
    genes <- cbind(genes, entrezid=rep(NA, nrow(genes)),
                   seq_coord_system=rep(NA, nrow(genes)))
    colnames(genes)[1:5] <- c("seq_name", "gene_seq_start", "gene_seq_end", "width",
                              "seq_strand")
    ## transforming seq_strand from +/- to +1, -1.
    strand <- rep(0L, nrow(genes))
    strand[as.character(genes$seq_strand) == "+"] <- 1L
    strand[as.character(genes$seq_strand) == "-"] <- -1L
    genes[ , "seq_strand"] <- strand
    ## rearranging data.frame...
    genes <- genes[ , c("gene_id", "gene_name", "entrezid", "gene_biotype",
                        "gene_seq_start", "gene_seq_end", "seq_name",
                        "seq_strand", "seq_coord_system")]
    dbWriteTable(con, name="gene", genes, overwrite=TRUE, row.names=FALSE)
    if(verbose){
        cat("OK\n")
    }
    ## ----------------------------
    ##
    ## process chromosomes
    if(verbose)
        cat("processing chromosomes...")
    ## problem is I don't have these available...
    chroms <- data.frame(seq_name=unique(as.character(genes$seq_name)))
    chroms <- cbind(chroms, seq_length=rep(NA, nrow(chroms)),
                    is_circular=rep(NA, nrow(chroms)))
    ## write the table.
    dbWriteTable(con, name="chromosome", chroms, overwrite=TRUE, row.names=FALSE)
    rm(genes)
    if(verbose)
        cat("OK\n")
    ## ----------------------------
    ##
    ## process transcripts
    if(verbose)
        cat("processing transcripts...")
    ## want to have: tx_id, tx_biotype, tx_seq_start, tx_seq_end, tx_cds_seq_start,
    ##               tx_cds_seq_end, gene_id
    reqCols <- c("transcript_id", "gene_id", txBiotypeCol)
    if(!any(reqCols %in% colnames(mcols(GTF))))
        stop(paste0("One or more required fields are not defined in the GTF! Need ",
                    paste(reqCols, collapse=","), " but got only ",
                    paste(reqCols[reqCols %in% colnames(mcols(GTF))], collapse=","),
                    "."))
    tx <- as.data.frame(GTF[GTF$type == "transcript" , reqCols])[ , -c(1, 4, 5)]
    ## process the CDS features to get the cds start and end of the transcript.
    CDS <- as.data.frame(GTF[GTF$type == "CDS", "transcript_id"])
    cdsStarts <- aggregate(CDS[, "start"], by=list(CDS$transcript_id), FUN=min, na.rm=TRUE)
    cdsEnds <- aggregate(CDS[, "end"], by=list(CDS$transcript_id), FUN=max, na.rm=TRUE)
    idx <- match(cdsStarts[, 1], tx$transcript_id)
    tx <- cbind(tx, tx_cds_seq_start=rep(NA, nrow(tx)), tx_cds_seq_end=rep(NA, nrow(tx)))
    tx[idx, "tx_cds_seq_start"] <- cdsStarts[, 2]
    tx[idx, "tx_cds_seq_end"] <- cdsEnds[, 2]
    colnames(tx) <- c("tx_seq_start", "tx_seq_end", "tx_id", "gene_id", "tx_biotype",
                      "tx_cds_seq_start", "tx_cds_seq_end")
    ## rearranging data.frame:
    tx <- tx[ , c("tx_id", "tx_biotype", "tx_seq_start", "tx_seq_end",
                  "tx_cds_seq_start", "tx_cds_seq_end", "gene_id")]
    ## write the table.
    dbWriteTable(con, name="tx", tx, overwrite=TRUE, row.names=FALSE)
    rm(tx)
    rm(CDS)
    rm(cdsStarts)
    rm(cdsEnds)
    if(verbose){
        cat("OK\n")
    }
    ## ----------------------------
    ##
    ## process exons
    if(verbose)
        cat("processing exons...")
    reqCols <- c("exon_id", "transcript_id", "exon_number")
    if(!any(reqCols %in% colnames(mcols(GTF))))
        stop(paste0("One or more required fields are not defined in the GTF! Need ",
                    paste(reqCols, collapse=","), " but got only ",
                    paste(reqCols[reqCols %in% colnames(mcols(GTF))], collapse=","),
                    "."))
    exons <- as.data.frame(GTF[GTF$type == "exon", reqCols])[, -c(1, 4, 5)]
    ## for table tx2exon we want to have:
    ##    tx_id, exon_id, exon_idx
    t2e <- unique(exons[ , c("transcript_id", "exon_id", "exon_number")])
    colnames(t2e) <- c("tx_id", "exon_id", "exon_idx")
    ## for table exons we want to have:
    ##    exon_id, exon_seq_start, exon_seq_end
    exons <- unique(exons[ , c("exon_id", "start", "end")])
    colnames(exons) <- c("exon_id", "exon_seq_start", "exon_seq_end")
    ## writing the tables.
    dbWriteTable(con, name="exon", exons, overwrite=TRUE, row.names=FALSE)
    dbWriteTable(con, name="tx2exon", t2e, overwrite=TRUE, row.names=FALSE)
    if(verbose)
        cat("OK\n")
    if(verbose)
        cat("generating index...")
    ## generating all indices...
    dbGetQuery(con, "create index seq_name_idx on chromosome (seq_name);")
    dbGetQuery(con, "create index gene_id_idx on gene (gene_id);")
    dbGetQuery(con, "create index tx_id_idx on tx (tx_id);")
    dbGetQuery(con, "create index exon_id_idx on exon (exon_id);")
    dbGetQuery(con, "create index t2e_tx_id_idx on tx2exon (tx_id);")
    dbGetQuery(con, "create index t2e_exon_id_idx on tx2exon (exon_id);")
    if(verbose)
        cat("OK\n")
    dbDisconnect(con)
    return(dbname)
}


## compare the contents of the EnsDb sqlite database generated from a GTF (file name submitted
## with x ) with the one provided by package "lib".
compareEnsDbs <- function(x, y){
    ## compare two EnsDbs...
    if(organism(x)!=organism(y))
        stop("Well, at least the organism should be the same for both databases!")
    Messages <- rep("OK", 5)
    names(Messages) <- c("metadata", "chromosome", "gene", "transcript", "exon")
    ## comparing metadata.
    metadataX <- metadata(x)
    metadataY <- metadata(y)
    rownames(metadataX) <- metadataX[, 1]
    rownames(metadataY) <- metadataY[, 1]
    metadataY <- metadataY[rownames(metadataX),]
    cat("\nComparing metadata:\n")
    idx <- which(metadataX[, "value"]!=metadataY[, "value"])
    if(length(idx)>0)
        Messages["metadata"] <- "NOTE"
    ## check ensembl version
    if(metadataX["ensembl_version", "value"] == metadataY["ensembl_version", "value"]){
        cat(" Ensembl versions match.\n")
    }else{
        cat(" WARNING: databases base on different Ensembl versions! Expect considerable differences!\n")
        Messages["metadata"] <- "WARN"
    }
    ## genome build
    if(metadataX["genome_build", "value"] == metadataY["genome_build", "value"]){
        cat(" Genome builds match.\n")
    }else{
        cat(" WARNING: databases base on different Genome builds! Expect considerable differences!\n")
        Messages["metadata"] <- "WARN"
    }
    if(length(idx)>0){
        cat(" All differences: <name>: <value x> != <value y>\n")
        for(i in idx){
            cat(paste("  - ", metadataX[i, "name"], ":", metadataX[i, "value"], " != ",
                      metadataY[i, "value"], "\n"))
        }
    }
    cat(paste0("Done. Result: ", Messages["metadata"],"\n"))
    ## now comparing chromosomes
    Messages["chromosome"] <- compareChromosomes(x, y)
    ## comparing genes
    Messages["gene"] <- compareGenes(x, y)
    ## comparing transcripts
    Messages["transcript"] <- compareTx(x, y)
    ## comparing exons
    Messages["exon"] <- compareExons(x, y)
    return(Messages)
}


compareChromosomes <- function(x, y){
    Ret <- "OK"
    cat("\nComparing chromosome data:\n")
    chromX <- as.data.frame(seqinfo(x))
    chromY <- as.data.frame(seqinfo(y))
    ## compare seqnames
    inboth <- rownames(chromX)[rownames(chromX) %in% rownames(chromY)]
    onlyX <- rownames(chromX)[!(rownames(chromX) %in% rownames(chromY))]
    onlyY <- rownames(chromY)[!(rownames(chromY) %in% rownames(chromX))]
    if(length(onlyX) > 0 | length(onlyY) > 0)
        Ret <- "WARN"
    cat(paste0( " Sequence names: (", length(inboth), ") common, (",
               length(onlyX), ") only in x, (", length(onlyY), ") only in y.\n" ))
    same <- length(which(chromX[inboth, "seqlengths"]==chromY[inboth, "seqlengths"]))
    different <- length(inboth) - same
    cat(paste0( " Sequence lengths: (",same, ") identical, (", different, ") different.\n" ))
    if(different > 0)
        Ret <- "WARN"
    cat(paste0("Done. Result: ", Ret,"\n"))
    return(Ret)
}

compareGenes <- function(x, y){
    cat("\nComparing gene data:\n")
    Ret <- "OK"
    genesX <- genes(x)
    genesY <- genes(y)
    inboth <- names(genesX)[names(genesX) %in% names(genesY)]
    onlyX <- names(genesX)[!(names(genesX) %in% names(genesY))]
    onlyY <- names(genesY)[!(names(genesY) %in% names(genesX))]
    if(length(onlyX) > 0 | length(onlyY) > 0)
        Ret <- "WARN"
    cat(paste0(" gene IDs: (", length(inboth), ") common, (",
               length(onlyX), ") only in x, (", length(onlyY), ") only in y.\n"))
    ## seq names
    same <- length(
        which(as.character(seqnames(genesX[inboth]))==as.character(seqnames(genesY[inboth])))
        )
    different <- length(inboth) - same
    if(different > 0)
        Ret <- "ERROR"
    cat(paste0( " Sequence names: (",same, ") identical, (", different, ") different.\n" ))
    ## start
    same <- length(
        which(start(genesX[inboth]) == start(genesY[inboth]))
        )
    different <- length(inboth) - same
    if(different > 0)
        Ret <- "ERROR"
    cat(paste0( " Gene start coordinates: (",same,
               ") identical, (", different, ") different.\n" ))
    ## end
    same <- length(
        which(end(genesX[inboth]) == end(genesY[inboth]))
        )
    different <- length(inboth) - same
    if(different > 0)
        Ret <- "ERROR"
    cat(paste0( " Gene end coordinates: (",same,
               ") identical, (", different, ") different.\n" ))
    ## strand
    same <- length(
        which(as.character(strand(genesX[inboth]))
              == as.character(strand(genesY[inboth])))
        )
    different <- length(inboth) - same
    if(different > 0)
        Ret <- "ERROR"
    cat(paste0( " Gene strand: (",same,
               ") identical, (", different, ") different.\n" ))
    ## name
    same <- length(
        which(genesX[inboth]$gene_name == genesY[inboth]$gene_name)
        )
    different <- length(inboth) - same
    if(different > 0 & Ret!="ERROR")
        Ret <- "WARN"
    cat(paste0( " Gene names: (",same,
               ") identical, (", different, ") different.\n" ))
    ## entrezid
    same <- length(
        which(genesX[inboth]$entrezid == genesY[inboth]$entrezid)
        )
    different <- length(inboth) - same
    if(different > 0 & Ret!="ERROR")
        Ret <- "WARN"
    cat(paste0( " Entrezgene IDs: (",same,
               ") identical, (", different, ") different.\n" ))
    ## gene biotype
    same <- length(
        which(genesX[inboth]$gene_biotype == genesY[inboth]$gene_biotype)
        )
    different <- length(inboth) - same
    if(different > 0 & Ret!="ERROR")
        Ret <- "WARN"
    cat(paste0( " Gene biotypes: (",same,
               ") identical, (", different, ") different.\n" ))
    cat(paste0("Done. Result: ", Ret,"\n"))
    return(Ret)
}

compareTx <- function(x, y){
    cat("\nComparing transcript data:\n")
    Ret <- "OK"
    txX <- transcripts(x)
    txY <- transcripts(y)
    inboth <- names(txX)[names(txX) %in% names(txY)]
    onlyX <- names(txX)[!(names(txX) %in% names(txY))]
    onlyY <- names(txY)[!(names(txY) %in% names(txX))]
    if(length(onlyX) > 0 | length(onlyY) > 0)
        Ret <- "WARN"
    cat(paste0(" transcript IDs: (", length(inboth), ") common, (",
               length(onlyX), ") only in x, (", length(onlyY), ") only in y.\n"))
    ## start
    same <- length(
        which(start(txX[inboth]) == start(txY[inboth]))
        )
    different <- length(inboth) - same
    if(different > 0)
        Ret <- "ERROR"
    cat(paste0( " Transcript start coordinates: (",same,
               ") identical, (", different, ") different.\n" ))
    ## end
    same <- length(
        which(end(txX[inboth]) == end(txY[inboth]))
        )
    different <- length(inboth) - same
    if(different > 0)
        Ret <- "ERROR"
    cat(paste0( " Transcript end coordinates: (",same,
               ") identical, (", different, ") different.\n" ))
    ## tx biotype
    same <- length(
        which(txX[inboth]$tx_biotype == txY[inboth]$tx_biotype)
        )
    different <- length(inboth) - same
    if(different > 0 & Ret!="ERROR")
        Ret <- "WARN"
    cat(paste0( " Transcript biotypes: (",same,
               ") identical, (", different, ") different.\n" ))
    ## gene id
    same <- length(
        which(txX[inboth]$gene_id == txY[inboth]$gene_id)
        )
    different <- length(inboth) - same
    if(different > 0)
        Ret <- "ERROR"
    cat(paste0( " Associated gene IDs: (",same,
               ") identical, (", different, ") different.\n" ))
    cat(paste0("Done. Result: ", Ret,"\n"))
    return(Ret)
}

compareExons <- function(x, y){
    cat("\nComparing exon data:\n")
    Ret <- "OK"
    exonX <- exons(x)
    exonY <- exons(y)
    inboth <- names(exonX)[names(exonX) %in% names(exonY)]
    onlyX <- names(exonX)[!(names(exonX) %in% names(exonY))]
    onlyY <- names(exonY)[!(names(exonY) %in% names(exonX))]
    if(length(onlyX) > 0 | length(onlyY) > 0)
        Ret <- "WARN"
    cat(paste0(" exon IDs: (", length(inboth), ") common, (",
               length(onlyX), ") only in x, (", length(onlyY), ") only in y.\n"))
    ## start
    same <- length(
        which(start(exonX[inboth]) == start(exonY[inboth]))
        )
    different <- length(inboth) - same
    if(different > 0)
        Ret <- "ERROR"
    cat(paste0( " Exon start coordinates: (",same,
               ") identical, (", different, ") different.\n" ))
    ## end
    same <- length(
        which(end(exonX[inboth]) == end(exonY[inboth]))
        )
    different <- length(inboth) - same
    if(different > 0)
        Ret <- "ERROR"
    cat(paste0( " Exon end coordinates: (",same,
               ") identical, (", different, ") different.\n" ))
    ## now getting also the exon index in tx:
    exonX <- exons(x, columns=c("exon_id", "tx_id", "exon_idx"),
                   return.type="DataFrame")
    rownames(exonX) <- paste(exonX$tx_id, exonX$exon_id, sep=":")
    exonY <- exons(y, columns=c("exon_id", "tx_id", "exon_idx"),
                   return.type="DataFrame")
    rownames(exonY) <- paste(exonY$tx_id, exonY$exon_id, sep=":")
    inboth <- rownames(exonX)[rownames(exonX) %in% rownames(exonY)]
    onlyX <- rownames(exonX)[!(rownames(exonX) %in% rownames(exonY))]
    onlyY <- rownames(exonY)[!(rownames(exonY) %in% rownames(exonX))]

    ## tx exon idx
    same <- length(
        which(exonX[inboth, ]$exon_idx == exonY[inboth, ]$exon_idx)
        )
    different <- length(inboth) - same
    if(different > 0 )
        Ret <- "ERROR"
    cat(paste0( " Exon index in transcript models: (",same,
               ") identical, (", different, ") different.\n" ))
    cat(paste0("Done. Result: ", Ret,"\n"))
    return(Ret)
}



organismFromGtfFileName <- function(x){
    tmp <- unlist(strsplit(x, split=.Platform$file.sep, fixed=TRUE))
    splitty <- unlist(strsplit(tmp[length(tmp)], split=".", fixed=TRUE))
    return(splitty[1])
}

ensemblVersionFromGtfFileName <- function(x){
    tmp <- unlist(strsplit(x, split=.Platform$file.sep, fixed=TRUE))
    splitty <- unlist(strsplit(tmp[length(tmp)], split=".", fixed=TRUE))
    return(splitty[(grep(splitty, pattern="gtf")-1)])
}

genomeVersionFromGtfFileName <- function(x){
    tmp <- unlist(strsplit(x, split=.Platform$file.sep, fixed=TRUE))
    splitty <- unlist(strsplit(tmp[length(tmp)], split=".", fixed=TRUE))
    return(splitty[2])
}

