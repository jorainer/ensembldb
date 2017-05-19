############################################################
## Functions related to the creation of EnsDb databases.

## Separate helper function for abbreviating the genus and species name strings
## this simply makes the first character uppercase
.organismName <- function(x){
    substring(x, 1, 1) <- toupper(substring(x, 1, 1))
    return(x)
}

.abbrevOrganismName <- function(organism){
  spc <- unlist(strsplit(organism, "_"))
  ## this assumes a binomial nomenclature has been maintained.
  return(paste0(substr(spc[[1]], 1, 1), spc[[2]]))
}

## x has to be the connection to the database.
.makePackageName <- function(x){
    species <- .getMetaDataValue(x, "Organism")
    ensembl_version <- .getMetaDataValue(x, "ensembl_version")
    pkgName <- paste0("EnsDb.",.abbrevOrganismName(.organismName(species)),
                      ".v", ensembl_version)
    return(pkgName)
}

.makeObjectName <- function(pkgName){
  strs <- unlist(strsplit(pkgName, "\\."))
  paste(c(strs[2:length(strs)],strs[1]), collapse="_")
}


## retrieve Ensembl data
## save all files to local folder.
## returns the path where files have been saved to.
fetchTablesFromEnsembl <- function(version, ensemblapi, user="anonymous",
                                   host="ensembldb.ensembl.org", pass="",
                                   port=5306, species="human"){
    if(missing(version))
        stop("The version of the Ensembl database has to be provided!")
    ## setting the stage for perl:
    fn <- system.file("perl", "get_gene_transcript_exon_tables.pl",
                      package="ensembldb")
    ## parameters: s, U, H, P, e
    ## replacing white spaces with _
    species <- gsub(species, pattern=" ", replacement="_")

    cmd <- paste0("perl ", fn, " -s ", species," -e ", version,
                  " -U ", user, " -H ", host, " -p ", port, " -P ", pass)
    if(!missing(ensemblapi)){
        Sys.setenv(ENS=ensemblapi)
    }
    system(cmd)
    if(!missing(ensemblapi)){
        Sys.unsetenv("ENS")
    }

    ## we should now have the files:
    in_files <- c("ens_gene.txt", "ens_tx.txt", "ens_exon.txt",
                  "ens_tx2exon.txt", "ens_chromosome.txt", "ens_metadata.txt")
    ## check if we have all files...
    all_files <- dir(pattern="txt")
    if(sum(in_files %in% all_files)!=length(in_files))
        stop("Something went wrong! I'm missing some of the txt files the perl script should have generated.")
}


####
##
## create a SQLite database containing the information defined in the txt files.
makeEnsemblSQLiteFromTables <- function(path=".", dbname){
    ## check if we have all files...
    in_files <- c("ens_gene.txt", "ens_tx.txt", "ens_exon.txt",
                  "ens_tx2exon.txt", "ens_chromosome.txt", "ens_metadata.txt")
    ## check if we have all files...
    all_files <- dir(path, pattern="txt")
    if(sum(in_files %in% all_files)!=length(in_files))
        stop("Something went wrong! I'm missing some of the txt files the",
             " perl script should have generated.")

    ## read information
    info <- read.table(paste0(path, .Platform$file.sep ,"ens_metadata.txt"),
                       sep="\t", as.is=TRUE, header=TRUE)
    species <- .organismName(info[ info$name=="Organism", "value" ])
    ##substring(species, 1, 1) <- toupper(substring(species, 1, 1))
    if(missing(dbname)){
        dbname <- paste0("EnsDb.",substring(species, 1, 1),
                         unlist(strsplit(species, split="_"))[ 2 ], ".v",
                         info[ info$name=="ensembl_version", "value" ], ".sqlite")
    }
    con <- dbConnect(dbDriver("SQLite"), dbname=dbname)

    ## write information table
    dbWriteTable(con, name="metadata", info, row.names=FALSE)

    ## process chromosome
    message("Processing 'chromosome' table ... ", appendLF = FALSE)
    tmp <- read.table(paste0(path, .Platform$file.sep ,"ens_chromosome.txt"),
                      sep="\t", as.is=TRUE, header=TRUE)
    tmp[, "seq_name"] <- as.character(tmp[, "seq_name"])
    dbWriteTable(con, name="chromosome", tmp, row.names=FALSE)
    rm(tmp)
    message("OK")

    message("Processing 'gene' table ... ", appendLF = FALSE)
    ## process genes: some gene names might have fancy names...
    tmp <- read.table(paste0(path, .Platform$file.sep, "ens_gene.txt"),
                      sep="\t", as.is=TRUE, header=TRUE,
                      quote="", comment.char="" )
    OK <- .checkIntegerCols(tmp)
    dbWriteTable(con, name="gene", tmp, row.names=FALSE)
    rm(tmp)
    message("OK")

    if (as.numeric(info[info$name == "DBSCHEMAVERSION", "value"]) > 1) {
        message("Processing 'entrezgene' table ... ", appendLF = FALSE)
        ## process genes: some gene names might have fancy names...
        tmp <- read.table(paste0(path, .Platform$file.sep, "ens_entrezgene.txt"),
                          sep="\t", as.is=TRUE, header=TRUE,
                          quote="", comment.char="" )
        dbWriteTable(con, name="entrezgene", tmp, row.names=FALSE)
        rm(tmp)
        message("OK")
    }
    
    message("Processing 'trancript' table ... ", appendLF = FALSE)
    ## process transcripts:
    tmp <- read.table(paste0(path, .Platform$file.sep, "ens_tx.txt"),
                      sep="\t", as.is=TRUE, header=TRUE)
    ## Fix the tx_cds_seq_start and tx_cds_seq_end columns: these should be integer!
    suppressWarnings(
        tmp[, "tx_cds_seq_start"] <- as.integer(tmp[, "tx_cds_seq_start"])
    )
    suppressWarnings(
        tmp[, "tx_cds_seq_end"] <- as.integer(tmp[, "tx_cds_seq_end"])
    )
    OK <- .checkIntegerCols(tmp)
    dbWriteTable(con, name="tx", tmp, row.names=FALSE)
    rm(tmp)
    message("OK")

    ## process exons:
    message("Processing 'exon' table ... ", appendLF = FALSE)
    tmp <- read.table(paste0(path, .Platform$file.sep, "ens_exon.txt"),
                      sep = "\t", as.is = TRUE, header = TRUE)
    OK <- .checkIntegerCols(tmp)
    dbWriteTable(con, name="exon", tmp, row.names=FALSE)
    rm(tmp)
    message("OK")
    message("Processing 'tx2exon' table ... ", appendLF = FALSE)
    tmp <- read.table(paste0(path, .Platform$file.sep, "ens_tx2exon.txt"),
                      sep = "\t", as.is = TRUE, header = TRUE)
    OK <- .checkIntegerCols(tmp)
    dbWriteTable(con, name="tx2exon", tmp, row.names = FALSE)
    rm(tmp)
    message("OK")

    ## process proteins; if available.
    prot_file <- paste0(path, .Platform$file.sep, "ens_protein.txt")
    if (file.exists(prot_file)) {
        message("Processing 'protein' table ... ", appendLF = FALSE)
        tmp <- read.table(prot_file, sep = "\t", as.is = TRUE, header = TRUE)
        OK <- .checkIntegerCols(tmp)
        dbWriteTable(con, name = "protein", tmp, row.names = FALSE)
        message("OK")
        message("Processing 'uniprot' table ... ", appendLF = FALSE)
        tmp <- read.table(paste0(path, .Platform$file.sep, "ens_uniprot.txt"),
                          sep = "\t", as.is = TRUE, header = TRUE)
        OK <- .checkIntegerCols(tmp)
        dbWriteTable(con, name = "uniprot", tmp, row.names = FALSE)
        message("OK")
        message("Processing 'protein_domain' table ... ", appendLF = FALSE)
        tmp <- read.table(paste0(path, .Platform$file.sep, "ens_protein_domain.txt"),
                          sep = "\t", as.is = TRUE, header = TRUE)
        OK <- .checkIntegerCols(tmp)
        dbWriteTable(con, name = "protein_domain", tmp, row.names = FALSE)
        message("OK")
    }

    ## Create indices
    message("Creating indices ... ", appendLF = FALSE)
    .createEnsDbIndices(con, proteins = file.exists(prot_file))
    message("OK")
    dbDisconnect(con)
    ## Check if the data could be loaded.
    message("Checking validity of the database ... ", appendLF = FALSE)
    msg <- validObject(EnsDb(dbname))
    if (!is.logical(msg))
        stop(msg)
    message("OK")
    ## done.
    return(dbname)
}

############################################################
## Simply checking that some columns are integer
.checkIntegerCols <- function(x, columns = c("gene_seq_start", "gene_seq_end",
                                             "tx_seq_start", "tx_seq_start",
                                             "exon_seq_start", "exon_seq_end",
                                             "exon_idx", "tx_cds_seq_start",
                                             "tx_cds_seq_end", "prot_dom_start",
                                             "prot_dom_end")) {
    cols <- columns[columns %in% colnames(x)]
    if(length(cols) > 0) {
        sapply(cols, function(z) {
            if(!is.integer(x[, z]))
                stop("Column '", z,"' is not of type integer!")
        })
    }
    return(TRUE)
}


####
## the function that creates the annotation package.
## ensdb should be a connection to an SQLite database, or a character string...
makeEnsembldbPackage <- function(ensdb,
                                 version,
                                 maintainer,
                                 author,
                                 destDir=".",
                                 license="Artistic-2.0"){
    if(class(ensdb)!="character")
        stop("ensdb has to be the name of the SQLite database!")
    ensdbfile <- ensdb
    ensdb <- EnsDb(x=ensdbfile)
    con <- dbconn(ensdb)
    pkgName <- .makePackageName(con)
    ensembl_version <- .getMetaDataValue(con, "ensembl_version")
    ## there should only be one template
    template_path <- system.file("pkg-template",package="ensembldb")
    ## We need to define some symbols in order to have the
    ## template filled out correctly.
    symvals <- list(
        PKGTITLE=paste("Ensembl based annotation package"),
        PKGDESCRIPTION="Exposes an annotation databases generated from Ensembl.",
        PKGVERSION=version,
        AUTHOR=author,
        MAINTAINER=maintainer,
        LIC=license,
        ORGANISM=.organismName(.getMetaDataValue(con ,'Organism')),
        SPECIES=.organismName(.getMetaDataValue(con,'Organism')),
        PROVIDER="Ensembl",
        PROVIDERVERSION=as.character(ensembl_version),
        RELEASEDATE= .getMetaDataValue(con ,'Creation time'),
        SOURCEURL= .getMetaDataValue(con ,'ensembl_host'),
        ORGANISMBIOCVIEW=gsub(" ","_",
                              .organismName(.getMetaDataValue(con ,'Organism'))),
        TXDBOBJNAME=pkgName ## .makeObjectName(pkgName)
       )
    ## Should never happen
    if (any(duplicated(names(symvals)))) {
        str(symvals)
        stop("'symvals' contains duplicated symbols")
    }
    createPackage(pkgname=pkgName,
                  destinationDir=destDir,
                  originDir=template_path,
                  symbolValues=symvals)
    ## then copy the contents of the database into the extdata dir
    sqlfilename <- unlist(strsplit(ensdbfile, split=.Platform$file.sep))
    sqlfilename <- sqlfilename[ length(sqlfilename) ]
    dir.create(paste(c(destDir, pkgName, "inst", "extdata"),
                     collapse=.Platform$file.sep),
               showWarnings=FALSE, recursive=TRUE)
    db_path <- file.path(destDir, pkgName, "inst", "extdata",
                         paste(pkgName,"sqlite",sep="."))
    file.copy(ensdbfile, to=db_path)
}


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
ensDbFromGtf <- function(gtf, outfile, path, organism, genomeVersion,
                         version, ...){
    options(useFancyQuotes=FALSE)
    message("Importing GTF file ... ", appendLF=FALSE)
    ## wanted.features <- c("gene", "transcript", "exon", "CDS")
    wanted.features <- c("exon")
    ## GTF <- import(con=gtf, format="gtf", feature.type=wanted.features)
    GTF <- import(con=gtf, format="gtf")
    message("OK")
    ## check what we've got...
    ## all wanted features?
    if(any(!(wanted.features %in% levels(GTF$type)))){
        stop(paste0("One or more required types are not in the gtf file. Need ",
                    paste(wanted.features, collapse=","), " but got only ",
                    paste(wanted.features[wanted.features %in% levels(GTF$type)],
                          collapse=","),
                    "."))
    }
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
    haveHeader <- FALSE
    if(length(tmp) > 0){
        ##message("GTF file has a header.")
        tmp <- gsub(tmp, pattern="^#", replacement="")
        tmp <- gsub(tmp, pattern="^!", replacement="")
        Header <- do.call(rbind, strsplit(tmp, split=" ", fixed=TRUE))
        colnames(Header) <- c("name", "value")
        haveHeader <- TRUE
    }
    ## Check parameters
    Parms <- .checkExtractVersions(gtf, organism, genomeVersion, version)
    ensemblVersion <- Parms["version"]
    organism <- Parms["organism"]
    genomeVersion <- Parms["genomeVersion"]

    if(haveHeader){
        if(genomeVersion!=Header[Header[, "name"] == "genome-version", "value"]){
            stop(paste0("The GTF file name is not as expected: <Organism>.",
                        "<genome version>.<Ensembl version>.gtf!",
                        " I've got genome version ", genomeVersion,
                        " but in the header of the GTF file ",
                        Header[Header[, "name"] == "genome-version", "value"],
                        " is specified!"))
        }
    }

    GTF <- fixCDStypeInEnsemblGTF(GTF)
    ## here on -> call ensDbFromGRanges.
    dbname <- ensDbFromGRanges(GTF, outfile = outfile, path = path,
                               organism = organism,
                               genomeVersion = genomeVersion,
                               version = ensemblVersion, ...)

    gtfFilename <- unlist(strsplit(gtf, split=.Platform$file.sep))
    gtfFilename <- gtfFilename[length(gtfFilename)]
    ## updating the Metadata information...
    lite <- dbDriver("SQLite")
    con <- dbConnect(lite, dbname = dbname )
    bla <- dbGetQuery(con, paste0("update metadata set value='",
                                  gtfFilename,
                                  "' where name='source_file';"))
    dbDisconnect(con)
    return(dbname)
}

####============================================================
##  fixCDStypeInEnsemblGTF
##
##  Takes an GRanges object as input and returns a GRanges object in
##  which the feature type stop_codon and start_codon is replaced by
##  feature type CDS. This is to fix a potential problem (bug?) in
##  GTF files from Ensembl, in which the stop_codon or start_codon for
##  some transcripts is outside of the CDS.
####------------------------------------------------------------
fixCDStypeInEnsemblGTF <- function(x){
    if(any(unique(x$type) %in% c("start_codon", "stop_codon"))){
        x$type[x$type %in% c("start_codon", "stop_codon")] <- "CDS"
    }
    return(x)
}

####============================================================
##  ensDbFromAH
##
##  Retrieve a GTF file from AnnotationHub and build a EnsDb object from that.
##
####------------------------------------------------------------
ensDbFromAH <- function(ah, outfile, path, organism, genomeVersion, version){
    options(useFancyQuotes=FALSE)
    ## Input checking...
    if(!is(ah, "AnnotationHub"))
        stop("Argument 'ah' has to be a (single) AnnotationHub object.")
    if(length(ah) != 1)
        stop("Argument 'ah' has to be a single AnnotationHub resource!")
    if(tolower(ah$dataprovider) != "ensembl")
        stop("Can only process GTF files provided by Ensembl!")
    if(tolower(ah$sourcetype) != "gtf")
        stop("Resource is not a GTF file!")
    ## Check parameters
    Parms <- .checkExtractVersions(ah$title, organism, genomeVersion, version)
    ensFromAH <- Parms["version"]
    orgFromAH <- Parms["organism"]
    genFromAH <- Parms["genomeVersion"]
    gtfFilename <- ah$title
    message("Fetching data ... ", appendLF=FALSE)
    suppressMessages(
        gff <- ah[[1]]
    )
    message("OK")
    message("  -------------")
    message("Proceeding to create the database.")

    gff <- fixCDStypeInEnsemblGTF(gff)
    ## Proceed.
    dbname <- ensDbFromGRanges(gff, outfile=outfile, path=path, organism=orgFromAH,
                               genomeVersion=genFromAH, version=ensFromAH)
    ## updating the Metadata information...
    lite <- dbDriver("SQLite")
    con <- dbConnect(lite, dbname = dbname )
    bla <- dbGetQuery(con, paste0("update metadata set value='",
                                  gtfFilename,
                                  "' where name='source_file';"))
    dbDisconnect(con)
    return(dbname)
}

.checkExtractVersions <- function(filename, organism, genomeVersion, version){
    if(isEnsemblFileName(filename)){
        ensFromFile <- ensemblVersionFromGtfFileName(filename)
        orgFromFile <- organismFromGtfFileName(filename)
        genFromFile <- genomeVersionFromGtfFileName(filename)
    }else{
        ensFromFile <- NA
        orgFromFile <- NA
        genFromFile <- NA
        if(missing(organism) | missing(genomeVersion) | missing(version))
            stop("The file name does not match the expected naming scheme",
                 " of Ensembl files hence I cannot extract any information",
                 " from it! Parameters 'organism', 'genomeVersion' and",
                 " 'version' are thus required!")
    }
    ## Do some more testing with versions provided from the user.
    if(!missing(organism)){
        if(!is.na(orgFromFile)){
            if(organism != orgFromFile){
                warning("User specified organism (", organism,
                        ") is different to the one extracted",
                        " from the file name (", orgFromFile,
                        ")! Using the one defined by the user.")
            }
        }
        orgFromFile <- organism
    }
    if(!missing(genomeVersion)){
        if(!is.na(genFromFile)){
            if(genomeVersion != genFromFile){
                warning("User specified genome version (", genomeVersion,
                        ") is different to the one extracted",
                        " from the file name (", genFromFile,
                        ")! Using the one defined by the user.")
            }
        }
        genFromFile <- genomeVersion
    }
    if(!missing(version)){
        if(!is.na(ensFromFile)){
            if(version != ensFromFile){
                warning("User specified Ensembl version (", version,
                        ") is different to the one extracted",
                        " from the file name (", ensFromFile,
                        ")! Using the one defined by the user.")
            }
        }
        ensFromFile <- version
    }
    res <- c(orgFromFile, genFromFile, ensFromFile)
    names(res) <- c("organism", "genomeVersion", "version")
    return(res)
}



####============================================================
##
##  ensDbFromGff
##
####------------------------------------------------------------
ensDbFromGff <- function(gff, outfile, path, organism, genomeVersion,
                         version, ...){
    options(useFancyQuotes=FALSE)

    ## Check parameters
    Parms <- .checkExtractVersions(gff, organism, genomeVersion, version)
    ensFromFile <- Parms["version"]
    orgFromFile <- Parms["organism"]
    genFromFile <- Parms["genomeVersion"]
    ## Reading some info from the header.
    tmp <- readLines(gff, n=500)
    if(length(grep(tmp[1], pattern="##gff-version")) == 0)
        stop("File ", gff, " does not seem to be a correct GFF file! ",
             "The ##gff-version line is missing!")
    gffVersion <- unlist(strsplit(tmp[1], split="[ ]+"))[2]
    if(gffVersion != "3")
        stop("This function supports only GFF version 3 files!")
    tmp <- tmp[grep(tmp, pattern="^#!")]
    if(length(tmp) > 0){
        ## Check if I can extract the genome-version
        idx <- grep(tmp, pattern = "^#!genome-version")
        if (length(idx) > 0) {
            genFromHeader <- sub(tmp[idx], pattern = "^#!genome-version",
                                 replacement = "")
            genFromHeader <- gsub(genFromHeader, pattern = " ",
                                  replacement = "", fixed = TRUE)
            if (genFromHeader != genFromFile) {
                warning("Genome version extracted from file name (",
                        genFromFile, ") does not match genome version",
                        " defined within the gff file (", genFromHeader,
                        "). Will use the version defined within the gff.")
            }
        }
    }

    message("Importing GFF ... ", appendLF=FALSE)
    suppressWarnings(
        theGff <- import(gff, format=paste0("gff", gffVersion))
    )
    message("OK")
    ## Works with Ensembl 83; eventually not for updated Ensembl gff files!

    ## what seems a little strange: exons have an ID of NA.
    ## Ensembl specific fields: gene_id, transcript_id, exon_id, rank, biotype.
    ## GFF3 fields: type, ID, Name, Parent
    ## check columns and subset...
    gffcols <- c("type", "ID", "Name", "Parent")
    if(!all(gffcols %in% colnames(mcols(theGff))))
        stop("Required columns/fields ",
             paste(gffcols[!(gffcols %in% colnames(mcols(theGff)))],
                   collapse=";"),
             " not present in the GFF file!")
    enscols <- c("gene_id", "transcript_id", "exon_id", "rank", "biotype")
    if(!all(enscols %in% colnames(mcols(theGff))))
        stop("Required columns/fields ",
             paste(enscols[!(enscols %in% colnames(mcols(theGff)))],
                   collapse=";"),
             " not present in the GFF file!")
    ## Subsetting to eventually speed up further processing.
    theGff <- theGff[, c(gffcols, enscols)]
    ## Renaming and fixing some columns:
    CN <- colnames(mcols(theGff))
    colnames(mcols(theGff))[CN == "Name"] <- "gene_name"
    colnames(mcols(theGff))[CN == "biotype"] <- "gene_biotype"
    colnames(mcols(theGff))[CN == "rank"] <- "exon_number"
    theGff$transcript_biotype <- theGff$gene_biotype

    ## Processing that stuff...
    ## Replace the ID format type:ID.
    ids <- strsplit(theGff$ID, split=":")
    message("Fixing IDs ... ", appendLF=FALSE)
    ## For those that have length > 1 use the second element.
    theGff$ID <- unlist(lapply(ids, function(z){
        if(length(z) > 1)
            return(z[2])
        return(z)
    }))
    message("OK")
    ## Process genes...
    message("Processing genes ... ", appendLF=FALSE)
    ## Bring the GFF into the correct format for EnsDb/ensDbFromGRanges.
    idx <- which(!is.na(theGff$gene_id))
    theGff$type[idx] <- "gene"
    message("OK")

    ## ## Can not use the lengths of chromosomes provided in the chromosome features!!!
    ## ## For whatever reasons chromosome Y length is incorrect!!!
    ## message("Processing seqinfo...", appendLF=FALSE)
    ## SI <- seqinfo(theGff)
    ## tmp <- theGff[theGff$ID %in% seqlevels(SI)]
    ## ## Check if we've got length for all.
    ## message("OK")

    ## Process transcripts...
    message("Processing transcripts ... ", appendLF=FALSE)
    idx <- which(!is.na(theGff$transcript_id))
    ## Check if I've got multiple parents...
    parentGenes <- theGff$Parent[idx]
    if(any(lengths(parentGenes) > 1))
        stop("Transcripts with multiple parents in GFF element 'Parent'",
             " not (yet) supported!")
    theGff$type[idx] <- "transcript"
    ## Setting the gene_id for these guys...
    theGff$gene_id[idx] <- unlist(sub(parentGenes, pattern="gene:",
                                      replacement="", fixed=TRUE))
    ## The CDS:
    idx <- which(theGff$type == "CDS")
    parentTx <- theGff$Parent[idx]
    if(any(lengths(parentTx) > 1))
        stop("CDS with multiple parent transcripts in GFF element 'Parent'",
             " not (yet) supported!")
    theGff$transcript_id[idx] <- unlist(sub(parentTx, pattern="transcript:",
                                            replacement="", fixed=TRUE))
    message("OK")

    message("Processing exons ... ", appendLF=FALSE)
    idx <- which(!is.na(theGff$exon_id))
    parentTx <- theGff$Parent[idx]
    if(any(lengths(parentTx) > 1))
        stop("Exons with multiple parent transcripts in GFF element 'Parent'",
             " not (yet) supported!")
    theGff$transcript_id[idx] <- unlist(sub(parentTx, pattern="transcript:",
                                            replacement="", fixed=TRUE))
    message("OK")

    theGff <- theGff[theGff$type %in% c("gene", "transcript", "exon", "CDS")]
    theGff <- keepSeqlevels(theGff, as.character(unique(seqnames(theGff))))
    ## Now we can proceed and pass that to the next function!

    message("  -------------")
    message("Proceeding to create the database.")

    ## Proceed.
    dbname <- ensDbFromGRanges(theGff, outfile = outfile, path = path,
                               organism = orgFromFile,
                               genomeVersion = genFromFile,
                               version = ensFromFile, ...)

    gtfFilename <- unlist(strsplit(gff, split=.Platform$file.sep))
    gtfFilename <- gtfFilename[length(gtfFilename)]
    ## updating the Metadata information...
    lite <- dbDriver("SQLite")
    con <- dbConnect(lite, dbname = dbname )
    bla <- dbGetQuery(con, paste0("update metadata set value='",
                                  gtfFilename,
                                  "' where name='source_file';"))
    dbDisconnect(con)
    return(dbname)
}



#### build a EnsDb SQLite database from the GRanges.
## we can however not get all of the information from the GRanges (yet), for example,
## the seqinfo might not be available in all GRanges objects. Also, there is no way
## we can guess the organism or the Ensembl version from the GRanges, thus, this
## information has to be provided by the user.
## x: the GRanges object or file name. If file name, the function tries to guess
##    the organism, genome build and ensembl version from the file name, if not
##    provided.
##
ensDbFromGRanges <- function(x, outfile, path, organism, genomeVersion,
                             version, ...){
    if(!is(x, "GRanges"))
        stop("This method can only be called on GRanges objects!")
    ## check for missing parameters
    if(missing(organism)){
        stop("The organism has to be specified (e.g. using",
             " organism=\"Homo_sapiens\")")
    }
    if(missing(version)){
        stop("The Ensembl version has to be specified!")
    }

    ## checking the seqinfo in the GRanges object...
    Seqinfo <- seqinfo(x)
    fetchSeqinfo <- FALSE
    ## check if we've got some information...
    if(any(is.na(seqlengths(Seqinfo)))){
        fetchSeqinfo <- TRUE   ## means we have to fetch the seqinfo ourselfs...
    }
    if(missing(genomeVersion)){
        ## is there a seqinfo in x that I could use???
        if(!fetchSeqinfo){
            genomeVersion <- unique(genome(Seqinfo))
            if(is.na(genomeVersion) | length(genomeVersion) > 1){
                stop(paste0("The genome version has to be specified as",
                            " it can not be extracted from the seqinfo!"))
            }
        }else{
            stop("The genome version has to be specified!")
        }
    }
    if(missing(outfile)){
        ## use the organism, genome version and ensembl version as the file name.
        outfile <- paste0(c(organism, genomeVersion, version, "sqlite"),
                          collapse=".")
        if(missing(path))
            path <- "."
        dbname <- paste0(path, .Platform$file.sep, outfile)
    }else{
        if(!missing(path))
            warning("outfile specified, thus I will discard the path argument.")
        dbname <- outfile
    }

    ## that's quite some hack
    ## transcript biotype?
    if(any(colnames(mcols(x))=="transcript_biotype")){
        txBiotypeCol <- "transcript_biotype"
    }else{
        ## that's a little weird, but it seems that certain gtf files from Ensembl
        ## provide the transcript biotype in the element "source"
        txBiotypeCol <- "source"
    }

    con <- dbConnect(dbDriver("SQLite"), dbname=dbname)
    on.exit(dbDisconnect(con))
    ## ----------------------------
    ## metadata table:
    message("Processing metadata ... ", appendLF=FALSE)
    Metadata <- buildMetadata(organism, version, host="unknown",
                              sourceFile="GRanges object",
                              genomeVersion=genomeVersion)
    dbWriteTable(con, name="metadata", Metadata, overwrite=TRUE,
                 row.names=FALSE)
    message("OK")
    ## Check if we've got column "type"
    if(!any(colnames(mcols(x)) == "type"))
        stop("The GRanges object lacks the required column 'type', sorry.")
    gotTypes <- as.character(unique(x$type))
    gotColumns <- colnames(mcols(x))
    ## ----------------------------
    ##
    ## process genes
    ## we're lacking NCBI Entrezids and also the coord system, but these are not
    ## required columns anyway...
    message("Processing genes ... ")
    ## want to have: gene_id, gene_name, entrezid, gene_biotype, gene_seq_start,
    ##               gene_seq_end, seq_name, seq_strand, seq_coord_system.
    wouldBeNice <- c("gene_id", "gene_name", "entrezid", "gene_biotype")
    dontHave <- wouldBeNice[!(wouldBeNice %in% gotColumns)]
    haveGot <- wouldBeNice[wouldBeNice %in% gotColumns]
    ## Just really require the gene_id...
    reqCols <- c("gene_id")
    if(length(dontHave) > 0){
        mess <- paste0(" I'm missing column(s): ", paste0(sQuote(dontHave),
                                                          collapse=","),
                       ".")
        warning(mess, " The corresponding database column(s) will be empty!")
    }
    message(" Attribute availability:", appendLF=TRUE)
    for(i in 1:length(wouldBeNice)){
        message("  o ", wouldBeNice[i], " ... ",
                ifelse(any(gotColumns == wouldBeNice[i]), yes="OK", no="Nope"))
    }
    if(!any(reqCols %in% haveGot))
        stop(paste0("One or more required fields are not defined in the",
                    " submitted GRanges object! Need ",
                    paste(sQuote(reqCols), collapse=","), " but got only ",
                    paste(reqCols[reqCols %in% gotColumns], collapse=","),
                    "."))
    ## Now gets tricky; special case Ensembl < 75: we've got NO gene type.
    if(any(gotTypes == "gene")){
        ## All is fine.
        genes <- as.data.frame(x[x$type == "gene", haveGot])
    }else{
        ## Well, have to split by gene_id and process...
        genes <- split(x[ , haveGot], x$gene_id)
        gnRanges <- unlist(range(genes))
        gnMcol <- as.data.frame(unique(mcols(unlist(genes))))
        genes <- as.data.frame(gnRanges)
        ## Adding mcols again.
        genes <- cbind(genes, gnMcol[match(rownames(genes), gnMcol$gene_id), ])
        rm(gnRanges)
        rm(gnMcol)
    }
    colnames(genes) <- c("seq_name", "gene_seq_start", "gene_seq_end", "width",
                         "seq_strand", haveGot)
    ## Add missing cols...
    if(length(dontHave) > 0){
        cn <- colnames(genes)
        for(i in 1:length(dontHave)){
            genes <- cbind(genes, rep(NA, nrow(genes)))
        }
        colnames(genes) <- c(cn, dontHave)
    }
    genes <- cbind(genes, seq_coord_system=rep(NA, nrow(genes)))

    ## transforming seq_strand from +/- to +1, -1.
    strand <- rep(0L, nrow(genes))
    strand[as.character(genes$seq_strand) == "+"] <- 1L
    strand[as.character(genes$seq_strand) == "-"] <- -1L
    genes[ , "seq_strand"] <- strand
    ## rearranging data.frame...
    genes <- genes[ , c("gene_id", "gene_name", "entrezid", "gene_biotype",
                        "gene_seq_start", "gene_seq_end", "seq_name",
                        "seq_strand", "seq_coord_system")]
    OK <- .checkIntegerCols(genes)
    dbWriteTable(con, name="gene", genes, overwrite=TRUE, row.names=FALSE)
    ## Done.

    message("OK")
    ## ----------------------------
    ##
    ## process transcripts
    message("Processing transcripts ... ", appendLF=TRUE)
    ## want to have: tx_id, tx_biotype, tx_seq_start, tx_seq_end, tx_cds_seq_start,
    ##               tx_cds_seq_end, gene_id
    wouldBeNice <- c("transcript_id", "gene_id", txBiotypeCol)
    dontHave <- wouldBeNice[!(wouldBeNice %in% gotColumns)]
    if(length(dontHave) > 0){
        mess <- paste0("I'm missing column(s): ", paste0(sQuote(dontHave),
                                                         collapse=","), ".")
        warning(mess, " The corresponding database columns will be empty!")
    }
    haveGot <- wouldBeNice[wouldBeNice %in% gotColumns]
    message(" Attribute availability:", appendLF=TRUE)
    for(i in 1:length(wouldBeNice)){
        message("  o ", wouldBeNice[i], " ... ",
                ifelse(any(gotColumns == wouldBeNice[i]), yes="OK", no="Nope"))
    }
    reqCols <- c("transcript_id", "gene_id")
    if(!any(reqCols %in% gotColumns))
        stop(paste0("One or more required fields are not defined in",
                    " the submitted GRanges object! Need ",
                    paste(reqCols, collapse=","), " but got only ",
                    paste(reqCols[reqCols %in% gotColumns], collapse=","),
                    "."))
    if(any(gotTypes == "transcript")){
        tx <- as.data.frame(x[x$type == "transcript" , haveGot])
    }else{
        tx <- split(x[, haveGot], x$transcript_id)
        txRanges <- unlist(range(tx))
        txMcol <- as.data.frame(unique(mcols(unlist(tx))))
        tx <- as.data.frame(txRanges)
        tx <- cbind(tx, txMcol[match(rownames(tx), txMcol$transcript_id), ])
        rm(txRanges)
        rm(txMcol)
    }
    ## Drop columns seqnames, width and strand
    tx <- tx[, -c(1, 4, 5)]
    ## Add empty columns, eventually
    if(length(dontHave) > 0){
        cn <- colnames(tx)
        for(i in 1:length(dontHave)){
            tx <- cbind(tx, rep(NA, nrow(tx)))
        }
        colnames(tx) <- c(cn, dontHave)
    }
    ## Add columns for UTR
    tx <- cbind(tx, tx_cds_seq_start=rep(NA, nrow(tx)),
                tx_cds_seq_end=rep(NA, nrow(tx)))
    ## Process CDS...
    if(any(gotTypes == "CDS")){
        ## Only do that if we've got type == "CDS"!
        ## process the CDS features to get the cds start and end of the transcript.
        CDS <- as.data.frame(x[x$type == "CDS", "transcript_id"])
        ##
        startByTx <- split(CDS$start, f=CDS$transcript_id)
        cdsStarts <- unlist(lapply(startByTx,
                                   function(z){return(min(z, na.rm=TRUE))}))
        endByTx <- split(CDS$end, f=CDS$transcript_id)
        cdsEnds <- unlist(lapply(endByTx,
                                 function(z){return(max(z, na.rm=TRUE))}))
        idx <- match(names(cdsStarts), tx$transcript_id)
        areNas <- is.na(idx)
        idx <- idx[!areNas]
        cdsStarts <- cdsStarts[!areNas]
        cdsEnds <- cdsEnds[!areNas]
        tx[idx, "tx_cds_seq_start"] <- cdsStarts
        tx[idx, "tx_cds_seq_end"] <- cdsEnds
    }else{
        mess <- paste0(" I can't find type=='CDS'! The resulting database",
                       " will lack CDS information!")
        message(mess, appendLF = TRUE)
        warning(mess)
    }
    colnames(tx) <- c("tx_seq_start", "tx_seq_end", "tx_id", "gene_id",
                      "tx_biotype", "tx_cds_seq_start", "tx_cds_seq_end")
    ## rearranging data.frame:
    tx <- tx[ , c("tx_id", "tx_biotype", "tx_seq_start", "tx_seq_end",
                  "tx_cds_seq_start", "tx_cds_seq_end", "gene_id")]
    ## write the table.
    OK <- .checkIntegerCols(tx)
    dbWriteTable(con, name="tx", tx, overwrite=TRUE, row.names=FALSE)
    rm(tx)
    rm(CDS)
    rm(cdsStarts)
    rm(cdsEnds)
    message("OK")
    ## ----------------------------
    ##
    ## process exons
    message("Processing exons ... ", appendLF=FALSE)
    reqCols <- c("exon_id", "transcript_id", "exon_number")
    if(!any(reqCols %in% gotColumns))
        stop(paste0("One or more required fields are not defined in",
                    " the submitted GRanges object! Need ",
                    paste(reqCols, collapse=","), " but got only ",
                    paste(reqCols[reqCols %in% gotColumns], collapse=","),
                    "."))
    exons <- as.data.frame(x[x$type == "exon", reqCols])[, -c(1, 4, 5)]
    ## for table tx2exon we want to have:
    ##    tx_id, exon_id, exon_idx
    t2e <- unique(exons[ , c("transcript_id", "exon_id", "exon_number")])
    colnames(t2e) <- c("tx_id", "exon_id", "exon_idx")
    ## Force exon_idx to be an integer!
    t2e[, "exon_idx"] <- as.integer(t2e[, "exon_idx"])
    ## Cross-check that we've got the corresponding tx_ids in the tx table!
    ## for table exons we want to have:
    ##    exon_id, exon_seq_start, exon_seq_end
    exons <- unique(exons[ , c("exon_id", "start", "end")])
    colnames(exons) <- c("exon_id", "exon_seq_start", "exon_seq_end")
    ## writing the tables.
    .checkIntegerCols(exons)
    .checkIntegerCols(t2e)
    dbWriteTable(con, name="exon", exons, overwrite=TRUE, row.names=FALSE)
    dbWriteTable(con, name="tx2exon", t2e, overwrite=TRUE, row.names=FALSE)
    message("OK")
    ## ----------------------------
    ##
    ## process chromosomes
    message("Processing chromosomes ... ", appendLF=FALSE)
    if (fetchSeqinfo) {
        ## problem is I don't have these available...
        chroms <- data.frame(seq_name = unique(as.character(genes$seq_name)))
        chroms <- cbind(chroms, seq_length = rep(NA, nrow(chroms)),
                        is_circular = rep(NA, nrow(chroms)))
        rownames(chroms) <- chroms$seq_name
        ## Try to get sequence lengths from Ensembl or Ensemblgenomes.
        sl <- tryGetSeqinfoFromEnsembl(organism, version,
                                       seqnames = chroms$seq_name)
        if (nrow(sl) > 0) {
            sl <- sl[sl[, "name"] %in% rownames(chroms), ]
            chroms[sl[, "name"], "seq_length"] <- sl[, "length"]
        }
    } else {
        ## have seqinfo available.
        chroms <- data.frame(seq_name = seqnames(Seqinfo),
                             seq_length = seqlengths(Seqinfo),
                             is_circular = isCircular(Seqinfo))
    }
    ## write the table.
    dbWriteTable(con, name="chromosome", chroms, overwrite=TRUE, row.names=FALSE)
    rm(genes)
    message("OK")
    message("Generating index ... ", appendLF=FALSE)
    ## generating all indices...
    .createEnsDbIndices(con)
    message("OK")
    message("  -------------")
    message("Verifying validity of the information in the database:")
    checkValidEnsDb(EnsDb(dbname))
    return(dbname)
}


## helper function that checks that the gene, transcript and exon data in the
## EnsDb database is correct (i.e. transcript within gene coordinates, exons within
## transcript coordinates, cds within transcript)
checkValidEnsDb <- function(x){
    message("Checking transcripts ... ", appendLF=FALSE)
    tx <- transcripts(x, columns=c("gene_id", "tx_id", "gene_seq_start",
                                   "gene_seq_end", "tx_seq_start",
                                   "tx_seq_end", "tx_cds_seq_start",
                                   "tx_cds_seq_end"), return.type="DataFrame")
    ## check if the tx are inside the genes...
    isInside <- tx$tx_seq_start >= tx$gene_seq_start &
        tx$tx_seq_end <= tx$gene_seq_end
    if(any(!isInside))
        stop("Start and end coordinates for ", sum(!isInside),
             "transcripts are not within the gene coordinates!")
    ## check cds coordinates
    notInside <- which(!(tx$tx_cds_seq_start >= tx$tx_seq_start &
                         tx$tx_cds_seq_end <= tx$tx_seq_end))
    if(length(notInside) > 0){
        stop("The CDS start and end coordinates for ", length(notInside),
             " transcripts are not within the transcript coordinates!")
    }
    rm(tx)
    message("OK\nChecking exons ... ", appendLF=FALSE)
    ex <- exons(x, columns=c("exon_id", "tx_id", "exon_seq_start", "exon_seq_end",
                             "tx_seq_start", "tx_seq_end", "seq_strand", "exon_idx"),
                return.type="data.frame")
    ## check if exons are within tx
    isInside <- ex$exon_seq_start >= ex$tx_seq_start &
        ex$exon_seq_end <= ex$tx_seq_end
    if(any(!isInside))
        stop("Start and end coordinates for ", sum(!isInside),
             " exons are not within the transcript coordinates!")
    ## checking the exon index...
    extmp <- ex[ex$seq_strand==1, c("exon_idx", "tx_id", "exon_seq_start")]
    extmp <- extmp[order(extmp$exon_seq_start), ]
    extmp.split <- split(extmp[ , c("exon_idx")], f=factor(extmp$tx_id))
    Different <- unlist(lapply(extmp.split, FUN=function(z){
                                   return(any(z != seq(1, length(z))))
                               }))
    if(any(Different)){
        stop(paste0("Provided exon index in transcript does not match with",
                    " ordering of the exons by chromosomal coordinates for",
                    sum(Different), "of the", length(Different),
                    "transcripts encoded on the + strand!"))
    }
    extmp <- ex[ex$seq_strand==-1, c("exon_idx", "tx_id", "exon_seq_end")]
    extmp <- extmp[order(extmp$exon_seq_end, decreasing=TRUE), ]
    extmp.split <- split(extmp[ , c("exon_idx")], f=factor(extmp$tx_id))
    Different <- unlist(lapply(extmp.split, FUN=function(z){
                                   return(any(z != seq(1, length(z))))
                               }))
    if(any(Different)){
        stop(paste0("Provided exon index in transcript does not match with",
                    " ordering of the exons by chromosomal coordinates for",
                    sum(Different), "of the", length(Different),
                    "transcripts encoded on the - strand!"))
    }
    message("OK")
    return(TRUE)
}


############################################################
##' Fetch chromosome sequence lengths from Ensembl.
##' @param organism The organism. Has to be in the form "Homo sapiens"
##' @param ensemblVersion The Ensembl version.
##' @param seqnames The names of the chromosomes/sequences; optional.
##' @return A matrix with two columns name and seq_length.
##' @noRd
tryGetSeqinfoFromEnsembl <- function(organism, ensemblVersion, seqnames,
                                     skip = FALSE){
    if (skip)
        return(matrix(nrow = 0, ncol = 2))
    message("Fetch seqlengths from ensembl ... ", appendLF=FALSE)
    tmp <- try(
        .getSeqlengthsFromMysqlFolder(organism = organism,
                                      ensembl = ensemblVersion,
                                      seqnames = seqnames)
    , silent = TRUE)
    if (is(tmp, "try-error") | is.null(tmp)) {
        message("FAIL")
        warning("Unable to retrieve sequence lengths from Ensembl.")
        return(matrix(nrow = 0, ncol = 2))
    }
    colnames(tmp) <- c("name", "length")
    return(tmp)
}


buildMetadata <- function(organism="", ensemblVersion="", genomeVersion="",
                          host="", sourceFile=""){
    MetaData <- data.frame(matrix(ncol=2, nrow=11))
    colnames(MetaData) <- c("name", "value")
    MetaData[1, ] <- c("Db type", "EnsDb")
    MetaData[2, ] <- c("Type of Gene ID", "Ensembl Gene ID")
    MetaData[3, ] <- c("Supporting package", "ensembldb")
    MetaData[4, ] <- c("Db created by", "ensembldb package from Bioconductor")
    MetaData[5, ] <- c("script_version", "0.0.1")
    MetaData[6, ] <- c("Creation time", date())
    MetaData[7, ] <- c("ensembl_version", ensemblVersion)
    MetaData[8, ] <- c("ensembl_host", host)
    MetaData[9, ] <- c("Organism", organism )
    MetaData[10, ] <- c("genome_build", genomeVersion)
    MetaData[11, ] <- c("DBSCHEMAVERSION", "1.0")
    MetaData[12, ] <- c("source_file", sourceFile)
    return(MetaData)
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
    if(metadataX["ensembl_version", "value"] ==
       metadataY["ensembl_version", "value"]){
        cat(" Ensembl versions match.\n")
    }else{
        cat(" WARNING: databases base on different Ensembl versions!",
            " Expect considerable differences!\n", sep = "")
        Messages["metadata"] <- "WARN"
    }
    ## genome build
    if(metadataX["genome_build", "value"] == metadataY["genome_build", "value"]){
        cat(" Genome builds match.\n")
    }else{
        cat(" WARNING: databases base on different Genome builds!",
            " Expect considerable differences!\n", sep = "")
        Messages["metadata"] <- "WARN"
    }
    if(length(idx)>0){
        cat(" All differences: <name>: <value x> != <value y>\n")
        for(i in idx){
            cat(paste("  - ", metadataX[i, "name"], ":", metadataX[i, "value"],
                      " != ", metadataY[i, "value"], "\n"))
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
    ## If we've got protein data in one of the two:
    if (hasProteinData(x) | hasProteinData(y)) {
        Messages <- c(Messages, protein = "OK")
        Messages["protein"] <- compareProteins(x, y)
    }
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
               length(onlyX), ") only in x, (", length(onlyY),
               ") only in y.\n" ))
    ## seqlengths:
    if (!all.equal(chromX[inboth, "seqlengths"],
                   chromY[inboth, "seqlengths"])) {
        same <- length(which(chromX[inboth, "seqlengths"] ==
                             chromY[inboth, "seqlengths"]))
        different <- length(inboth) - same
    } else {
        same <- length(inboth)
        different <- 0
    }
    cat(paste0( " Sequence lengths: (",same, ") identical, (",
               different, ") different.\n" ))
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
        which(as.character(seqnames(genesX[inboth])) ==
              as.character(seqnames(genesY[inboth])))
        )
    different <- length(inboth) - same
    if(different > 0)
        Ret <- "ERROR"
    cat(paste0( " Sequence names: (",same, ") identical, (",
               different, ") different.\n" ))
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
               length(onlyX), ") only in x, (", length(onlyY),
               ") only in y.\n"))
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
    ## cds start
    ## Makes sense to just compare for those that have the same tx!
    txXSub <- txX[inboth]
    txYSub <- txY[inboth]
    txCdsX <- names(txXSub)[!is.na(txXSub$tx_cds_seq_start)]
    txCdsY <- names(txYSub)[!is.na(txYSub$tx_cds_seq_start)]
    cdsInBoth <- txCdsX[txCdsX %in% txCdsY]
    cdsOnlyX <- txCdsX[!(txCdsX %in% txCdsY)]
    cdsOnlyY <- txCdsY[!(txCdsY %in% txCdsX)]
    if((length(cdsOnlyX) > 0 | length(cdsOnlyY)) & Ret!="ERROR")
        Ret <- "ERROR"
    cat(paste0(" Common transcripts with defined CDS: (", length(cdsInBoth),
               ") common, (", length(cdsOnlyX), ") only in x, (",
               length(cdsOnlyY), ") only in y.\n"))
    same <- length(
        which(txX[cdsInBoth]$tx_cds_seq_start == txY[cdsInBoth]$tx_cds_seq_start)
    )
    different <- length(cdsInBoth) - same
    if(different > 0 & Ret!="ERROR")
        Ret <- "ERROR"
    cat(paste0( " CDS start coordinates: (",same,
               ") identical, (", different, ") different.\n" ))
    ## cds end
    same <- length(
        which(txX[cdsInBoth]$tx_cds_seq_end == txY[cdsInBoth]$tx_cds_seq_end)
    )
    different <- length(cdsInBoth) - same
    if(different > 0 & Ret!="ERROR")
        Ret <- "ERROR"
    cat(paste0( " CDS end coordinates: (",same,
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

compareProteins <- function(x, y){
    cat("\nComparing protein data:\n")
    Ret <- "OK"
    if (!hasProteinData(x) | !hasProteinData(y)) {
        Ret <- "WARN"
        cat(paste0("No protein data available for one or both EnsDbs."))
        return(Ret)
    }
    X <- proteins(x)
    Y <- proteins(y)
    inboth <- X$protein_id[X$protein_id %in% Y$protein_id]
    onlyX <- X$protein_id[!(X$protein_id %in% Y$protein_id)]
    onlyY <- Y$protein_id[!(Y$protein_id %in% X$protein_id)]
    if(length(onlyX) > 0 | length(onlyY) > 0)
        Ret <- "WARN"
    cat(paste0(" protein IDs: (", length(inboth), ") common, (",
               length(onlyX), ") only in x, (", length(onlyY),
               ") only in y.\n"))
    X <- X[X$protein_id %in% inboth, ]
    Y <- Y[Y$protein_id %in% inboth, ]
    ## sorting both by protein_id should be enough.
    X <- X[order(X$protein_id), ]
    Y <- Y[order(Y$protein_id), ]

    ## tx_id
    same <- length(which(X$tx_id == Y$tx_id))
    different <- length(inboth) - same
    if(different > 0)
        Ret <- "ERROR"
    cat(paste0( " Transcript IDs: (",same,
               ") identical, (", different, ") different.\n" ))
    ## sequence
    same <- length(which(X$protein_sequence == Y$protein_sequence))
    different <- length(inboth) - same
    if(different > 0)
        Ret <- "ERROR"
    cat(paste0( " Protein sequence: (",same,
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

####============================================================
##  isEnsemblFileName
##
##  evaluate whether the file name is "most likely" corresponding
##  to a file name from Ensembl, i.e. following the convention
##  <organism>.<genome version>.<ensembl version>.[chr].gff/gtf.gz
##  The problem is that the genome version can also be . separated.
####------------------------------------------------------------
isEnsemblFileName <- function(x){
    x <- basename(x)
    ## If we split by ., do we get at least 4 elements?
    els <- unlist(strsplit(x, split=".", fixed=TRUE))
    if(length(els) < 4)
        return(FALSE)
    ## Can we get an Ensembl version?
    ensVer <- ensemblVersionFromGtfFileName(x)
    if(is.na(ensVer))
        return(FALSE)
    ## If we got one, do we still have enough fields left of the version?
    idx <- which(els == ensVer)
    idx <- idx[length(idx)]
    if(idx < 3){
        ## No way, we're missing the organism and the genome build field!
        return(FALSE)
    }
    ## Well, can not think of any other torture... let's assume it's OK.
    return(TRUE)
}
organismFromGtfFileName <- function(x){
    return(elementFromEnsemblFilename(x, 1))
}
####============================================================
##  ensemblVersionFromGtfFileName
##
##  Tries to extract the Ensembl version from the file name. If it
##  finds a numeric value it returns it, otherwise it returns NA.
####------------------------------------------------------------
ensemblVersionFromGtfFileName <- function(x){
    x <- basename(x)
    els <- unlist(strsplit(x, split=".", fixed=TRUE))
    ## Ensembl version is the last numeric value in the file name.
    for(elm in rev(els)){
        suppressWarnings(
            if(!is.na(as.numeric(elm))){
                return(elm)
            }
        )
    }
    return(NA)
}
####============================================================
##  genomeVersionFromGtfFileName
##
## the genome build can also contain .! thus, I return everything which is not
## the first element (i.e. organism), or the ensembl version, that is one left of
## the gtf.
genomeVersionFromGtfFileName <- function(x){
    x <- basename(x)
    els <- unlist(strsplit(x, split=".", fixed=TRUE))
    ensVer <- ensemblVersionFromGtfFileName(x)
    if(is.na(ensVer)){
        stop("Can not extract the genome version from the file name!",
             " The file name does not follow the expected naming convention from Ensembl!")
    }
    idx <- which(els == ensVer)
    idx <- idx[length(idx)]
    if(idx < 3)
        stop("Can not extract the genome version from the file name!",
             " The file name does not follow the expected naming convention from Ensembl!")
    return(paste(els[2:(idx-1)], collapse="."))
}

## Returns NULL if there was a problem.
elementFromEnsemblFilename <- function(x, which=1){
    tmp <- unlist(strsplit(x, split=.Platform$file.sep, fixed=TRUE))
    splitty <- unlist(strsplit(tmp[length(tmp)], split=".", fixed=TRUE))
    if(length(splitty) < which){
        warning("File ", x, " does not conform to the Ensembl file naming convention.")
        return(NULL)
    }
    return(splitty[which])
}

############################################################
## Utilities to fetch sequence lengths from Ensembl's ftp server, more
## specifically from the MySQL tables there.
## These replace the (unexported) functions from GenomicFeatures used thus far.
.ENSEMBL_URL <- "ftp://ftp.ensembl.org/pub/"
.ENSEMBLGENOMES_URL <- "ftp://ftp.ensemblgenomes.org/pub/"
##' Get the base url containing the mysql database for the specified host,
##' orgnism and Ensembl version.
##' @details The function will first build an approximate database name (without
##' the trailing <_genome version number> as this is not easy to guess).
##' Next all directories in the base MySQL folder will be scanned for the best
##' matching folder.
##' @param type Either "ensembl" or "ensemblgenomes"
##' @param organism Character specifying the organism. Has to be the full name,
##' i.e "homo_sapiens" or "Homo sapiens".
##' @param ensembl The Ensembl version.
##' @param genone The Genome version.
##' @noRd
.getEnsemblMysqlUrl <- function(type = "ensembl", organism, ensembl, genome) {
    type <- match.arg(type, c("ensembl", "ensemblgenomes"))
    if (type == "ensembl") {
        my_url <- paste0(.ENSEMBL_URL, "release-", ensembl, "/mysql/")
        db_name <- .guessDatabaseName(organism, ensembl)
        ## List folders; GenomicFeatures does it without 'dirlistonly',
        ## eventually that's what breaks on Windows?
        ## res <- getURL(my_url, dirlistonly = TRUE)
        res <- readLines(curl(my_url))
        res <- gsub(res, pattern = "\r", replacement = "", fixed = TRUE)
        if (length(res) > 0) {
            ## dirs <- unlist(strsplit(res, split = "\n"))
            ## ## Remove the \r on Windows.
            ## dirs <- sub(dirs, pattern = "\r", replacement = "", fixed = TRUE)
            ## idx <- grep(dirs, pattern = db_name)
            idx <- grep(res, pattern = db_name)
            if (length(idx) > 1)
                stop("Found more than one database matching '", db_name,
                     "' in Ensembl's ftp server!")
            if (length(idx) == 0)
                stop("No database matching '", db_name, "' found in Ensembl's",
                     " ftp server.")
            db_dir <- unlist(strsplit(res[idx], split = " ", fixed = TRUE))
            db_dir <- db_dir[length(db_dir)]
            return(paste0(my_url, db_dir))
        }
    } else {
        ## That's tricky! Have to find out whether the species is in plants,
        ## fungi, bacteria etc.
        ## List dirs of bacteria, fungi, metazoa, plants, protists
        sub_folders <- c("bacteria", "fungi", "metazoa", "plants", "protists")
        db_name <- .guessDatabaseName(organism, ensembl)
        for (fold in sub_folders) {
            my_url <- paste0(.ENSEMBLGENOMES_URL, "release-", ensembl, "/",
                             fold, "/mysql/")
            res <- try(readLines(curl(my_url)))
            if (is(res, "try-error")| length(res) == 0)
                next
            if (length(res) > 0) {
                res <- gsub(res, pattern = "\r", replacement = "", fixed = TRUE)
                idx <- grep(res, pattern = db_name)
                if (length(idx) > 1)
                    stop("Found more than one database matching '", db_name,
                         "' in Ensemblgenomes' ftp server!")
                if (length(idx) == 1) {
                    db_dir <- unlist(strsplit(res[idx], split = " ",
                                              fixed = TRUE))
                    db_dir <- db_dir[length(db_dir)]
                    return(paste0(my_url, db_dir))
                }
            }
            ## ## Catch eventual errors
            ## res <- try(getURL(my_url, dirlistonly = TRUE), silent = TRUE)
            ## if (is(res, "try-error") | length(res) == 0)
            ##     next
            ## if (length(res) > 0) {
            ##     dirs <- unlist(strsplit(res, split = "\n"))
            ##     ## Remove the \r on Windows.
            ##     dirs <- sub(dirs, pattern = "\r", replacement = "", fixed = TRUE)
            ##     idx <- grep(dirs, pattern = db_name)
            ##     if (length(idx) == 1)
            ##         return(paste0(my_url, dirs[idx]))
            ##     if (length(idx) > 1)
            ##         stop("Found more than one database matching '", db_name,
            ##              "' in Ensembl's ftp server!")
            ##     ## Well, then let's go to the next one.
            ## }
        }
        stop("No database matching '", db_name, "' found in Ensembl's",
             " ftp server")
    }
}

############################################################
## .guessDatabaseName
##' build the database name from species, ensembl version and genome version.
##' The latter is specifically difficult, as it is not quite clear how Ensembl
##' defines the Genome version number.
##' @param organism Character specifying the organism. Has to be in the format
##' "homo_sapiens" or "Homo sapiens", i.e. the full name.
##' @param ensembl The Ensembl version number.
##' @param genome The Genome version, e.g. GRCh38 (optional!).
##' @return A character representing the guessed database name in Ensembl.
##' @noRd
.guessDatabaseName <- function(organism, ensembl, genome) {
    if (missing(organism) & missing(ensembl))
        stop("'organism' and 'ensembl' are required!")
    ## Organism: all lower case, replace . with _
    organism <- tolower(gsub(organism, pattern = ".", replacement = "_",
                             fixed = TRUE))
    organism <- gsub(organism, pattern = " ", replacement = "_",
                     fixed = TRUE)
    dbname <- paste0(organism, "_core_", ensembl)
    ## Genome: remove all letters and keep just the numbers.
    if (!missing(genome)) {
        genome <- gsub(genome, pattern = "[a-zA-Z]", replacement = "")
        ## replace .0 at the end
        genome <- gsub(genome, pattern = ".0$", replacement = "")
        genome <- gsub(genome, pattern = ".", replacement = "", fixed = TRUE)
        genome <- gsub(genome, pattern = "_", replacement = "", fixed = TRUE)
        dbname <- paste0(dbname, "_", genome)
    }
    return(dbname)
}

############################################################
## .getSeqlengthsFromMysqlFolder
##' Fetch the coord_system.txt.gz and seq_region.txt.gz and extract the
##' seqlengths from there.
##' @noRd
.getSeqlengthsFromMysqlFolder <- function(organism, ensembl, seqnames) {
    ## Test whether we have the database in ensembl
    mysql_url <- try(.getEnsemblMysqlUrl(type = "ensembl", organism = organism,
                                         ensembl = ensembl), silent = TRUE)
    if (is(mysql_url, "try-error")) {
        mysql_url <- try(.getEnsemblMysqlUrl(type = "ensemblgenomes",
                                             organism = organism,
                                             ensembl = ensembl), silent = TRUE)
    }
    if (is(mysql_url, "try-error")) {
        warning("Can not get the sequence lengths from Ensembl or",
                " Ensemblgenomes. Seqinfo will lack the sequence lengths.")
        return(NULL)
    }
    ## Get the coord_system table
    coord_syst <- .getReadMysqlTable(mysql_url, "coord_system.txt.gz",
                                     colnames = c("coord_system_id",
                                                  "species_id", "name",
                                                  "version", "rank", "attrib"))
    ## Subset to the ones with "default" in "attrib"
    coord_syst <- coord_syst[grep(coord_syst$attrib, pattern = "default"),
                           , drop = FALSE]
    rownames(coord_syst) <- as.character(coord_syst$coord_system_id)
    ## Get the seq_region table
    seq_region <- .getReadMysqlTable(mysql_url, "seq_region.txt.gz",
                                     colnames = c("seq_region_id", "name",
                                                  "coord_system_id", "length"))
    ## Sub-set to the ones matching the coord_syst_ids and from these, select
    ## the one entry with the smallest rank, if more than one present.
    seq_region <- seq_region[seq_region$coord_system_id %in%
                             coord_syst$coord_system_id,
                           , drop = FALSE]
    seq_region <- cbind(seq_region,
                        rank = coord_syst[as.character(seq_region$coord_system_id),
                                          "rank"])
    ## Sub-set to the seqlevels we've got.
    if (!missing(seqnames)) {
        seq_region <- seq_region[seq_region$name %in% seqnames, , drop = FALSE]
        if (!all(seqnames %in% seq_region$name))
            warning("Could not determine length for all seqnames.")
    }
    sr <- split(seq_region, f = seq_region$name)
    if (length(sr) == 0)
        return(NULL)
    sr <- lapply(sr, function(z) {
        if (nrow(z) == 1)
            return(z)
        z <- z[order(z$rank), ]
        return(z[1, , drop = FALSE])
    })
    sr <- do.call(rbind, sr)
    rownames(sr) <- sr$name
    return(sr[, c("name", "length")])
}

##' Download and read a table from Ensembl
##' @param base_url the base url to the mysql folder on the server.
##' @param file_name the file name of the table.
##' @param colnames the column names.
##' @return A data.frame with the table's content.
##' @noRd
.getReadMysqlTable <- function(base_url, file_name, colnames) {
    tmp_file <- tempfile()
    download.file(url = paste0(base_url, "/", file_name), destfile = tmp_file,
                  quiet = TRUE)
    tmp <- read.table(tmp_file, sep = "\t", quote = "", comment.char = "",
                      as.is = TRUE)
    colnames(tmp) <- colnames
    return(tmp)
}

