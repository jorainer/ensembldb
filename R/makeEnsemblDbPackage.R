## part of this code is from GenomicFeatures makeTxDbPackage.R
## So to make a package we need a couple things:
## 1) we need a method called makeTxDbPackage (that will take a txdb object)
## 2) we will need a package template to use



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
    fn <- system.file("perl", "get_gene_transcript_exon_tables.pl", package="ensembldb")
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
        stop("Something went wrong! I'm missing some of the txt files the perl script should have generated.")

    ## read information
    info <- read.table(paste0(path, .Platform$file.sep ,"ens_metadata.txt"), sep="\t",
                       as.is=TRUE, header=TRUE)
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
    tmp <- read.table(paste0(path, .Platform$file.sep ,"ens_chromosome.txt"), sep="\t", as.is=TRUE, header=TRUE)
    dbWriteTable(con, name="chromosome", tmp, row.names=FALSE)
    rm(tmp)
    ## make index
    dbGetQuery(con, "create index seq_name_idx on chromosome (seq_name);")

    ## process genes: some gene names might have fancy names...
    tmp <- read.table(paste0(path, .Platform$file.sep, "ens_gene.txt"), sep="\t", as.is=TRUE, header=TRUE,
                     quote="", comment.char="" )
    dbWriteTable(con, name="gene", tmp, row.names=FALSE)
    rm(tmp)
    ## make index
    dbGetQuery(con, "create index gene_id_idx on gene (gene_id);")

    ## process transcripts:
    tmp <- read.table(paste0(path, .Platform$file.sep, "ens_tx.txt"), sep="\t", as.is=TRUE, header=TRUE)
    dbWriteTable(con, name="tx", tmp, row.names=FALSE)
    rm(tmp)
    ## make index
    dbGetQuery(con, "create index tx_id_idx on tx (tx_id);")

    ## process exons:
    tmp <- read.table(paste0(path, .Platform$file.sep, "ens_exon.txt"), sep="\t", as.is=TRUE, header=TRUE)
    dbWriteTable(con, name="exon", tmp, row.names=FALSE)
    rm(tmp)
    ## make index
    dbGetQuery(con, "create index exon_id_idx on exon (exon_id);")
    tmp <- read.table(paste0(path, .Platform$file.sep, "ens_tx2exon.txt"), sep="\t", as.is=TRUE, header=TRUE)
    dbWriteTable(con, name="tx2exon", tmp, row.names=FALSE)
    rm(tmp)
    ## make index
    dbGetQuery(con, "create index t2e_tx_id_idx on tx2exon (tx_id);")
    dbGetQuery(con, "create index t2e_exon_id_idx on tx2exon (exon_id);")
    dbDisconnect(con)
    ## done.
    return(dbname)
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
        PKGDESCRIPTION=paste("Exposes an annotation databases generated from Ensembl."),
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
        ORGANISMBIOCVIEW=gsub(" ","_",.organismName(.getMetaDataValue(con ,'Organism'))),
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
                      collapse=.Platform$file.sep), showWarnings=FALSE, recursive=TRUE)
    db_path <- file.path(destDir, pkgName, "inst", "extdata",
                         paste(pkgName,"sqlite",sep="."))
    file.copy(ensdbfile, to=db_path)
}

