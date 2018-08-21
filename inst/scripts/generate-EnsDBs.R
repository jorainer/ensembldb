## Functions related to create EnsDbs by downloading and installing MySQL
## databases from Ensembl.
library(RCurl)
library(RMariaDB)
library(ensembldb)

#' @description Get core database names from the specified folder.
#' 
#' @param ftp_folder The ftp url to the per-species mysql folders.
#' 
#' @author Johannes Rainer
#' 
#' @noRd
listCoreDbsInFolder <- function(ftp_folder) {
    if (missing(ftp_folder))
        stop("Argument 'ftp_folder' missing!")
    folders <- unlist(strsplit(getURL(ftp_folder,
                                      dirlistonly = TRUE), split = "\n"))
    res <- t(sapply(folders, function(z) {
        tmp <- unlist(strsplit(z, split = "_"))
        return(c(folder = z,
                 organism = paste0(tmp[1:2], collapse = "_"),
                 type = tmp[3],
                 version = paste0(tmp[4:length(tmp)], collapse = "_")))
    }))
    return(res[which(res[, "type"] == "core"), ])
}

#' @description Creates an EnsDb for the specified species by first downloading
#'     the corresponding MySQL database from Ensembl, installing it and
#'     subsequently creating the EnsDb database from it.
#'
#' @param ftp_folder The ftp url to the per-species mysql folders. If not
#'     provided it will use the default Ensembl ftp:
#'     \code{ftp://ftp.ensembl.org/pub/release-<ens_version>/mysql/}.
#' 
#' @param ens_version The Ensembl version (version of the Ensembl Perl API).
#' 
#' @param species The name of the species (e.g. "homo_sapiens").
#' 
#' @param user The user name for the MySQL/MariaDB database (write access).
#' 
#' @param host The host on which the MySQL/MariaDB database is running.
#' 
#' @param pass The password for the MySQL/MariaDB database.
#' 
#' @param port The port of the MySQL/MariaDB database.
#' 
#' @param local_tmp Local directory that will be used to temporarily store the
#'     downloaded MySQL/MariaDB database files.
#' 
#' @param dropDb Whether the Ensembl core database should be deleted once the
#'     EnsDb has been created.
#'
#' @author Johannes Rainer
#' 
#' @examples
#'
#' ## For Ensemblgenomes:
#' ftp_folder <- "ftp://ftp.ensemblgenomes.org/pub/release-33/fungi/mysql/"
#' @noRd
createEnsDbForSpecies <- function(ftp_folder,
                                  ens_version = 86, species, user, host, pass,
                                  port = 3306, local_tmp = tempdir(),
                                  sub_dir = "",
                                  dropDb = TRUE) {
    ## if ftp_folder is missing use the default one:
    base_url = "ftp://ftp.ensembl.org/pub"
    ## (1) Get all directories from Ensembl
    if (missing(ftp_folder))
        ftp_folder <- paste0(base_url, "/release-", ens_version, "/mysql/")
    res <- listCoreDbsInFolder(ftp_folder)

    folders <- unlist(strsplit(getURL(ftp_folder,
                                      dirlistonly = TRUE), split = "\n"))
    res <- t(sapply(folders, function(z) {
        tmp <- unlist(strsplit(z, split = "_"))
        return(c(folder = z,
                 organism = paste0(tmp[1:2], collapse = "_"),
                 type = tmp[3],
                 version = paste0(tmp[4:length(tmp)], collapse = "_")))
    }))
    res <- res[which(res[, "type"] == "core"), ]
    if (nrow(res) == 0)
        stop("No directories found!")
    if (missing(species))
        species <- res[, "organism"]
    rownames(res) <- res[, "organism"]
    ##     Check if we've got the species available
    got_specs <- species %in% rownames(res)
    if (!all(got_specs))
        warning("No core database for species ",
                paste0(species[!got_specs], collapse = ", "), " found.")
    species <- species[got_specs]
    res <- res[species, , drop = FALSE]
    if (length(species) == 0)
        stop("No database for any provided species found!")
    ## (2) Process each species
    message("Going to process ", nrow(res), " species.")
    for (i in 1:nrow(res)) {
        message("Processing species: ", res[i, "organism"], " (", i, " of ",
                nrow(res), ")")
        processOneSpecies(ftp_folder = paste0(ftp_folder, res[i, "folder"]),
                          ens_version = ens_version,
                          species = species[i], user = user, host = host,
                          pass = pass, port = port, local_tmp = local_tmp,
                          dropDb = dropDb)
        message("Done with species: ", res[i, "organism"], ", ",
                nrow(res) - i, " left.")
    }
}

#' @description This function performs the actual tasks of downloading the
#'     database files, installing them, deleting the download, creating the
#'     EnsDb and deleting the database.
#'
#' @details While the location of the downloaded temporary MySQL database file
#'     can be specified, the final SQLite file as well as all intermediate files
#'     will be placed in the current working directory.
#'
#' @param ftp_folder The folder on Ensembl's ftp server containing the mysql
#'     database files. Has to be the full path to these files.
#' 
#' @param ens_version The Ensembl version (version of the Ensembl Perl API).
#' 
#' @param species The name of the species (e.g. "homo_sapiens").
#' 
#' @param user The user name for the MySQL/MariaDB database (write access).
#' 
#' @param host The host on which the MySQL/MariaDB database is running.
#' 
#' @param pass The password for the MySQL/MariaDB database.
#' 
#' @param port The port of the MySQL/MariaDB database.
#' 
#' @param local_tmp Local directory that will be used to temporarily store the
#'     downloaded MySQL/MariaDB database files.
#' 
#' @param dropDb Whether the Ensembl core database should be deleted once the
#'     EnsDb has been created.
#'
#' @author Johannes Rainer
#' 
#' @noRd
processOneSpecies <- function(ftp_folder, ens_version = 86, species, user,
                              host = "localhost",
                              pass, port = 3306, local_tmp = tempdir(),
                              dropDb = TRUE) {
    if (missing(ftp_folder))
        stop("'ftp_folder' has to be specified!")
    if (missing(user))
        stop("'user' has to be specified!")
    if (missing(species))
        stop("'species' has to be specified!")
    ## (1) Download database files.
    res <- downloadFilesFromFtpFolder(url = ftp_folder, dest = local_tmp)
    ## (2) Install database.
    db_name <- basename(ftp_folder)
    res <- installEnsemblDb(dir = local_tmp, host = host, dbname = db_name,
                            user = user, pass = pass, port = port)
    ## (3) Delete the downloads.
    fls <- dir(local_tmp, full.names = TRUE)
    res <- sapply(fls, unlink)
    ## (4) Create the EnsDb (requires the correct Ensembl API)
    ##     They are created in the local directory.
    fetchTablesFromEnsembl(version = ens_version, species = species,
                           user = user, host = host, pass = pass, port = port)
    DBFile <- makeEnsemblSQLiteFromTables()
    unlink("*.txt")
    ## (5) Delete the database.
    if (dropDb) {
        con <- dbConnect(MariaDB(), host = host, user = user, pass = pass,
                         port = port, dbname = "mysql")
        dbSendStatement(con, paste("drop database ", db_name))
        dbDisconnect(con)
    }
}


#'
#' @description Download all files from an ftp directory to a local directory.
#'
#' @param url A character string specifying the url of the directory.
#'
#' @param dest A character string specifying the local directory.
#'
#' @return A character string with the path of the local directory.
#'
#' @author Johannes Rainer
#' 
#' @noRd
#'
#' @examples
#'
#' ftp_dir <- "ftp://ftp.ensembl.org/pub/release-88/mysql/homo_sapiens_core_88_38"
#' local_dir <- downloadFilesFromFtpFolder(ftp_dir)
downloadFilesFromFtpFolder <- function(url, dest = tempdir()) {
    fls <- getURL(paste0(url, "/"), dirlistonly = TRUE)
    fls <- unlist(strsplit(fls, split = "\n"))
    message("Downloading ", length(fls), " files ... ", appendLF = FALSE)
    for (i in 1:length(fls)) {
        download.file(url = paste0(url, "/", fls[i]),
                      destfile = paste0(dest, "/", fls[i]), quiet = TRUE)
    }
    message("OK")
    return(dest)
}

#' @description Install an Ensembl MySQL database downloaded from the Ensembl
#'     ftp server (e.g. using \link{downloadFilesFromFtpFolder}).
#'
#' @note The local directory is expected to correspond to the name of the
#'     database, i.e. \code{basename(dir)} will be used as the database name if
#'     argument \code{dbname} is missing.
#'
#' @param dir The path to the local directory where the database files are.
#' 
#' @param host The host running the MySQL database.
#' 
#' @param dbname The name of the database. If not provided the name of the
#'     provided directory will be used instead.
#' 
#' @param user The user name for the MySQL database (rw access).
#' 
#' @param pass The password for the MySQL database.
#' 
#' @param port The port of the MySQL database.
#' 
#' @author Johannes Rainer
#' 
#' @noRd
#'
#' @examples
#' user <- "user"
#' pass <- "pass"
#' dbname <- "homo_sapiens_core_88_38"
#' ## set to directory returned by the downloadFilesFromFtpFolder
#' dir <- local_dir
#' 
#' installEnsemblDb(dir = dir, dbname = dbname, user = user, pass = pass)
installEnsemblDb <- function(dir, host = "localhost", dbname, user, pass,
                              port = 3306) {
    if (missing(dir))
        stop("Argument 'dir' missing!")
    if (missing(dbname))
        dbname <- basename(dir)
    if (missing(user))
        stop("Argument 'user' missing!")
    ## Eventually unzip the files...
    tmp <- system(paste0("gunzip ", dir, "/*.gz"))
    ## Create the database
    con <- dbConnect(MariaDB(), host = host, user = user, pass = pass, port = port,
                     dbname = "mysql")
    dbSendStatement(con, paste0("create database ", dbname))
    dbDisconnect(con)
    ## Now create the tables and stuff.
    tmp <- system(paste0("mysql -h ", host, " -u ", user, " --password=", pass,
                         " -P ", port, " ", dbname, " < ", dir, "/", dbname,
                         ".sql"))
    ## Importing the data.
    cmd <- paste0("mysqlimport -h ", host, " -u ", user,
                  " --password=", pass, " -P ", port,
                  " ", dbname, " -L ", dir, "/*.txt")
    tmp <- system(cmd)
}

#' @description Creates EnsDb packages from all sqlite database files found in
#' the directory specified with parameter \code{dir}.
#' @param dir The path to the directory where the SQLite files can be found.
#' @param author The author of the package.
#' @param maintainer The maintainer of the package.
#' @param version The version of the package.
#' @noRd
createPackagesFromSQLite <- function(dir = ".", author, maintainer, version) {
    if (missing(author) | missing(maintainer) | missing(version))
        stop("Parameter 'author', 'maintainer' and 'version' are required!")
    edbs <- dir(dir, full.names = TRUE, pattern = ".sqlite")
    if (length(edbs) == 0)
        stop("Found no SQLite database files in the specified directory!")
    message("Processing ", length(edbs), " packages.")
    for (i in 1:length(edbs)) {
        message("Processing ", basename(edbs[i]), " (", i, " of ",
                length(edbs), ")", appendLF = FALSE)
        makeEnsembldbPackage(ensdb = edbs[i], version = version,
                             author = author, maintainer = maintainer)
        message("OK")
    }
}


## ftpf <- paste0("ftp://ftp.ensembl.org/pub/release-86/mysql/",
##                "anas_platyrhynchos_core_86_1")
## local_dir <- tempdir()
## downloadFilesFromFtpFolder(ftpf, dest = local_dir)
## installEnsemblDb(dir = local_dir, host = "localhost", user = "jo",
##                  pass = "jo123", dbname = "anas_platyrhynchos_core_86_1")
## fls <- dir(local_dir, full.names = TRUE)
## res <- sapply(fls, unlink)

## fetchTablesFromEnsembl(86, species = "anas_platyrhynchos", user = "jo",
##                        host = "localhost", pass = "jo123", port = 3306)
## DBFile <- makeEnsemblSQLiteFromTables()
## unlink("*.txt")

## system.time(fetchTablesFromEnsembl(86, species = "anas_platyrhynchos"))


## ftpf <- paste0("ftp://ftp.ensembl.org/pub/release-86/mysql/",
##                "homo_sapiens_core_86_38")
## local_dir <- tempdir()
## processOneSpecies(ftp_folder = ftpf, version = 86,
##                   species = "homo_sapiens", user = "jo",
##                   host = "localhost",
##                   pass = "jo123", port = 3306, local_tmp = local_dir,
##                   dropDb = FALSE)


## Add an issue:
## + Fix problem of non-defined sequence type "chromosome" in anas platyrhynchos
##   database. -> update to the perl script.
## + Compare Hsapiens EnsDb created with new script and the "original" one.
