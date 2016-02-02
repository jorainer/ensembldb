####============================================================
##  Methods and functions to allow usage of EnsDb objects also
##  with genomic resources that do not use Ensembl based
##  seqnames
##  We're storing the seqname style into the .properties slot
##  of the EnsDb object.
####------------------------------------------------------------
.ENSOPT.SEQNOTFOUND="ensembldb.seqnameNotFound"
####============================================================
##  formatSeqnamesForQuery
##
##  Formating/renamaing the seqname(s) according to the specified
##  style.
##  x is an EnsDb,
##  sn the seqnames to convert...
##  If a seqname can not be mapped NA will be returned.
####------------------------------------------------------------
setMethod("formatSeqnamesForQuery", "EnsDb", function(x, sn, ifNotFound){
    return(.formatSeqnameByStyleForQuery(x, sn, ifNotFound))
})
## Little helper function that returns eventually the argument.
## Returns MISSING if the argument was not set.
.getSeqnameNotFoundOption <- function(){
    notFound <- "MISSING"
    if(any(names(options()) == .ENSOPT.SEQNOTFOUND)){
        notFound <- getOption(.ENSOPT.SEQNOTFOUND)
        ## Do some sanity checks?
    }
    return(notFound)
}
.formatSeqnameByStyleForQuery <- function(x, sn, ifNotFound){
    ## Fixing ifNotFound, allowing that this can be set using options.
    if(missing(ifNotFound)){
        ifNotFound <- .getSeqnameNotFoundOption()
    }
    ## Map whatever to Ensembl seqnames, such that we can perform queries.
    ## Use mapSeqlevels, or rather genomeStyles and do it hand-crafted!
    sst <- seqlevelsStyle(x)
    dbSst <- dbSeqlevelsStyle(x)
    if(sst == dbSst)
        return(sn)
    ## Don't like that the genomeStyles is reading the stuff form file.
    map <- getProperty(x, "genomeStyle")
    if(!is(map, "data.frame"))
        map <- genomeStyles(organism(x))
    ## sn are supposed to be in sst style, map them to dbSst
    idx <- match(sn, map[, sst])
    mapped <- map[idx, dbSst]
    if(any(is.na(mapped))){
        noMap <- which(is.na(mapped))
        seqNoMap <- unique(sn[noMap])
        if(length(seqNoMap) > 5){
            theMess <- paste0("More than 5 seqnames could not be mapped to ",
                              "the seqlevels style of the database (", dbSst, ")!")
        }else{
            theMess <- paste0("Seqnames: ", paste(seqNoMap, collapse=", "),
                              " could not be mapped to ",
                              " the seqlevels style of the database (", dbSst, ")!")
        }
        if(is.na(ifNotFound) | is.null(ifNotFound)){
            ## Replacing the missing seqname mappings with ifNotFound.
            mapped[noMap] <- ifNotFound
            warnMess <- paste0(" Returning ", ifNotFound, " for these.")
        }else{
            ## If MISSING -> STOP
            if(ifNotFound == "MISSING"){
                stop(theMess)
            }else{
                ## Next special case: use the original names, i.e. don't map at all.
                if(ifNotFound == "ORIGINAL"){
                    mapped[noMap] <- sn[noMap]
                    warnMess <- "Returning the orginal seqnames for these."
                }else{
                    mapped[noMap] <- ifNotFound
                    warnMess <- paste0(" Returning ", ifNotFound, " for these.")
                }
            }
        }
        warning(theMess, warnMess)
    }
    return(mapped)
}
setMethod("formatSeqnamesFromQuery", "EnsDb", function(x, sn, ifNotFound){
    return(.formatSeqnameByStyleFromQuery(x, sn, ifNotFound))
})
.formatSeqnameByStyleFromQuery <- function(x, sn, ifNotFound){
    ## Fixing ifNotFound, allowing that this can be set using options.
    if(missing(ifNotFound)){
        ifNotFound <- .getSeqnameNotFoundOption()
    }
    ## Map Ensembl seqnames resulting form queries to the seqlevel style by
    ## seqlevelsStyle.
    sst <- seqlevelsStyle(x)
    dbSst <- dbSeqlevelsStyle(x)
    if(sst == dbSst)
        return(sn)
    ## Otherwise...
    map <- getProperty(x, "genomeStyle")
    if(!is(map, "data.frame"))
        map <- genomeStyles(organism(x))
    ## sn are supposed to be in sst style, map them to dbSst
    idx <- match(sn, map[, dbSst])
    mapped <- map[idx, sst]
    if(any(is.na(mapped))){
        noMap <- which(is.na(mapped))
        seqNoMap <- unique(sn[noMap])
        if(length(seqNoMap) > 5){
            theMess <- paste0("More than 5 seqnames with seqlevels style of the database (",
                              dbSst, ") could not be mapped to the seqlevels style: ",
                              sst, "!)")
        }else{
            theMess <- paste0("Seqnames: ", paste(seqNoMap, collapse=", "),
                              " with seqlevels style of the database (", dbSst,
                              ") could not be mapped to seqlevels style: ", sst,
                              "!")
        }

        if(is.na(ifNotFound) | is.null(ifNotFound)){
            ## Replacing the missing seqname mappings with ifNotFound.
            mapped[noMap] <- ifNotFound
            warnMess <- paste0(" Returning ", ifNotFound, " for these.")
        }else{
            ## If MISSING -> STOP
            if(ifNotFound == "MISSING"){
                stop(theMess)
            }else{
                ## Next special case: use the original names, i.e. don't map at all.
                if(ifNotFound == "ORIGINAL"){
                    mapped[noMap] <- sn[noMap]
                    warnMess <- " Returning the orginal seqnames for these."
                }else{
                    mapped[noMap] <- ifNotFound
                    warnMess <- paste0(" Returning ", ifNotFound, " for these.")
                }
            }
        }
        warning(theMess, warnMess)
    }
    return(mapped)
}


####============================================================
##  dbSeqlevelsStyle
##
##  Returns the seqname style used by the database. Defaults to
##  Ensembl and reads the property: dbSeqlevelsStyle.
####------------------------------------------------------------
setMethod("dbSeqlevelsStyle", "EnsDb", function(x){
    stl <- getProperty(x, "dbSeqlevelsStyle")
    if(is.na(stl))
        stl <- "Ensembl"
    return(stl)
})
####============================================================
##  seqlevelStyle
##
##  Get or set the seqlevel style. If we can't find the stype in
##  GenomeInfoDb throw and error.
####------------------------------------------------------------
setMethod("seqlevelsStyle", "EnsDb", function(x){
    st <- getProperty(x, "seqlevelsStyle")
    if(is.na(st))
        st <- "Ensembl"
    return(st)
})
setReplaceMethod("seqlevelsStyle", "EnsDb", function(x, value){
    if(value == dbSeqlevelsStyle(x)){
        ## Not much to do; that's absolutely fine.
        x <- setProperty(x, seqlevelsStyle=value)
    }else{
        ## Have to check whether I have the mapping available in GenomeInfoDb, if not
        ## -> throw an error.
        dbStyle <- dbSeqlevelsStyle(x)
        ## Note that both, the db seqlevel style and the style have to be available!
        ## Check if we could use the mapping provided by GenomeInfoDb.
        genSt <- try(genomeStyles(organism(x)), silent=TRUE)
        if(is(genSt, "try-error")){
            stop("No mapping of seqlevel styles available in GenomeInfoDb for",
                 " species ", organism(x), "! Please refer to the Vignette of the",
                 " GenomeInfoDb package if you would like to provide this mapping.")
        }
        if(!any(colnames(genSt) == value)){
            stop("The provided seqlevels style is not known to GenomeInfoDb!")
        }
        if(!any(colnames(genSt) == dbStyle)){
            stop("The seqlevels style of the database (", dbStyle,
                 ") is not known to GenomeInfoDb!")
        }
        ## If we got that far it should be OK
        x <- setProperty(x, seqlevelsStyle=value)
        x <- setProperty(x, genomeStyle=list(genSt))
    }
    return(x)
})

####============================================================
##  supportedSeqlevelsStyles
##
##  Get all supported seqlevels styles for the species of the EnsDb
####------------------------------------------------------------
setMethod("supportedSeqlevelsStyles", "EnsDb", function(x){
    map <- genomeStyles(organism(x))
    cn <- colnames(map)
    cn <- cn[!(cn %in% c("circular", "auto", "sex"))]
    return(colnames(cn))
})


####==================== OLD STUFF BELOW ====================

###==============================================================
##  Prefix chromosome names with "chr" if ucscChromosomeNames option
##  is set, otherwise, use chromosome names "as is".
##  This function should be used in functions that return results from
##  EnsDbs.
###--------------------------------------------------------------
prefixChromName <- function(x, prefix="chr"){
    ucsc <- getOption("ucscChromosomeNames", default=FALSE)
    if(ucsc){
        ## TODO fix also the mitochondrial chromosome name.
        mapping <- ucscToEnsMapping()
        for(i in 1:length(mapping)){
            x <- sub(x, pattern=names(mapping)[i], replacement=mapping[i],
                     fixed=TRUE)
        }
        ## Replace chr if it's already there
        x <- gsub(x, pattern="^chr", replacement="", ignore.case=TRUE)
        x <- paste0(prefix, x)
    }
    return(x)
}

###==============================================================
##  Remove leading "chr" to fit Ensembl based chromosome names.
##  This function should be called in functions that fetch data from
##  EnsDbs.
###--------------------------------------------------------------
ucscToEns <- function(x){
    ## TODO rename all additional chromosome names.
    mapping <- ucscToEnsMapping()
    for(i in 1:length(mapping)){
        x <- sub(x, pattern=mapping[i], replacement=names(mapping)[i],
                 fixed=TRUE)
    }
    x <- gsub(x, pattern="^chr", replacement="", ignore.case=TRUE)
    return(x)
}
###============================================================
##  Returns a character vector, elements representing UCSC chromosome
##  names with their names corresponding to the respective Ensembl
##  chromosome names.
###------------------------------------------------------------
ucscToEnsMapping <- function(){
    theMap <- c(MT="chrM")
    return(theMap)
}

