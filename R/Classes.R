##***********************************************************************
##
##     EnsBb classes
##
##     Main class providing access and functionality for the database.
##
##***********************************************************************
setClass("EnsDb",
         representation(ensdb="DBIConnection", tables="list", .properties="list"),
         prototype=list(ensdb=NULL, tables=list(), .properties=list())
        )


##***********************************************************************
##
##     BasicFilter classes
##
##     Allow to filter the results fetched from the database.
##
##     gene:
##     - GeneidFilter
##     - GenebiotypeFilter
##     - GenenameFilter
##     - EntrezidFilter
##
##     transcript:
##     - TxidFilter
##     - TxbiotypeFilter
##
##     exon:
##     - ExonidFilter
##
##     chrom position (using info from exon):
##     - SeqnameFilter
##     - SeqstartFilter
##     - SeqendFilter
##     - SeqstrandFilter
##     alternative: GRangesFilter. See below.
##
##***********************************************************************
setClass("BasicFilter",
         representation(
             "VIRTUAL",
             condition="character",
             value="character",
             .valueIsCharacter="logical"
            ),
         prototype=list(
             condition="=",
             value="",
             .valueIsCharacter=TRUE
            )
        )

## Table gene
## filter for gene_id
setClass("GeneidFilter", contains="BasicFilter",
         prototype=list(
             condition="=",
             value="",
             .valueIsCharacter=TRUE
            )
        )
GeneidFilter <- function(value, condition="="){
    if(missing(value)){
        stop("A filter without a value makes no sense!")
    }
    if(length(value) > 1){
        if(condition=="=")
            condition="in"
        if(condition=="!=")
            condition="not in"
    }
    return(new("GeneidFilter", condition=condition, value=as.character(value)))
}
## filter for gene_biotype
setClass("GenebiotypeFilter", contains="BasicFilter",
         prototype=list(
             condition="=",
             value="",
             .valueIsCharacter=TRUE
            )
        )
GenebiotypeFilter <- function(value, condition="="){
    if(missing(value)){
        stop("A filter without a value makes no sense!")
    }
    if(length(value) > 1){
        if(condition=="=")
            condition="in"
        if(condition=="!=")
            condition="not in"
    }
    return(new("GenebiotypeFilter", condition=condition, value=as.character(value)))
}
## filter for gene_name
setClass("GenenameFilter", contains="BasicFilter",
         prototype=list(
             condition="=",
             value="",
             .valueIsCharacter=TRUE
            )
        )
GenenameFilter <- function(value, condition="="){
    if(missing(value)){
        stop("A filter without a value makes no sense!")
    }
    if(length(value) > 1){
        if(condition=="=")
            condition="in"
        if(condition=="!=")
            condition="not in"
    }
    return(new("GenenameFilter", condition=condition, value=as.character(value)))
}
## filter for entrezid
setClass("EntrezidFilter", contains="BasicFilter",
         prototype=list(
             condition="=",
             value="",
             .valueIsCharacter=TRUE
            )
        )
EntrezidFilter <- function(value, condition="="){
    if(missing(value)){
        stop("A filter without a value makes no sense!")
    }
    if(length(value) > 1){
        if(condition=="=")
            condition="in"
        if(condition=="!=")
            condition="not in"
    }
    return(new("EntrezidFilter", condition=condition, value=as.character(value)))
}


## Table transcript
## filter for tx_id
setClass("TxidFilter", contains="BasicFilter",
         prototype=list(
             condition="=",
             value="",
             .valueIsCharacter=TRUE
            )
        )
TxidFilter <- function(value, condition="="){
    if(missing(value)){
        stop("A filter without a value makes no sense!")
    }
    if(length(value) > 1){
        if(condition=="=")
            condition="in"
        if(condition=="!=")
            condition="not in"
    }
    return(new("TxidFilter", condition=condition, value=as.character(value)))
}
## filter for gene_biotype
setClass("TxbiotypeFilter", contains="BasicFilter",
         prototype=list(
             condition="=",
             value="",
             .valueIsCharacter=TRUE
            )
        )
TxbiotypeFilter <- function(value, condition="="){
    if(missing(value)){
        stop("A filter without a value makes no sense!")
    }
    if(length(value) > 1){
        if(condition=="=")
            condition="in"
        if(condition=="!=")
            condition="not in"
    }
    return(new("TxbiotypeFilter", condition=condition, value=as.character(value)))
}

## Table exon
## filter for exon_id
setClass("ExonidFilter", contains="BasicFilter",
         prototype=list(
             condition="=",
             value="",
             .valueIsCharacter=TRUE
            )
        )
ExonidFilter <- function(value, condition="="){
    if(missing(value)){
        stop("A filter without a value makes no sense!")
    }
    if(length(value) > 1){
        if(condition=="=")
            condition="in"
        if(condition=="!=")
            condition="not in"
    }
    return(new("ExonidFilter", condition=condition, value=as.character(value)))
}

## Table tx2exon
## filter for exon_idx
setClass("ExonrankFilter", contains="BasicFilter",
         prototype=list(
             condition="=",
             value="",
             .valueIsCharacter=FALSE
            )
        )
ExonrankFilter <- function(value, condition="="){
    if(missing(value)){
        stop("A filter without a value makes no sense!")
    }
    if(any(is.na(as.numeric(value))))
        stop("Argument 'value' has to be numeric!")
    if(length(value) > 1){
        if(condition=="=")
            condition="in"
        if(condition=="!=")
            condition="not in"
    }
    return(new("ExonrankFilter", condition=condition, value=as.character(value)))
}


## chromosome positions
## basic chromosome/seqname filter.
setClass("SeqnameFilter", contains="BasicFilter",
         prototype=list(
             condition="=",
             value="",
             .valueIsCharacter=TRUE
            )
        )
## builder...
SeqnameFilter <- function(value, condition="="){
    if(missing(value)){
        stop("A filter without a value makes no sense!")
    }
    if(length(value) > 1){
        if(condition=="=")
            condition="in"
        if(condition=="!=")
            condition="not in"
    }
    return(new("SeqnameFilter", condition=condition, value=as.character(value)))
}

## basic chromosome strand filter.
setClass("SeqstrandFilter", contains="BasicFilter",
         prototype=list(
             condition="=",
             value="",
             .valueIsCharacter=FALSE
            )
        )
## builder...
SeqstrandFilter <- function(value, condition="="){
    if(missing(value)){
        stop("A filter without a value makes no sense!")
    }
    ## checking value: should be +, -, will however be translated to -1, 1
    if(class(value)=="character"){
        value <- match.arg(value, c("1", "-1", "+1", "-", "+"))
        if(value=="-")
            value <- "-1"
        if(value=="+")
            value <- "+1"
        ## OK, now transforming to number
        value <- as.numeric(value)
    }
    if(!(value==1 | value==-1))
        stop("The strand has to be either 1 or -1 (or \"+\" or \"-\")")
    return(new("SeqstrandFilter", condition=condition, value=as.character(value)))
}

## chromstart filter
setClass("SeqstartFilter", contains="BasicFilter",
         representation(
             feature="character"
            ),
         prototype=list(
             condition=">",
             value="",
             .valueIsCharacter=FALSE,
             feature="gene"
            )
        )
SeqstartFilter <- function(value, condition="=", feature="gene"){
    if(missing(value)){
        stop("A filter without a value makes no sense!")
    }
    if(length(value) > 1){
        value <- value[ 1 ]
        warning("Multiple values provided, but only the first (", value,") will be considered")
    }
    return(new("SeqstartFilter", condition=condition, value=as.character(value),
                feature=feature))
}

## chromend filter
setClass("SeqendFilter", contains="BasicFilter",
         representation(
             feature="character"
            ),
         prototype=list(
             condition="<",
             value="",
             .valueIsCharacter=FALSE,
             feature="gene"
            )
        )
SeqendFilter <- function(value, condition="=", feature="gene"){
    if(missing(value)){
        stop("A filter without a value makes no sense!")
    }
    if(length(value) > 1){
        value <- value[ 1 ]
        warning("Multiple values provided, but only the first (", value,") will be considered")
    }
    return(new("SeqendFilter", condition=condition, value=as.character(value),
                feature=feature))
}


###============================================================
##  GRangesFilter
##  adding new arguments since we can not overwrite the data type
##  of the BasicFilter class... unfortunately.
##  + grange <- value
##  + location <- condition
###------------------------------------------------------------
setClass("GRangesFilter", contains="BasicFilter",
         representation(grange="GRanges",
                        feature="character",
                        location="character"),
         prototype=list(
             grange=GRanges(),
             .valueIsCharacter=FALSE,
             condition="=",
             location="within",
             feature="gene",
             value=""
         ))
## Constructor
GRangesFilter <- function(value, condition="within", feature="gene"){
    if(missing(value))
        stop("No value provided for the filter!")
    if(!is(value, "GRanges"))
        stop("'value' has to be a GRanges object!")
    if(length(value) == 0)
        stop("No value provided for the filter!")
    ## if(length(value) > 1){
    ##     warning(paste0("GRanges in 'value' has length ", length(value),
    ##                    "! Using only the first element!"))
    ##     value <- value[1]
    ## }
    grf <- new("GRangesFilter", grange=value, location=condition,
               feature=feature)
    ##validObject(grf)
    return(grf)
}
###------------------------------------------------------------


