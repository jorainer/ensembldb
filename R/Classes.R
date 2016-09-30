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


###============================================================
##  SymbolFilter
###------------------------------------------------------------
setClass("SymbolFilter", contains = "BasicFilter",
         prototype = list(
             condition = "=",
             value = "",
             .valueIsCharacter = TRUE
         )
         )
SymbolFilter <- function(value, condition = "=") {
    if(missing(value)){
        stop("A filter without a value makes no sense!")
    }
    if(length(value) > 1) {
        if(condition == "=")
            condition = "in"
        if(condition == "!=")
            condition = "not in"
    }
    return(new("SymbolFilter", condition = condition,
               value = as.character(value)))
}

############################################################
## OnlyCodingTx
##
## That's a special case filter that just returns transcripts
## that have tx_cds_seq_start defined (i.e. not NULL).
setClass("OnlyCodingTx", contains = "BasicFilter",
         prototype = list(
             condition = "=",
             value = "",
             .valueIsCharacter = TRUE
         ))
OnlyCodingTx <- function() {
    return(new("OnlyCodingTx"))
}

#' Protein annotation related filters
#'
#' @description These filters allow to query specific information from an
#' \code{\linkS4class{EnsDb}} database. Filters should be created using the
#' dedicated constructor functions \code{ProteinidFilter},
#' \code{ProtdomidFilter} and \code{UniprotidFilter}. Protein annotation-based
#' filters are:
#'
#' @slot condition The condition to be used in the filter.
#' @slot value The value(s) for the filter.
#' @note These filters can only be used if the \code{\linkS4class{EnsDb}}
#' database contains protein annotations, i.e. if \code{\link{hasProteinData}}
#' is \code{TRUE}. Also, only protein coding transcripts will have protein
#' annotations available, thus, non-coding transcripts/genes will not be
#' returned by the queries.
#' @name ProteinFilters
#' @seealso \link{RNA-DNA-filters} for filters using RNA/DNA annotations.
#' @author Johannes Rainer
#' @examples
#' library(EnsDb.Hsapiens.v75)
#' edb <- EnsDb.Hsapiens.v75
#' ## Create a ProteinidFilter:
#' pf <- ProteinidFilter("ENSP00000362111")
#' pf
#' ## Using this filter would retrieve all database entries matching the
#' ## condition:
#' where(pf)
#' ## i.e. all entries that are associated with a protein with the ID
#' ## "ENSP00000362111"
#' if (hasProteinData(edb)) {
#'     res <- genes(edb, filter = pf)
#'     res
#' }
#'
#' ## To get all entries except those associated with that protein:
#' condition(pf) <- "!="
#' where(pf)
#'
#' ## UniprotidFilter:
#' uf <- UniprotidFilter("O60762")
#' ## Get the transcripts encoding that protein:
#' if (hasProteinData(edb)) {
#'     transcripts(edb, filter = uf)
#'     ## The mapping Ensembl protein ID to Uniprot ID can however be 1:n:
#'     transcripts(edb, filter = TxidFilter("ENST00000371588"),
#'         columns = c("protein_id", "uniprot_id"))
#' }
#'
#' ## ProtdomidFilter:
#' pdf <- ProtdomidFilter("PF00335")
#' ## Also here we could get all transcripts related to that protein domain
#' if (hasProteinData(edb)) {
#'     transcripts(edb, filter = pdf, columns = "protein_id")
#' }
NULL
#> NULL

############################################################
## ProteinidFilter
##' @description The \code{ProteinidFilter} allows to filter based on the
##' (Ensembl) protein ID of the (coding) transcripts' proteins.
##' @rdname ProteinFilters
setClass("ProteinidFilter", contains = "BasicFilter",
         prototype = list(
             condition = "=",
             value = "",
             .valueIsCharacter = TRUE
         ))
##' @param value A character vector of length 1 or larger with the value(s)
##' for the filter.
##' @param condition A character of length 1 specifying the condition of the
##' filter. Can be one of \code{"="}, \code{"!="}, \code{"like"}, or, for
##' \code{value} of length larger 1 \code{"in"} and \code{"not in"}.
##' @return For \code{ProteinidFilter}: A \code{ProteinidFilter} object.
##' @rdname ProteinFilters
ProteinidFilter <- function(value, condition = "=") {
    if (missing(value)) {
        stop("A filter without a value makes no sense!")
    }
    if (length(value) > 1) {
        if(condition == "=")
            condition = "in"
        if(condition == "!=")
            condition = "not in"
    }
    return(new("ProteinidFilter", condition = condition,
               value = as.character(value)))
}

############################################################
## UniprotFilter
##' @description The \code{UniprotidFilter} allows to retrieve annotations
##' filtering by the provided Uniprot ID associated with the transcript's
##' protein.
##' @rdname ProteinFilters
setClass("UniprotidFilter", contains = "BasicFilter",
         prototype = list(
             condition = "=",
             value = "",
             .valueIsCharacter = TRUE
         ))
##' @return For \code{UniprotidFilter}: A \code{UniprotidFilter} object.
##' @rdname ProteinFilters
UniprotidFilter <- function(value, condition = "=") {
    if (missing(value)) {
        stop("A filter without a value makes no sense!")
    }
    if (length(value) > 1) {
        if(condition == "=")
            condition = "in"
        if(condition == "!=")
            condition = "not in"
    }
    return(new("UniprotidFilter", condition = condition,
               value = as.character(value)))
}

############################################################
## ProtdomidFilter
##' @description The \code{ProtdomidFilter} allows to retrieve entries from
##' the database matching the provided filter criteria based on their protein
##' domain ID (\emph{protein_domain_id}).
##' @rdname ProteinFilters
setClass("ProtdomidFilter", contains = "BasicFilter",
         prototype = list(
             condition = "=",
             value = "",
             .valueIsCharacter = TRUE
         ))
##' @return For \code{ProtdomidFilter}: A \code{ProtdomidFilter} object.
##' @rdname ProteinFilters
ProtdomidFilter <- function(value, condition = "=") {
    if (missing(value)) {
        stop("A filter without a value makes no sense!")
    }
    if (length(value) > 1) {
        if(condition == "=")
            condition = "in"
        if(condition == "!=")
            condition = "not in"
    }
    return(new("ProtdomidFilter", condition = condition,
               value = as.character(value)))
}
