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
## setClass("BasicFilter",
##          representation(
##              "VIRTUAL",
##              condition="character",
##              value="character",
##              .valueIsCharacter="logical"
##             ),
##          prototype=list(
##              condition="=",
##              value="",
##              .valueIsCharacter=TRUE
##             )
##         )

## Table gene
## filter for gene_id
## setClass("GeneidFilter", contains="AnnotationFilter",
##          prototype = list(
##              condition = "==",
##              value = character(),
##              .valueIsCharacter = TRUE
##          )
##          )
## ## filter for gene_biotype
## setClass("GenebiotypeFilter", contains="AnnotationFilter",
##          prototype = list(
##              condition = "==",
##              value = "",
##              .valueIsCharacter = TRUE
##          )
##          )
## filter for gene_name
## setClass("GenenameFilter", contains = "AnnotationFilter",
##          prototype = list(
##              condition = "==",
##              value = "",
##              .valueIsCharacter = TRUE
##          )
##          )
## GenenameFilter <- function(value, condition = "=="){
##     if(missing(value))
##         stop("A filter without a value makes no sense!")
##     return(new("GenenameFilter", condition = condition,
##                value = as.character(value)))
## }
## filter for entrezid
## setClass("EntrezidFilter", contains = "AnnotationFilter",
##          prototype = list(
##              condition = "==",
##              value = "",
##              .valueIsCharacter = TRUE
##             )
##         )


## Table transcript
## filter for tx_id
## setClass("TxidFilter", contains="AnnotationFilter",
##          prototype = list(
##              condition = "==",
##              value = "",
##              .valueIsCharacter = TRUE
##          )
##          )
## ## filter for gene_biotype
## setClass("TxbiotypeFilter", contains="AnnotationFilter",
##          prototype=list(
##              condition="==",
##              value="",
##              .valueIsCharacter=TRUE
##             )
##         )

## Table exon
## filter for exon_id
## setClass("ExonidFilter", contains="AnnotationFilter",
##          prototype=list(
##              condition="==",
##              value="",
##              .valueIsCharacter=TRUE
##             )
##         )

## Table tx2exon
## filter for exon_idx
## setClass("ExonrankFilter", contains="AnnotationFilter",
##          prototype=list(
##              condition="==",
##              value=integer(),
##              .valueIsCharacter=FALSE
##             )
##         )


## chromosome positions
## basic chromosome/seqname filter.
## setClass("SeqnameFilter", contains="AnnotationFilter",
##          prototype=list(
##              condition="==",
##              value="",
##              .valueIsCharacter=TRUE
##             )
##         )
## builder...

## ## basic chromosome strand filter.
## setClass("SeqstrandFilter", contains="AnnotationFilter",
##          prototype=list(
##              condition="==",
##              value=1L,
##              .valueIsCharacter=FALSE
##             )
##         )
## builder...

## ## chromstart filter
## setClass("SeqstartFilter", contains="AnnotationFilter",
##          representation(
##              feature="character"
##             ),
##          prototype=list(
##              condition=">",
##              value=0L,
##              .valueIsCharacter=FALSE,
##              feature="gene"
##             )
##         )

## ## chromend filter
## setClass("SeqendFilter", contains="AnnotationFilter",
##          representation(
##              feature="character"
##             ),
##          prototype=list(
##              condition="<",
##              value=0L,
##              .valueIsCharacter=FALSE,
##              feature="gene"
##             )
##         )


###============================================================
##  GRangesFilter
##  adding new arguments since we can not overwrite the data type
##  of the BasicFilter class... unfortunately.
##  + grange <- value
##  + location <- condition
###------------------------------------------------------------
setClass("GRangesFilter", contains="AnnotationFilter",
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
## setClass("SymbolFilter", contains = "BasicFilter",
##          prototype = list(
##              condition = "=",
##              value = "",
##              .valueIsCharacter = TRUE
##          )
##          )
## SymbolFilter <- function(value, condition = "=") {
##     if(missing(value)){
##         stop("A filter without a value makes no sense!")
##     }
##     if(length(value) > 1) {
##         if(condition == "=")
##             condition = "in"
##         if(condition == "!=")
##             condition = "not in"
##     }
##     return(new("SymbolFilter", condition = condition,
##                value = as.character(value)))
## }

############################################################
## OnlyCodingTxFilter
##
## That's a special case filter that just returns transcripts
## that have tx_cds_seq_start defined (i.e. not NULL).
setClass("OnlyCodingTxFilter", contains = "AnnotationFilter",
         prototype = list(
             condition = "=",
             value = character(),
             .valueIsCharacter = TRUE
         ))
OnlyCodingTxFilter <- function() {
    return(new("OnlyCodingTxFilter"))
}

#' Protein annotation related filters
#'
#' @description These filters allow to query specific information from an
#' \code{\linkS4class{EnsDb}} database. Filters should be created using the
#' dedicated constructor functions \code{ProteinidFilter},
#' \code{ProtdomidFilter}, \code{UniprotidFilter}, \code{UniprotdbFilter} or
#' \code{UniprotmappingtypeFilter}. Protein annotation-based
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
#' \code{\link{listUniprotDbs}} and \code{\link{listUniprotMappingTypes}} to
#' list all Uniprot database names respectively mapping method types from the
#' database.
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
#'
NULL
#> NULL

## ############################################################
## ## ProteinidFilter
## ##' @description The \code{ProteinidFilter} allows to filter based on the
## ##' (Ensembl) protein ID of the (coding) transcripts' proteins.
## ##' @rdname ProteinFilters
## setClass("ProteinidFilter", contains = "AnnotationFilter",
##          prototype = list(
##              condition = "==",
##              value = "",
##              .valueIsCharacter = TRUE
##          ))

## ##' @param value A character vector of length 1 or larger with the value(s)
## ##' for the filter.
## ##' @param condition A character of length 1 specifying the condition of the
## ##' filter. Can be one of \code{"="}, \code{"!="}, \code{"like"}, or, for
## ##' \code{value} of length larger 1 \code{"in"} and \code{"not in"}.
## ##' @return For \code{ProteinidFilter}: A \code{ProteinidFilter} object.
## ##' @rdname ProteinFilters
## ProteinidFilter <- function(value, condition = "==") {
##     .Deprecated("ProteinIdFilter")
##     if (missing(value))
##         stop("A filter without a value makes no sense!")
##     return(new("ProteinIdFilter", condition = condition,
##                value = as.character(value)))
## }

## ############################################################
## ## UniprotFilter
## ##' @description The \code{UniprotidFilter} allows to retrieve annotations
## ##' filtering by the provided Uniprot ID associated with the transcript's
## ##' protein.
## ##' @rdname ProteinFilters
## setClass("UniprotidFilter", contains = "AnnotationFilter",
##          prototype = list(
##              condition = "==",
##              value = "",
##              .valueIsCharacter = TRUE
##          ))
## ##' @return For \code{UniprotidFilter}: A \code{UniprotidFilter} object.
## ##' @rdname ProteinFilters
## UniprotidFilter <- function(value, condition = "==") {
##     .Deprecated("UniprotFilter")
##     if (missing(value))
##         stop("A filter without a value makes no sense!")
##     return(new("UniprotFilter", condition = condition,
##                value = as.character(value)))
## }

############################################################
## ProtdomidFilter
##' @description The \code{ProtDomIdFilter} allows to retrieve entries from
##' the database matching the provided filter criteria based on their protein
##' domain ID (\emph{protein_domain_id}).
##' @rdname ProteinFilters
setClass("ProtDomIdFilter", contains = "AnnotationFilter",
         prototype = list(
             condition = "==",
             value = "",
             .valueIsCharacter = TRUE
         ))
##' @return For \code{ProtDomIdFilter}: A \code{ProtDomIdFilter} object.
##' @rdname ProteinFilters
ProtDomIdFilter <- function(value, condition = "==") {
    return(new("ProtDomIdFilter", condition = condition,
               value = as.character(value)))
}

############################################################
## UniprotdbFilter
##' @description The \code{UniprotDbFilter} allows to filter results based on
##' the specified Uniprot database names.
##' @rdname ProteinFilters
setClass("UniprotDbFilter", contains = "AnnotationFilter",
         prototype = list(
             condition = "==",
             values = "",
             .valueIsCharacter = TRUE
         ))
##' @return For \code{UniprotDbFilter}: A \code{UniprotDbFilter} object.
##' @rdname ProteinFilters
UniprotDbFilter <- function(value, condition = "==") {
    return(new("UniprotDbFilter", condition = condition,
               value = as.character(value)))
}

############################################################
## UniprotmappingtypeFilter
##' @description The \code{UniprotMappingTypeFilter} allows to filter results
##' based on the mapping method/type that was used to assign Uniprot IDs to
##' Ensembl protein IDs.
##' @rdname ProteinFilters
setClass("UniprotMappingTypeFilter", contains = "AnnotationFilter",
         prototype = list(
             condition = "==",
             values = "",
             .valueIsCharacter = TRUE
         ))
##' @return For \code{UniprotMappingTypeFilter}: A
##' \code{UniprotMappingTypeFilter} object.
##' @rdname ProteinFilters
UniprotMappingTypeFilter <- function(value, condition = "==") {
    return(new("UniprotMappingTypeFilter", condition = condition,
               value = as.character(value)))
}

