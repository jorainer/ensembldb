## Deprecated functions.

#' @aliases ensembldb-deprecated
#' 
#' @title Deprecated functionality
#'
#' @description All functions, methods and classes listed on this page are
#' deprecated and might be removed in future releases.
#'
#' @param value The value for the filter.
#' @param condition The condition for the filter.
#' 
#' @name Deprecated 
NULL
#> NULL

#' @description \code{GeneidFilter} creates a \code{GeneIdFilter}. Use
#' \code{GeneIdFilter} from the \code{AnnotationFilter} package instead.
#' 
#' @rdname Deprecated
GeneidFilter <- function(value, condition = "==") {
    .Deprecated("GeneIdFilter")
    if (missing(value))
        stop("A filter without a value makes no sense!")
    ## if(length(value) > 1){
    ##     if(condition=="=")
    ##         condition="in"
    ##     if(condition=="!=")
    ##         condition="not in"
    ## }
    return(new("GeneIdFilter", condition = condition,
               value = as.character(value), field = "gene_id"))
}

#' @description \code{GenebiotypeFilter} creates a \code{GeneBiotypeFilter}. Use
#' \code{GeneBiotypeFilter} from the \code{AnnotationFilter} package instead.
#' 
#' @rdname Deprecated
GenebiotypeFilter <- function(value, condition = "=="){
    .Deprecated("GeneBiotypeFilter")
    if(missing(value))
        stop("A filter without a value makes no sense!")
    return(new("GeneBiotypeFilter", condition = condition,
               value = as.character(value), field = "gene_biotype"))
}

#' @description \code{EntrezidFilter} creates a \code{EntrezFilter}. Use
#' \code{EntrezFilter} from the \code{AnnotationFilter} package instead.
#' 
#' @rdname Deprecated
EntrezidFilter <- function(value, condition = "=="){
    .Deprecated("EntrezFilter")
    if(missing(value))
        stop("A filter without a value makes no sense!")
    return(new("EntrezFilter", condition = condition,
               value = as.character(value), field = "entrez"))
}

#' @description \code{TxidFilter} creates a \code{TxIdFilter}. Use
#' \code{TxIdFilter} from the \code{AnnotationFilter} package instead.
#' 
#' @rdname Deprecated
TxidFilter <- function(value, condition = "=="){
    .Deprecated("TxIdFilter")
    if(missing(value))
        stop("A filter without a value makes no sense!")
    return(new("TxIdFilter", condition = condition,
               value = as.character(value), field = "tx_id"))
}

#' @description \code{TxbiotypeFilter} creates a \code{TxBiotypeFilter}. Use
#' \code{TxBiotypeFilter} from the \code{AnnotationFilter} package instead.
#' 
#' @rdname Deprecated
TxbiotypeFilter <- function(value, condition="=="){
    .Deprecated("TxBiotypeFilter")
    if(missing(value))
        stop("A filter without a value makes no sense!")
    return(new("TxBiotypeFilter", condition=condition,
               value=as.character(value), field = "tx_biotype"))
}

#' @description \code{ExonidFilter} creates a \code{ExonIdFilter}. Use
#' \code{ExonIdFilter} from the \code{AnnotationFilter} package instead.
#' 
#' @rdname Deprecated
ExonidFilter <- function(value, condition="=="){
    .Deprecated("ExonIdFilter")
    if(missing(value))
        stop("A filter without a value makes no sense!")
    return(new("ExonIdFilter", condition=condition,
               value=as.character(value), field = "exon_id"))
}

#' @description \code{ExonrankFilter} creates a \code{ExonRankFilter}. Use
#' \code{ExonRankFilter} from the \code{AnnotationFilter} package instead.
#' 
#' @rdname Deprecated
ExonrankFilter <- function(value, condition="=="){
    .Deprecated("ExonRankFilter")
    if(missing(value))
        stop("A filter without a value makes no sense!")
    return(new("ExonRankFilter", condition=condition, value=as.integer(value),
               field = "exon_rank"))
}

#' @description \code{SeqNameFilter} creates a \code{SeqNameFilter}. Use
#' \code{SeqNameFilter} from the \code{AnnotationFilter} package instead.
#' 
#' @rdname Deprecated
SeqnameFilter <- function(value, condition="=="){
    .Deprecated("SeqNameFilter")
    if(missing(value))
        stop("A filter without a value makes no sense!")
    return(new("SeqNameFilter", condition=condition, value=as.character(value),
               field = "seq_name"))
}

#' @description \code{SeqstrandFilter} creates a \code{SeqStrandFilter}. Use
#' \code{SeqStrandFilter} from the \code{AnnotationFilter} instead.
#' 
#' @rdname Deprecated
SeqstrandFilter <- function(value, condition="=="){
    .Deprecated("SeqStrandFilter")
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
    return(new("SeqStrandFilter", condition=condition, value=as.character(value),
               field = "seq_strand"))
}

#' @description
#'
#' \code{SeqstartFilter} creates a \code{GeneStartFilter}, \code{TxStartFilter}
#' or \code{ExonStartFilter} depending on the value of the parameter
#' \code{feature}. Use \code{GeneStartFilter}, \code{TxStartFilter} and
#' \code{ExonStartFilter} instead.
#'
#' @param feature For \code{SeqstartFilter} and \code{SeqendFilter}: on what type
#' of feature should the filter be applied? Supported are \code{"gene"},
#' \code{"tx"} and \code{"exon"}.
#' 
#' @rdname Deprecated
SeqstartFilter <- function(value, condition=">", feature="gene"){
    .Deprecated(msg = paste0("The use of 'SeqstartFilter' is deprecated. Use ",
                             "one of 'GeneStartFilter', 'ExonStartFilter'",
                             " or 'TxStartFilter' instead."))
    feature <- match.arg(feature, c("gene", "exon", "tx"))
    return(new(paste0(sub("^([[:alpha:]])", "\\U\\1", feature, perl=TRUE),
                      "StartFilter"), value = as.integer(value),
               condition = condition,
               field = paste0(feature, "_start")))
}

#' @description
#'
#' \code{SeqendFilter} creates a \code{GeneEndFilter}, \code{TxEndFilter}
#' or \code{ExonEndFilter} depending on the value of the parameter
#' \code{feature}. Use \code{GeneEndFilter}, \code{TxEndFilter} and
#' \code{ExonEndFilter} instead.
#'
#' @rdname Deprecated
SeqendFilter <- function(value, condition="<", feature="gene"){
    .Deprecated(msg = paste0("The use of 'SeqendFilter' is deprecated. Use ",
                             "one of 'GeneEndFilter', 'ExonEndFilter'",
                             " or 'TxEndFilter' instead."))
    feature <- match.arg(feature, c("gene", "exon", "tx"))
    return(new(paste0(sub("^([[:alpha:]])", "\\U\\1", feature, perl=TRUE)),
               "EndFilter"), value = as.integer(value),
           condition = condition,
           field = paste0(feature, "_end"))
}
