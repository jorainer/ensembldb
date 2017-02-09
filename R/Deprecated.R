## Deprecated functions.


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
               value = as.character(value)))
}

GenebiotypeFilter <- function(value, condition = "=="){
    .Deprecated("GeneBiotypeFilter")
    if(missing(value))
        stop("A filter without a value makes no sense!")
    return(new("GeneBiotypeFilter", condition = condition,
               value = as.character(value)))
}

EntrezidFilter <- function(value, condition = "=="){
    .Deprecated("EntrezFilter")
    if(missing(value))
        stop("A filter without a value makes no sense!")
    return(new("EntrezFilter", condition = condition,
               value = as.character(value)))
}

TxidFilter <- function(value, condition = "=="){
    .Deprecated("TxIdFilter")
    if(missing(value))
        stop("A filter without a value makes no sense!")
    return(new("TxIdFilter", condition = condition,
               value = as.character(value)))
}

TxbiotypeFilter <- function(value, condition="=="){
    .Deprecated("TxBiotypeFilter")
    if(missing(value))
        stop("A filter without a value makes no sense!")
    return(new("TxBiotypeFilter", condition=condition,
               value=as.character(value)))
}

ExonidFilter <- function(value, condition="=="){
    .Deprecated("ExonIdFilter")
    if(missing(value))
        stop("A filter without a value makes no sense!")
    return(new("ExonIdFilter", condition=condition,
               value=as.character(value)))
}

ExonrankFilter <- function(value, condition="=="){
    .Deprecated("ExonRankFilter")
    if(missing(value))
        stop("A filter without a value makes no sense!")
    return(new("ExonRankFilter", condition=condition, value=as.integer(value)))
}

SeqnameFilter <- function(value, condition="=="){
    .Deprecated("SeqNameFilter")
    if(missing(value))
        stop("A filter without a value makes no sense!")
    return(new("SeqNameFilter", condition=condition, value=as.character(value)))
}

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
    return(new("SeqStrandFilter", condition=condition, value=as.character(value)))
}

SeqstartFilter <- function(value, condition=">", feature="gene"){
    .Deprecated(msg = paste0("The use of 'SeqstartFilter' is deprecated. Use ",
                             "one of 'GeneStartFilter', 'ExonStartFilter'",
                             " or 'TxStartFilter' instead."))
    feature <- match.arg(feature, c("gene", "exon", "tx"))
    return(new(paste0(sub("^([[:alpha:]])", "\\U\\1", feature, perl=TRUE)),
               "StartFilter", value = as.integer(value),
               condition = "condition"))
    ## if(missing(value)){
    ##     stop("A filter without a value makes no sense!")
    ## }
    ## if(length(value) > 1){
    ##     value <- value[ 1 ]
    ##     warning("Multiple values provided, but only the first (", value,
    ##             ") will be considered")
    ## }
    ## return(new("SeqstartFilter", condition=condition, value=as.integer(value),
    ##             feature=feature))
}

SeqendFilter <- function(value, condition="<", feature="gene"){
    .Deprecated(msg = paste0("The use of 'SeqendFilter' is deprecated. Use ",
                             "one of 'GeneEndFilter', 'ExonEndFilter'",
                             " or 'TxEndFilter' instead."))
    feature <- match.arg(feature, c("gene", "exon", "tx"))
    return(new(paste0(sub("^([[:alpha:]])", "\\U\\1", feature, perl=TRUE)),
               "StartFilter", value = as.integer(value),
               condition = "condition"))
    ## if(missing(value)){
    ##     stop("A filter without a value makes no sense!")
    ## }
    ## if(length(value) > 1){
    ##     value <- value[ 1 ]
    ##     warning("Multiple values provided, but only the first (", value,") will be considered")
    ## }
    ## return(new("SeqendFilter", condition=condition, value=as.integer(value),
    ##             feature=feature))
}
