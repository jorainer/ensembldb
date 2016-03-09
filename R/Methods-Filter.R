##***********************************************************************
##
##     Methods for BasicFilter classes.
##
##***********************************************************************
validateConditionFilter <- function(object){
    if(object@.valueIsCharacter){
        ## condition has to be either = or in
        if(!any(c("=", "in", "not in", "like", "!=")==object@condition)){
            return(paste("only \"=\", \"!=\", \"in\" , \"not in\" and \"like\"",
                         "allowed for condition",
                         ", I've got", object@condition))
        }
    }else{
        ## condition has to be = < > >= <=
        if(!any(c("=", ">", "<", ">=", "<=", "in", "not in")==object@condition)){
            return(paste("only \"=\", \">\", \"<\", \">=\", \"<=\" , \"in\" and \"not in\"",
                         " are allowed for condition, I've got", object@condition))
        }
    }
    if(length(object@value) > 1){
        if(any(!object@condition %in% c("in", "not in")))
            return(paste("only \"in\" and \"not in\" are allowed if value",
                         "is a vector with more than one value!"))
    }
    if(!object@.valueIsCharacter){
        ## value has to be numeric!!!
        if(object@value!=""){
            suppressWarnings(
                if(any(is.na(is.numeric(object@value))))
                    return(paste("value has to be numeric!!!"))
               )
        }
    }
    return(TRUE)
}
setValidity("BasicFilter", validateConditionFilter)
setMethod("initialize", "BasicFilter", function(.Object, ...){
    OK <- validateConditionFilter(.Object)
    if(class(OK)=="character"){
        stop(OK)
    }
    callNextMethod(.Object, ...)
})

.where <- function(object, db=NULL){
    if(is.null(db)){
        Vals <- value(object)
    }else{
        Vals <- value(object, db)
    }
    ## if not a number we have to single quote!
    if(object@.valueIsCharacter){
        Vals <- sQuote(gsub(unique(Vals),pattern="'",replacement="''"))
    }else{
        Vals <- unique(Vals)
    }
    ## check, if there are more than one, concatenate in that case on put () aroung
    if(length(Vals) > 1){
        Vals <- paste0("(", paste(Vals, collapse=",") ,")")
    }
    return(paste(condition(object), Vals))
}
setMethod("where", signature(object="BasicFilter", db="missing", with.tables="missing"),
          function(object, db, with.tables, ...){
    return(.where(object))
})
setMethod("where", signature(object="BasicFilter", db="EnsDb", with.tables="missing"),
          function(object, db, with.tables, ...){
    return(.where(object, db=db))
})
setMethod("where", signature(object="BasicFilter", db="EnsDb", with.tables="character"),
          function(object, db, with.tables, ...){
    return(.where(object, db=db))
})
setMethod("condition", "BasicFilter", function(x, ...){
    if(length(unique(value(x))) > 1){
        if(x@condition=="in" | x@condition=="not in")
            return(x@condition)
        if(x@condition=="!="){
            return("not in")
        }else if(x@condition=="="){
            return("in")
        }else{
            stop("With more than 1 value only conditions \"=\" and \"!=\" are allowed!")
        }
    }else{
        ## check first if we do have "in" or "not in" and if
        ## cast it to a = and != respectively
        if(x@condition=="in")
            return("=")
        if(x@condition=="not in")
            return("!=")
        return(x@condition)
    }
})
setReplaceMethod("condition", "BasicFilter", function(x, value){
    if(x@.valueIsCharacter){
        allowed <- c("=", "!=", "in", "not in", "like")
        if(!any(allowed == value)){
            stop("Only ", paste(allowed, collapse=", "), " are allowed if the value from",
                 " the filter is of type character!")
        }
        if(value == "=" & length(x@value) > 1)
            value <- "in"
        if(value == "!=" & length(x@value) > 1)
            value <- "not in"
        if(value == "in" & length(x@value) == 1)
            value <- "="
        if(value == "not in" & length(x@value) == 1)
            value <- "!="
    }else{
        allowed <- c("=", ">", "<", ">=", "<=")
        if(!any(allowed == value)){
            stop("Only ", paste(allowed, collapse=", "), " are allowed if the value from",
                 " the filter is numeric!")
        }
    }
    x@condition <- value
    validObject(x)
    return(x)
})
setMethod("value", signature(x="BasicFilter", db="missing"),
          function(x, db, ...){
              return(x@value)
          })
setMethod("value", signature(x="BasicFilter", db="EnsDb"),
          function(x, db, ...){
              return(x@value)
          })
setReplaceMethod("value", "BasicFilter", function(x, value){
    if(is.numeric(value)){
        x@.valueIsCharacter <- FALSE
    }else{
        x@.valueIsCharacter <- TRUE
    }
    x@value <- as.character(value)
    ## Checking if condition matches the value.
    if(length(value) > 1){
        if(x@condition == "=")
            x@condition <- "in"
        if(x@condition == "!=")
            x@condition <- "not in"
    }else{
        if(x@condition == "in")
            x@condition <- "="
        if(x@condition == "not in")
            x@condition <- "!="
    }
    ## Test validity
    validObject(x)
    return(x)
})
## setMethod("requireTable", "EnsFilter", function(object, ...){
##     return(object@required.table)
## })
setMethod("print", "BasicFilter", function(x, ...){
    show(x)
})
setMethod("show", "BasicFilter", function(object){
    cat("| ", class(object), "\n")
    cat("| condition: ", object@condition, "\n")
    cat("| value: ", value(object), "\n")
})

##***********************************************************************
##
##     where for a list.
##
##***********************************************************************
setMethod("where", signature(object="list",db="missing", with.tables="missing"),
          function(object, db, with.tables, ...){
              wherequery <- paste(" where", paste(unlist(lapply(object, where)),
                                                  collapse=" and "))
              return(wherequery)
          })
setMethod("where", signature(object="list",db="EnsDb", with.tables="missing"),
          function(object, db, with.tables, ...){
              wherequery <- paste(" where", paste(unlist(lapply(object, where, db)),
                                                  collapse=" and "))
              return(wherequery)
          })
setMethod("where", signature(object="list",db="EnsDb", with.tables="character"),
          function(object, db, with.tables, ...){
              wherequery <- paste(" where", paste(unlist(lapply(object, where, db,
                                                                with.tables=with.tables)),
                                                  collapse=" and "))
              return(wherequery)
          })



##***********************************************************************
##
##     Methods for GeneidFilter classes.
##
##***********************************************************************
setMethod("where", signature(object="GeneidFilter", db="missing", with.tables="missing"),
          function(object, db, with.tables, ...){
              suff <- callNextMethod()
              return(paste(column(object), suff))
          })
setMethod("column", signature(object="GeneidFilter", db="missing", with.tables="missing"),
          function(object, db, with.tables, ...){
              return("gene_id")
          })
setMethod("where", signature(object="GeneidFilter", db="EnsDb", with.tables="missing"),
          function(object, db, with.tables, ...){
              tn <- names(listTables(db))
              return(where(object, db, with.tables=tn))
          })
setMethod("column", signature(object="GeneidFilter", db="EnsDb", with.tables="missing"),
          function(object, db, with.tables, ...){
              tn <- names(listTables(db))
              return(column(object, db, with.tables=tn))
          })
setMethod("where", signature(object="GeneidFilter", db="EnsDb", with.tables="character"),
          function(object, db, with.tables, ...){
              suff <- callNextMethod()
              return(paste(column(object, db, with.tables), suff))
          })
setMethod("column", signature("GeneidFilter", db="EnsDb", with.tables="character"),
          function(object, db, with.tables, ...){
              return(unlist(prefixColumns(db, column(object), with.tables=with.tables),
                            use.names=FALSE))
          })



##***********************************************************************
##
##     Methods for EntrezidFilter classes.
##
##***********************************************************************
setMethod("where", signature(object="EntrezidFilter", db="missing", with.tables="missing"),
          function(object, db, with.tables, ...){
              suff <- callNextMethod()
              return(paste(column(object), suff))
          })
setMethod("column", signature(object="EntrezidFilter", db="missing", with.tables="missing"),
          function(object, db, with.tables, ...){
              return("entrezid")
          })
setMethod("where", signature(object="EntrezidFilter", db="EnsDb", with.tables="missing"),
          function(object, db, with.tables, ...){
              tn <- names(listTables(db))
              return(where(object, db, with.tables=tn))
          })
setMethod("column", signature(object="EntrezidFilter", db="EnsDb", with.tables="missing"),
          function(object, db, with.tables, ...){
              tn <- names(listTables(db))
              return(column(object, db, with.tables=tn))
          })
setMethod("where", signature(object="EntrezidFilter", db="EnsDb", with.tables="character"),
          function(object, db, with.tables, ...){
              suff <- callNextMethod()
              return(paste(column(object, db, with.tables=with.tables), suff))
          })
setMethod("column", signature("EntrezidFilter", db="EnsDb", with.tables="character"),
          function(object, db, with.tables, ...){
              return(unlist(prefixColumns(db, column(object), with.tables=with.tables),
                            use.names=FALSE))
          })


##***********************************************************************
##
##     Methods for GenebiotypeFilter classes.
##
##***********************************************************************
setMethod("where", signature(object="GenebiotypeFilter", db="missing", with.tables="missing"),
          function(object, db, with.tables, ...){
              suff <- callNextMethod()
              return(paste(column(object), suff))
          })
setMethod("column", signature(object="GenebiotypeFilter", db="missing", with.tables="missing"),
          function(object, db, with.tables, ...){
              return("gene_biotype")
          })
setMethod("where", signature(object="GenebiotypeFilter", db="EnsDb", with.tables="missing"),
          function(object, db, with.tables, ...){
              tn <- names(listTables(db))
              return(where(object, db, with.tables=tn))
          })
setMethod("column", signature(object="GenebiotypeFilter", db="EnsDb", with.tables="missing"),
          function(object, db, with.tables, ...){
              tn <- names(listTables(db))
              return(column(object, db, with.tables=tn))
          })
setMethod("where", signature(object="GenebiotypeFilter", db="EnsDb", with.tables="character"),
          function(object, db, with.tables, ...){
              suff <- callNextMethod()
              return(paste(column(object, db, with.tables=with.tables), suff))
          })
setMethod("column", signature(object="GenebiotypeFilter", db="EnsDb", with.tables="character"),
          function(object, db, with.tables, ...){
              return(unlist(prefixColumns(db, column(object), with.tables=with.tables),
                            use.names=FALSE))
          })



##***********************************************************************
##
##     Methods for GenenameFilter classes.
##
##***********************************************************************
setMethod("where", signature(object="GenenameFilter", db="missing", with.tables="missing"),
          function(object, db, with.tables, ...){
              suff <- callNextMethod()
              return(paste(column(object), suff))
          })
setMethod("column", signature(object="GenenameFilter", db="missing", with.tables="missing"),
          function(object, db, with.tables, ...){
              return("gene_name")
          })
setMethod("where", signature(object="GenenameFilter", db="EnsDb", with.tables="missing"),
          function(object, db, with.tables, ...){
              tn <- names(listTables(db))
              return(where(object, db, with.tables=tn))
          })
setMethod("column", signature(object="GenenameFilter", db="EnsDb", with.tables="missing"),
          function(object, db, with.tables, ...){
              tn <- names(listTables(db))
              return(column(object, db, with.tables=tn))
          })
setMethod("where", signature(object="GenenameFilter", db="EnsDb", with.tables="character"),
          function(object, db, with.tables="character", ...){
              suff <- callNextMethod()
              return(paste(column(object, db, with.tables=with.tables), suff))
          })
setMethod("column", signature(object="GenenameFilter", db="EnsDb", with.tables="character"),
          function(object, db, with.tables, ...){
              return(unlist(prefixColumns(db, column(object), with.tables=with.tables),
                            use.names=FALSE))
          })




##***********************************************************************
##
##     Methods for TxidFilter classes.
##
##***********************************************************************
setMethod("where", signature(object="TxidFilter", db="missing", with.tables="missing"),
          function(object, db, with.tables, ...){
              suff <- callNextMethod()
              return(paste(column(object), suff))
          })
setMethod("column", signature(object="TxidFilter", db="missing", with.tables="missing"),
          function(object, db, with.tables, ...){
              return("tx_id")
          })
setMethod("where", signature(object="TxidFilter", db="EnsDb", with.tables="missing"),
          function(object, db, with.tables, ...){
              tn <- names(listTables(db))
              return(where(object, db, with.tables=tn))
          })
setMethod("column", signature(object="TxidFilter", db="EnsDb", with.tables="missing"),
          function(object, db, with.tables, ...){
              tn <- names(listTables(db))
              return(column(object, db, with.tables=tn))
          })
setMethod("where", signature(object="TxidFilter", db="EnsDb", with.tables="character"),
          function(object, db, with.tables, ...){
              suff <- callNextMethod()
              return(paste(column(object, db, with.tables=with.tables), suff))
          })
setMethod("column", signature(object="TxidFilter", db="EnsDb", with.tables="character"),
          function(object, db, with.tables, ...){
              return(unlist(prefixColumns(db, column(object), with.tables=with.tables),
                            use.names=FALSE))
          })




##***********************************************************************
##
##     Methods for TxbiotypeFilter classes.
##
##***********************************************************************
setMethod("where", signature(object="TxbiotypeFilter", db="missing", with.tables="missing"),
          function(object, db, with.tables, ...){
              suff <- callNextMethod()
              return(paste(column(object), suff))
          })
setMethod("column", signature(object="TxbiotypeFilter", db="missing", with.tables="missing"),
          function(object, db, with.tables...){
              return("tx_biotype")
          })
setMethod("where", signature(object="TxbiotypeFilter", db="EnsDb", with.tables="missing"),
          function(object, db, with.tables, ...){
              tn <- names(listTables(db))
              return(where(object, db, with.tables=tn))
          })
setMethod("column", signature(object="TxbiotypeFilter", db="EnsDb", with.tables="missing"),
          function(object, db, with.tables, ...){
              tn <- names(listTables(db))
              return(column(object, db, with.tables=tn))
          })
setMethod("where", signature(object="TxbiotypeFilter", db="EnsDb", with.tables="character"),
          function(object, db, with.tables, ...){
              suff <- callNextMethod()
              return(paste(column(object, db, with.tables=with.tables), suff))
          })
setMethod("column", signature(object="TxbiotypeFilter", db="EnsDb", with.tables="character"),
          function(object, db, with.tables, ...){
              return(unlist(prefixColumns(db, column(object), with.tables=with.tables),
                            use.names=FALSE))
          })




##***********************************************************************
##
##     Methods for ExonidFilter classes.
##
##***********************************************************************
setMethod("where", signature(object="ExonidFilter", db="missing", with.tables="missing"),
          function(object, db, with.tables, ...){
              suff <- callNextMethod()
              return(paste(column(object), suff))
          })
setMethod("column", signature(object="ExonidFilter", db="missing", with.tables="missing"),
          function(object, db, with.tables, ...){
              return("exon_id")
          })
setMethod("where", signature(object="ExonidFilter", db="EnsDb", with.tables="missing"),
          function(object, db, with.tables, ...){
              tn <- names(listTables(db))
              return(where(object, db, with.tables=tn))
          })
setMethod("column", signature(object="ExonidFilter", db="EnsDb", with.tables="missing"),
          function(object, db, with.tables, ...){
              tn <- names(listTables(db))
              return(column(object, db, with.tables=tn))
          })
setMethod("where", signature(object="ExonidFilter", db="EnsDb", with.tables="character"),
          function(object, db, with.tables, ...){
              suff <- callNextMethod()
              return(paste(column(object, db, with.tables=with.tables), suff))
          })
setMethod("column", signature(object="ExonidFilter", db="EnsDb", with.tables="character"),
          function(object, db, with.tables, ...){
              return(unlist(prefixColumns(db, column(object), with.tables=with.tables),
                            use.names=FALSE))
          })


##***********************************************************************
##
##     Methods for ExonrankFilter classes.
##
##***********************************************************************
setMethod("where", signature(object="ExonrankFilter", db="missing", with.tables="missing"),
          function(object, db, with.tables, ...){
              suff <- callNextMethod()
              return(paste(column(object), suff))
          })
setMethod("column", signature(object="ExonrankFilter", db="missing", with.tables="missing"),
          function(object, db, with.tables, ...){
              return("exon_idx")
          })
setMethod("where", signature(object="ExonrankFilter", db="EnsDb", with.tables="missing"),
          function(object, db, with.tables, ...){
              tn <- names(listTables(db))
              return(where(object, db, with.tables=tn))
          })
setMethod("column", signature(object="ExonrankFilter", db="EnsDb", with.tables="missing"),
          function(object, db, with.tables, ...){
              tn <- names(listTables(db))
              return(column(object, db, with.tables=tn))
          })
setMethod("where", signature(object="ExonrankFilter", db="EnsDb", with.tables="character"),
          function(object, db, with.tables, ...){
              suff <- callNextMethod()
              return(paste(column(object, db, with.tables=with.tables), suff))
          })
setMethod("column", signature(object="ExonrankFilter", db="EnsDb", with.tables="character"),
          function(object, db, with.tables, ...){
              return(unlist(prefixColumns(db, column(object), with.tables=with.tables),
                            use.names=FALSE))
          })
setReplaceMethod("value", "ExonrankFilter", function(x, value){
    if(any(is.na(as.numeric(value))))
        stop("Argument 'value' has to be numeric!")
    x@value <- value
    validObject(x)
    return(x)
})


##***********************************************************************
##
##     Methods for SeqnameFilter classes.
##
##***********************************************************************
setMethod("where", signature(object="SeqnameFilter", db="missing", with.tables="missing"),
          function(object, db, with.tables, ...){
              suff <- callNextMethod()
              return(paste(column(object), suff))
          })
setMethod("column", signature(object="SeqnameFilter", db="missing", with.tables="missing"),
          function(object, db, with.tables, ...){
              return("seq_name")
          })
setMethod("where", signature(object="SeqnameFilter", db="EnsDb", with.tables="missing"),
          function(object, db, with.tables, ...){
              tn <- names(listTables(db))
              return(where(object, db, with.tables=tn))
          })
setMethod("column", signature(object="SeqnameFilter", db="EnsDb", with.tables="missing"),
          function(object, db, with.tables, ...){
              tn <- names(listTables(db))
              return(column(object, db, with.tables=tn))
          })
setMethod("where", signature(object="SeqnameFilter", db="EnsDb", with.tables="character"),
          function(object, db, with.tables, ...){
              suff <- callNextMethod()
              return(paste(column(object, db, with.tables=with.tables), suff))
          })
setMethod("column", signature(object="SeqnameFilter", db="EnsDb", with.tables="character"),
          function(object, db, with.tables, ...){
              return(unlist(prefixColumns(db, column(object), with.tables=with.tables),
                            use.names=FALSE))
          })
## Overwriting the value method allows us to fix chromosome names (e.g. with prefix chr)
## to be usable for EnsDb and Ensembl based chromosome names (i.e. without chr).
setMethod("value", signature(x="SeqnameFilter", db="EnsDb"),
          function(x, db, ...){
              val <- formatSeqnamesForQuery(db, value(x))
              if(any(is.na(val))){
                  stop("A value of <NA> is not allowed for a SeqnameFilter!")
              }
              return(val)
              ##return(ucscToEns(value(x)))
          })


##***********************************************************************
##
##     Methods for SeqstrandFilter classes.
##
##***********************************************************************
setMethod("where", signature(object="SeqstrandFilter", db="missing", with.tables="missing"),
          function(object, db, with.tables, ...){
              suff <- callNextMethod()
              return(paste(column(object), suff))
          })
setMethod("column", signature(object="SeqstrandFilter", db="missing", with.tables="missing"),
          function(object, db, with.tables, ...){
              return("seq_strand")
          })
setMethod("where", signature(object="SeqstrandFilter", db="EnsDb", with.tables="missing"),
          function(object, db, with.tables, ...){
              tn <- names(listTables(db))
              return(where(object, db, with.tables=tn))
          })
setMethod("column", signature(object="SeqstrandFilter", db="EnsDb", with.tables="missing"),
          function(object, db, with.tables, ...){
              tn <- names(listTables(db))
              return(column(object, db, with.tables=tn))
          })
setMethod("where", signature(object="SeqstrandFilter", db="EnsDb", with.tables="character"),
          function(object, db, with.tables, ...){
              suff <- callNextMethod()
              return(paste(column(object, db, with.tables=with.tables), suff))
          })
setMethod("column", signature(object="SeqstrandFilter", db="EnsDb", with.tables="character"),
          function(object, db, with.tables, ...){
              return(unlist(prefixColumns(db, column(object), with.tables=with.tables),
                            use.names=FALSE))
          })



##***********************************************************************
##
##     Methods for SeqstartFilter classes.
##
##***********************************************************************
setMethod("where", signature(object="SeqstartFilter", db="missing", with.tables="missing"),
          function(object, db, with.tables, ...){
              suff <- callNextMethod()
              return(paste(column(object), suff))
          })
setMethod("column", signature(object="SeqstartFilter", db="missing", with.tables="missing"),
          function(object, db, with.tables, ...){
              ## assuming that we follow the naming convention:
              ## <feature>_seq_end for the naming of the database columns.
              feature <- object@feature
              feature <- match.arg(feature, c("gene", "transcript", "exon", "tx"))
              if(object@feature=="transcript")
                  feature <- "tx"
              return(paste0(feature, "_seq_start"))
          })
setMethod("where", signature(object="SeqstartFilter", db="EnsDb", with.tables="missing"),
          function(object, db, with.tables, ...){
              tn <- names(listTables(db))
              return(where(object, db, with.tables=tn))
          })
setMethod("column", signature(object="SeqstartFilter", db="EnsDb", with.tables="missing"),
          function(object, db, with.tables, ...){
              tn <- names(listTables(db))
              return(column(object, db, with.tables=tn))
          })
setMethod("where", signature(object="SeqstartFilter", db="EnsDb", with.tables="character"),
          function(object, db, with.tables, ...){
              suff <- callNextMethod()
              return(paste(column(object, db, with.tables=with.tables), suff))
          })
setMethod("column", signature(object="SeqstartFilter", db="EnsDb", with.tables="character"),
          function(object, db, with.tables, ...){
              return(unlist(prefixColumns(db, column(object), with.tables=with.tables),
                            use.names=FALSE))
          })



##***********************************************************************
##
##     Methods for SeqendFilter classes.
##
##***********************************************************************
setMethod("where", signature(object="SeqendFilter", db="missing", with.tables="missing"),
          function(object, db, with.tables, ...){
              suff <- callNextMethod()
              return(paste(column(object), suff))
          })
setMethod("column", signature(object="SeqendFilter", db="missing", with.tables="missing"),
          function(object, db, with.tables, ...){
              ## assuming that we follow the naming convention:
              ## <feature>_seq_end for the naming of the database columns.
              feature <- object@feature
              feature <- match.arg(feature, c("gene", "transcript", "exon", "tx"))
              if(object@feature=="transcript")
                  feature <- "tx"
              return(paste0(feature, "_seq_end"))
          })
setMethod("where", signature(object="SeqendFilter", db="EnsDb", with.tables="missing"),
          function(object, db, with.tables, ...){
              tn <- names(listTables(db))
              return(where(object, db, with.tables=tn))
          })
setMethod("column", signature(object="SeqendFilter", db="EnsDb", with.tables="missing"),
          function(object, db, with.tables, ...){
              tn <- names(listTables(db))
              return(column(object, db, with.tables=tn))
          })
setMethod("where", signature(object="SeqendFilter", db="EnsDb", with.tables="character"),
          function(object, db, with.tables, ...){
              suff <- callNextMethod()
              return(paste(column(object, db, with.tables=with.tables), suff))
          })
setMethod("column", signature(object="SeqendFilter", db="EnsDb", with.tables="character"),
          function(object, db, with.tables, ...){
              return(unlist(prefixColumns(db, column(object), with.tables=with.tables),
                            use.names=FALSE))
          })


###============================================================
##    Methods for GRangesFilter
##    + show
##    + condition
##    + value
##    + where
##    + column
##    + start
##    + end
##    + seqnames
##    + strand
###------------------------------------------------------------
## Overwrite the validation method.
setValidity("GRangesFilter", function(object){
    if(!any(object@location == c("within", "overlapping") )){
        return(paste0("Argument condition should be either 'within' or 'overlapping'! Got ",
                      object@location, "!"))
    }
    ## GRanges has to have valid values for start, end and seqnames!
    if(length(start(object)) == 0)
        return("start coordinate of the range is missing!")
    if(length(end(object)) == 0)
        return("end coordinate of the range is missing!")
    if(length(seqnames(object)) == 0)
        return("A valid seqname is required from the submitted GRanges!")
    return(TRUE)
})
setMethod("show", "GRangesFilter", function(object){
    cat("|", class(object), "\n")
    cat("| region:\n")
    cat("| + start: ", paste0(start(object), collapse=", "), "\n")
    cat("| + end:   ", paste0(end(object), collapse=", "), "\n")
    cat("| + seqname: ", paste0(seqnames(object), collapse=", "), "\n")
    cat("| + strand: ", paste0(strand(object), collapse=", "), "\n")
    cat("| condition: ", condition(object), "\n")
})
setMethod("condition", "GRangesFilter", function(x, ...){
    return(x@location)
})
setReplaceMethod("condition", "GRangesFilter", function(x, value){
    value <- match.arg(value, c("within", "overlapping"))
    x@location <- value
    validObject(x)
    return(x)
})
setMethod("value", signature(x="GRangesFilter", db="missing"),
          function(x, db, ...){
              return(x@grange)
          })
setMethod("value", signature(x="GRangesFilter", db="EnsDb"),
          function(x, db, ...){
              return(x@grange)
          })
setMethod("start", signature(x="GRangesFilter"),
          function(x, ...){
              return(start(value(x)))
          })
setMethod("end", signature(x="GRangesFilter"),
          function(x, ...){
              return(end(value(x)))
          })
setMethod("strand", signature(x="GRangesFilter"),
          function(x, ...){
              strnd <- as.character(strand(value(x)))
              return(strnd)
          })
setMethod("seqnames", signature(x="GRangesFilter"),
          function(x){
              return(as.character(seqnames(value(x))))
          })
setMethod("seqlevels", signature(x="GRangesFilter"),
          function(x){
              return(seqlevels(value(x)))
          })
## The column method for GRangesFilter returns all columns required for the query, i.e.
## the _seq_start, _seq_end for the feature, seq_name and seq_strand.
## Note: this method has to return a named vector!
setMethod("column", signature(object="GRangesFilter", db="missing", with.tables="missing"),
          function(object, db, with.tables, ...){
              ## assuming that we follow the naming convention:
              ## <feature>_seq_end for the naming of the database columns.
              feature <- object@feature
              feature <- match.arg(feature, c("gene", "transcript", "exon", "tx"))
              if(object@feature=="transcript")
                  feature <- "tx"
              cols <- c(start=paste0(feature, "_seq_start"),
                        end=paste0(feature, "_seq_end"),
                        seqname="seq_name",
                        strand="seq_strand")
              return(cols)
          })
setMethod("column", signature(object="GRangesFilter", db="EnsDb", with.tables="missing"),
          function(object, db, with.tables, ...){
              tn <- names(listTables(db))
              return(column(object, db, with.tables=tn))
          })
## Providing also the columns.
setMethod("column", signature(object="GRangesFilter", db="EnsDb", with.tables="character"),
          function(object, db, with.tables, ...){
              cols <- unlist(prefixColumns(db, column(object), with.tables=with.tables),
                             use.names=FALSE)
              ## We have to give the vector the required names!
              names(cols) <- 1:length(cols)
              names(cols)[grep(cols, pattern="seq_name")] <- "seqname"
              names(cols)[grep(cols, pattern="seq_strand")] <- "strand"
              names(cols)[grep(cols, pattern="seq_start")] <- "start"
              names(cols)[grep(cols, pattern="seq_end")] <- "end"
              return(cols[c("start", "end", "seqname", "strand")])
          })
## Where for GRangesFilter only.
setMethod("where", signature(object="GRangesFilter", db="missing", with.tables="missing"),
          function(object, db, with.tables, ...){
              ## Get the names of the columns we're going to query.
              cols <- column(object)
              query <- buildWhereForGRanges(object, cols)
              return(query)
          })
setMethod("where", signature(object="GRangesFilter", db="EnsDb", with.tables="missing"),
          function(object, db, with.tables, ...){
              tn <- names(listTables(db))
              return(where(object, db, with.tables=tn))
          })
setMethod("where", signature(object="GRangesFilter", db="EnsDb", with.tables="character"),
          function(object, db, with.tables, ...){
              cols <- column(object, db, with.tables)
              query <- buildWhereForGRanges(object, cols, db=db)
              return(query)
          })


## grf: GRangesFilter
buildWhereForGRanges <- function(grf, columns, db=NULL){
    condition <- condition(grf)
    if(!any(condition == c("within", "overlapping")))
        stop(paste0("'condition' for GRangesFilter should either be ",
                    "'within' or 'overlapping', got ", condition, "."))
    if(is.null(names(columns))){
        stop(paste0("The vector with the required column names for the",
                    " GRangesFilter query has to have names!"))
    }
    if(!all(c("start", "end", "seqname", "strand") %in% names(columns)))
        stop(paste0("'columns' has to be a named vector with names being ",
                    "'start', 'end', 'seqname', 'strand'!"))
    ## Build the query to fetch all features that are located within the range
    quers <- sapply(value(grf), function(z){
        if(!is.null(db)){
            seqn <- formatSeqnamesForQuery(db, as.character(seqnames(z)))
        }else{
            seqn <- as.character(seqnames(z))
        }
        if(condition == "within"){
            query <- paste0(columns["start"], " >= ", start(z), " and ",
                            columns["end"], " <= ", end(z), " and ",
                            columns["seqname"], " == '", seqn, "'")
        }
        ## Build the query to fetch all features (partially) overlapping the range. This
        ## includes also all features (genes or transcripts) that have an intron at that
        ## position.
        if(condition == "overlapping"){
            query <- paste0(columns["start"], " <= ", end(z), " and ",
                            columns["end"], " >= ", start(z), " and ",
                            columns["seqname"], " = '", seqn, "'")
        }
        ## Include the strand, if it's not "*"
        if(as.character(strand(z)) != "*"){
            query <- paste0(query, " and ", columns["strand"], " = ",
                            strand2num(as.character(strand(z))))
        }
        return(query)
    })
    if(length(quers) > 1)
        quers <- paste0("(", quers, ")")
    query <- paste0(quers, collapse=" or ")
    ## Collapse now the queries.
    return(query)
}



## map chromosome strand...
strand2num <- function(x){
    if(x == "+" | x == "-"){
        return(as.numeric(paste0(x, 1)))
    }else{
        stop("Only '+' and '-' supported!")
    }
}
num2strand <- function(x){
    if(x < 0){
        return("-")
    }else{
        return("+")
    }
}

