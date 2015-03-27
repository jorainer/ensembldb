##***********************************************************************
##
##     Methods for BasicFilter classes.
##
##***********************************************************************
validateConditionFilter <- function(object){
    if(object@.valueIsCharacter){
        ## condition has to be either = or in
        if(!any(c("=", "in", "not in", "like", "!=")==object@condition)){
            return(paste("only \"=\", \"!=\", \"in\" , \"not in\" and \"like\" allowed for condition, I've got", object@condition))
        }
    }else{
        ## condition has to be = < > >= <=
        if(!any(c("=", ">", "<", ">=", "<=")==object@condition)){
            return(paste("only \"=\", \">\", \"<\", \">=\" and \"<=\" are allowed for condition, I've got", object@condition))
        }
    }
    if(length(object@value) > 1){
        if(any(!object@condition %in% c("in", "not in")))
            return(paste("only \"in\" and \"not in\" are allowed if value is a vector with more than one value!"))
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

.where <- function(object){
    ## if not a number we have to single quote!
    if(object@.valueIsCharacter){
        Vals <- sQuote(gsub(unique(object@value),pattern="'",replacement="''"))
    }else{
        Vals <- unique(object@value)
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
    return(.where(object))
})
setMethod("where", signature(object="BasicFilter", db="EnsDb", with.tables="character"),
          function(object, db, with.tables, ...){
    return(.where(object))
})
setMethod("condition", "BasicFilter", function(x, ...){
    if(length(unique(x@value)) > 1){
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

## setMethod("requireTable", "EnsFilter", function(object, ...){
##     return(object@required.table)
## })
setMethod("print", "BasicFilter", function(x, ...){
    show(x)
})
setMethod("show", "BasicFilter", function(object){
    cat("| condition: ", object@condition, "\n")
    cat("| value: ", object@value, "\n")
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

