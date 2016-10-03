############################################################
## Utility functions

############################################################
## orderDataFrameBy
##
## Simply orders the data.frame x based on the columns specified
## with by.
orderDataFrameBy <- function(x, by = "", decreasing = FALSE) {
    if (all(by == "") | all(is.null(by)))
        return(x)
    return(x[do.call(order,
                     args = c(list(method = "radix",
                                   decreasing = decreasing),
                              as.list(x[, by, drop = FALSE]))), ])
}

############################################################
## checkOrderBy
##
## Check the orderBy argument.
## o orderBy can be a character vector or a , separated list.
## o Ensure that the columns are valid by comparing with 'supported'.
## Returns a character vector, each element representing a column
## on which sorting should be performed.
checkOrderBy <- function(orderBy, supported = character()) {
    if (is.null(orderBy) | all(orderBy == "")) {
        return(orderBy)
    }
    if (length(orderBy) == 1 & length(grep(orderBy, pattern = ",")) > 0) {
        orderBy <- unlist(strsplit(orderBy, split = ","), use.names = FALSE)
        orderBy <- gsub(orderBy, pattern = " ", replacement = "", fixed = TRUE)
    }
    not_supported <- !(orderBy %in% supported)
    if (any(not_supported)) {
        warning("Columns in 'order.by' (",
                paste(orderBy[not_supported], collapse = ", "),
                ") are not in 'columns' and were thus removed.")
        orderBy <- orderBy[!not_supported]
        if (length(orderBy) == 0)
            orderBy <- ""
    }
    return(orderBy)
}

############################################################
## addFilterColumns
##
## This function checks the filter objects and adds, depending on the
## returnFilterColumns setting of the EnsDb, also columns for each of the
## filters, ensuring that:
## a) "Symlink" filters are added correctly (the column returned by the
##    column call without db are added).
## b) GRangesFilter: the feature is set based on the specified feature parameter
## Args:
addFilterColumns <- function(cols, filter = list(), edb) {
    if (missing(cols))
        cols <- NULL
    gimmeAll <- returnFilterColumns(edb)
    if (!missing(filter)) {
        if(!is.list(filter))
            filter <- list(filter)
    } else {
        return(cols)
    }
    if (!gimmeAll)
        return(cols)
    ## Or alternatively process the filters and add columns.
    symFilts <- c("SymbolFilter")
    addC <- unlist(lapply(filter, function(z) {
        if(class(z) %in% symFilts)
            return(column(z))
        return(column(z))
    }))
    return(unique(c(cols, addC)))
}

############################################################
## SQLiteName2MySQL
##
## Convert the SQLite database name (file name) to the corresponding
## MySQL database name.
SQLiteName2MySQL <- function(x) {
    return(tolower(gsub(x, pattern = ".", replacement = "_", fixed = TRUE)))
}


## running the shiny web app.
runEnsDbApp <- function(...){
    if(requireNamespace("shiny", quietly=TRUE)){
        message("Starting the EnsDb shiny web app. Use Ctrl-C to stop.")
        shiny::runApp(appDir=system.file("shinyHappyPeople", package="ensembldb"), ...)
    }else{
        stop("Package shiny not installed!")
    }
}




