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

