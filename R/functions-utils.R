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
addFilterColumns <- function(cols, filter = AnnotationFilterList(), edb) {
    if (missing(cols))
        cols <- NULL
    gimmeAll <- returnFilterColumns(edb)
    if (!gimmeAll)
        return(cols)
    ## Put filter into an AnnotationFilterList if it's not already one
    if (is(filter, "AnnotationFilter"))
        filter <- AnnotationFilterList(filter)
    ## Or alternatively process the filters and add columns.
    symFilts <- c("SymbolFilter")
    addC <- unlist(lapply(filter, function(z) {
        if(class(z) %in% symFilts)
            return(z@field)
        if (is(z, "AnnotationFilterList"))
            return(addFilterColumns(cols = cols, filter = z, edb))
        return(ensDbColumn(z))
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
        shiny::runApp(appDir=system.file("shinyHappyPeople",
                                         package="ensembldb"), ...)
    }else{
        stop("Package shiny not installed!")
    }
}

############################################################
## anyProteinColumns
##
## Check if any of 'x' are protein columns.
anyProteinColumns <- function(x){
    return(any(x %in% unlist(.ensdb_protein_tables(), use.names = FALSE)))
}

############################################################
## listProteinColumns
##
#' @description The \code{listProteinColumns} function allows to conveniently
#'     extract all database columns containing protein annotations from
#'     an \code{\linkS4class{EnsDb}} database.
#' 
#' @return The \code{listProteinColumns} function returns a character vector
#'     with the column names containing protein annotations or throws an error
#'     if no such annotations are available.
#' 
#' @rdname ProteinFunctionality
#' 
#' @examples
#'
#' ## List all columns containing protein annotations
#' library(EnsDb.Hsapiens.v75)
#' edb <- EnsDb.Hsapiens.v75
#' if (hasProteinData(edb))
#'     listProteinColumns(edb)
listProteinColumns <- function(object) {
    if (missing(object))
        stop("'object' is missing with no default.")
    if (!is(object, "EnsDb"))
        stop("'object' has to be an instance of an 'EnsDb' object.")
    if (!hasProteinData(object))
        stop("The provided EnsDb database does not contain protein annotations!")
    return(listColumns(object, c("protein", "uniprot", "protein_domain")))
}

############################################################
## .ProteinsFromDataframe
#' @param x \code{EnsDb} object.
#' 
#' @param data \code{data.frame} with the results from a call to the
#'     \code{proteins} method; has to have required columns \code{"protein_id"}
#'     and \code{"protein_sequence"}.
#' 
#' @noRd
.ProteinsFromDataframe <- function(x, data) {
    if (!all(c("protein_id", "protein_sequence") %in% colnames(data)))
        stop("Reguired columns 'protein_id' and 'protein_sequence' not in 'data'!")
    ## Get the column names for uniprot and protein_domain
    uniprot_cols <- listColumns(x, "uniprot")
    uniprot_cols <- uniprot_cols[uniprot_cols != "protein_id"]
    uniprot_cols <- uniprot_cols[uniprot_cols %in% colnames(data)]
    if (length(uniprot_cols) > 0)
        warning("Don't know yet how to handle the 1:n mapping between",
                " protein_id and uniprot_id!")

    prot_dom_cols <- listColumns(x, "protein_domain")
    prot_dom_cols <- prot_dom_cols[prot_dom_cols != "protein_id"]
    prot_dom_cols <- prot_dom_cols[prot_dom_cols %in% colnames(data)]

    ## Create the protein part of the object, i.e. the AAStringSet.
    ## Use all columns other than protein_id, protein_sequence
    prot_cols <- colnames(data)
    prot_cols <- prot_cols[!(prot_cols %in% c(uniprot_cols, prot_dom_cols))]
    protein_sub <- unique(data[, prot_cols, drop = FALSE])
    aass <- AAStringSet(protein_sub$protein_sequence)
    names(aass) <- protein_sub$protein_id
    prot_cols <- prot_cols[!(prot_cols %in% c("protein_id", "protein_sequence"))]
    if (length(prot_cols) > 0) {
        mcols(aass) <- DataFrame(protein_sub[, prot_cols, drop = FALSE])
        ## drop these columns from data to eventually speed up splits
        data <- data[, !(colnames(data) %in% prot_cols), drop = FALSE]
    }

    ## How to process the Uniprot here??? have a 1:n mapping!

    ## Create the protein domain part
    if (length(prot_dom_cols) > 0) {
        message("Processing protein domains not yet implemented!")
        ## Split the dataframe by protein_id
        ## process this list to create the IRangesList.
        ## pranges should have the same order and the same names
    } else {
        pranges <- IRangesList(replicate(length(aass), IRanges()))
        names(pranges) <- names(aass)
    }
    metadata <- list(created = date())

    ##return(new("Proteins", aa = aass, pranges = pranges, metadata = metadata))
}

## map chromosome strand...
strand2num <- function(x){
    if (is.numeric(x)) {
        if (x >= 0) return(1)
        else return(-1)
    }
    xm <- x
    if(xm == "+" | xm == "-")
        xm <- paste0(xm, 1)
    xm <- as.numeric(xm)
    if (is.na(xm))
        stop("'", x, "' can not be converted to a strand!")
    return(xm)
}

num2strand <- function(x){
    if(x < 0){
        return("-")
    }else{
        return("+")
    }
}

#' @description Collapses entries in the \code{"entrezid"} column of a
#'     \code{data.frame} or \code{DataFrame} making the rest of \code{x} unique.
#'
#' @param x Either a \code{data.frame} or a \code{DataFrame}.
#'
#' @param by \code{character(1)} defining the column by which the
#'     \code{"entrezid"} column should be splitted.
#' 
#' @author Johannes Rainer
#' 
#' @noRd
.collapseEntrezidInTable <- function(x, by = "gene_id") {
    ## Slow version: use unique call.
    eg_idx <- which(colnames(x) == "entrezid")
    if (length(eg_idx)) {
        ## Avoid an additional lapply unique call.
        tmp <- unique(x[, c(by, "entrezid")])
        egs <- split(tmp[, "entrezid"],
                     f = factor(tmp[, by], levels = unique(tmp[, by])))
        ## Use a unique call.
        ## x_sub <- x[match(names(egs), x[, by]), , drop = FALSE] would be much
        ## faster but does not work e.g. for exons or transcripts.
        x_sub <- unique(x[, -eg_idx, drop = FALSE])
        x_sub$entrezid <- egs[x_sub[, by]]
        return(x_sub)
    }
    x
}
