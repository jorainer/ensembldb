##***********************************************************************
##
##     Generic methods
##
##***********************************************************************
## A

## B
setGeneric("buildQuery", function(x, ...)
    standardGeneric("buildQuery"))

## C
setGeneric("cleanColumns", function(x, columns, ...)
    starndardGeneric("cleanColumns"))
setGeneric("column", function(object, db, with.tables, ...)
    standardGeneric("column"))
setGeneric("condition", function(x, ...)
    standardGeneric("condition"))
setGeneric("condition<-", function(x, value)
    standardGeneric("condition<-"))

## D
setGeneric("dbSeqlevelsStyle", function(x, ...)
    standardGeneric("dbSeqlevelsStyle"))

## E
setGeneric("ensemblVersion", function(x)
    standardGeneric("ensemblVersion"))
## setGeneric("exons", function(x, ...)
##     standardGeneric("exons"))
## setGeneric("exonsBy", function(x, ...)
##     standardGeneric("exonsBy"))

## F
setGeneric("formatSeqnamesForQuery", function(x, sn, ...)
    standardGeneric("formatSeqnamesForQuery"))
setGeneric("formatSeqnamesFromQuery", function(x, sn, ...)
    standardGeneric("formatSeqnamesFromQuery"))

## G
## setGeneric("genes", function(x, ...)
##     standardGeneric("genes"))
setGeneric("getGeneRegionTrackForGviz", function(x, ...)
    standardGeneric("getGeneRegionTrackForGviz"))
setGeneric("getGenomeFaFile", function(x, ...)
    standardGeneric("getGenomeFaFile"))
setGeneric("getGenomeTwoBitFile", function(x, ...)
    standardGeneric("getGenomeTwoBitFile"))
setGeneric("getMetadataValue", function(x, name)
    standardGeneric("getMetadataValue"))
setGeneric("getProperty", function(x, name=NULL, ...)
    standardGeneric("getProperty"))
setGeneric("getWhat", function(x, ...)
    standardGeneric("getWhat"))

## H
setGeneric("hasProteinData", function(x)
    standardGeneric("hasProteinData"))

## I

## J

## K

## L
setGeneric("lengthOf", function(x, ...)
    standardGeneric("lengthOf"))
setGeneric("listColumns", function(x, ...)
    standardGeneric("listColumns"))
setGeneric("listGenebiotypes", function(x, ...)
    standardGeneric("listGenebiotypes"))
setGeneric("listTables", function(x, ...)
    standardGeneric("listTables"))
setGeneric("listTxbiotypes", function(x, ...)
    standardGeneric("listTxbiotypes"))
setGeneric("listUniprotDbs", function(object, ...)
    standardGeneric("listUniprotDbs"))
setGeneric("listUniprotMappingTypes", function(object, ...)
    standardGeneric("listUniprotMappingTypes"))

## M

## N

## O
setGeneric("orderResultsInR", function(x)
    standardGeneric("orderResultsInR"))
setGeneric("orderResultsInR<-", function(x, value)
    standardGeneric("orderResultsInR<-"))

## P
setGeneric("print", function(x, ...)
    standardGeneric("print"))
setGeneric("properties", function(x, ...)
    standardGeneric("properties"))

## Q

## R
setGeneric("requireTable", function(x, db, ...)
    standardGeneric("requireTable"))
setGeneric("returnFilterColumns", function(x)
    standardGeneric("returnFilterColumns"))
setGeneric("returnFilterColumns<-", function(x, value)
    standardGeneric("returnFilterColumns<-"))

## S
setGeneric("seqinfo", function(x)
    standardGeneric("seqinfo"))
setGeneric("setProperty", function(x, value=NULL, ...)
    standardGeneric("setProperty"))
## setGeneric("show", function(object, ...)
##     standardGeneric("show"))
setGeneric("supportedSeqlevelsStyles", function(x)
    standardGeneric("supportedSeqlevelsStyles"))

## T
setGeneric("tablesByDegree", function(x, ...)
    standardGeneric("tablesByDegree"))
setGeneric("tablesForColumns", function(x, attributes, ...)
    standardGeneric("tablesForColumns"))
setGeneric("toSAF", function(x, ...)
    standardGeneric("toSAF"))
setGeneric("transcriptLengths", function(x, with.cds_len=FALSE,
                                         with.utr5_len=FALSE,
                                         with.utr3_len=FALSE, ...)
    standardGeneric("transcriptLengths"))
## setGeneric("transcripts", function(x, ...)
##     standardGeneric("transcripts"))
## setGeneric("transcriptsBy", function(x, ...)
##     standardGeneric("transcriptsBy"))

## U
setGeneric("updateEnsDb", function(x, ...)
    standardGeneric("updateEnsDb"))
setGeneric("useMySQL", function(x, host = "localhost", port = 3306, user, pass)
    standardGeneric("useMySQL"))

## V
setGeneric("value", function(x, db, ...)
    standardGeneric("value"))
setGeneric("value<-", function(x, value)
    standardGeneric("value<-"))

## W
setGeneric("where", function(object, db, with.tables, ...)
    standardGeneric("where"))

