##***********************************************************************
##
##     Generic methods
##
##***********************************************************************
if(!isGeneric("column"))
    setGeneric("column", function(object, db, with.tables, ...)
        standardGeneric("column"))
if(!isGeneric("buildQuery"))
    setGeneric("buildQuery", function(x, ...)
        standardGeneric("buildQuery"))
if(!isGeneric("cleanColumns"))
    setGeneric("cleanColumns", function(x, columns, ...)
        starndardGeneric("cleanColumns"))
if(!isGeneric("condition"))
    setGeneric("condition", function(x, ...)
        standardGeneric("condition"))
setGeneric("condition<-", function(x, value)
        standardGeneric("condition<-"))
setGeneric("dbSeqlevelsStyle", function(x, ...)
    standardGeneric("dbSeqlevelsStyle"))

if(!isGeneric("genes"))
    setGeneric("genes", function(x, ...)
        standardGeneric("genes"))
if(!isGeneric("getWhat"))
    setGeneric("getWhat", function(x, ...)
        standardGeneric("getWhat"))
if(!isGeneric("ensemblVersion"))
    setGeneric("ensemblVersion", function(x)
        standardGeneric("ensemblVersion"))
if(!isGeneric("exons"))
    setGeneric("exons", function(x, ...)
        standardGeneric("exons"))
if(!isGeneric("exonsBy"))
    setGeneric("exonsBy", function(x, ...)
        standardGeneric("exonsBy"))

setGeneric("getGeneRegionTrackForGviz", function(x, ...)
    standardGeneric("getGeneRegionTrackForGviz"))

if(!isGeneric("getGenomeFaFile"))
    setGeneric("getGenomeFaFile", function(x, ...)
        standardGeneric("getGenomeFaFile"))
if(!isGeneric("getGenomeTwoBitFile"))
    setGeneric("getGenomeTwoBitFile", function(x, ...)
        standardGeneric("getGenomeTwoBitFile"))
if(!isGeneric("getMetadataValue"))
    setGeneric("getMetadataValue", function(x, name)
        standardGeneric("getMetadataValue"))
if(!isGeneric("listColumns")){
    setGeneric("listColumns", function(x, ...)
        standardGeneric("listColumns"))
}
if(!isGeneric("listGenebiotypes")){
    setGeneric("listGenebiotypes", function(x, ...)
        standardGeneric("listGenebiotypes"))
}
if(!isGeneric("listTxbiotypes")){
    setGeneric("listTxbiotypes", function(x, ...)
        standardGeneric("listTxbiotypes"))
}
if(!isGeneric("lengthOf"))
    setGeneric("lengthOf", function(x, ...)
        standardGeneric("lengthOf"))
if(!isGeneric("print"))
    setGeneric("print", function(x, ...)
        standardGeneric("print"))
if(!isGeneric("requireTable"))
    setGeneric("requireTable", function(x, db, ...)
        standardGeneric("requireTable"))

setGeneric("supportedSeqlevelsStyles", function(x)
           standardGeneric("supportedSeqlevelsStyles"))

if(!isGeneric("seqinfo"))
    setGeneric("seqinfo", function(x)
        standardGeneric("seqinfo"))
if(!isGeneric("show"))
    setGeneric("show", function(object, ...)
        standardGeneric("show"))
if(!isGeneric("toSAF"))
    setGeneric("toSAF", function(x, ...)
        standardGeneric("toSAF"))
if(!isGeneric("listTables")){
    setGeneric("listTables", function(x, ...)
        standardGeneric("listTables"))
}
if(!isGeneric("tablesByDegree")){
    setGeneric("tablesByDegree", function(x, ...)
        standardGeneric("tablesByDegree"))
}
if(!isGeneric("tablesForColumns"))
    setGeneric("tablesForColumns", function(x, attributes, ...)
        standardGeneric("tablesForColumns"))

if(!isGeneric("transcriptLengths"))
    setGeneric("transcriptLengths", function(x, with.cds_len=FALSE,
                                             with.utr5_len=FALSE,
                                             with.utr3_len=FALSE, ...)
        standardGeneric("transcriptLengths"))

if(!isGeneric("transcripts"))
    setGeneric("transcripts", function(x, ...)
        standardGeneric("transcripts"))
if(!isGeneric("transcriptsBy"))
    setGeneric("transcriptsBy", function(x, ...)
        standardGeneric("transcriptsBy"))
setGeneric("updateEnsDb", function(x, ...)
    standardGeneric("updateEnsDb"))
##if(!isGeneric("value"))
    setGeneric("value", function(x, db, ...)
        standardGeneric("value"))
setGeneric("value<-", function(x, value)
    standardGeneric("value<-"))
if(!isGeneric("where"))
    setGeneric("where", function(object, db, with.tables, ...)
        standardGeneric("where"))

####============================================================
##  Private methods
##
####------------------------------------------------------------
setGeneric("properties", function(x, ...)
    standardGeneric("properties"))
## setGeneric("properties<-", function(x, name, value, ...)
##             standardGeneric("properties<-"))
setGeneric("getProperty", function(x, name=NULL, ...)
    standardGeneric("getProperty"))
setGeneric("setProperty", function(x, value=NULL, ...)
    standardGeneric("setProperty"))
setGeneric("formatSeqnamesForQuery", function(x, sn, ...)
    standardGeneric("formatSeqnamesForQuery"))
setGeneric("formatSeqnamesFromQuery", function(x, sn, ...)
    standardGeneric("formatSeqnamesFromQuery"))
