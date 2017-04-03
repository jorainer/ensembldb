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

#' @title Filters supported by ensembldb
#'
#' @description \code{ensembldb} supports most of the filters from the
#'     \code{\link{AnnotationFilter}} package to retrieve specific content from
#'     \code{\linkS4class{EnsDb}} databases.
#'
#' @note For users of \code{ensembldb} version < 2.0: in the
#'     \code{\link[AnnotationFilter]{GRangesFilter}} from the
#'     \code{AnnotationFilter} package the \code{condition} parameter was
#'     renamed to \code{type} (to be consistent with the \code{IRanges} package)
#'     . In addition, the \code{condition = "overlapping"} is no longer
#'     recognized. To retrieve all features overlapping the range
#'     \code{type = "any"} has to be used.
#'     
#' @details \code{ensembldb} supports the following filters from the
#' \code{AnnotationFilter} package:
#' 
#' \describe{
#' 
#' \item{GeneIdFilter}{
#'     filter based on the Ensembl gene ID.
#' }
#'
#' \item{GenenameFilter}{
#'     filter based on the name of the gene as provided by Ensembl. In most cases
#'     this will correspond to the official gene symbol.
#' }
#'
#' \item{SymbolFilter}{
#'     filter based on the gene names. \code{\linkS4class{EnsDb}} objects don't
#'     have a dedicated \emph{symbol} column, the filtering is hence based on the
#'     gene names.
#' }
#'
#' \item{GeneBiotype}{
#'     filter based on the biotype of genes (e.g. \code{"protein_coding"}).
#' }
#'
#' \item{GeneStartFilter}{
#'     filter based on the genomic start coordinate of genes.
#' }
#' 
#' \item{GeneEndFilter}{
#'     filter based on the genomic end coordinate of genes.
#' }
#' 
#' \item{EntrezidFilter}{
#'     filter based on the genes' NCBI Entrezgene ID.
#' }
#' 
#' \item{TxIdFilter}{
#'     filter based on the Ensembld transcript ID.
#' }
#' 
#' \item{TxNameFilter}{
#'     filter based on the Ensembld transcript ID; no transcript names are
#'     provided in \code{\linkS4class{EnsDb}} databases.
#' }
#' 
#' \item{TxBiotypeFilter}{
#'     filter based on the transcripts' biotype.
#' }
#' 
#' \item{TxStartFilter}{
#'     filter based on the genomic start coordinate of the transcripts.
#' }
#' 
#' \item{TxEndFilter}{
#'     filter based on the genonic end coordinates of the transcripts.
#' }
#' 
#' \item{ExonIdFilter}{
#'     filter based on Ensembl exon IDs.
#' }
#'
#' \item{ExonRankFilter}{
#'     filter based on the index/rank of the exon within the transcrips.
#' }
#'
#' \item{ExonStartFilter}{
#'     filter based on the genomic start coordinates of the exons.
#' }
#' 
#' \item{ExonEndFilter}{
#'     filter based on the genomic end coordinates of the exons.
#' }
#'
#' \item{GRangesFilter}{
#'     Allows to fetch features within or overlapping specified genomic region(s)/
#'     range(s). This filter takes a \code{\link[GenomicRanges]{GRanges}} object
#'     as input and, if \code{type = "any"} (the default) will restrict
#'     results to features (genes, transcripts or exons) that are partially
#'     overlapping the region. Alternatively, by specifying
#'     \code{condition = "within"} it will return features located within the
#'     range. In addition, the \code{\link[AnnotationFilter]{GRangesFilter}}
#'     supports \code{condition = "start"}, \code{condition = "end"} and
#'     \code{condition = "equal"} filtering for features with the same start or
#'     end coordinate or that are equal to the \code{GRanges}.
#'
#'     Note that the type of feature on which the filter is applied depends on
#'     the method that is called, i.e. \code{\link{genes}} will filter on the
#'     genomic coordinates of genes, \code{\link{transcripts}} on those of
#'     transcripts and \code{\link{exons}} on exon coordinates.
#'
#'     Calls to the methods \code{\link{exonsBy}}, \code{\link{cdsBy}} and
#'     \code{\link{transcriptsBy}} use the start and end coordinates of the
#'     feature type specified with argument \code{by} (i.e. \code{"gene"},
#'     \code{"transcript"} or \code{"exon"}) for the filtering.
#'
#'     If the specified \code{GRanges} object defines multiple regions, all
#'     features within (or overlapping) any of these regions are returned.
#'
#'     Chromosome names/seqnames can be provided in UCSC format (e.g.
#'     \code{"chrX"}) or Ensembl format (e.g. \code{"X"}); see
#'     \code{\link{seqlevelsStyle}} for more information. 
#' }
#'
#' \item{SeqNameFilter}{
#'     filter based on chromosome names.
#' }
#'
#' \item{SeqStrandFilter}{
#'     filter based on the chromosome strand. The strand can be specified with
#'     \code{value = "+"}, \code{value = "-"}, \code{value = -1} or
#'     \code{value = 1}.
#' }
#' 
#' \item{ProteinIdFilter}{
#'     filter based on Ensembl protein IDs. This filter is only supported if the
#'     \code{\linkS4class{EnsDb}} provides protein annotations; use the
#'     \code{\link{hasProteinData}} method to evaluate.
#' }
#'
#' \item{UniprotFilter}{
#'     filter based on Uniprot IDs. This filter is only supported if the
#'     \code{\linkS4class{EnsDb}} provides protein annotations; use the
#'     \code{\link{hasProteinData}} method to evaluate.
#' }
#'
#' }
#'
#' In addition, the following filters are defined by \code{ensembldb}:
#' \describe{
#' 
#' \item{UniprotDbFilter}{
#'     allows to filter results based on the specified Uniprot database name(s).
#' }
#' 
#' \item{UniprotMappingTypeFilter}{
#'     allows to filter results based on the mapping method/type that was used
#'     to assign Uniprot IDs to Ensembl protein IDs.
#' }
#'
#' \item{ProtDomIdFilter}{
#'     allows to retrieve entries from the database matching the provided filter
#'     criteria based on their protein  domain ID (\emph{protein_domain_id}).
#' }
#'
#' \item{OnlyCodingTxFilter}{
#'     allows to retrieve entries only for protein coding transcripts, i.e.
#'     transcripts with a CDS. This filter does not take any input arguments.
#' }
#' 
#' }
#'
#' @param condition \code{character(1)} specifying the \emph{condition} of the
#'     filter. For \code{character}-based filters (such as
#'     \code{\link[AnnotationFilter]{GeneIdFilter}}) \code{"=="}, \code{"!="},
#'     \code{"startsWith"} and \code{"endsWith"} are supported. Allowed values
#'     for \code{integer}-based filters (such as
#'     \code{\link[AnnotationFilter]{GeneStartFilter}}) are \code{"=="},
#'     \code{"!="}, \code{"<"}. \code{"<="}, \code{">"} and \code{">="}.
#' 
#' @param value The value(s) for the filter. For
#'     \code{\link[AnnotationFilter]{GRangesFilter}} it has to be a
#'     \code{\link[GenomicRanges]{GRanges}} object.
#' 
#' @note Protein annotation based filters can only be used if the
#'     \code{\linkS4class{EnsDb}} database contains protein annotations, i.e.
#'     if \code{\link{hasProteinData}} is \code{TRUE}. Also, only protein coding
#'     transcripts will have protein annotations available, thus, non-coding
#'     transcripts/genes will not be returned by the queries using protein
#'     annotation filters.
#' 
#' @name Filter-classes
#' @seealso
#' \code{\link{supportedFilters}} to list all filters supported for \code{EnsDb}
#'     objects.
#'     \code{\link{listUniprotDbs}} and \code{\link{listUniprotMappingTypes}} to
#'     list all Uniprot database names respectively mapping method types from
#'     the database.
#'
#'     \code{\link[AnnotationFilter]{GeneIdFilter}} for more details on the
#'     filter objects.
#'
#'     \code{\link{genes}}, \code{\link{transcripts}}, \code{\link{exons}},
#'     \code{\link{listGenebiotypes}}, \code{\link{listTxbiotypes}}.
#' 
#' @author Johannes Rainer
#' @examples
#'
#' ## Create a filter that could be used to retrieve all informations for
#' ## the respective gene.
#' gif <- GeneIdFilter("ENSG00000012817")
#' gif
#' 
#' ## Create a filter for a chromosomal end position of a gene
#' sef <- GeneEndFilter(10000, condition = ">")
#' sef
#' 
#' ## For additional examples see the help page of "genes".
#' 
#' 
#' ## Example for GRangesFilter:
#' ## retrieve all genes overlapping the specified region
#' grf <- GRangesFilter(GRanges("11", ranges = IRanges(114000000, 114000050),
#'                              strand = "+"), type = "any")
#' library(EnsDb.Hsapiens.v75)
#' edb <- EnsDb.Hsapiens.v75
#' genes(edb, filter = grf)
#' 
#' ## Get also all transcripts overlapping that region.
#' transcripts(edb, filter = grf)
#' 
#' ## Retrieve all transcripts for the above gene
#' gn <- genes(edb, filter = grf)
#' txs <- transcripts(edb, filter = GenenameFilter(gn$gene_name))
#' ## Next we simply plot their start and end coordinates.
#' plot(3, 3, pch=NA, xlim=c(start(gn), end(gn)), ylim=c(0, length(txs)),
#' yaxt="n", ylab="")
#' ## Highlight the GRangesFilter region
#' rect(xleft=start(grf), xright=end(grf), ybottom=0, ytop=length(txs),
#' col="red", border="red")
#' for(i in 1:length(txs)){
#'     current <- txs[i]
#'     rect(xleft=start(current), xright=end(current), ybottom=i-0.975, ytop=i-0.125, border="grey")
#'     text(start(current), y=i-0.5,pos=4, cex=0.75, labels=current$tx_id)
#' }
#' ## Thus, we can see that only 4 transcripts of that gene are indeed
#' ## overlapping the region.
#' 
#' 
#' ## No exon is overlapping that region, thus we're not getting anything
#' exons(edb, filter = grf)
#' 
#' 
#' ## Example for ExonRankFilter
#' ## Extract all exons 1 and (if present) 2 for all genes encoded on the
#' ## Y chromosome
#' exons(edb, columns = c("tx_id", "exon_idx"),
#'       filter=list(SeqNameFilter("Y"),
#'                   ExonRankFilter(3, condition = "<")))
#' 
#' 
#' ## Get all transcripts for the gene SKA2
#' transcripts(edb, filter = GenenameFilter("SKA2"))
#' 
#' ## Which is the same as using a SymbolFilter
#' transcripts(edb, filter = SymbolFilter("SKA2"))
#' 
#' 
#' ## Create a ProteinIdFilter:
#' pf <- ProteinIdFilter("ENSP00000362111")
#' pf
#' ## Using this filter would retrieve all database entries that are associated
#' ## with a protein with the ID "ENSP00000362111"
#' if (hasProteinData(edb)) {
#'     res <- genes(edb, filter = pf)
#'     res
#' }
#'
#' ## UniprotFilter:
#' uf <- UniprotFilter("O60762")
#' ## Get the transcripts encoding that protein:
#' if (hasProteinData(edb)) {
#'     transcripts(edb, filter = uf)
#'     ## The mapping Ensembl protein ID to Uniprot ID can however be 1:n:
#'     transcripts(edb, filter = TxIdFilter("ENST00000371588"),
#'         columns = c("protein_id", "uniprot_id"))
#' }
#'
#' ## ProtDomIdFilter:
#' pdf <- ProtDomIdFilter("PF00335")
#' ## Also here we could get all transcripts related to that protein domain
#' if (hasProteinData(edb)) {
#'     transcripts(edb, filter = pdf, columns = "protein_id")
#' }
#'
NULL

############################################################
## OnlyCodingTxFilter
##
## That's a special case filter that just returns transcripts
## that have tx_cds_seq_start defined (i.e. not NULL).
#' @rdname Filter-classes
setClass("OnlyCodingTxFilter", contains = "CharacterFilter",
         prototype = list(
             condition = "==",
             value = character(),
             field = "empty"
         ))
#' @rdname Filter-classes
OnlyCodingTxFilter <- function() {
    new("OnlyCodingTxFilter")
}

############################################################
## ProtDomIdFilter
#' @rdname Filter-classes
setClass("ProtDomIdFilter", contains = "CharacterFilter",
         prototype = list(
             condition = "==",
             value = "",
             field = "prot_dom_id"
         ))
#' @return For \code{ProtDomIdFilter}: A \code{ProtDomIdFilter} object.
#' @rdname Filter-classes
ProtDomIdFilter <- function(value, condition = "==") {
    new("ProtDomIdFilter", condition = condition,
        value = as.character(value))
}

############################################################
## UniprotDbFilter
#' @rdname Filter-classes
setClass("UniprotDbFilter", contains = "CharacterFilter",
         prototype = list(
             condition = "==",
             values = "",
             field = "uniprot_db"
         ))
#' @return For \code{UniprotDbFilter}: A \code{UniprotDbFilter} object.
#' @rdname Filter-classes
UniprotDbFilter <- function(value, condition = "==") {
    new("UniprotDbFilter", condition = condition,
        value = as.character(value))
}

############################################################
## UniprotMappingTypeFilter
#' @rdname Filter-classes
setClass("UniprotMappingTypeFilter", contains = "CharacterFilter",
         prototype = list(
             condition = "==",
             values = "",
             field = "uniprot_mapping_type"
         ))
#' @return For \code{UniprotMappingTypeFilter}: A
#' \code{UniprotMappingTypeFilter} object.
#' @rdname Filter-classes
UniprotMappingTypeFilter <- function(value, condition = "==") {
    new("UniprotMappingTypeFilter", condition = condition,
        value = as.character(value))
}

