##***********************************************************************
##
##     EnsBb classes
##
##     Main class providing access and functionality for the database.
##
##***********************************************************************
setClass("EnsDb",
         slots = c(ensdb = "DBIConnection",
                   tables = "list",
                   .properties = "list"),
         prototype = list(ensdb = NULL,
                          tables = list(),
                          .properties = list())
         )

#' @title Filters supported by ensembldb
#'
#' @description
#'
#' `ensembldb` supports most of the filters from the [AnnotationFilter]
#' package to retrieve specific content from [EnsDb] databases. These filters
#' can be passed to the methods such as [genes()] with the `filter` parameter
#' or can be added as a *global* filter to an `EnsDb` object (see
#' [addFilter()] for more details). Use [supportedFilters()] to get an
#' overview of all filters supported by `EnsDb` object.
#'
#' @note
#'
#' For users of `ensembldb` version < 2.0: in the `GRangesFilter` from the
#' `AnnotationFilter` package the `condition` parameter was renamed to `type`
#' (to be consistent with the `IRanges` package). In addition,
#' `condition = "overlapping"` is no longer recognized. To retrieve all
#' features overlapping the range `type = "any"` has to be used.
#'
#' @details
#'
#' `ensembldb` supports the following filters from the `AnnotationFilter`
#' package:
#'
#' - `GeneIdFilter`: filter based on the Ensembl gene ID.
#' - `GeneNameFilter`: filter based on the name of the gene as provided
#'   Ensembl. In most cases this will correspond to the official gene symbol.
#' - `SymbolFilter` filter based on the gene names. `EnsDb` objects don't
#'   have a dedicated *symbol* column, the filtering is hence based on the
#'   gene names.
#' - `GeneBiotype`: filter based on the biotype of genes (e.g.
#'   `"protein_coding"`).
#' - `GeneStartFilter`: filter based on the genomic start coordinate of genes.
#' - `GeneEndFilter`: filter based on the genomic end coordinate of genes.
#' - `EntrezidFilter`: filter based on the genes' NCBI Entrezgene ID.
#' - `TxIdFilter`: filter based on the Ensembld transcript ID.
#' - `TxNameFilter`: filter based on the Ensembld transcript ID; no transcript
#'   names are provided in `EnsDb` databases.
#' - `TxBiotypeFilter`: filter based on the transcripts' biotype.
#' - `TxStartFilter`: filter based on the genomic start coordinate of the
#'   transcripts.
#' - `TxEndFilter`: filter based on the genonic end coordinates of the
#'   transcripts.
#' - `ExonIdFilter`: filter based on Ensembl exon IDs.
#' - `ExonRankFilter`: filter based on the index/rank of the exon within the
#'   transcrips.
#' - `ExonStartFilter`: filter based on the genomic start coordinates of the
#'   exons.
#' - `ExonEndFilter`: filter based on the genomic end coordinates of the exons.
#'
#' - `GRangesFilter`: Allows to fetch features within or overlapping specified
#'   genomic region(s)/range(s). This filter takes a `GRanges` object
#'   as input and, if `type = "any"` (the default) will restrict results to
#'   features (genes, transcripts or exons) that are partially overlapping the
#'   region. Alternatively, by specifying `condition = "within"` it will
#'   return features located within the range. In addition, the `GRangesFilter`
#'   `condition = "start"`, `condition = "end"` and `condition = "equal"`
#'   filtering for features with the same start or end coordinate or that are
#'   equal to the `GRanges`.
#'
#'   Note that the type of feature on which the filter is applied depends on
#'   the method that is called, i.e. [genes()] will filter on the
#'   genomic coordinates of genes, [transcripts()] on those of
#'   transcripts and [exons()] on exon coordinates.
#'
#'   Calls to the methods [exonsBy()], [cdsBy()] and
#'   [transcriptsBy()] use the start and end coordinates of the
#'   feature type specified with argument `by` (i.e. `"gene"`,
#'   `"transcript"` or `"exon"`) for the filtering.
#'
#'   If the specified `GRanges` object defines multiple regions, all
#'   features within (or overlapping) any of these regions are returned.
#'
#'   Chromosome names/seqnames can be provided in UCSC format (e.g.
#'   `"chrX"`) or Ensembl format (e.g. `"X"`); see [seqlevelsStyle()] for
#'   more information.
#'
#' - `SeqNameFilter`: filter based on chromosome names.
#' - `SeqStrandFilter`: filter based on the chromosome strand. The strand can
#'   be specified with `value = "+"`, `value = "-"`, `value = -1` or
#'   `value = 1`.
#' - `ProteinIdFilter`: filter based on Ensembl protein IDs. This filter is
#'   only supported if the `EnsDb` provides protein annotations; use the
#'   [hasProteinData()] method to check.
#' - `UniprotFilter`: filter based on Uniprot IDs. This filter is only
#'   supported if the `EnsDb` provides protein annotations; use the
#'   [hasProteinData()] method to check.
#'
#' In addition, the following filters are defined by `ensembldb`:
#'
#' - `TxSupportLevel`: allows to filter results using the provided transcript
#'   support level. Support levels for transcripts are defined by Ensembl
#'   based on the available evidences for a transcript with 1 being the
#'   highest evidence grade and 5 the lowest level. This filter is only
#'   supported on `EnsDb` databases with a db schema version higher 2.1.
#' - `UniprotDbFilter`: allows to filter results based on the specified Uniprot
#'   database name(s).
#' - `UniprotMappingTypeFilter`: allows to filter results based on the mapping
#'   method/type that was used to assign Uniprot IDs to Ensembl protein IDs.
#' - `ProtDomIdFilter`, `ProteinDomainIdFilter`: allows to retrieve entries
#'   from the database matching the provided filter criteria based on their
#'   protein domain ID (*protein_domain_id*).
#' - `ProteinDomainSourceFilter`: filter results based on the source
#'   (database/method) defining the protein domain (e.g. `"pfam"`).
#' - `OnlyCodingTxFilter`: allows to retrieve entries only for protein coding
#'   transcripts, i.e. transcripts with a CDS. This filter does not take any
#'   input arguments.
#'
#' @param condition `character(1)` specifying the *condition* of the
#'     filter. For `character`-based filters (such as
#'     `GeneIdFilter`) `"=="`, `"!="`, `"startsWith"` and `"endsWith"` are
#'     supported. Allowed values for `integer`-based filters (such as
#'     `GeneStartFilter`) are `"=="`, `"!="`, `"<"`. `"<="`, `">"` and `">="`.
#'
#' @param value The value(s) for the filter. For `GRangesFilter` it has to be a
#'     `GRanges` object.
#'
#' @note Protein annotation based filters can only be used if the
#'     `EnsDb` database contains protein annotations, i.e. if `hasProteinData`
#'     is `TRUE`. Also, only protein coding transcripts will have protein
#'     annotations available, thus, non-coding transcripts/genes will not be
#'     returned by the queries using protein annotation filters.
#'
#' @name Filter-classes
#'
#' @md
#'
#' @seealso
#'
#' [supportedFilters()] to list all filters supported for `EnsDb` objects.
#'
#' [listUniprotDbs()] and [listUniprotMappingTypes()] to list all Uniprot
#' database names respectively mapping method types from the database.
#'
#' [GeneIdFilter()] in the `AnnotationFilter` package for more details on the
#' filter objects.
#'
#' [genes()], [transcripts()], [exons()], [listGenebiotypes()],
#' [listTxbiotypes()].
#'
#' [addFilter()] and [filter()] for globally adding filters to an `EnsDb`.
#'
#' @author Johannes Rainer
#'
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
#' grf <- GRangesFilter(GRanges("11", ranges = IRanges(114129278, 114129328),
#'                              strand = "+"), type = "any")
#' library(EnsDb.Hsapiens.v86)
#' edb <- EnsDb.Hsapiens.v86
#' genes(edb, filter = grf)
#'
#' ## Get also all transcripts overlapping that region.
#' transcripts(edb, filter = grf)
#'
#' ## Retrieve all transcripts for the above gene
#' gn <- genes(edb, filter = grf)
#' txs <- transcripts(edb, filter = GeneNameFilter(gn$gene_name))
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
#' transcripts(edb, filter = GeneNameFilter("SKA2"))
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
setClass("OnlyCodingTxFilter",
         contains = "CharacterFilter",
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
setClass("ProtDomIdFilter",
         contains = "CharacterFilter",
         prototype = list(
             condition = "==",
             value = "",
             field = "prot_dom_id"
         ))
#' @return For `ProtDomIdFilter`: A `ProtDomIdFilter` object.
#'
#' @md
#'
#' @rdname Filter-classes
ProtDomIdFilter <- function(value, condition = "==") {
    new("ProtDomIdFilter", condition = condition,
        value = as.character(value))
}
#' @rdname Filter-classes
setClass("ProteinDomainIdFilter",
         contains = "CharacterFilter",
         prototype = list(
             condition = "==",
             value = "",
             field = "protein_domain_id"
         ))
#' @return For `ProteinDomainIdFilter`: A `ProteinDomainIdFilter` object.
#'
#' @md
#'
#' @rdname Filter-classes
ProteinDomainIdFilter <- function(value, condition = "==") {
    new("ProteinDomainIdFilter", condition = condition,
        value = as.character(value))
}

#' @rdname Filter-classes
setClass("ProteinDomainSourceFilter",
         contains = "CharacterFilter",
         prototype = list(
             condition = "==",
             value = "",
             field = "protein_domain_source"
         ))
#' @return For `ProteinDomainSourceFilter`: A `ProteinDomainSourceFilter`
#'     object.
#'
#' @md
#'
#' @rdname Filter-classes
ProteinDomainSourceFilter <- function(value, condition = "==") {
    new("ProteinDomainSourceFilter", condition = condition,
        value = as.character(value))
}

############################################################
## UniprotDbFilter
#' @rdname Filter-classes
setClass("UniprotDbFilter",
         contains = "CharacterFilter",
         prototype = list(
             condition = "==",
             value = "",
             field = "uniprot_db"
         ))
#' @return For `UniprotDbFilter`: A `UniprotDbFilter` object.
#'
#' @md
#'
#' @rdname Filter-classes
UniprotDbFilter <- function(value, condition = "==") {
    new("UniprotDbFilter", condition = condition,
        value = as.character(value))
}

############################################################
## UniprotMappingTypeFilter
#' @rdname Filter-classes
setClass("UniprotMappingTypeFilter",
         contains = "CharacterFilter",
         prototype = list(
             condition = "==",
             value = "",
             field = "uniprot_mapping_type"
         ))
#' @return For `UniprotMappingTypeFilter`: A `UniprotMappingTypeFilter` object.
#'
#' @md
#'
#' @rdname Filter-classes
UniprotMappingTypeFilter <- function(value, condition = "==") {
    new("UniprotMappingTypeFilter", condition = condition,
        value = as.character(value))
}

#' @rdname Filter-classes
setClass("TxSupportLevelFilter",
         contains = "IntegerFilter",
         prototype = list(
             condition = "==",
             value = 0L,
             field = "tx_support_level"
         ))
#' @return For `TxSupportLevel`: A `TxSupportLevel` object.
#'
#' @md
#'
#' @rdname Filter-classes
TxSupportLevelFilter <- function(value, condition = "==") {
    if (!is.numeric(value))
        stop("Parameter 'value' has to be numeric")
    new("TxSupportLevelFilter", condition = condition,
        value = as.integer(value))
}
