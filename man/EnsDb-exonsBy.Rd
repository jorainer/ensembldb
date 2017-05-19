\name{exonsBy}
\Rdversion{1.1}
\docType{class}
\alias{disjointExons,EnsDb-method}
\alias{cdsBy}
\alias{cdsBy,EnsDb-method}
\alias{fiveUTRsByTranscript,EnsDb-method}
\alias{threeUTRsByTranscript,EnsDb-method}
\alias{exons}
\alias{exons,EnsDb-method}
\alias{exonsBy}
\alias{exonsBy,EnsDb-method}
\alias{exonsByOverlaps,EnsDb-method}
\alias{genes}
\alias{genes,EnsDb-method}
\alias{toSAF}
\alias{toSAF,GRangesList-method}
\alias{transcripts}
\alias{transcripts,EnsDb-method}
\alias{transcriptsBy}
\alias{transcriptsBy,EnsDb-method}
\alias{transcriptsByOverlaps,EnsDb-method}
\alias{promoters}
\alias{promoters,EnsDb-method}

\title{Retrieve annotation data from an Ensembl based package}
\description{
  Retrieve gene/transcript/exons annotations stored in an Ensembl based
  database package generated with the \code{\link{makeEnsembldbPackage}}
  function. Parameter \code{filter} enables to define filters to
  retrieve only specific data.
}
\usage{

\S4method{exons}{EnsDb}(x, columns = listColumns(x,"exon"),
        filter = AnnotationFilterList(), order.by,
        order.type = "asc", return.type = "GRanges")

\S4method{exonsBy}{EnsDb}(x, by = c("tx", "gene"),
        columns = listColumns(x, "exon"), filter =
        AnnotationFilterList(), use.names = FALSE)

\S4method{exonsByOverlaps}{EnsDb}(x, ranges, maxgap = 0L, minoverlap = 1L,
        type = c("any", "start", "end"), columns = listColumns(x, "exon"),
        filter = AnnotationFilterList())

\S4method{transcripts}{EnsDb}(x, columns = listColumns(x, "tx"),
        filter = AnnotationFilterList(), order.by, order.type = "asc",
        return.type = "GRanges")

\S4method{transcriptsBy}{EnsDb}(x, by = c("gene", "exon"),
        columns = listColumns(x, "tx"), filter = AnnotationFilterList())

\S4method{transcriptsByOverlaps}{EnsDb}(x, ranges, maxgap = 0L,
        minoverlap = 1L, type = c("any", "start", "end"),
        columns = listColumns(x, "tx"), filter = AnnotationFilterList())

\S4method{promoters}{EnsDb}(x, upstream = 2000, downstream = 200, ...)

\S4method{genes}{EnsDb}(x, columns = c(listColumns(x, "gene"), "entrezid"),
        filter = AnnotationFilterList(), order.by, order.type = "asc",
        return.type = "GRanges")

\S4method{disjointExons}{EnsDb}(x, aggregateGenes = FALSE,
        includeTranscripts = TRUE, filter = AnnotationFilterList(), ...)

\S4method{cdsBy}{EnsDb}(x, by = c("tx", "gene"), columns = NULL,
        filter = AnnotationFilterList(), use.names = FALSE)

\S4method{fiveUTRsByTranscript}{EnsDb}(x, columns = NULL,
        filter = AnnotationFilterList())

\S4method{threeUTRsByTranscript}{EnsDb}(x, columns = NULL,
        filter = AnnotationFilterList())

\S4method{toSAF}{GRangesList}(x, ...)

}
\arguments{

  (In alphabetic order)

  \item{...}{
    For \code{promoters}: additional arguments to be passed to the
    \code{transcripts} method.
  }

  \item{aggregateGenes}{
    For \code{disjointExons}: When \code{FALSE} (default) exon fragments
    that overlap multiple genes are dropped. When \code{TRUE}, all
    fragments are kept and the \code{gene_id} metadata column includes
    all gene IDs that overlap the exon fragment.
  }

  \item{by}{
    For \code{exonsBy}: wheter exons sould be fetched by genes
    or by transcripts; as in the corresponding function of the
    \code{GenomicFeatures} package.
    For \code{transcriptsBy}: whether
    transcripts should be fetched by genes or by exons; fetching
    transcripts by cds as supported by the
    \code{\link[GenomicFeatures]{transcriptsBy}} method in the
    \code{GenomicFeatures} package is currently not implemented.
    For \code{cdsBy}: whether cds should be fetched by transcript of by
    gene.
  }

  \item{columns}{
    Columns to be retrieved from the database tables.

    Default values for \code{genes}  are all columns from the \code{gene}
    database table, for \code{exons} and \code{exonsBy} the column names of
    the \code{exon} database table table and for \code{transcript} and
    \code{transcriptBy} the columns of the \code{tx} data base table
    (see details below for more information).

    Note that any of the column names of the database tables can be
    submitted to any of the methods (use \code{\link{listTables}} or
    \code{\link{listColumns}} methods for a complete list of allowed
    column names).

    For \code{cdsBy}: this argument is only supported for for
    \code{by="tx"}.
  }

  \item{downstream}{
    For method \code{promoters}: the number of nucleotides downstream of
    the transcription start site that should be included in the promoter region.
  }

  \item{filter}{
    A filter describing which results to retrieve from the database. Can
    be a single object extending
    \code{\link[AnnotationFilter]{AnnotationFilter}}, an
    \code{\link[AnnotationFilter]{AnnotationFilterList}} object
    combining several such objects or a \code{formula} representing a
    filter expression (see examples below or
    \code{\link[AnnotationFilter]{AnnotationFilter}} for more details).
  }

  \item{includeTranscripts}{
    For \code{disjointExons}: When \code{TRUE} (default) a
    \code{tx_name} metadata column is included that lists all transcript
    IDs that overlap the exon fragment. Note: this is different to the
    \code{\link[GenomicFeatures]{disjointExons}} function in the
    \code{GenomicFeatures} package, that lists the transcript names, not
    IDs.
  }

  \item{maxgap}{
    For \code{exonsByOverlaps} and \code{transcriptsByOverlaps}: see
    \code{\link[GenomicFeatures]{exonsByOverlaps}} help page in the
    \code{GenomicFeatures} package.
  }

  \item{minoverlap}{
    For \code{exonsByOverlaps} and \code{transcriptsByOverlaps}: see
    \code{\link[GenomicFeatures]{exonsByOverlaps}} help page in the
    \code{GenomicFeatures} package.
  }

  \item{order.by}{
    Character vector specifying the column(s) by which the result should
    be ordered. This can be either in the form of
    \code{"gene_id, seq_name"} or \code{c("gene_id", "seq_name")}.
  }

  \item{order.type}{
    If the results should be ordered ascending
    (\code{asc}, default) or descending (\code{desc}).
  }

  \item{ranges}{
    For \code{exonsByOverlaps} and \code{transcriptsByOverlaps}: a
    \code{GRanges} object specifying the genomic regions.
  }

  \item{return.type}{
    Type of the returned object. Can be either
    \code{"data.frame"}, \code{"DataFrame"} or \code{"GRanges"}. In the
    latter case the return object will be a \code{GRanges} object with
    the GRanges specifying the chromosomal start and end coordinates of
    the feature (gene, transcript or exon, depending whether \code{genes},
    \code{transcripts} or \code{exons} was called). All additional
    columns are added as metadata columns to the GRanges object.
  }

  \item{type}{
    For \code{exonsByOverlaps} and \code{transcriptsByOverlaps}: see
    \code{\link[GenomicFeatures]{exonsByOverlaps}} help page in the
    \code{GenomicFeatures} package.
  }

  \item{upstream}{
    For method \code{promoters}: the number of nucleotides upstream of
    the transcription start site that should be included in the promoter region.
  }

  \item{use.names}{
    For \code{cdsBy} and \code{exonsBy}: only for \code{by="gene"}: use the names of the
    genes instead of their IDs as names of the resulting
    \code{GRangesList}.
  }

  \item{x}{
    For \code{toSAF} a \code{GRangesList} object.
    For all other methods an \code{EnsDb} instance.
  }

}
\section{Methods and Functions}{
  \describe{

    \item{exons}{
      Retrieve exon information from the database. Additional
      columns from transcripts or genes associated with the exons can be specified
      and are added to the respective exon annotation.
    }

    \item{exonsBy}{
      Retrieve exons grouped by transcript or by gene. This
      function returns a \code{GRangesList} as does the analogous function
      in the \code{GenomicFeatures} package. Using the \code{columns}
      parameter it is possible to determine which additional values should
      be retrieved from the database. These will be included in the
      \code{GRanges} object for the exons as metadata columns.
      The exons in the inner \code{GRanges} are ordered by the exon
      index within the transcript (if \code{by="tx"}), or increasingly by the
      chromosomal start position of the exon or decreasingly by the chromosomal end
      position of the exon depending whether the gene is encoded on the
      + or - strand (for \code{by="gene"}).
      The \code{GRanges} in the \code{GRangesList} will be ordered by
      the name of the gene or transcript.
    }

    \item{exonsByOverlaps}{
      Retrieve exons overlapping specified genomic ranges. For
      more information see
      \code{\link[GenomicFeatures]{exonsByOverlaps}} method in the
      \code{GenomicFeatures} package. The functionality is to some
      extent similar and redundant to the \code{exons} method in
      combination with \code{\link{GRangesFilter}} filter.
    }

    \item{transcripts}{
      Retrieve transcript information from the database. Additional
      columns from genes or exons associated with the transcripts can be specified
      and are added to the respective transcript annotation.
    }

    \item{transcriptsBy}{
      Retrieve transcripts grouped by gene or exon. This
      function returns a \code{GRangesList} as does the analogous function
      in the \code{GenomicFeatures} package. Using the \code{columns}
      parameter it is possible to determine which additional values should
      be retrieved from the database. These will be included in the
      \code{GRanges} object for the transcripts as metadata columns.
      The transcripts in the inner \code{GRanges} are ordered  increasingly by the
      chromosomal start position of the transcript for genes encoded on
      the + strand and in a decreasing manner by the chromosomal end
      position of the transcript for genes encoded on the - strand.
      The \code{GRanges} in the \code{GRangesList} will be ordered by
      the name of the gene or exon.
    }

    \item{transcriptsByOverlaps}{
      Retrieve transcripts overlapping specified genomic ranges. For
      more information see
      \code{\link[GenomicFeatures]{transcriptsByOverlaps}} method in the
      \code{GenomicFeatures} package. The functionality is to some
      extent similar and redundant to the \code{transcripts} method in
      combination with \code{\link{GRangesFilter}} filter.
    }

    \item{promoters}{
      Retrieve promoter information from the database. Additional
      columns from genes or exons associated with the promoters can be specified
      and are added to the respective promoter annotation.
    }
    \item{genes}{
      Retrieve gene information from the database. Additional columns
      from transcripts or exons associated with the genes can be
      specified and are added to the respective gene annotation. Note
      that column \code{"entrezid"} is a \code{list} of Entrezgene
      identifiers to accomodate the potential 1:n mapping between
      Ensembl genes and Entrezgene IDs.
    }

    \item{disjointExons}{
      This method is identical to
      \code{\link[GenomicFeatures]{disjointExons}} defined in the
      \code{GenomicFeatures} package. It creates a \code{GRanges} of
      non-overlapping exon parts with metadata columns of \code{gene_id}
      and \code{exonic_part}. Exon parts that overlap more than one gene
      can be dropped with \code{aggregateGenes=FALSE}.
    }

    \item{cdsBy}{
      Returns the coding region grouped either by transcript or by
      gene. Each element in the \code{GRangesList} represents the cds
      for one transcript or gene, with the individual ranges
      corresponding to the coding part of its exons.
      For \code{by="tx"} additional annotation columns can be added to
      the individual \code{GRanges} (in addition to the default columns
      \code{exon_id} and \code{exon_rank}).
      Note that the \code{GRangesList} is sorted by its names.
    }

    \item{fiveUTRsByTranscript}{
      Returns the 5' untranslated region for protein coding
      transcripts.
    }

    \item{threeUTRsByTranscript}{
      Returns the 3' untranslated region for protein coding
      transcripts.
    }

    \item{toSAF}{
      Reformats a \code{GRangesList} object into a
      \code{data.frame} corresponding to a standard SAF (Simplified
      Annotation Format) file (i.e. with column names \code{"GeneID"},
      \code{"Chr"}, \code{"Start"}, \code{"End"} and
      \code{"Strand"}). Note: this method makes only sense on a
      \code{GRangesList} that groups features (exons, transcripts) by gene.
    }

  }
}
\details{
  A detailed description of all database tables and the associated
  attributes/column names is also given in the vignette of this package.
  An overview of the columns is given below:
  \describe{
    \item{gene_id}{the Ensembl gene ID of the gene.}
    \item{gene_name}{the name of the gene (in most cases its official symbol).}
    \item{entrezid}{the NCBI Entrezgene ID of the gene. Note that this
      column contains a \code{list} of Entrezgene identifiers to
      accommodate the potential 1:n mapping between Ensembl genes and
      Entrezgene IDs.}
    \item{gene_biotype}{the biotype of the gene.}
    \item{gene_seq_start}{the start coordinate of the gene on the
      sequence (usually a chromosome).}
    \item{gene_seq_end}{the end coordinate of the gene.}
    \item{seq_name}{the name of the sequence the gene is encoded
      (usually a chromosome).}
    \item{seq_strand}{the strand on which the gene is encoded}
    \item{seq_coord_system}{the coordinate system of the sequence.}
    \item{tx_id}{the Ensembl transcript ID.}
    \item{tx_biotype}{the biotype of the transcript.}
    \item{tx_seq_start}{the chromosomal start coordinate of the transcript.}
    \item{tx_seq_end}{the chromosomal end coordinate of the transcript.}
    \item{tx_cds_seq_start}{the start coordinate of the coding region of
      the transcript (NULL for non-coding transcripts).}
    \item{tx_cds_seq_end}{the end coordinate of the coding region.}
    \item{exon_id}{the ID of the exon. In Ensembl, each exon specified
      by a unique chromosomal start and end position has its own
      ID. Thus, the same exon might be part of several transcripts.}
    \item{exon_seq_start}{the chromosomal start coordinate of the exon.}
    \item{exon_seq_end}{the chromosomal end coordinate of the exon.}
    \item{exon_idx}{the index of the exon in the transcript model. As
      noted above, an exon can be part of several transcripts and thus
      its position inside these transcript might differ.}
  }

  Many \code{EnsDb} databases provide also protein related
  annotations. See \code{\link{listProteinColumns}} for more information.
}
\note{
  Ensembl defines genes not only on standard chromosomes, but also on
  patched chromosomes and chromosome variants. Thus it might be
  advisable to restrict the queries to just those chromosomes of
  interest (e.g. by specifying a \code{SeqNameFilter(c(1:22, "X", "Y"))}).
  In addition, also so called LRG genes (Locus Reference Genomic) are defined in
  Ensembl. Their gene id starts with LRG instead of ENS for Ensembl
  genes, thus, a filter can be applied to specifically select those
  genes or exclude those genes (see examples below).

  Depending on the value of the global option
  \code{"ucscChromosomeNames"} (use
  \code{getOption(ucscChromosomeNames, FALSE)} to get its value or
  \code{option(ucscChromosomeNames=TRUE)} to change its value)
  the sequence/chromosome names of the returned \code{GRanges} objects
  or provided in the returned \code{data.frame} or \code{DataFrame}
  correspond to Ensembl chromosome names (if value is \code{FALSE}) or
  UCSC chromosome names (if \code{TRUE}). This ensures a better
  integration with the \code{Gviz} package, in which this option is set
  by default to \code{TRUE}.
}

\value{
  For \code{exons}, \code{transcripts} and \code{genes},
  a \code{data.frame}, \code{DataFrame}
  or a \code{GRanges}, depending on the value of the
  \code{return.type} parameter. The result is ordered as specified by
  the parameter \code{order.by} or, if not provided, by \code{seq_name}
  and chromosomal start coordinate, but NOT by any ordering of values in
  eventually submitted filter objects.

  For \code{exonsBy}, \code{transcriptsBy}:
  a \code{GRangesList}, depending on the value of the
  \code{return.type} parameter. The results are ordered by the value of the
  \code{by} parameter.

  For \code{exonsByOverlaps} and \code{transcriptsByOverlaps}: a
  \code{GRanges} with the exons or transcripts overlapping the specified
  regions.

  For \code{toSAF}: a \code{data.frame} with column names
  \code{"GeneID"} (the group name from the \code{GRangesList}, i.e. the
  ID by which the \code{GRanges} are split), \code{"Chr"} (the seqnames
  from the \code{GRanges}), \code{"Start"} (the start coordinate),
  \code{"End"} (the end coordinate) and \code{"Strand"} (the strand).

  For \code{disjointExons}: a \code{GRanges} of non-overlapping exon
  parts.

  For \code{cdsBy}: a \code{GRangesList} with \code{GRanges} per either
  transcript or exon specifying the start and end coordinates of the
  coding region of the transcript or gene.

  For \code{fiveUTRsByTranscript}: a \code{GRangesList} with
  \code{GRanges} for each protein coding transcript representing the
  start and end coordinates of full or partial exons that constitute the
  5' untranslated region of the transcript.

  For \code{threeUTRsByTranscript}: a \code{GRangesList} with
  \code{GRanges} for each protein coding transcript representing the
  start and end coordinates of full or partial exons that constitute the
  3' untranslated region of the transcript.

}
\note{
  While it is possible to request values from a column \code{"tx_name"}
  (with the \code{columns} argument), no such column is present in the
  database. The returned values correspond to the ID of the transcripts.
}
\author{
  Johannes Rainer, Tim Triche
}
\seealso{
  \code{\link{supportedFilters}} to get an overview of supported filters.
  \code{\link{makeEnsembldbPackage}},
  \code{\link{listColumns}}, \code{\link{lengthOf}}
}
\examples{

library(EnsDb.Hsapiens.v75)
edb <- EnsDb.Hsapiens.v75

######   genes
##
## Get all genes encoded on chromosome Y
AllY <- genes(edb, filter = SeqNameFilter("Y"))
AllY

## Return the result as a DataFrame; also, we use a filter expression here
## to define which features to extract from the database.
AllY.granges <- genes(edb,
                      filter = ~ seq_name == "Y",
                      return.type="DataFrame")
AllY.granges

## Include all transcripts of the gene and their chromosomal
## coordinates, sort by chrom start of transcripts and return as
## GRanges.
AllY.granges.tx <- genes(edb,
                         filter = SeqNameFilter("Y"),
                         columns = c("gene_id", "seq_name",
                                     "seq_strand", "tx_id", "tx_biotype",
                                     "tx_seq_start", "tx_seq_end"),
                         order.by = "tx_seq_start")
AllY.granges.tx



######   transcripts
##
## Get all transcripts of a gene
Tx <- transcripts(edb,
                  filter = GeneIdFilter("ENSG00000184895"),
                  order.by = "tx_seq_start")
Tx

## Get all transcripts of two genes along with some information on the
## gene and transcript
Tx <- transcripts(edb,
                  filter = GeneIdFilter(c("ENSG00000184895",
                                          "ENSG00000092377")),
                  columns = c("gene_id", "gene_seq_start", "gene_seq_end",
                              "gene_biotype", "tx_biotype"))
Tx

######   promoters
##
## Get the bona-fide promoters (2k up- to 200nt downstream of TSS)
promoters(edb, filter = GeneIdFilter(c("ENSG00000184895",
                                       "ENSG00000092377")))

######   exons
##
## Get all exons of protein coding transcript for the gene ENSG00000184895
Exon <- exons(edb,
              filter = ~ gene_id == "ENSG00000184895" &
                  tx_biotype == "protein_coding",
              columns = c("gene_id", "gene_seq_start", "gene_seq_end",
                          "tx_biotype", "gene_biotype"))
Exon



#####    exonsBy
##
## Get all exons for transcripts encoded on chromosomes X and Y.
ETx <- exonsBy(edb, by = "tx",
               filter = SeqNameFilter(c("X", "Y")))
ETx
## Get all exons for genes encoded on chromosome 1 to 22, X and Y and
## include additional annotation columns in the result
EGenes <- exonsBy(edb, by = "gene",
                  filter = SeqNameFilter(c("X", "Y")),
                  columns = c("gene_biotype", "gene_name"))
EGenes

## Note that this might also contain "LRG" genes.
length(grep(names(EGenes), pattern="LRG"))

## to fetch just Ensemblgenes, use an GeneIdFilter with value
## "ENS%" and condition "like"
eg <- exonsBy(edb, by = "gene",
              filter = AnnotationFilterList(SeqNameFilter(c("X", "Y")),
                                            GeneIdFilter("ENS", "startsWith")),
              columns = c("gene_biotype", "gene_name"))
eg
length(grep(names(eg), pattern="LRG"))

#####    transcriptsBy
##
TGenes <- transcriptsBy(edb, by = "gene",
                        filter = SeqNameFilter(c("X", "Y")))
TGenes

## convert this to a SAF formatted data.frame that can be used by the
## featureCounts function from the Rsubreader package.
head(toSAF(TGenes))


#####   transcriptsByOverlaps
##
ir <- IRanges(start = c(2654890, 2709520, 28111770),
              end = c(2654900, 2709550, 28111790))
gr <- GRanges(rep("Y", length(ir)), ir)

## Retrieve all transcripts overlapping any of the regions.
txs <- transcriptsByOverlaps(edb, gr)
txs

## Alternatively, use a GRangesFilter
grf <- GRangesFilter(gr, type = "any")
txs <- transcripts(edb, filter = grf)
txs


####    cdsBy
## Get the coding region for all transcripts on chromosome Y.
## Specifying also additional annotation columns (in addition to the default
## exon_id and exon_rank).
cds <- cdsBy(edb, by = "tx", filter = SeqNameFilter("Y"),
             columns = c("tx_biotype", "gene_name"))

####    the 5' untranslated regions:
fUTRs <- fiveUTRsByTranscript(edb, filter = SeqNameFilter("Y"))

####    the 3' untranslated regions with additional column gene_name.
tUTRs <- threeUTRsByTranscript(edb, filter = SeqNameFilter("Y"),
                               columns = "gene_name")


}
\keyword{classes}





