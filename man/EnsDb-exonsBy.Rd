\name{exonsBy}
\Rdversion{1.1}
\docType{class}
\alias{disjointExons,EnsDb-method}
\alias{exons}
\alias{exons,EnsDb-method}
\alias{exonsBy}
\alias{exonsBy,EnsDb-method}
\alias{genes}
\alias{genes,EnsDb-method}
\alias{getGenomeFaFile}
\alias{getGenomeFaFile,EnsDb-method}
\alias{toSAF}
\alias{toSAF,GRangesList-method}
\alias{transcripts}
\alias{transcripts,EnsDb-method}
\alias{transcriptsBy}
\alias{transcriptsBy,EnsDb-method}
\alias{promoters}
\alias{promoters,EnsDb-method}

\title{Retrieve annotation data from an Ensembl based package}
\description{
  Retrieve gene/transcript/exons annotations stored in an Ensembl based
  database package generated with the \code{\link{makeEnsembldbPackage}}
  function.
}
\usage{

\S4method{exons}{EnsDb}(x, columns=listColumns(x,"exon"),
                        filter, order.by, order.type="asc",
                        return.type="GRanges")

\S4method{exonsBy}{EnsDb}(x, by=c("tx", "gene"),
                          columns=listColumns(x, "exon"), filter)

\S4method{transcripts}{EnsDb}(x, columns=listColumns(x, "tx"),
                              filter, order.by, order.type="asc",
                              return.type="GRanges")

\S4method{transcriptsBy}{EnsDb}(x, by=c("gene", "exon"),
                                columns=listColumns(x, "tx"), filter)

\S4method{promoters}{EnsDb}(x, upstream=2000, downstream=200, ...)

\S4method{genes}{EnsDb}(x, columns=listColumns(x, "gene"), filter,
                        order.by, order.type="asc",
                        return.type="GRanges")

\S4method{disjointExons}{EnsDb}(x, aggregateGenes=FALSE,
                                includeTranscripts=TRUE, filter, ...)

\S4method{toSAF}{GRangesList}(x, ...)

\S4method{getGenomeFaFile}{EnsDb}(x, pattern="dna.toplevel.fa")

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
  }

  \item{downstream}{
    For method \code{promoters}: the number of nucleotides downstream of
    the transcription start site that should be included in the promoter region.
  }

  \item{filter}{
    A filter object extending \code{\linkS4class{BasicFilter}} or a list
    of such object(s) to select specific entries from the database (see examples below).
  }

  \item{includeTranscripts}{
    For \code{disjointExons}: When \code{TRUE} (default) a
    \code{tx_name} metadata column is included that lists all transcript
    IDs that overlap the exon fragment. Note: this is different to the
    \code{\link[GenomicFeatures]{disjointExons}} function in the
    \code{GenomicFeatures} package, that lists the transcript names, not
    IDs.
  }

  \item{order.by}{
    Name of one of the columns above on which the
    results should be sorted.
  }

  \item{order.type}{
    If the results should be ordered ascending
    (\code{asc}, default) or descending (\code{desc}).
  }

  \item{pattern}{
    For method \code{getGenomeFaFile}: the pattern to be used to
    identify the fasta file representing genomic DNA sequence.
  }

  \item{return.type}{
    Type of the returned object. Can be either
    \code{"data.frame"}, \code{"DataFrame"} or \code{"GRanges"}. In the latter case the return
    object will be a \code{GRanges} object with the GRanges specifying the
    chromosomal start and end coordinates of the feature (gene,
    transcript or exon, depending whether \code{genes},
    \code{transcripts} or \code{exons} was called). All additional
    columns are added as metadata columns to the GRanges object.
  }

  \item{upstream}{
    For method \code{promoters}: the number of nucleotides upstream of
    the transcription start site that should be included in the promoter region.
  }

  \item{x}{
    For \code{toSAF} a \code{GRangesList} object. For all other
    methods an \code{EnsDb} instance.
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

    \item{promoters}{
      Retrieve promoter information from the database. Additional
      columns from genes or exons associated with the promoters can be specified
      and are added to the respective promoter annotation.
    }
    \item{genes}{
      Retrieve gene information from the database. Additional
      columns from transcripts or exons associated with the genes can be specified
      and are added to the respective gene annotation.
    }

    \item{disjointExons}{
      This method is identical to
      \code{\link[GenomicFeatures]{disjointExons}} defined in the
      \code{GenomicFeatures} package. It creates a \code{GRanges} of
      non-overlapping exon parts with metadata columns of \code{gene_id}
      and \code{exonic_part}. Exon parts that overlap more than one gene
      can be dropped with \code{aggregateGenes=FALSE}.
    }

    \item{toSAF}{
      Reformats a \code{GRangesList} object into a
      \code{data.frame} corresponding to a standard SAF (Simplified
      Annotation Format) file (i.e. with column names \code{"GeneID"},
      \code{"Chr"}, \code{"Start"}, \code{"End"} and
      \code{"Strand"}). Note: this method makes only sense on a
      \code{GRangesList} that groups features (exons, transcripts) by gene.
    }

    \item{getGenomeFaFile}{
      Returns a \code{\link[Rsamtools]{FaFile-class}} (defined in
      \code{Rsamtools}) with the genomic sequence of the genome build
      matching the Ensembl version of the \code{EnsDb} object. The file
      is retrieved using the \code{AnnotationHub} package. See the
      vignette for an example to work with such files.
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
    \item{entrezid}{the NCBI Entrezgene ID of the gene; note that this
      can also be a \code{";"} separated list of IDs for Ensembl genes
      mapped to more than one Entrezgene.}
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

  Also, the vignette provides examples on how to retrieve sequences for
  genes/transcripts/exons.
}
\note{
  Ensembl defines genes not only on standard chromosomes, but also on
  patched chromosomes and chromosome variants. Thus it might be
  advisable to restrict the queries to just those chromosomes of
  interest (e.g. by specifying a \code{SeqnameFilter(c(1:22, "X", "Y"))}).
  In addition, also so called LRG genes (Locus Reference Genomic) are defined in
  Ensembl. Their gene id starts with LRG instead of ENS for Ensembl
  genes, thus, a filter can be applied to specifically select those
  genes or exclude those genes (see examples below).
}

\value{
  For \code{exons}, \code{transcripts} and \code{genes},
  a \code{data.frame}, \code{DataFrame}
  or a \code{GRanges}, depending on the value of the
  \code{return.type} parameter. The result
  is ordered as specified by the parameter \code{order.by}, NOT by any
  ordering of values in eventually submitted filter objects.

  For \code{exonsBy}, \code{transcriptsBy}:
  a \code{GRangesList}, depending on the value of the
  \code{return.type} parameter. The results are ordered by the value of the
  \code{by} parameter.

  For \code{toSAF}: a \code{data.frame} with column names
  \code{"GeneID"} (the group name from the \code{GRangesList}, i.e. the
  ID by which the \code{GRanges} are split), \code{"Chr"} (the seqnames
  from the \code{GRanges}), \code{"Start"} (the start coordinate),
  \code{"End"} (the end coordinate) and \code{"Strand"} (the strand).

  For \code{disjointExons}: a \code{GRanges} of non-overlapping exon
  parts.

  For \code{getGenomeFaFile}: a \code{\link[Rsamtools]{FaFile-class}}
  object with the genomic DNA sequence.
}
\author{
  Johannes Rainer, Tim Triche
}
\seealso{
  \code{\link{makeEnsembldbPackage}}, \code{\linkS4class{BasicFilter}},
      \code{\link{listColumns}}, \code{\link{lengthOf}}
}
\examples{

library(EnsDb.Hsapiens.v75)

######   genes
##
## get all genes endcoded on chromosome Y
AllY <- genes(EnsDb.Hsapiens.v75, filter=SeqnameFilter("Y"))
AllY

## return result as DataFrame.
AllY.granges <- genes(EnsDb.Hsapiens.v75,
                      filter=SeqnameFilter("Y"),
                      return.type="DataFrame")
AllY.granges

## include all transcripts of the gene and their chromosomal
## coordinates, sort by chrom start of transcripts and return as
## GRanges.
AllY.granges.tx <- genes(EnsDb.Hsapiens.v75,
                         filter=SeqnameFilter("Y"),
                         columns=c("gene_id", "seq_name",
                             "seq_strand", "tx_id", "tx_biotype",
                             "tx_seq_start", "tx_seq_end"),
                         order.by="tx_seq_start")
AllY.granges.tx



######   transcripts
##
## get all transcripts of a gene
Tx <- transcripts(EnsDb.Hsapiens.v75,
                  filter=GeneidFilter("ENSG00000184895"),
                  order.by="tx_seq_start")
Tx

## get all transcripts of two genes along with some information on the
## gene and transcript
Tx <- transcripts(EnsDb.Hsapiens.v75,
                  filter=GeneidFilter(c("ENSG00000184895",
                      "ENSG00000092377")),
                      columns=c("gene_id", "gene_seq_start",
                          "gene_seq_end", "gene_biotype", "tx_biotype"))
Tx

######   promoters
##
## get the bona-fide promoters (2k up- to 200nt downstream of TSS)
promoters(EnsDb.Hsapiens.v75, filter=GeneidFilter(c("ENSG00000184895",
                                                    "ENSG00000092377")))

######   exons
##
## get all exons of the provided genes
Exon <- exons(EnsDb.Hsapiens.v75,
              filter=GeneidFilter(c("ENSG00000184895",
                  "ENSG00000092377")),
              order.by="exon_seq_start",
              columns=c( "gene_id", "gene_seq_start",
                  "gene_seq_end", "gene_biotype"))
Exon



#####    exonsBy
##
## get all exons for transcripts encoded on chromosomes X and Y.
ETx <- exonsBy(EnsDb.Hsapiens.v75, by="tx",
               filter=SeqnameFilter(c("X", "Y")))
ETx
## get all exons for genes encoded on chromosome 1 to 22, X and Y and
## include additional annotation columns in the result
EGenes <- exonsBy(EnsDb.Hsapiens.v75, by="gene",
                  filter=SeqnameFilter(c("X", "Y")),
                  columns=c("gene_biotype", "gene_name"))
EGenes

## Note that this might also contain "LRG" genes.
length(grep(names(EGenes), pattern="LRG"))

## to fetch just Ensemblgenes, use an GeneidFilter with value
## "ENS%" and condition "like"


#####    transcriptsBy
##
TGenes <- transcriptsBy(EnsDb.Hsapiens.v75, by="gene",
                        filter=SeqnameFilter(c("X", "Y")))
TGenes

## convert this to a SAF formatted data.frame that can be used by the
## featureCounts function from the Rsubreader package.
head(toSAF(TGenes))

}
\keyword{classes}

