\name{lengthOf}
\Rdversion{1.1}
\alias{lengthOf}
\alias{lengthOf,GRangesList-method}
\alias{lengthOf,EnsDb-method}
%\alias{transcriptLengths}
%\alias{transcriptLengths,EnsDb-method}
%\alias{transcriptLengths,TxDb-method}

\title{Calculating lengths of features}
\description{
  These methods allow to calculate the lengths of features (transcripts, genes,
  CDS, 3' or 5' UTRs) defined in an \code{EnsDb} object or database.
}
\usage{

\S4method{lengthOf}{EnsDb}(x, of="gene", filter = AnnotationFilterList())

}
\arguments{

  (In alphabetic order)

  \item{filter}{
    A filter describing which results to retrieve from the database. Can
    be a single object extending
    \code{\link[AnnotationFilter]{AnnotationFilter}}, an
    \code{\link[AnnotationFilter]{AnnotationFilterList}} object
    combining several such objects or a \code{formula} representing a
    filter expression (see examples below or
    \code{\link[AnnotationFilter]{AnnotationFilter}} for more details).
  }

  \item{of}{
    for \code{lengthOf}: whether the length of genes or
    transcripts should be retrieved from the database.
  }

  \item{x}{
    For \code{lengthOf}: either an \code{EnsDb} or a
    \code{GRangesList} object. For all other methods an \code{EnsDb}
    instance.
  }

}
\section{Methods and Functions}{
  \describe{

    \item{lengthOf}{
      Retrieve the length of genes or transcripts from the
      database. The length is the sum of the lengths of all exons of a
      transcript or a gene. In the latter case the exons are first reduced
      so that the length corresponds to the part of the genomic sequence covered by
      the exons.

      Note: in addition to this method, also the
      \code{\link[GenomicFeatures]{transcriptLengths}} function in the
      \code{GenomicFeatures} package can be used.
    }

  }
}

\value{
  For \code{lengthOf}: see method description above.
}
\author{
  Johannes Rainer
}
\seealso{
  \code{\link{exonsBy}}
  \code{\link{transcripts}}
  \code{\link[GenomicFeatures]{transcriptLengths}}
}
\examples{

library(EnsDb.Hsapiens.v75)
edb <- EnsDb.Hsapiens.v75

#####    lengthOf
##
## length of a specific gene.
lengthOf(edb, filter = GeneIdFilter("ENSG00000000003"))

## length of a transcript
lengthOf(edb, of = "tx", filter = TxIdFilter("ENST00000494424"))

## Average length of all protein coding genes encoded on chromosomes X
mean(lengthOf(edb, of = "gene",
              filter = ~ gene_biotype == "protein_coding" &
                  seq_name == "X"))

## Average length of all snoRNAs
mean(lengthOf(edb, of = "gene",
              filter = ~ gene_biotype == "snoRNA" &
                  seq_name == "X"))

##### transcriptLengths
##
## Calculate the length of transcripts encoded on chromosome Y, including
## length of the CDS, 5' and 3' UTR.
len <- transcriptLengths(edb, with.cds_len = TRUE, with.utr5_len = TRUE,
                         with.utr3_len = TRUE, filter = SeqNameFilter("Y"))
head(len)

}
\keyword{classes}





