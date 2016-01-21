\name{lengthOf}
\Rdversion{1.1}
\alias{lengthOf}
\alias{lengthOf,GRangesList-method}
\alias{lengthOf,EnsDb-method}
\alias{transcriptLengths}
\alias{transcriptLengths,EnsDb-method}
\alias{transcriptLengths,TxDb-method}

\title{Calculating lengths of features}
\description{
  These methods allow to calculate the lengths of features (transcripts, genes,
  CDS, 3' or 5' UTRs) defined in an \code{EnsDb} object or database.
}
\usage{

\S4method{lengthOf}{EnsDb}(x, of="gene", filter=list())

\S4method{transcriptLengths}{EnsDb}(x, with.cds_len=FALSE,
                                    with.utr5_len=FALSE,
                                    with.utr3_len=FALSE, filter=list())

}
\arguments{

  (In alphabetic order)

  \item{filter}{
    list of \code{\linkS4class{BasicFilter}} instance(s) to
    select specific entries from the database (see examples below).
  }

  \item{of}{
    for \code{lengthOf}: whether the length of genes or
    transcripts should be retrieved from the database.
  }

  \item{with.cds_len}{
    If a column with the length of the CDS should be included in the
    returned \code{data.frame}.
  }

  \item{with.utr3_len}{
    If a column with the length of the 3' UTR should be included in the
    returned \code{data.frame}.
  }

  \item{with.utr5_len}{
    If a column with the length of the 5' UTR should be included in the
    returned \code{data.frame}.
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
    }

    \item{transcriptLengths}{
      Similar to the \code{\link[GenomicFeatures]{transcriptLengths}}
      function in the \code{GenomicFeatures} package, this function
      calculates the length of transcripts and other summary statistics
      (number of exons and, depending on the input parameter, also
      length of the coding region or the 5' or 3' UTR). The ordering of
      the transcripts returned by this method is the same as the
      ordering of transcripts returned by the method
      \code{\link{transcripts}}.

      The method returns a \code{data.frame} with columns
      \code{"tx_id"}, \code{"gene_id"}, \code{"nexon"}, \code{"tx_len"}
      with the transcript id, the id of the gene, the number of exons and
      the length of the transcript. Depending on the input arguments the
      \code{data.frame} can also have columns \code{"cds_len"},
      \code{"utr5_len"} and \code{"utr3_len"} with the lengths of the
      CDS, 5' UTR or 3' UTR (or \code{NA} for non protein coding
      transcripts).

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
}
\examples{

library(EnsDb.Hsapiens.v75)
edb <- EnsDb.Hsapiens.v75

#####    lengthOf
##
## length of a specific gene.
lengthOf(edb,
         filter=list(GeneidFilter("ENSG00000000003")))

## length of a transcript
lengthOf(edb, of="tx",
         filter=list(TxidFilter("ENST00000494424")))

## average length of all protein coding genes encoded on chromosomes X
## and Y
mean(lengthOf(edb, of="gene",
              filter=list(GenebiotypeFilter("protein_coding"),
                  SeqnameFilter(c("X", "Y")))))

## average length of all snoRNAs
mean(lengthOf(edb, of="gene",
              filter=list(GenebiotypeFilter("snoRNA"),
                  SeqnameFilter(c("X", "Y")))))

##### transcriptLengths
##
## Calculate the length of transcripts encoded on chromosome Y, including
## length of the CDS, 5' and 3' UTR.
len <- transcriptLengths(edb, with.cds_len=TRUE, with.utr5_len=TRUE,
                         with.utr3_len=TRUE, filter=SeqnameFilter("Y"))
head(len)

}
\keyword{classes}





