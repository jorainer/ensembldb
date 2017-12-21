\name{getGeneRegionTrackForGviz}
\Rdversion{1.1}
\alias{getGeneRegionTrackForGviz}
\alias{getGeneRegionTrackForGviz,EnsDb-method}

\title{Utility functions}
\description{
  Utility functions integrating \code{EnsDb} objects with other
  Bioconductor packages.
}
\usage{

\S4method{getGeneRegionTrackForGviz}{EnsDb}(x,
        filter = AnnotationFilterList(), chromosome = NULL,
        start = NULL, end = NULL, featureIs = "gene_biotype")
}
\arguments{

  (In alphabetic order)

  \item{chromosome}{
    For \code{getGeneRegionTrackForGviz}: optional chromosome name to
    restrict the returned entry to a specific chromosome.
  }

  \item{end}{
    For \code{getGeneRegionTrackForGviz}: optional chromosomal end
    coordinate specifying, together with \code{start}, the chromosomal
    region from which features should be retrieved.
  }

  \item{featureIs}{
    For \code{getGeneRegionTrackForGviz}: whether the gene
    (\code{"gene_biotype"}) or the transcript biotype
    (\code{"tx_biotype"}) should be returned in column \code{"feature"}.
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

  \item{start}{
    For \code{getGeneRegionTrackForGviz}: optional chromosomal start
    coordinate specifying, together with \code{end}, the chromosomal
    region from which features should be retrieved.
  }

  \item{x}{
    For \code{toSAF} a \code{GRangesList} object. For all other
    methods an \code{EnsDb} instance.
  }

}
\section{Methods and Functions}{
  \describe{

    \item{getGeneRegionTrackForGviz}{
      Retrieve a \code{GRanges} object with transcript features from the
      \code{EnsDb} that can be used directly in the \code{Gviz} package
      to create a \code{GeneRegionTrack}. Using the \code{filter},
      \code{chromosome}, \code{start} and \code{end} arguments it is
      possible to fetch specific features (e.g. lincRNAs) from the
      database.

      If \code{chromosome}, \code{start} and \code{end} is provided the
      function internally first retrieves all transcripts that have an
      exon or an intron in the specified chromosomal region and
      subsequently fetch all of these transcripts. This ensures that all
      transcripts of the region are returned, even those that have
      \emph{only} an intron in the region.

      The function returns a \code{GRanges} object with additional
      annotation columns \code{"feature"}, \code{"gene"}, \code{"exon"},
      \code{"exon_rank"}, \code{"trancript"}, \code{"symbol"} specifying
      the feature type (either gene or transcript biotype), the
      (Ensembl) gene ID, the exon ID, the rank/index of the exon in the
      transcript, the transcript ID and the gene symbol/name.
    }
  }
}

\value{
  For \code{getGeneRegionTrackForGviz}: see method description above.
}
\author{
  Johannes Rainer
}
\seealso{
  \code{\link{transcripts}}
}
\examples{

library(EnsDb.Hsapiens.v86)
edb <- EnsDb.Hsapiens.v86
######   getGeneRegionTrackForGviz
##
## Get all genes encoded on chromosome Y in the specifyed region.
AllY <- getGeneRegionTrackForGviz(edb, chromosome = "Y", start = 5131959,
                                  end = 7131959)
## We could plot this now using plotTracks(GeneRegionTrack(AllY))


}
\keyword{classes}





