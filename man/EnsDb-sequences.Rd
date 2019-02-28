\name{getGenomeFaFile}
\Rdversion{1.1}
%\alias{extractTranscriptSeqs}
%\alias{extractTranscriptSeqs,ANY-method}
%\alias{extractTranscriptSeqs,ANY,ANY}
%\alias{extractTranscriptSeqs,ANY,EnsDb-method}
\alias{getGenomeFaFile}
\alias{getGenomeFaFile,EnsDb-method}
\alias{getGenomeTwoBitFile}
\alias{getGenomeTwoBitFile,EnsDb-method}

\title{Functionality related to DNA/RNA sequences}
\description{
  Utility functions related to RNA/DNA sequences, such as extracting
  RNA/DNA sequences for features defined in \code{Ensb}.
}
\usage{

\S4method{getGenomeFaFile}{EnsDb}(x, pattern="dna.toplevel.fa")

\S4method{getGenomeTwoBitFile}{EnsDb}(x)

%\S4method{extractTranscriptSeqs}{ANY,EnsDb}(x, transcripts, filter)

}
\arguments{

  (In alphabetic order)

  \item{pattern}{
    For method \code{getGenomeFaFile}: the pattern to be used to
    identify the fasta file representing genomic DNA sequence.
  }

  \item{x}{
    An \code{EnsDb} instance.
  }

}
\section{Methods and Functions}{
  \describe{

    \item{getGenomeFaFile}{
      Returns a \code{\link[Rsamtools]{FaFile-class}} (defined in
      \code{Rsamtools}) with the genomic sequence of the genome build
      matching the Ensembl version of the \code{EnsDb} object.
      The file is retrieved using the \code{AnnotationHub} package,
      thus, at least for the first invocation, an internet connection is
      required to locate and download the file; subsequent calls will
      load the cached file instead.
      If no fasta file for the actual Ensembl version is available the
      function tries to identify a file matching the species and genome
      build version of the closest Ensembl release and returns that
      instead.
      See the vignette for an example to work with such files.
    }

    \item{getGenomeTwoBitFile}{
      Returns a \code{\link{TwoBitFile-class}} (defined in the
      \code{rtracklayer} package) with the genomeic sequence of the
      genome build matching the Ensembl version of the \code{EnsDb}
      object. The file is retrieved from \code{AnnotationHub} and hence
      requires (at least for the first query) an active internet
      connection to download the respective resource. If no DNA sequence
      matching the Ensembl version of \code{x} is available, the
      function tries to find the genomic sequence of the best matching
      genome build (closest Ensembl release) and returns that.

      See the \code{ensembldb} vignette for details.
    }
  }
}

\value{
  For \code{getGenomeFaFile}: a \code{\link[Rsamtools]{FaFile-class}}
  object with the genomic DNA sequence.

  For \code{getGenomeTwoBitFile}: a \code{\link{TwoBitFile-class}}
  object with the genome sequence.
}
\author{
  Johannes Rainer
}
\seealso{
  \code{\link{transcripts}}
  \code{\link{exonsBy}}
}
\examples{

## Loading an EnsDb for Ensembl version 86 (genome GRCh38):
library(EnsDb.Hsapiens.v86)
edb <- EnsDb.Hsapiens.v86

\dontrun{
    ## Retrieve a TwoBitFile with the gneomic DNA sequence matching the organism,
    ## genome release version and, if possible, the Ensembl version of the
    ## EnsDb object.
    Dna <- getGenomeTwoBitFile(edb)
    ## Extract the transcript sequence for all transcripts encoded on chromosome
    ## Y.
    ##extractTranscriptSeqs(Dna, edb, filter=SeqNameFilter("Y"))

}

}
\keyword{classes}
