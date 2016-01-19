\name{getGenomeFaFile}
\Rdversion{1.1}
%\alias{extractTranscriptSeqs}
%\alias{extractTranscriptSeqs,ANY-method}
%\alias{extractTranscriptSeqs,ANY,ANY}
%\alias{extractTranscriptSeqs,ANY,EnsDb-method}
\alias{getGenomeFaFile}
\alias{getGenomeFaFile,EnsDb-method}

\title{Functionality related to DNA/RNA sequences}
\description{
  Utility functions related to RNA/DNA sequences, such as extracting
  RNA/DNA sequences for features defined in \code{Ensb}.
}
\usage{

\S4method{getGenomeFaFile}{EnsDb}(x, pattern="dna.toplevel.fa")

%\S4method{extractTranscriptSeqs}{ANY,EnsDb}(x, transcripts, filter)

}
\arguments{

  (In alphabetic order)

  %% \item{filter}{
  %%   A filter object extending \code{\linkS4class{BasicFilter}} or a list
  %%   of such object(s) to select specific entries from the database (see
  %%   examples below).
  %% }

  \item{pattern}{
    For method \code{getGenomeFaFile}: the pattern to be used to
    identify the fasta file representing genomic DNA sequence.
  }

  %% \item{transcripts}{
  %%   For \code{extractTranscriptSeqs}: the \code{EnsDb} object from which
  %%   the transcript definitions should be extracted.
  %% }

  \item{x}{
    %% For \code{extractTranscriptSeqs}: An object representing a single
    %% chromosome or a collection of chromosomes. Refer to the help of the
    %% \code{\link[GenomicFeatures]{extractTranscriptSeqs}} method in
    %% \code{GenomicFeatures} package for more information.
    For all other methods an \code{EnsDb} instance.
  }

}
\section{Methods and Functions}{
  \describe{
    %% \item{extractTranscriptSeqs}{
    %%   Extract transcript sequences. This method adapts the
    %%   \code{\link[GenomicFeatures]{extractTranscriptSeqs}} from the
    %%   \code{GenomicFeatures} package to allow the usage of filters to
    %%   specify the transcripts from which the sequence should be
    %%   extracted.
    %% }

    \item{getGenomeFaFile}{
      Returns a \code{\link[Rsamtools]{FaFile-class}} (defined in
      \code{Rsamtools}) with the genomic sequence of the genome build
      matching the Ensembl version of the \code{EnsDb} object.
      The file is retrieved using the \code{AnnotationHub} package,
      thus, at least for the first invocation, an internet connection is
      required to locate and download the file; subsequent calls will
      load the cached file instead.
      If no fasta file for the actual Ensembl version is available the
      function tries to identify a file matchint the species and genome
      build version of the closest Ensembl release and returns that
      instead.
      See the vignette for an example to work with such files.
    }

  }
}

\value{
  For \code{getGenomeFaFile}: a \code{\link[Rsamtools]{FaFile-class}}
  object with the genomic DNA sequence.

  %% For \code{extractTranscriptSeqs}: A \code{DNAStringSet} object
  %% parallel to \code{transcripts} (i.e. the i-th element in it is the
  %% sequence of the i-th transcript in \code{transcripts}).
}
\author{
  Johannes Rainer
}
\seealso{
  \code{\linkS4class{BasicFilter}}
  \code{\link{transcripts}}
  \code{\link{exonsBy}}
}
\examples{

## Loading an EnsDb for Ensembl version 75 (genome GRCh37):
library(EnsDb.Hsapiens.v75)
edb <- EnsDb.Hsapiens.v75

\dontrun{
    ## Retrieve a FaFile with the gneomic DNA sequence matching the organism,
    ## genome release version and, if possible, the Ensembl version of the
    ## EnsDb object.
    Dna <- getGenomeFaFile(edb)
    ## Extract the transcript sequence for all transcripts encoded on chromosome
    ## Y.
    ##extractTranscriptSeqs(Dna, edb, filter=SeqnameFilter("Y"))

}

}
\keyword{classes}





