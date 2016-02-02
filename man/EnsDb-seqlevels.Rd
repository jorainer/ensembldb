\name{seqlevelsStyle}
\Rdversion{1.1}
\alias{seqlevelsStyle}
\alias{seqlevelsStyle,EnsDb-method}
\alias{seqlevelsStyle<-}
\alias{seqlevelsStyle<-,EnsDb-method}
\alias{supportedSeqlevelsStyles}
\alias{supportedSeqlevelsStyles,EnsDb-method}

\title{Support for other than Ensembl seqlevel style}
\description{
  The methods and functions on this help page allow to integrate
  \code{EnsDb} objects and the annotations they provide with other
  Bioconductor annotation packages that base on chromosome names
  (seqlevels) that are different from those defined by Ensembl.
}
\usage{

\S4method{seqlevelsStyle}{EnsDb}(x)

\S4method{seqlevelsStyle}{EnsDb}(x) <- value

\S4method{supportedSeqlevelsStyles}{EnsDb}(x)

}
\arguments{

  (In alphabetic order)

  \item{value}{
    For \code{seqlevelsStyle<-}: a character string specifying the
    seqlevels style that should be set. Use the
    \code{supportedSeqlevelsStyle} to list all available and supported
    seqlevel styles.
  }

  \item{x}{
    An \code{EnsDb} instance.
  }

}
\section{Methods and Functions}{
  \describe{

    \item{seqlevelsStyle}{
      Get the style of the seqlevels in which results returned from the
      \code{EnsDb} object are encoded. By default, and internally,
      seqnames as provided by Ensembl are used.

      The method returns a character string specifying the currently used
      seqlevelstyle.
    }

    \item{seqlevelsStyle<-}{
      Change the style of the seqlevels in which results returned from
      the \code{EnsDb} object are encoded. Changing the seqlevels helps
      integrating annotations from \code{EnsDb} objects e.g. with
      annotations from packages that base on UCSC annotations.
    }

    \item{supportedSeqlevelsStyles}{
      Lists all seqlevel styles for which mappings between seqlevel
      styles are available in the \code{GenomeInfoDb} package.

      The method returns a character vector with supported seqlevel
      styles for the organism of the \code{EnsDb} object.
    }

  }
}

\note{
  The mapping between different seqname styles is performed based on
  data provided by the \code{GenomeInfoDb} package. Note that in most
  instances no mapping is provided for seqnames other than for primary
  chromosomes. By default functions from the \code{ensembldb} package
  return the \emph{original} seqname is in such cases. This behaviour
  can be changed with the \code{ensembldb.seqnameNotFound} global
  option. For the special keyword \code{"ORIGINAL"} (the default), the
  original seqnames are returned, for \code{"MISSING"} an error is
  thrown if a seqname can not be mapped. In all other cases, the value
  of the option is returned as seqname if no mapping is available
  (e.g. setting \code{options(ensembldb.seqnameNotFound=NA)} returns an
  \code{NA} if the seqname is not mappable).
}

\value{
  For \code{seqlevelsStyle}: see method description above.

  For \code{supportedSeqlevelsStyles}: see method description above.
}
\author{
  Johannes Rainer
}
\seealso{
  \code{\linkS4class{EnsDb}}
  \code{\link{transcripts}}
}
\examples{

library(EnsDb.Hsapiens.v75)
edb <- EnsDb.Hsapiens.v75

## Get the internal, default seqlevel style.
seqlevelsStyle(edb)

## Get the seqlevels from the database.
seqlevels(edb)

## Get all supported mappings for the organism of the EnsDb.
supportedSeqlevelsStyles(edb)

## Change the seqlevels to UCSC style.
seqlevelsStyle(edb) <- "UCSC"
seqlevels(edb)

## Change the option ensembldb.seqnameNotFound to return NA in case
## the seqname can not be mapped form Ensembl to UCSC.
options(ensembldb.seqnameNotFound=NA)

seqlevels(edb)

## Restoring the original setting.
options(ensembldb.seqnameNotFound="ORIGINAL")


## Integrate Ensembl based annotations with a BSgenome package that is based on
## UCSC style seqnames.
library(BSgenome.Hsapiens.UCSC.hg19)
bsg <- BSgenome.Hsapiens.UCSC.hg19

## Get the genome version
unique(genome(bsg))
unique(genome(edb))
## Although differently named, both represent genome build GRCh37.

## Extract the full transcript sequences of all lincRNAs encoded on chromsome Y.
yTxSeqs <- extractTranscriptSeqs(bsg, exonsBy(edb, "tx",
                                              filter=list(SeqnameFilter("chrY"),
                                                          GenebiotypeFilter("lincRNA"))))
yTxSeqs

}
\keyword{classes}





