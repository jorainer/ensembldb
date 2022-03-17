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
    seqlevel styles. As an alternative, it is also possible to submit a
    \code{data.frame} with custom mapping. This needs to have two
    columns, one of them being called \code{"Ensembl"} with the original
    chromosome names and another column with the new names.
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
      annotations from packages that base on UCSC annotations. The
      function also supports using/defining custom mappings by
      submitting a mapping \code{data.frame} (see examples below).
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

library(EnsDb.Hsapiens.v86)
edb <- EnsDb.Hsapiens.v86

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
options(ensembldb.seqnameNotFound = NA)

seqlevels(edb)

## Defining custom mapping for chromosome names. The `data.frame` should have
## one column named `"Ensembl"` with the original name and an additional column
## with the new names
mymap <- data.frame(Ensembl = c(4, 7, 9, 10),
                    myway = c("a", "b", "c", "d"))
seqlevelsStyle(edb) <- mymap
seqlevels(edb)

## This allows us also to rename individual chromosomes but keeping all
## original names for the others.
options(ensembldb.seqnameNotFound = "ORIGINAL")
seqlevels(edb)

}
\keyword{classes}
