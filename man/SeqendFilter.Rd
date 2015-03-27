\name{SeqendFilter}
\alias{EntrezidFilter}
\alias{GeneidFilter}
\alias{GenenameFilter}
\alias{GenebiotypeFilter}
\alias{TxidFilter}
\alias{TxbiotypeFilter}
\alias{ExonidFilter}
\alias{SeqnameFilter}
\alias{SeqstrandFilter}
\alias{SeqstartFilter}
\alias{SeqendFilter}
\title{
  Constructor functions for filter objects
}
\description{
  These functions allow to create filter objects that can be used to
  retrieve specific elements from the annotation database.
}
\usage{
EntrezidFilter(value, condition = "=")
GeneidFilter(value, condition = "=")
GenenameFilter(value, condition = "=")
GenebiotypeFilter(value, condition = "=")
TxidFilter(value, condition = "=")
TxbiotypeFilter(value, condition = "=")
ExonidFilter(value, condition = "=")
SeqnameFilter(value, condition = "=")
SeqstrandFilter(value, condition = "=")
SeqstartFilter(value, condition = "=", feature = "gene")
SeqendFilter(value, condition = "=", feature = "gene")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{value}{
    The filter value, e.g., for \code{GeneidFilter} the id of the gene
    for which the data should be retrieved. For character values (all
    filters except \code{SeqstartFilter} and \code{SeqendFilter}) also a
    character vector of values is allowed. Allowed values for
    \code{SeqstrandFilter} are: \code{"+"}, \code{"-"}, \code{"1"} or \code{"-1"}.
  }
  \item{condition}{
    The condition to be used in the comparison. For character values
    \code{"="}, \code{"in"} and \code{"like"} are allowed, for numeric values
    (\code{SeqstartFilter} and \code{SeqendFilter}) \code{"="},
    \code{">"}, \code{">="}, \code{"<"} and \code{"<="}. Note that for
    \code{"like"} \code{value} should be a SQL pattern (e.g. \code{"ENS\%"}).
  }
  \item{feature}{
    For \code{SeqstartFilter} and \code{SeqendFilter}: the chromosomal
    position of which features should be used in the filter (either
    \code{"gene"}, \code{"transcript"} or \code{"exon"})?
  }
}
\details{
  \describe{
    \item{\code{EntrezidFilter}}{
      Filter results based on the NCBI Entrezgene ID of the genes.
    }
    \item{\code{GeneidFilter}}{
      Filter results based on Ensembl gene IDs.
    }
    \item{\code{GenenameFilter}}{
      Filter results based on gene names (gene symbols).
    }
    \item{\code{GenebiotypeFilter}}{
      Filter results based on the biotype of the genes. For a complete
      list of available gene biotypes use the
      \code{\link{listGenebiotypes}} method.
    }
    \item{\code{TxidFilter}}{
      Filter results based on the Ensembl transcript IDs.
    }
    \item{\code{TxbiotypeFilter}}{
      Filter results based on the biotype of the transcripts. For a
      complete list of available transcript biotypes use the
      \code{\link{listTxbiotypes}} method.
    }
    \item{\code{ExonidFilter}}{
      Filter based on the Ensembl exon ID.
    }
    \item{\code{SeqnameFilter}}{
      Filter results based on the name of the sequence the features are encoded.
    }
    \item{\code{SeqstrandFilter}}{
      Filter results based on the strand on which the features are encoded.
    }
    \item{\code{SeqstartFilter}}{
      Filter results based on the (chromosomal) start coordinate of the
      features (exons, genes or transcripts).
    }
    \item{\code{SeqendFilter}}{
      Filter results based on the (chromosomal) end coordinates.
    }
  }
}
\value{
  Depending on the function called an instance of:
  \code{\linkS4class{EntrezidFilter}},
  \code{\linkS4class{GeneidFilter}},
  \code{\linkS4class{GenenameFilter}},
  \code{\linkS4class{GenebiotypeFilter}},
  \code{\linkS4class{TxidFilter}},
  \code{\linkS4class{TxbiotypeFilter}},
  \code{\linkS4class{ExonidFilter}},
  \code{\linkS4class{SeqnameFilter}},
  \code{\linkS4class{SeqstrandFilter}},
  \code{\linkS4class{SeqstartFilter}},
  \code{\linkS4class{SeqendFilter}}
}
\author{
  Johannes Rainer
}
\seealso{
  \code{\linkS4class{EntrezidFilter}},
  \code{\linkS4class{GeneidFilter}},
  \code{\linkS4class{GenenameFilter}},
  \code{\linkS4class{GenebiotypeFilter}},
  \code{\linkS4class{TxidFilter}},
  \code{\linkS4class{TxbiotypeFilter}},
  \code{\linkS4class{ExonidFilter}},
  \code{\linkS4class{SeqnameFilter}},
  \code{\linkS4class{SeqstrandFilter}},
  \code{\linkS4class{SeqstartFilter}},
  \code{\linkS4class{SeqendFilter}}
}
\examples{

## create a filter that could be used to retrieve all informations for
## the respective gene.
Gif <- GeneidFilter("ENSG00000012817")
Gif
## returns the where clause of the SQL querys
where(Gif)

## create a filter for a chromosomal end position of a gene
Sef <- SeqendFilter(100000, condition="<", "gene")
Sef

## To find genes within a certain chromosomal position filters should be
## combined:
Ssf <- SeqstartFilter(10000, condition=">", "gene")
Snf <- SeqnameFilter("2")
## combine the filters
Filter <- list(Ssf, Sef, Snf)

Filter

## generate the where SQL call for these filters:
where(Filter)

}
\keyword{data}
