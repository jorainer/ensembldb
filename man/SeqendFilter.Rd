\name{SeqendFilter}
\alias{EntrezidFilter}
\alias{GeneidFilter}
\alias{GenenameFilter}
\alias{GenebiotypeFilter}
\alias{TxidFilter}
\alias{TxbiotypeFilter}
\alias{ExonidFilter}
\alias{ExonrankFilter}
\alias{SeqnameFilter}
\alias{SeqstrandFilter}
\alias{SeqstartFilter}
\alias{SeqendFilter}
\alias{GRangesFilter}
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

GRangesFilter(value, condition="within", feature="gene")

TxidFilter(value, condition = "=")

TxbiotypeFilter(value, condition = "=")

ExonidFilter(value, condition = "=")

ExonrankFilter(value, condition = "=")

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
    \code{SeqstrandFilter} are: \code{"+"}, \code{"-"}, \code{"1"} or
    \code{"-1"}.

    For \code{GRangeFilter} this has to be a \code{GRanges} object.
  }
  \item{condition}{
    The condition to be used in the comparison. For character values
    \code{"="}, \code{"in"} and \code{"like"} are allowed, for numeric values
    (\code{SeqstartFilter} and \code{SeqendFilter}) \code{"="},
    \code{">"}, \code{">="}, \code{"<"} and \code{"<="}. Note that for
    \code{"like"} \code{value} should be a SQL pattern
    (e.g. \code{"ENS\%"}).

    For \code{GRangesFilter}, \code{"within"} and \code{"overlapping"}
    are allowed. See below for details.
  }
  \item{feature}{
    For \code{SeqstartFilter} and \code{SeqendFilter}: the chromosomal
    position of which features should be used in the filter (either
    \code{"gene"}, \code{"transcript"} or \code{"exon"}).

    For \code{GRangesFilter}: the submitted value is overwritten
    internally depending on the called method, i.e. calling \code{genes}
    will set feature to \code{"gene"}, \code{transcripts} to \code{"tx"}
    and \code{exons} to \code{"exon"}.

  }
}
\details{
  \describe{
    \item{EntrezidFilter}{
      Filter results based on the NCBI Entrezgene ID of the genes.
    }
    \item{GeneidFilter}{
      Filter results based on Ensembl gene IDs.
    }
    \item{GenenameFilter}{
      Filter results based on gene names (gene symbols).
    }
    \item{GenebiotypeFilter}{
      Filter results based on the biotype of the genes. For a complete
      list of available gene biotypes use the
      \code{\link{listGenebiotypes}} method.
    }
    \item{GRangesFilter}{
      Allows to fetch features within or overlapping the specified genomic
      region(s)/range(s). This filter takes a \code{GRanges} object as input
      and, if \code{condition="within"} (the default) will restrict
      results to features (genes, transcripts or exons) that are
      completely within the region. Alternatively, by specifying
      \code{condition="overlapping"} it will return all features that
      are partially overlapping with the region, i.e. which start
      coordinate is smaller than the end coordinate of the region and
      which end coordinate is larger than the start coordinate of the
      region. Thus, genes and transcripts that have an intron
      overlapping the region will also be returned.

      Note: if the specified \code{GRanges} object defines multiple
      region, all features within (or overlapping) any of these regions
      are returned.

      See \code{\linkS4class{GRangesFilter}} for more details.
    }
    \item{TxidFilter}{
      Filter results based on the Ensembl transcript IDs.
    }
    \item{TxbiotypeFilter}{
      Filter results based on the biotype of the transcripts. For a
      complete list of available transcript biotypes use the
      \code{\link{listTxbiotypes}} method.
    }
    \item{ExonidFilter}{
      Filter based on the Ensembl exon ID.
    }
    \item{ExonrankFilter}{
      Filter results based on exon ranks (indices) of exons within
      transcripts.
    }
    \item{SeqnameFilter}{
      Filter results based on the name of the sequence the features are
      encoded.
    }
    \item{SeqstrandFilter}{
      Filter results based on the strand on which the features are encoded.
    }
    \item{SeqstartFilter}{
      Filter results based on the (chromosomal) start coordinate of the
      features (exons, genes or transcripts).
    }
    \item{SeqendFilter}{
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
  \code{\linkS4class{GRangesFilter}},
  \code{\linkS4class{TxidFilter}},
  \code{\linkS4class{TxbiotypeFilter}},
  \code{\linkS4class{ExonidFilter}},
  \code{\linkS4class{ExonrankFilter}},
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
  \code{\linkS4class{GRangesFilter}},
  \code{\linkS4class{TxidFilter}},
  \code{\linkS4class{TxbiotypeFilter}},
  \code{\linkS4class{ExonidFilter}},
  \code{\linkS4class{ExonrankFilter}},
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


## Create a GRangesFilter
GRangesFilter(GRanges("X", IRanges(123, 5454)))

## Create a GRangesFilter with multiple ranges
grf <- GRangesFilter(GRanges(c("X", "Y"),
                             IRanges(start=c(123, 900),
                                     end=c(5454, 910))))
## Evaluate the "where" SQL clause that would be applied.
where(grf)
## Change the "condition" of the filter and evaluate the
## where clause again.
condition(grf) <- "overlapping"
where(grf)

}
\keyword{data}
