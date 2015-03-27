\name{GeneidFilter-class}
\Rdversion{1.1}
\docType{class}
\alias{BasicFilter-class}
\alias{EntrezidFilter-class}
\alias{GeneidFilter-class}
\alias{GenebiotypeFilter-class}
\alias{GenenameFilter-class}
\alias{TxidFilter-class}
\alias{TxbiotypeFilter-class}
\alias{ExonidFilter-class}
\alias{SeqnameFilter-class}
\alias{SeqstrandFilter-class}
\alias{SeqstartFilter-class}
\alias{SeqendFilter-class}
\alias{column,EntrezidFilter,missing,missing-method}
\alias{column,GeneidFilter,missing,missing-method}
\alias{column,GenenameFilter,missing,missing-method}
\alias{column,GenebiotypeFilter,missing,missing-method}
\alias{column,TxidFilter,missing,missing-method}
\alias{column,TxbiotypeFilter,missing,missing-method}
\alias{column,ExonidFilter,missing,missing-method}
\alias{column,SeqnameFilter,missing,missing-method}
\alias{column,SeqstrandFilter,missing,missing-method}
\alias{column,SeqstartFilter,missing,missing-method}
\alias{column,SeqendFilter,missing,missing-method}
\alias{where,EntrezidFilter,missing,missing-method}
\alias{where,GeneidFilter,missing,missing-method}
\alias{where,GenenameFilter,missing,missing-method}
\alias{where,GenebiotypeFilter,missing,missing-method}
\alias{where,TxidFilter,missing,missing-method}
\alias{where,TxbiotypeFilter,missing,missing-method}
\alias{where,ExonidFilter,missing,missing-method}
\alias{where,SeqnameFilter,missing,missing-method}
\alias{where,SeqstrandFilter,missing,missing-method}
\alias{where,SeqstartFilter,missing,missing-method}
\alias{where,SeqendFilter,missing,missing-method}
% EnsDb, missing
\alias{column,EntrezidFilter,EnsDb,missing-method}
\alias{column,GeneidFilter,EnsDb,missing-method}
\alias{column,GenenameFilter,EnsDb,missing-method}
\alias{column,GenebiotypeFilter,EnsDb,missing-method}
\alias{column,TxidFilter,EnsDb,missing-method}
\alias{column,TxbiotypeFilter,EnsDb,missing-method}
\alias{column,ExonidFilter,EnsDb,missing-method}
\alias{column,SeqnameFilter,EnsDb,missing-method}
\alias{column,SeqstrandFilter,EnsDb,missing-method}
\alias{column,SeqstartFilter,EnsDb,missing-method}
\alias{column,SeqendFilter,EnsDb,missing-method}
\alias{where,EntrezidFilter,EnsDb,missing-method}
\alias{where,GeneidFilter,EnsDb,missing-method}
\alias{where,GenenameFilter,EnsDb,missing-method}
\alias{where,GenebiotypeFilter,EnsDb,missing-method}
\alias{where,TxidFilter,EnsDb,missing-method}
\alias{where,TxbiotypeFilter,EnsDb,missing-method}
\alias{where,ExonidFilter,EnsDb,missing-method}
\alias{where,SeqnameFilter,EnsDb,missing-method}
\alias{where,SeqstrandFilter,EnsDb,missing-method}
\alias{where,SeqstartFilter,EnsDb,missing-method}
\alias{where,SeqendFilter,EnsDb,missing-method}
% EnsDb, character
\alias{column,EntrezidFilter,EnsDb,character-method}
\alias{column,GeneidFilter,EnsDb,character-method}
\alias{column,GenenameFilter,EnsDb,character-method}
\alias{column,GenebiotypeFilter,EnsDb,character-method}
\alias{column,TxidFilter,EnsDb,character-method}
\alias{column,TxbiotypeFilter,EnsDb,character-method}
\alias{column,ExonidFilter,EnsDb,character-method}
\alias{column,SeqnameFilter,EnsDb,character-method}
\alias{column,SeqstrandFilter,EnsDb,character-method}
\alias{column,SeqstartFilter,EnsDb,character-method}
\alias{column,SeqendFilter,EnsDb,character-method}
\alias{where,EntrezidFilter,EnsDb,character-method}
\alias{where,GeneidFilter,EnsDb,character-method}
\alias{where,GenenameFilter,EnsDb,character-method}
\alias{where,GenebiotypeFilter,EnsDb,character-method}
\alias{where,TxidFilter,EnsDb,character-method}
\alias{where,TxbiotypeFilter,EnsDb,character-method}
\alias{where,ExonidFilter,EnsDb,character-method}
\alias{where,SeqnameFilter,EnsDb,character-method}
\alias{where,SeqstrandFilter,EnsDb,character-method}
\alias{where,SeqstartFilter,EnsDb,character-method}
\alias{where,SeqendFilter,EnsDb,character-method}
%
\alias{condition,BasicFilter-method}
\alias{show,BasicFilter-method}
\alias{print,BasicFilter-method}
\alias{where,BasicFilter,missing,missing-method}
\alias{where,BasicFilter,EnsDb,missing-method}
\alias{where,BasicFilter,EnsDb,character-method}
\alias{where,list,EnsDb,character-method}
\alias{where,list,EnsDb,missing-method}
\alias{where,list,missing,missing-method}
\alias{condition}
\alias{column}
\alias{where}

\title{Filter results fetched from the Ensembl database}
\description{
  These classes allow to specify which entries (i.e. genes, transcripts
  or exons) should be retrieved from the database.
}
\section{Objects from the Class}{
  While objects can be created by calls e.g. of the form
  \code{new("GeneidFilter", ...)} users are strongly encouraged to use the
  specific functions: \code{\link{GeneidFilter}}, \code{\link{EntrezidFilter}},
  \code{\link{GenenameFilter}}, \code{\link{GenebiotypeFilter}},
  \code{\link{TxidFilter}}, \code{\link{TxbiotypeFilter}},
  \code{\link{ExonidFilter}}, \code{\link{SeqnameFilter}},
  \code{\link{SeqstrandFilter}}, \code{\link{SeqstartFilter}} and \code{\link{SeqendFilter}}.
  See examples below for usage.
}
\section{Slots}{
  \describe{
    \item{\code{condition}:}{
      Object of class \code{"character"}: can be
      either \code{"="}, \code{"in"} or \code{"like"} to filter on character values
      (e.g. gene id, gene biotype, seqname etc), or \code{"="}, \code{">"}
      or \code{"<"} for numerical values (chromosome/seq
      coordinates). Note that for \code{"like"} \code{value} should be a
      SQL pattern (e.g. \code{"ENS\%"}).
    }

    \item{\code{value}:}{
      Object of class \code{"character"}: the value
      to be used for filtering.
    }

  }
}
\section{Extends}{
  Class \code{\linkS4class{BasicFilter}}, directly.
}
\section{Methods}{
  \describe{
    Note: these methods are applicable to all classes extending the
    \code{BasicFilter} class.

    \item{column}{\code{signature(object = "GeneidFilter", db="EnsDb",
	with.tables="character")}:
      returns the column (attribute name) to be used for the
      filtering. Submitting the \code{db} parameter ensures that
      returned column is valid in the corresponding database schema. The
      optional argument \code{with.tables} allows to specify which in
      which database table the function should look for the
      attribute/column name. By default the method will check all
      database tables.
    }

    \item{column}{\code{signature(object = "GeneidFilter", db="EnsDb",
	with.tables="missing")}:
      returns the column (attribute name) to be used for the
      filtering. Submitting the \code{db} parameter ensures that
      returned column is valid in the corresponding database schema.
    }

    \item{column}{\code{signature(object = "GeneidFilter", db="missing",
	with.tables="missing")}:
      returns the column (table column name) to be used for the
      filtering.
    }

    \item{condition}{\code{signature(x="BasicFilter")}: returns
      the value for the \code{condition} slot.
    }

    \item{where}{\code{signature(object = "GeneidFilter", db="EnsDb",
	with.tables="character")}:
      returns the where clause for the SQL call. Submitting also the
      \code{db} parameter ensures that
      the columns are valid in the corresponding database schema. The
      optional argument \code{with.tables} allows to specify which in
      which database table the function should look for the
      attribute/column name. By default the method will check all
      database tables.
    }

    \item{where}{\code{signature(object = "GeneidFilter", db="EnsDb",
	with.tables="missing")}:
      returns the
      where clause for the SQL call. Submitting also the \code{db}
      parameter ensures that
      the columns are valid in the corresponding database schema.
    }

    \item{where}{\code{signature(object = "GeneidFilter", db="missing",
	with.tables="missing")}:
      returns the where clause for the SQL call.
    }
  }
}
\details{
  \describe{
    \item{\code{ExonidFilter}}{
      Allows to filter based on the (Ensembl) exon identifier.
    }

    \item{\code{EntrezidFilter}}{
      Filter results based on the NCBI Entrezgene identifierts of the
      genes. Use the \code{\link{listGenebiotypes}} method to get a
      complete list of all available gene biotypes.
    }

    \item{\code{GenebiotypeFilter}}{
      Filter results based on the gene biotype as defined in the Ensembl
      database.

    }

    \item{\code{GeneidFilter}}{
      Filter results based on the Ensembl gene identifiers.
    }

    \item{\code{GenenameFilter}}{
      Allows to filter on the gene names (symbols) of the genes.
    }

    \item{\code{SeqendFilter}}{
      Filter based on the chromosomal end coordinate of the exons,
      transcripts or genes.
    }

    \item{\code{SeqnameFilter}}{
      Filter on the sequence name on which the features are encoded
      (mostly the chromosome names).
    }

    \item{\code{SeqstartFilter}}{
      Filter based on the chromosomal start coordinates of the exons,
      transcripts or genes.
    }

    \item{\code{SeqstrandFilter}}{
      Filter based on the strand on which the features are encoded.
    }

    \item{\code{TxbiotypeFilter}}{
      Filter on the transcript biotype defined in Ensembl. Use the
      \code{\link{listTxbiotypes}} method to get a complete list of all
      available transcript biotypes.
    }

    \item{\code{TxidFilter}}{
      Filter on the Ensembl transcript identifiers.
    }
  }
}
\note{
  The \code{column} and \code{where} methods should be always called
  along with the \code{EnsDb} object, as this ensures that the
  returned column names are valid for the database schema. The optional
  argument \code{with.tables} should on the other hand only be used
  rarely as it is more intended for internal use.

  Note that the database column \code{"entrezid"} queried for
  \code{EntrezidFilter} classes can contain multiple, \code{";"}
  separated, Entrezgene IDs, thus, using this filter at present might
  not return all entries from the database.
}
\author{
  Johannes Rainer
}
\seealso{
  \code{\link{genes}}, \code{\link{transcripts}}, \code{\link{exons}},
  \code{\link{listGenebiotypes}}, \code{\link{listTxbiotypes}}
}
\examples{

## create a filter that could be used to retrieve all informations for
## the respective gene.
Gif <- GeneidFilter("ENSG00000012817")
Gif
## returns the where clause of the SQL querys
where(Gif)

## create a filter for a chromosomal end position of a gene
Sef <- SeqendFilter(10000, condition=">", "gene")
Sef

## for additional examples see the help page of "genes"


}
\keyword{classes}

