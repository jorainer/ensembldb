\name{EnsDb-class}
\Rdversion{1.1}
\docType{class}
\alias{EnsDb-class}
\alias{EnsDb}
\alias{buildQuery}
\alias{buildQuery,EnsDb-method}
\alias{dbconn}
\alias{dbconn,EnsDb-method}
\alias{ensemblVersion}
\alias{ensemblVersion,EnsDb-method}
\alias{listColumns}
\alias{listColumns,EnsDb-method}
\alias{metadata}
\alias{metadata,EnsDb-method}
\alias{seqinfo}
\alias{seqinfo,EnsDb-method}
\alias{seqlevels}
\alias{seqlevels,EnsDb-method}
\alias{organism}
\alias{organism,EnsDb-method}
\alias{show}
\alias{show,EnsDb-method}
\alias{listGenebiotypes}
\alias{listGenebiotypes,EnsDb-method}
\alias{listTxbiotypes}
\alias{listTxbiotypes,EnsDb-method}
\alias{listTables}
\alias{listTables,EnsDb-method}
\alias{updateEnsDb}
\alias{updateEnsDb,EnsDb-method}

\title{Basic usage of an Ensembl based annotation database}
\description{
  Get some basic information from an Ensembl based annotation package
  generated with \code{\link{makeEnsembldbPackage}}.

}
\section{Objects from the Class}{
  A connection to the respective annotation database is created upon
  loading of an annotation package created with the
  \code{\link{makeEnsembldbPackage}} function. In addition, the
  \code{EnsDb} constructor specifying the SQLite database file can be
  called to generate an instance of the object (see
  \code{\link{makeEnsemblSQLiteFromTables}} for an example).
}
\usage{

\S4method{buildQuery}{EnsDb}(x, columns=c("gene_id", "gene_biotype",
                                    "gene_name"), filter=list(), order.by,
                             order.type="asc", skip.order.check=FALSE)

\S4method{dbconn}{EnsDb}(x)

EnsDb(x)

\S4method{ensemblVersion}{EnsDb}(x)

\S4method{listColumns}{EnsDb}(x, table, skip.keys=TRUE, ...)

\S4method{listGenebiotypes}{EnsDb}(x, ...)

\S4method{listTxbiotypes}{EnsDb}(x, ...)

\S4method{listTables}{EnsDb}(x, ...)

\S4method{metadata}{EnsDb}(x, ...)

\S4method{organism}{EnsDb}(object)

\S4method{seqinfo}{EnsDb}(x)

\S4method{seqlevels}{EnsDb}(x)

\S4method{updateEnsDb}{EnsDb}(x, ...)

}
\arguments{

  (in alphabetic order)

  \item{...}{Additional arguments.
    Not used.
  }

  \item{columns}{
    Columns (attributes) to be retrieved from the database tables. Use the
    \code{listColumns} or \code{listTables} method for a list of
    supported columns.
  }

  \item{filter}{
    list of \code{\linkS4class{BasicFilter}} instance(s) to
    select specific entries from the database (see examples below).
  }

  \item{object}{
    For \code{organism}: an \code{EnsDb} instance.
  }

  \item{order.by}{name of one of the columns above on which the
    results should be sorted.
  }

  \item{order.type}{if the results should be ordered ascending
    (\code{asc}, default) or descending (\code{desc}).
  }

  \item{skip.keys}{
    for \code{listColumns}: whether primary and foreign keys (not
    being e.g. \code{"gene_id"} or alike) should be returned or not. By
    default these will not be returned.
  }

  \item{skip.order.check}{
    if paramter \code{order.by} should be checked for allowed column
    names. If \code{TRUE} the function checks if the provided order
    criteria orders on columns present in the database tables.
  }

  \item{table}{
    For \code{listColumns}: optionally specify the table name for
    which the columns should be returned.
  }

  \item{x}{
    For \code{EnsDb}: the file name of the SQLite database.
  }

}
\section{Slots}{
  \describe{
    \item{ensdb}{
      Object of class \code{"DBIConnection"}: the
      connection to the database.
    }

    \item{tables}{
      Named list of database table columns with the names being the
      database table names. The tables are ordered by their degree,
      i.e. the number of other tables they can be joined with.
    }

    \item{.properties}{
      Internal list storing user-defined properties. Should not be
      directly accessed.
    }
  }
}
\section{Methods and Functions}{
  \describe{

    \item{buildQuery}{
      Helper function building the SQL query to be used to retrieve the
      wanted information. Usually there is no need to call this method.
    }

    \item{dbconn}{
      Returns the connection to the internal SQL database.
    }

    \item{ensemblVersion}{
      Returns the Ensembl version on which the package was built.
    }

    \item{listColumns}{
      Lists all columns of all tables in the database, or, if
      \code{table} is specified, of the respective table.
    }

    \item{listGenebiotypes}{
      Lists all gene biotypes defined in the database.
    }

    \item{listTxbiotypes}{
      Lists all transcript biotypes defined in the database.
    }

    \item{listTables}{
      Returns a named list of database table columns (names of the
      list being the database table names).
    }

    \item{metadata}{
      Returns a \code{data.frame} with the metadata information from the
      database, i.e. informations about the Ensembl version or Genome
      build the database was build upon.
    }

    \item{organism}{
      Returns the organism name (e.g. \code{"homo_sapiens"}).
    }

    \item{seqinfo}{
      Returns the sequence/chromosome information from the database.
    }

    \item{seqlevels}{
      Returns the chromosome/sequence names that are available in the
      database.
    }

    \item{show}{
      Displays some informations from the database.
    }

    \item{updateEnsDb}{
      Updates the \code{EnsDb} object to the most recent implementation.
    }

  }
}
\value{
  \describe{
    \item{For \code{buildQuery}}{
      A character string with the SQL query.
    }

    \item{For \code{connection}}{
      The SQL connection to the RSQLite database.
    }

    \item{For \code{EnsDb}}{
      An \code{EnsDb} instance.
    }

    \item{For \code{lengthOf}}{
      A named integer vector with the length of the genes or transcripts.
    }

    \item{For \code{listColumns}}{
      A character vector with the column names.
    }

    \item{For \code{listGenebiotypes}}{
      A character vector with the biotypes of the genes in the database.
    }

    \item{For \code{listTxbiotypes}}{
      A character vector with the biotypes of the transcripts in the database.
    }

    \item{For \code{listTables}}{
      A list with the names corresponding to the database table names
      and the elements being the attribute (column) names of the table.
    }

    \item{For \code{metadata}}{
      A \code{data.frame}.
    }

    \item{For \code{organism}}{
      A character string.
    }

    \item{For \code{seqinfo}}{
      A \code{Seqinfo} class.
    }

    \item{For \code{updateEnsDb}}{
      A \code{EnsDb} object.
    }
  }
}
\author{
  Johannes Rainer
}
\seealso{
  \code{\link{makeEnsembldbPackage}}, \code{\linkS4class{BasicFilter}},
      \code{\link{exonsBy}}, \code{\link{genes}},
      \code{\link{transcripts}},
      \code{\link{makeEnsemblSQLiteFromTables}}
}
\examples{

library(EnsDb.Hsapiens.v75)

## Display some information:
EnsDb.Hsapiens.v75

## Show the tables along with its columns
listTables(EnsDb.Hsapiens.v75)

## For what species is this database?
organism(EnsDb.Hsapiens.v75)

## What Ensembl version if the database based on?
ensemblVersion(EnsDb.Hsapiens.v75)

## Get some more information from the database
metadata(EnsDb.Hsapiens.v75)

## Get all the sequence names.
seqlevels(EnsDb.Hsapiens.v75)

######    buildQuery
##
## Join tables gene and transcript and return gene_id and tx_id
buildQuery(EnsDb.Hsapiens.v75, columns=c("gene_id", "tx_id"))


## Get all exon_ids and transcript ids of genes encoded on chromosome Y.
buildQuery(EnsDb.Hsapiens.v75, columns=c("exon_id", "tx_id"),
           filter=list(SeqnameFilter( "Y")))

## List all available gene biotypes from the database:
listGenebiotypes(EnsDb.Hsapiens.v75)

## List all available transcript biotypes:
listTxbiotypes(EnsDb.Hsapiens.v75)

## Update the EnsDb; this is in most instances not necessary at all.
updateEnsDb(EnsDb.Hsapiens.v75)

}
\keyword{classes}

