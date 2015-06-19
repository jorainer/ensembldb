\name{runEnsDbApp}
\alias{runEnsDbApp}
\title{
  Search annotations interactively
}
\description{
	This function starts the interactive \code{EnsDb} shiny web application that
	allows to look up gene/transcript/exon annotations from an \code{EnsDb}
	annotation package installed locally.
}
\usage{

	runEnsDbApp(...)

}
\arguments{

	\item{...}{
		Additional arguments passed to the \code{\link[shiny]{runApp}} function
		from the \code{shiny} package.
	}

}
\details{
	The \code{shiny} based web application allows to look up any annotation
	available in any of the locally installed \code{EnsDb} annotation packages.
}
\value{
	Nothing is returned by the function.
}
\seealso{
  \code{\link{EnsDb}}, \code{\link{genes}}
}
\author{
Johannes Rainer
}
\keyword{data}
\keyword{shiny}

