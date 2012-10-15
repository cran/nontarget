\name{plotall}
\alias{plotall}

\title{
RT vs. m/z scatterplot marking isotope pattern and adduct group peaks.
}
\usage{
plotall(pattern, adduct)
}

\description{
RT vs. m/z scatterplot marking isotope pattern and adduct group peaks.
}

\arguments{
  \item{pattern}{
	List of type pattern produced by \code{\link[nontarget]{pattern.search}}.
	}
  \item{adduct}{
	List of type adduct produced by \code{\link[nontarget]{adduct.search}}.
	}
}
\author{
Martin Loos
}
\seealso{
	\code{\link[nontarget]{pattern.search}}
	\code{\link[nontarget]{adduct.search}}
}
\examples{
\donttest{
plotall(pattern, adduct)
}
}
