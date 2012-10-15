\name{make.isos}
\alias{make.isos}

\title{
Deriving list of m/z isotope differences for input into \code{\link[nontarget]{pattern.search}}.
}

\description{
Deriving list of m/z isotope differences for input into \code{\link[nontarget]{pattern.search}}.
}

\usage{
make.isos(iso_list, elements = c("C", "N", "Cl", "S", "Br", "K"), charges = c(1:2), 
rm.isotopes = c("33S", "36S", "36S", "40K", "40K", "15N"), rm.charges = c(2, 1, 2, 1, 2, 2))
}

\arguments{
  \item{iso_list}{
	List of type \code{\link[nontarget]{iso_list}}
}
  \item{elements}{
	Vector of character strings, elements to be used. Corresponds to character strings in column one of \code{\link[nontarget]{iso_list}}.
}
  \item{charges}{
	Numeric vector, charge states to be used (e.g. c(1,2,3) for positive ionization or c(-1,-2,-3) for negative ionization).
}
  \item{rm.isotopes}{
	Vector of character strings. Individual isotopes of elements to be excluded from output list, corresponds to second column entries in \code{\link[nontarget]{iso_list}}.
}
  \item{rm.charges}{
	Numeric vector of that length of rm.isotopes, charges at which isotopes of elements are to be excluded from output list. 
}
}
\details{
\code{rm.isotopes} and \code{rm.charges} allow to exclude those isotopes of low abundance and/or those having strong m/z overlap with other more abundant isotopes.

Charge state and z should not be confused.
}
\value{

List of type iso with 5 entries

\item{iso[[1]]}{\code{list of isotopes}.}
\item{iso[[2]]}{\code{list of isotope masses}.}
\item{iso[[3]]}{\code{charges}.}
\item{iso[[4]]}{\code{number of isotope m/z}.}
\item{iso[[5]]}{\code{elements}.}

}

\author{
Martin Loos
}

\seealso{
	\code{\link[nontarget]{iso_list}}
	\code{\link[nontarget]{pattern.search}}
}
\examples{
data(iso_list)
iso<-make.isos(iso_list,charges=c(1),
	elements=c("C","Cl","Br","N","S","K"),
	rm.isotopes=c("33S","36S","40K"),rm.charges=c(1,1,1))
}

