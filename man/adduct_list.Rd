\name{adduct_list}
\alias{adduct_list}
\docType{data}
\title{
Adduct list
}
\description{
List of common adducts observed for ESI-MS measurements in soft positive and negative ionization modes.
}
\usage{data(adduct_list)}
\format{
  A data frame with 47 observations on the following 6 variables.
  \describe{
    \item{\code{Name}}{Adduct name}
    \item{\code{calc}}{Equation for calculating adduct m/z from uncharged non-adduct molecular mass M (m/z = M/z + X)}
    \item{\code{Charge}}{z}
    \item{\code{Mult}}{1/z}
    \item{\code{Mass}}{X}
    \item{\code{Ion_mode}}{Ionization mode (positive or negative)}
  }
}

\seealso{
	\code{\link[nontarget]{adduct.search}}	
}

\source{
\url{http://fiehnlab.ucdavis.edu/staff/kind/Metabolomics/MS-Adduct-Calculator/}
}

\references{
Huang N., Siegel M.M., Kruppa G.H., Laukien F.H., J. Am. Soc. Mass. Spectrom. 1999, 10.
Automation of a Fourier transform ion cyclotron resonance mass spectrometer for acquisition, analysis, 
and e-mailing of high-resolution exact-mass electrospray ionization mass spectral data
}

\examples{
data(adduct_list);
adduct_list
}
\keyword{datasets}
