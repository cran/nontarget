\name{iso_list}
\alias{iso_list}
\docType{data}
\title{
Isotope list
}
\description{
Listing of isotopes, their masses, abundances and their ratio to C.
}
\usage{data(iso_list)}
\format{
  A data frame with 59 observations on the following 5 variables.
  \describe{
    \item{\code{element}}{Element name}
    \item{\code{isotope}}{Isotope name of an element}
    \item{\code{mass}}{Isotope mass}
    \item{\code{abundance}}{Natural isotope abundance}
    \item{\code{ratioC}}{Isotope ratio to C in 99.99 \% of cases}
  }
}
\details{
The \code{ratioC}-value stems from a database survey conducted by Kind&Fiehn (2007); to disable, set value to 0.
The list serves as input into several package nontarget-functions. The first column of the data frame also 
contains names of specific isotopes used for labeled compounds. 
}
\source{
\url{http://iupac.org/publications/pac/2003/pdf/7506x0683.pdf} 
}
\references{
Kind, T. and Fiehn, O., 2007. Seven golden rules for heuristic filtering  of molecular formulas obtained by 
accurate mass spectrometry. BMC Bioinformatics, 8:105.
}
\keyword{datasets}
