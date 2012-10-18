\name{homol.search}
\alias{homol.search}

\title{
Detecting homologue series in HRMS data.
}
\description{
An algorithm for unsupervised detection of homologue series in HRMS data.
}
\usage{
homol.search(peaklist, iso_list, minmz = 3, maxmz = 60, minrt = 0.05, maxrt = 2, 
ppm = TRUE, mztol = 3.5, rttol = 0.4, minlen = 5, plotit = TRUE, restr = FALSE)
}

\arguments{
  \item{peaklist}{
	List of HRMS peaks with three entries for (a) m/z, (b) intensity and (c) retention time, such as \code{\link[nontarget]{peaklist}}.
}
  \item{iso_list}{
	List of type \code{\link[nontarget]{iso_list}}
}
  \item{minmz}{
	Defines the lower limit of the m/z window to look for homologue series peaks in, relative to the m/z of the one peak to search from. Given in absolute m/z values [u]. 
}
  \item{maxmz}{
	Defines the upper limit of the m/z window to look for homologue series peaks in, relative to the m/z of the one peak to search from. Given in absolute m/z values [u]. 
}
  \item{minrt}{
	Defines the lower limit of the retention time (RT) window to look for other homologue peaks in, relative to the RT of the one peak to search from. 
}
  \item{maxrt}{
	Defines the upper limit of the retention time (RT) window to look for other homologue peaks in, relative to the RT of the one peak to search from. 
}
  \item{ppm}{
	Should \code{mztol} be set in ppm (\code{TRUE}) or in absolute m/z [u] (\code{FALSE})?
}
  \item{mztol}{
	m/z tolerance setting: value by which the m/z of a peak may vary from its expected value. If parameter \code{ppm=TRUE} (see below) given in ppm, otherwise, if \code{ppm=FALSE}, in absolute m/z [u]. 
}
  \item{rttol}{
	Retention time (RT) tolerance by which differences in RT between two adjacent peaks of a homologue series are allowed to differ. Units as given in column 3 of peaklist argument, e.g. [min].
}
  \item{minlen}{
	Minimum number of peaks in a homologue series.
}
  \item{plotit}{
	Plot progress report? The algorithm might need several minutes to complete.
}
  \item{restr}{
	Vector of numerics. Can be used to filter output for homologue series with certain differences in m/z along the series, given the tolerances in \code{mztol}.
}
}
\note{
Input \code{iso_list} is needed to limit intermediate numbers of m/z differences to screen over, based on changes in mass defect. 
Similarly, intermediate numbers are limited by the retention time and m/z windows defined by the \code{minrt / maxrt} and  \code{minmz / maxmz} arguments. The latter are always set relative of
the individual RT and m/z values of the peaks to be searched from. 

Input \code{peaklist} is internally sorted and saved in the lists returned by (a) increasing retention time and (b) m/z by all \code{\link[nontarget]{pattern.search}}, \code{\link[nontarget]{adduct.search}} 
and \code{\link[nontarget]{homol.search}}. Peak IDs refer to this very order - in contrast to group IDs. Different group IDs exist for adduct groups, isotope pattern groups and homologue series peak 
groups. Moreover, and at the highest level, yet other IDs exist for the individual components (see note section of \code{\link[nontarget]{combine}}).

Here, homologue series group IDs are given both in the function output \code{homol[[1]]} and \code{homol[[3]]}, with one homologue series stating one group of interrelated peaks.
}

\section{Warning}{
The \code{rttol} argument of \code{\link[nontarget]{homol.search}} must not be mixed with that of \code{\link[nontarget]{pattern.search}}.
}
\value{

List of type homol with 5 entries 

\item{homol[[1]]}{\code{Homologue Series}. Dataframe with peaks (\code{mass},\code{intensity},\code{rt},\code{peak ID}) and their homologue series 
relations (\code{to ID},\code{m/z increment},\code{RT increment}) within
groups of homologue series (\code{group ID},\code{series level}).}
\item{homol[[2]]}{\code{Parameters}. Parameters used.}
\item{homol[[3]]}{\code{Peaks in homologue series}. Dataframe listing all peaks (\code{peak IDs}) per homologue series (\code{group ID}), the 
underlying mean m/z & RT increments (\code{m/z increments} and \code{RT increments}) 
and the minimum and maximum RT changes between individual peaks of the series.}
\item{homol[[4]]}{\code{m/z restrictions used}. See function argument \code{restr}.}
\item{homol[[5]]}{\code{Peaks per level}. List of peak IDs per level in the individual series.}

}

\author{
Martin Loos
}
\seealso{
	\code{\link[nontarget]{rm.sat}}
	\code{\link[nontarget]{iso_list}}
	\code{\link[nontarget]{peaklist}}
	\code{\link[nontarget]{plothomol}}
}


\examples{
\donttest{
data(peaklist);
data(iso_list)
homol<-homol.search(
  peaklist,
  iso_list,
  minmz=6,
  maxmz=60,
  minrt=0.5,
  maxrt=1.8,
  ppm=TRUE,
  mztol=3.5,
  rttol=0.4,
  minlen=5,
  plotit=TRUE,
  restr=FALSE
);
plothomol(homol,xlim=FALSE,ylim=FALSE,plotlegend=FALSE);
}
}