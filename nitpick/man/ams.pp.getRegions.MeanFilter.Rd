\name{ams.pp.getRegions.MeanFilter}
\alias{ams.pp.getRegions.MeanFilter}
\title{Determine relevant regions}
\description{
  Reimplementation of Golotvin's / Kaur's method for relevant region detection.
}
\usage{
ams.pp.getRegions.MeanFilter(spectrum, mz.nIntervals, mz.neighborhood, noise.threshold.factor)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{spectrum}{The spectrum in which one searches for relevant regions}
  \item{mz.nIntervals}{The number of intervals the spectrum is divided in}
  \item{mz.neighborhood}{the size of the m/z neighborhood}
  \item{noise.threshold.factor}{the noise threshold factor}
}
\details{
 All details are available in the supplementary material of Kirchner et al. Bioinformatics (2007)
}
\value{
  Relevant regions vector: FALSE == relevant.
}
%\references{ ~put references to the literature/web site here ~ }
\author{ Marc Kirchner \email{marc.kirchner@iwr.uni-heidelberg.de} }
\examples{}
\keyword{multivariate}
