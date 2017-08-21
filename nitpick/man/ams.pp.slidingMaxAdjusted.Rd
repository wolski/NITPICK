\name{ams.pp.slidingMaxAdjusted}
\alias{ams.pp.slidingMaxAdjusted}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Apply a local maximum postprocessing filter to a peak list}
\description{
  Implements a local maximum postprocessing filter for peak lists.
}
\usage{
ams.pp.slidingMaxAdjusted(peaks, charge.state.wise = F, n = 3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{peaks}{peak list (as returned from \code{ams.pp.strap.getPeakList}}
  \item{charge.state.wise}{charge state-wise application of the filter?}
  \item{n}{size of the neighborhood (also called G)}
}
\details{
  FILL ME OUT.
}
\value{
	A peak list, consisting of all local maxima, according to the 
	given parameterization.
}
%\references{ ~put references to the literature/web site here ~ }
\author{Marc Kirchner \email{marc.kirchner@iwr.uni-heidleberg.de}}
\examples{}
\keyword{multivariate}
