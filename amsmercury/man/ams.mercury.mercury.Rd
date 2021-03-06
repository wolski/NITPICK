% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ams.mercury.mercury.R
\name{ams.mercury.mercury}
\alias{ams.mercury.mercury}
\title{calculate accurate isotopic masses and abundances}
\usage{
ams.mercury.mercury(composition, charge, limit)
}
\arguments{
\item{composition}{a vector of integers; needs to have \code{length==MAXLENGTH} of the libmercury++ library.
This is a current shortcoming that is going to be resolvedin the near future.}

\item{charge}{a single integer giving the desired charge}

\item{limit}{a pruning limit, recommended values are between 1e-25 and 1e-30}
}
\value{
Returns a two-column matrix. The first columns holds the m/z values, the second the corrsponding abundances.
}
\description{
\code{ams.mercury.mercury} uses the mercury7 algortihm for
the quick and accurate calculation of isotopic masses and
abundances
}
\details{
See the C++ source
}
\examples{
	# CH4
x <- ams.mercury.mercury(c(4,4,0,0,0), 1, 0);
plot(x[,1], x[,2], type="h")
}
\references{
Rockwood AL, Haimi P, Efficient calculation of accurate masses of isotopic peaks, J Am Soc Mass Spectrom, 17, 415-419, (2006).
}
\author{
Marc Kirchner
}
