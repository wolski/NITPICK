\name{ams.mercury.mercury}
\alias{ams.mercury.mercury}
\title{ calculate accurate isotopic masses and abundances }
\description{
	\code{ams.mercury.mercury} uses the mercury7 algortihm for
	the quick and accurate calculation of isotopic masses and
	abundances
}
\usage{
ams.mercury.mercury(composition, charge, limit)
}
\arguments{
  \item{composition}{ 
  	a vector of integers; needs to have \code{length==MAXLENGTH} of 
	the libmercury++ library. This is a current shortcoming that is
	going to be resolvedin the near future.
  }
  \item{charge}{ a single integer giving the desired charge }
  \item{limit}{ a pruning limit, recommended values are between 1e-25 and 1e-30 }
}
\details{
  See the C++ source
}
\value{
	Returns a two-column matrix. The first columns holds the m/z values,
	the second the corrsponding abundances.
}
\references{ 
	Rockwood AL, Haimi P, Efficient calculation of accurate masses of
	isotopic peaks, J Am Soc Mass Spectrom, 17, 415-419, (2006).
}
\author{ Marc Kirchner \email{marc.kirchner@iwr.uni-heidelberg.de} }
\examples{
	# CH4
	x <- ams.mercury.mercury(c(4,4,0,0,0), 1, 0);
}
\keyword{ datagen }
