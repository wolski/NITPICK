\name{ams.averagine}
\alias{ams.averagine}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Averagine model calculation}
\description{
  This calculates the number of H, C, N, O and S atoms in a peptide of
  specified mass, according to the averagine model of (Senko et al., 1995).
}
\usage{
ams.averagine(avgmass)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{avgmass}{the mass for which the element distribution is desired}
}
\details{
  NA.
}
\value{
  \item{model}{
  	list of the number of atoms for each element, e.g.
	\code{list(H=9, C=4, N=1, O=3, S=0)}.
  }
  \item{masserror}{
  	the mass approximation error \emph{after} hydrogen padding
  }
  \item{hydrogencorrection}{
  	the number of hydrogen atoms used to pad to approximate 
	the specified average mass
  }
}
\references{
	Senko MW, Beu SC, McLafferty FW, Determination of Monoisotopic Masses
	and Ion Populations for Large Biomolecules from Resolved Isotopic
	Distributions; J. Am. Soc. Mass Spec. 1995, 6, 229-233.
}
\author{Marc Kirchner \email{marc.kirchner@iwr.uni-heidelberg.de}}
\keyword{multivariate}
