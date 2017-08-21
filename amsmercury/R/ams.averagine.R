# $Id: ams.averagine.R 124 2006-10-20 10:02:25Z mkirchner $
#' Averagine model calculation
#' 
#' This calculates the number of H, C, N, O and S atoms in a peptide of specified mass, according to the averagine model of (Senko et al., 1995).
#' @param avgmass the mass for which the element distribution is desired
#' @references Senko MW, Beu SC, McLafferty FW, Determination of Monoisotopic Masses and Ion Populations for Large Biomolecules from Resolved Isotopic Distributions; J. Am. Soc. Mass Spec. 1995, 6, 229-233.
#' @author Marc Kirchner
#' @export
#' @return
#'   \item{model}{
#'  list of the number of atoms for each element, e.g.
#'  \code{list(H=9, C=4, N=1, O=3, S=0)}.
#'  }
#' \item{masserror}{
#'  the mass approximation error \emph{after} hydrogen padding
#' }
#' \item{hydrogencorrection}{
#'  the number of hydrogen atoms used to pad to approximate 
#'  the specified average mass
#'
#' }
#' @examples 
#' ams.averagine(100)
#' ams.averagine(1000)
#' ams.averagine(2000)
ams.averagine <-  function(avgmass) 
{
	# Senko et al., 1995
	model.freq <- list(H=7.7583, C=4.9384, N=1.3577, O=1.4773, S=0.0417);
	# avg. atom masses from GPMAW (http://welcome.to/gpmaw)
	elem.avgmass <- list(H=1.00794, C=12.011, N=14.00670, O=15.99940, S=32.06600);
	model.mass.avg <- 111.1254;

	# number of averagine AAs
	n <- avgmass / model.mass.avg;
	# number of atoms
	freq <- round(sapply(model.freq, '*', n));
	# pad with hydrogen
	modelmass <- sum(freq * as.numeric(elem.avgmass));
	deviance.mass <- avgmass - modelmass; # may be negative for round towards larger ints
	deviance.H <- deviance.mass %/% elem.avgmass$H;
	freq["H"] <- freq["H"] + deviance.H;
	# calculate error
	err.mass <- deviance.mass - deviance.H*elem.avgmass$H;
	
	return(list(model=as.list(freq), masserror=err.mass, hydrogencorrection=deviance.H));	
}
