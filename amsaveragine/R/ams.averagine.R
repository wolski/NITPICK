# $Id: ams.averagine.R 124 2006-10-20 10:02:25Z mkirchner $

"ams.averagine" <- 
function(avgmass) 
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
