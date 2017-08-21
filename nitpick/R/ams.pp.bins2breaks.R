"ams.pp.bins2breaks" <-
function(bins) 
{
	# returns the n+1 breaks enclosing n bins
	n <- length(bins);
	return(
		c(
		bins[1]-0.5*diff(bins[1:2]), 
		bins[1:(n-1)]+0.5*diff(bins), 
		bins[n]+0.5*diff(bins[(n-1):n])
		)
	);
}

