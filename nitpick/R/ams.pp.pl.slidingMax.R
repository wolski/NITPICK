"ams.pp.pl.slidingMax" <-
function(peaks, charge.state.wise=F, n=3)
{
	discardPeaks <- function(p, n)
	{
		keep <- NULL;
		current.row <- 1;
		while(current.row <= nrow(p)) {
			# p is sorted; so it suffices to look
			# (n-1)/2 bins ahead 
			neighborhood <- which( 
				p[current.row:nrow(p),2] 
				<= p[current.row,2]+(n-1)/2 
			) + current.row -1;
			# find max
			is.max <- which.max(p[neighborhood, 6]);
			keep <- c(keep, neighborhood[is.max]);
			current.row <- max(neighborhood)+1;
		}
		return(p[keep,,drop=F]);
	}

	# make sure peak list is sorted by mass
	if (is.unsorted(peaks[,1])) {
		peaks <- peaks[sort(peaks[,1], index.return=T)$ix,,drop=F];
	}

	if (charge.state.wise) {
		# split data charge state wise
		ret <- NULL;
		for (i in sort(unique(peaks[,5]))) {
			ret <- rbind(ret, discardPeaks(peaks[peaks[,5]==i,,drop=F], n));
		}
		# re-sort by mass
		ret <- ret[sort(ret[,1], index.return=T)$ix,,drop=F];
	} else {
		ret <- discardPeaks(peaks, n)
	}
	return(ret);
}

