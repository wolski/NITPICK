ams.pp.pl.slidingMaxAdjusted <-
function(result, charge.state.wise=F, n=3)
{peaks<-matrix(NA,nrow=length(result[,1]),ncol=7)
 peaks[,1]<-result[,1]
 peaks[,2]<-result[,1]
 peaks[,5]<-result[,2]
 peaks[,6]<-result[,3]
 peaks[,7]<-result[,4]
 print(result)
 print(dim(peaks))
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
 			sul<-sum(p[(p[neighborhood, 7]==1),6])
			nonsul<-sum(p[(p[neighborhood, 7]==2),6])
			midsul<-sum(p[(p[neighborhood, 7]==0),6])
			if (nonsul+sul+midsul==0)
				{sulCont<-0.0417  #from Averagine
				}
			else	{sulCont<-sul/(nonsul+sul+midsul)
				}
			#print(cbind(sulCont,nonsul,sul))
			is.max <- which.max(p[neighborhood, 6]);
			keep <- c(keep, neighborhood[is.max]);
			current.row <- max(neighborhood)+1;
			p[neighborhood[is.max],7]<-sulCont
			
		}
		
		return(p[keep,]);
	}

	# make sure peak list is sorted by mass
	if (is.unsorted(peaks[,1])) {
		peaks <- peaks[sort(peaks[,1], index.return=T)$ix,];
	}

	if (charge.state.wise) {
		# split data charge state wise
		ret <- NULL;
		for (i in sort(unique(peaks[,5]))) {
			ret <- rbind(ret, discardPeaks(peaks[peaks[,5]==i,], n));
		}
		# re-sort by mass
		ret <- ret[sort(ret[,1], index.return=T)$ix,];
	} else {
		ret <- discardPeaks(peaks, n)
	}
 	ret<-ret[,c(1,5,6,7)]
	print(ret)
	return(ret);
}
