ams.pp.runSingleRegression<-function(data,bins,charges=c(2,3,4,5),thresholdR2,width=0.0001,sulLevels=3,thresholdBetaHat=NULL)
{ flaggedLower<-list()
	flaggedUpper<-list()
	flaggedOff<-list()
 	resultSinglePos<-list()
	resultSingleInt<-list()
	resultSingleCha<-list()
	resultSingleLower<-list()
	resultSingleUpper<-list()
	resultSingleMax<-list()
	resultSingleS<-list()
	for (i in 1:(length(bins[,1])))
          {if (bins[i,1]==0)
             {  	if(is.null(thresholdBetaHat))
                  	{
                  	   thresholdBetaHat<-median(data[bins[i,2]:bins[i,3],2])
                  	}

                distance<-bins[i,3]-bins[i,2]+1
                exactness<-1/((width)/distance)#1/((width/charge)/distance)
                dataBin<-data[bins[i,2]:bins[i,3],2]
                comb<-2*length(charges)
                fftModel<-matrix(ncol=comb,nrow=length(data[bins[i,2]:bins[i,3],1]))
                fftModelPart2<-matrix(ncol=comb,nrow=1)
                filterRes<-matrix(ncol=comb,nrow=length(data[bins[i,2]:bins[i,3],1]))
                resThreshold<-matrix(ncol=comb,nrow=length(data[bins[i,2]:bins[i,3],1]))
                for(j in 1:length(charges))
                   {fftModel[,2*j]<-buildAveragineModelShort((data[bins[i,2],1]+data[bins[i,3],1])/2,data[bins[i,2]:bins[i,3],1],charges[j],width=width,S="Y")[[1]][1:length(data[bins[i,2]:bins[i,3],1])]
             		    fftModel[,2*j-1]<-buildAveragineModelShort((data[bins[i,2],1]+data[bins[i,3],1])/2,data[bins[i,2]:bins[i,3],1],charges[j],width=width,S="N")[[1]][1:length(data[bins[i,2]:bins[i,3],1])]
                   }
                dataBin<-t(cbind(t(dataBin),t(rep(0,max(length(fftModel[,1])-length(dataBin),0)))))

                for(j in 1:comb)
                   {filterRes[,j]<-fft(dataBin)*fftModel[,j]
                    resThreshold[,j]<-Re(fft(filterRes[,j],inverse=T))
                   }

                unsolved<-c(1:(length(resThreshold[,1])))

    	 	         resultsr<-"none"
                while (resultsr!="small"&!is.null(unsolved)&(length(unsolved)>1))
                          { maxVec<-apply(resThreshold[unsolved,],2,max)
			    maxValue<-max(resThreshold[unsolved,])
			    combVec<-c(1:comb)
			    SInd<-2-min(combVec[maxVec==maxValue])%%2
			    S<-"stop"
			if(SInd==1)
			    {S<-"N"
			     }
			if(SInd==2)
			    {S<-"Y"
			     }
	        maxCharge<-ceiling(min(combVec[maxVec==maxValue])/2)+min(charges)-1
	        maxChargePos<-(min(combVec[maxVec==maxValue]))
			    posCharge<-which(combVec==maxChargePos)

			    filteredSignal<-resThreshold[,posCharge]
          maxXPeak<-findMaxPeakAndPositionUns(cbind(c(1:(length(resThreshold[,1]))),data[bins[i,2]:bins[i,3],2]),unsolved)

          xMovedGeneral<-buildAveragineModelShort(data[bins[i,2]+maxXPeak[2]-1,1],data[bins[i,2]:bins[i,3],1],maxCharge,width=width,S=S)
			    xMoved<-Re(fft(xMovedGeneral[[1]],inverse=T))
			    maxPeak<-findMaxPeakAndPosition(cbind(c(1:(length(xMoved))),xMoved))
          xRot<-t(cbind(t(xMoved[min((maxPeak[2]+1),(length(xMoved))):(length(xMoved))]),t(xMoved[0:(maxPeak[2])])))
          offset<-max(2*max(c(1:(length(xRot)/2))[xRot[1:(length(xRot)/2)]>.1*maxPeak[1]])+2,2)

			    if((maxPeak[2]+offset>length(xMoved))||(maxPeak[2]-offset<1))
			     { 	if ((maxPeak[2]+offset>length(xMoved))&&(maxPeak[2]-offset<1))
 					{x<-c(rep(0,offset-maxPeak[2]+1),xMoved,rep(0,offset+maxPeak[2]-length(xMoved)))
					}
 				 else{if (maxPeak[2]-offset<1)
 					{x<-c(rep(0,offset-maxPeak[2]+1),xMoved[1:(maxPeak[2]+offset)])
					}
				     if (maxPeak[2]+offset>length(xMoved))
 					{x<-c(xMoved[(maxPeak[2]-offset):length(xMoved)],rep(0,offset+maxPeak[2]-length(xMoved)))
					}
				    }
			     }
			     else
 			      {
 			        x<-xMoved[(maxPeak[2]-offset):(maxPeak[2]+offset)]
			      }
			    if((bins[i,2]+maxXPeak[2]-offset-1<1)||(bins[i,2]+maxXPeak[2]+offset-1>length(data[,2])))
			     {
				if (bins[i,2]+maxXPeak[2]-offset-1<1)
 				 	{ y<-c(rep(0,-bins[i,2]-maxXPeak[2]+offset+2),data[max(bins[i,2]+maxXPeak[2]-offset-1,1):min(bins[i,2]+maxXPeak[2]+offset-1,length(data[,2])),2])
 					}
 			        if (bins[i,2]+maxXPeak[2]+offset>length(data[,2]))
 					{y<-c(data[max(bins[i,2]+maxXPeak[2]-offset,1):min(bins[i,2]+maxXPeak[2]+offset-1,length(data[,2])),2],rep(0,bins[i,2]+maxXPeak[2]+offset-length(data[,2])))
					}
			     }
			     else
 			      {
   				     y<-data[max(bins[i,2]+maxXPeak[2]-offset-1,1):min(bins[i,2]+maxXPeak[2]+offset-1,length(data[,2])),2]
 			      }
      regress<-lm(y~0+x)
      if(regress$coefficients[1]>thresholdBetaHat)
       {expMax<-maxXPeak[2]
    		lower<-expMax-(offset/2)
      			while(!(lower %in% unsolved) & (lower< expMax))
           				{lower<-lower+1
           				}
         			upper<-expMax+(offset/2)
        			while(!(upper %in% unsolved) & upper> expMax)
           				{upper<-upper-1
           				}
            unsolved<-unsolved[!unsolved %in% (max(expMax-offset/2,1)):min(expMax+offset/2,length(resThreshold[,1]))]

        			if(lower==upper)
          			{lower<-min(unsolved)
           			 upper<-max(unsolved)
           			 unsolved<-NULL
           			}

       			Rsquared<-summary(regress)$r.squared#1-mse/var(y)
       			
       			if(is.na(Rsquared))
     				{  resultsr<-"pos"
         	  }
 	    			
            if(!is.na(Rsquared)&&(Rsquared<thresholdR2))
     				{  resultsr<-"pos"
         			 flaggedLower<-c(flaggedLower,max(bins[i,2],bins[i,2]+maxXPeak[2]-offset/2-1))
				    	 flaggedUpper<-c(flaggedUpper,min(bins[i,3],bins[i,2]+maxXPeak[2]+offset/2-1))
					     flaggedOff<-c(flaggedOff,offset)
 			      }
 	    			else
           {
					  resultsr<-"none"
   					resultSinglePos<-c(resultSinglePos,bins[i,2]+maxXPeak[2]-1-maxPeak[2]+xMovedGeneral[[2]])
	  				resultSingleInt<-c(resultSingleInt,regress$coefficients[1])
		  			resultSingleCha<-c(resultSingleCha,maxCharge)
			  		resultSingleLower<-c(resultSingleLower,max(bins[i,2],bins[i,2]+maxXPeak[2]-offset/2-1))
				  	resultSingleUpper<-c(resultSingleUpper,min(bins[i,3],bins[i,2]+maxXPeak[2]+offset/2-1))
				  	resultSingleMax<-c(resultSingleMax,bins[i,2]+maxXPeak[2]-1)
				  	resultSingleS<-c(resultSingleS,SInd)
				   }
		     }
    		 else
       			{resultsr<-"small"
        		}
		    }
		}
	}

singleResult<-matrix(nrow=length(resultSinglePos),ncol=7)
flagged<-matrix(nrow=length(flaggedLower),ncol=3)

if(nrow(singleResult)>0)
{for (i in 1:length(resultSinglePos))
	{singleResult[i,1]<-resultSinglePos[[i]]
	 singleResult[i,3]<-resultSingleInt[[i]]
	 singleResult[i,2]<-resultSingleCha[[i]]
	 singleResult[i,4]<-resultSingleLower[[i]]
	 singleResult[i,5]<-resultSingleUpper[[i]]
	 singleResult[i,6]<-resultSingleMax[[i]]
	 singleResult[i,7]<-resultSingleS[[i]]
	}
}
if(nrow(flagged)>0)
{for (i in 1:length(flaggedLower))
	{flagged[i,1]<-flaggedLower[[i]]
	 flagged[i,2]<-flaggedUpper[[i]]
	 flagged[i,3]<-flaggedOff[[i]]
	}
}
return(list(0,singleResult,flagged))
}
