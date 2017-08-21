ams.pp.runHierarchical<-function(data,bins,charges,thresholdBetaHat=0,thresholdR2, limit,width,EFA=EFA,gdf=FALSE,bic.steps=0,maxRangeSize=10)
     {print("started running Simple Peak Detection")
      resultSingle<-ams.pp.runSingleRegression(data=data, bins=bins, charges=charges, thresholdBetaHat=thresholdBetaHat,thresholdR2=thresholdR2,width=width)
      nextN<-1
      flaggedAreas<-as.matrix(resultSingle[[3]])
      order<-order(flaggedAreas[,2])
      flaggedAreas<-cbind(rep(0,length(order)),flaggedAreas[order,])
      lengthFlaggedAreasOld<-0
      while(length(flaggedAreas)!=lengthFlaggedAreasOld)
        {lengthFlaggedAreasOld<-length(flaggedAreas)
         flaggedAreas<-joinAdjacent(flaggedAreas)
        }
       flaggedAreas<-flaggedAreas[,2:4] 
 #filter out those flaggedAreas which have been already explained by simple regression
	solvedAreas<-as.matrix(resultSingle[[2]][,4:5])
if(length(solvedAreas)==2)
{solvedAreas<-matrix(solvedAreas,nrow=1,ncol=2)
}
if(length(solvedAreas)>0)
{
 for (i in 1:nrow(solvedAreas))
	{ for (j in 1: nrow(flaggedAreas))
		{
		 if((solvedAreas[i,1]>=flaggedAreas[j,1])&(solvedAreas[i,1]<=flaggedAreas[j,2]))
      {
			 oldEnd<-flaggedAreas[j,2]
			 flaggedAreas[j,2]<-solvedAreas[i,1]-1
			 if(oldEnd>solvedAreas[i,2])
				{rbind(flaggedAreas,c(solvedAreas[i,2]+1,oldEnd,flaggedAreas[j,3]))
				}
			}
		if((solvedAreas[i,2]>=flaggedAreas[j,1])&(solvedAreas[i,2]<=flaggedAreas[j,2]))
			{
			 flaggedAreas[j,1]<-solvedAreas[i,2]+1

			}
		}
	}

}
#remove only one elemental flagged areas
flaggedAreas<-flaggedAreas[flaggedAreas[,1]!=flaggedAreas[,2],]


print("started running Multiple Peak Detection")

      resultMult<-matrix(nrow=0,ncol=4)
      for (k in 1:(length(flaggedAreas[,1])))
             {    kong<-ams.pp.buildModelMatrix(flaggedAreas[k,1],flaggedAreas[k,2],100,data=data,EFA=EFA,maxRangeSize=maxRangeSize,charge=charges,width=width,limit=limit)
                  fooResult<-predict.ams.pp.strap.model(kong,EFA=EFA,data=data,includingZero=TRUE,gdf=gdf,bic.steps=bic.steps)
                  fight<-ams.pp.listPeaks(fooResult,kong$iCollection)
                  resultMult=rbind(resultMult,fight)

             }

      #join result with Result from Single
      final<-rbind(resultMult,resultSingle[[2]][,c(1:3,7)])
            return(final)
     }
