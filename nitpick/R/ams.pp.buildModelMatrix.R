ams.pp.buildModelMatrix<-function(begin,end,offset,data, EFA=TRUE,maxRangeSize=1,thresholdInt=0,limit=1e-26, charge=c(2,3,4,5), width=width)
{   if (begin<=offset)
      { offsetLower<-begin-1  #set lower offset zero if begin <offset
      }
    else
      { offsetLower<-offset
      }
    if(end+offset>length(data[,1]))
      { offsetUpper<-length(data[,1])-end  #set upper offset zero if begin <offset
      }
    else
      { offsetUpper<-offset
      }


    leng<-end-begin+offsetLower+offsetUpper+1


    lastModelBuilt<-0
    xCombined<-NULL
    iCollection<-matrix(nrow=0,ncol=5)
    preselection<-NULL
  
    for (i in (begin+1):(end-1)) #find all local maxima which are above a certain threshold
      { if ((data[i,2]>data[i-1,2]) &(data[i,2]>data[i+1,2])&(data[i,2]>thresholdInt))
         {  if((i-lastModelBuilt>=maxRangeSize)||lastModelBuilt==0)
            { lastModelBuilt<-i
              if(EFA)
                { modelsPerBin<-length(charge)*2
                  xChargeMatrix<-matrix(nrow=leng,ncol=modelsPerBin)
                  maxXPeak<-vector(length=modelsPerBin)
                  monoisotopic<-vector(length=modelsPerBin)

                 for(m in 1:modelsPerBin)
                   {  if((floor(m/2)-m/2)==0)
                        {  SInd<-"Y"
                        }
                      else
                       {  SInd<-"N"
                       }
                     xChargefftGeneral<-buildAveragineModelShort(data[i,1], data[(begin-offsetLower):(end+offsetUpper),1],charge[floor((m+1)/2)],S=SInd,width=width)
                     xChargefft<-xChargefftGeneral[[1]]  #can be reused!!
                     x2<-Re(fft(xChargefft,inv=T))[1:leng]
                     maxXPeak[m]<-findMaxPeakAndPosition(cbind(c(1:(length(x2))),x2))[2]
                     monoisotopic[m]<-xChargefftGeneral[[2]]
                     xChargeMatrix[,m]<-x2
                  }

                }
              else
               {  modelsPerBin<-length(charge)
                  xChargeMatrix<-matrix(nrow=leng,ncol=modelsPerBin)
                  maxXPeak<-vector(length=modelsPerBin)
                  monoisotopic<-vector(length=modelsPerBin)

                  for(m in 1:modelsPerBin)
                    { SInd<-"Ave"
                       xChargefftGeneral<-buildAveragineModelShort(data[i,1], data[(begin-offsetLower):(end+offsetUpper),1],charge[m],S=SInd,width=width)
                       xChargefft<-xChargefftGeneral[[1]]  #can be reused!!
                       x2<-Re(fft(xChargefft,inv=T))[1:leng]
                       maxXPeak[m]<-findMaxPeakAndPosition(cbind(c(1:(length(x2))),x2))[2]
                       monoisotopic[m]<-xChargefftGeneral[[2]]
                       xChargeMatrix[,m]<-x2
                    }
                }
             }

            if(((i-1)>max(iCollection[,1],0)))
              { xAdd<-NULL
                for(m in 1:modelsPerBin)
                 {  xAdd<-cbind(xAdd,xMoved(i-begin,xChargeMatrix[,m],maxXPeak[m],offsetLower,offsetUpper))
                    iCollection<-rbind(iCollection,c(i-1,charge[(floor((m+1)/2))*EFA+m*(1-EFA)],2-2*(m/2-floor(m/2)),monoisotopic[m],maxXPeak[m]))
                 }
                xCombined<-cbind(xCombined,xAdd)
              }
            xAdd<-NULL
            for(m in 1:modelsPerBin)
              { xAdd<-cbind(xAdd,xMoved(i-begin+1,xChargeMatrix[,m],maxXPeak[m],offsetLower,offsetUpper))
                 iCollection<-rbind(iCollection,c(i,charge[(floor((m+1)/2))*EFA+m*(1-EFA)],2-2*(m/2-floor(m/2)),monoisotopic[m],maxXPeak[m]))
                   }
            xCombined<-cbind(xCombined,xAdd)
            if(!EFA)
              preselection<-c(preselection,(dim(xCombined)[2]-(modelsPerBin-1)):dim(xCombined)[2])
            if(EFA)
              preselection<-c(preselection,((dim(xCombined)[2]-(modelsPerBin-1)):dim(xCombined)[2])[seq(2,modelsPerBin,2)])
              
             xAdd<-NULL
            for(m in 1:modelsPerBin)
              {  xAdd<-cbind(xAdd,xMoved(i-begin+2,xChargeMatrix[,m],maxXPeak[m],offsetLower,offsetUpper))
                 iCollection<-rbind(iCollection,c(i+1,charge[(floor((m+1)/2))*EFA+m*(1-EFA)],2-2*(m/2-floor(m/2)),monoisotopic[m],maxXPeak[m]))
              }
            xCombined<-cbind(xCombined,xAdd)
 
        }
     }
  beginWithOffset<-(begin-offsetLower)
  endWithOffset<-(end+offsetUpper)
  ret <- list(
    xCombined=xCombined,
    preselection=preselection,
    iCollection=iCollection,
    beginWithOffset=beginWithOffset,
    endWithOffset=endWithOffset
  )
  class(ret) <- "ams.pp.strap.model"
  return(ret)
}


xMoved<-function(i,x,maxXPeak,offsetLower,offsetUpper)
{distance<-i-maxXPeak+offsetLower #Rotating the xVector such that its peak corresponds to the peak of the filteredSignal
   lengthX<-length(x)
   xMoved<-c( rep(0,max(distance,0)),
              x[max(-distance+1,0):min(lengthX,lengthX-distance)],
              rep(0,max(-distance,0))
            )
   return(xMoved)
}

