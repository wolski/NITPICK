runPostProcReFit<-function(width,quant,g,pathIn,pathOut,mz.bin=mz.bins, spectrum,chargestates=c(1,2,3,4,5,6),database=FALSE,recalc=FALSE,databasePath=NULL,threshCounts=0){
load(pathIn)
pp_resultList<-matrix(nrow=0,ncol=4)
charges<-chargestates
#source("~/projects/sulveragine/ams.pp.slidingMaxAdjustedSul.R")
##recalc the breakpoints
result<-ams.pp.pl.slidingMaxAdjusted(resultList,FALSE,g)

##if database comparison wanted
if(database)
{
##recalc the breakpoints
if(recalc)
{breakCountMatrix<-matrix(ncol=max(chargestates),nrow=length(mz.bin)+1)
 for (i in charges)
 {breakCountMatrix[,i]<-compareSwissprot(i,mz.bin,databasePath)
 }
 save(breakCountMatrix,file=paste(path,"/breakCountMatrix.RDA",sep=""))
}
else
{load(paste(path,"/breakCountMatrix.RDA",sep=""))
}
resultListSample<-NULL
for (k in chargestates)
{resultListSampleCharge<-result[result[,2]==k,,drop=F]
 resultListSample<- rbind(resultListSampleCharge[breakCountMatrix[resultListSampleCharge[,1],k]>threshCounts, ],resultListSample)
}
result<-resultListSample
}	 
#print(result)
result<-result[result[,1]>0,]
#print(i)
   
  
xCombined<-NULL
for(j in 1:length(result[,1]))
{
 if(result[j,4]>.5)
 {model<-buildAveragineModelShort(mz.bin[result[j,1]],mz.bin, result[j,2],S="Y",width=width)
  xNew<-Re(fft(model[[1]],inverse=TRUE))
  xNew<-c(xNew[max((model[[2]]-result[j,1]),0):length(xNew)],rep(0,max(0,(model[[2]]-result[j,1]))))
 }
 else
  {model<-buildAveragineModelShort(mz.bin[result[j,1]],mz.bin, result[j,2],S="N",width=width)
   xNew<-Re(fft(model[[1]],inverse=TRUE))
	 xNew<-c(xNew[max((model[[2]]-result[j,1]),0):length(xNew)],rep(0,max(0,(model[[2]]-result[j,1]))))
	}
 xCombined<-cbind(xCombined,xNew)
}
print(dim(xCombined))
result[,3]<-nnls.fit(xCombined[1:length(mz.bin),],spectrum[,2])
 
    
 xCombined<-xCombined[,result[,3]>.Machine$double.eps]
 result<-result[result[,3]>.Machine$double.eps,]
 pp_resultList<-result
 save(xCombined,file=paste(pathOut,"xCombined_","_",g,"_",width*1000,".RDA",sep=""))
save(pp_resultList,file=paste(pathOut,"pp_resultList_","_",g,"_",width*1000,".RDA",sep=""))
}



postProc<-function(result){
chargeOrd<-order(result[,2],result[,1])
result<-result[chargeOrd,]
#get rid of similar charge states next to each other, only using max
ppResult<-matrix(0,nrow=length(result)/3,ncol=3)
 nextN<-1
     for (i in 2:(length(result)/3)){
         if (( (result[i,1]-1==result[i-1,1]) && (result[i,2]==result[i-1,2])) |( (i>=3) && (result[i,1]-1==result[i-2,1]) && (result[i,2]==result[i-2,2]))){
           ppResult[nextN,]<-(result[i,3]>=ppResult[nextN,3])*result[i,]+(result[i,3]<ppResult[nextN,3])*ppResult[nextN,]
           }
          else{
            nextN<-nextN+1
            ppResult[nextN,]<-result[i,]
          }
        }
if(ppResult[1,1]<=0){
      ppResult<-ppResult[2:nextN,]
      nextN<-nextN-1
    }
else{
     ppResult<-ppResult[1:nextN,]
    }
#get rid of peaks at same position for two different charge states
positionOrd<-order(ppResult[,1])
ppResult<-ppResult[positionOrd,]
if(ppResult[1,1]<=0){
      ppResult<-ppResult[2:nextN,]
    }

nextN<-1
for (j in 2:(length(ppResult)/3))
{ if(ppResult[j,1]==ppResult[nextN,1]){
    ppResult[nextN,]<-(ppResult[j,3]>=ppResult[nextN,3])*ppResult[j,]+(ppResult[j,3]<ppResult[nextN,3])*ppResult[nextN,]
    #print(cbind(ppResult[j,1],nextN,ppResult[nextN,1]))
  }
  else{
     nextN<-nextN+1
     ppResult[nextN,]<-ppResult[j,]
  }
}
return(ppResult[1:nextN,])
}




compareSwissprot<-function(chargestate,bins,databasePath)
{
breaks<-ams.pp.bins2breaks(bins)
load(databasePath)
monoMass<-vector(length=length(db$H))
for (i in 1:length(db$H))
{
 monoMass[i]<-ams.mercury.mercury(c(db$H[i],db$C[i],db$N[i],db$O[i],db$S[i]),chargestate,1e-26)[1,1]
}

monoBins<-findInterval(monoMass,breaks)

breakCount<-vector(length=length(breaks))
for(j in 1:length(breaks))
{breakCount[j]<-sum(monoBins==j)
}

return(breakCount)
}
