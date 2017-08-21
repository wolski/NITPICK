# Publication script 


library(nitpick)


##real data

##load spectrum 2xN matrix of 
spectrum<-read.table("TOF-MS-yylBSAstd-sample4-23.817-29.278-rebinned.txt")

#adjust spectrum for relevant region detection to not have zero entries (and thus increase robustness)
minSpec<-min(spectrum[spectrum[,2]>0,2])
spectrum[,2]<-spectrum[,2]+minSpec

#relevant region detection
rr<-ams.pp.getRegions.MeanFilter(spectrum,30,15,5)


resultList<-matrix(nrow=0,ncol=4)

##bin data according to relevant regions
bins<-matrix(nrow=1,ncol=3)
bins[1,1]<-rr[1]
bins[1,2]<-1
nextn<-1
addBins<-matrix(nrow=1,ncol=3)

for (j in 2:length(rr))
{  if(rr[j]!=rr[j-1])
  {bins[nextn,3]<-j-1
   nextn<-nextn+1
   addBins[1,1]<-rr[j]
   addBins[1,2]<-j
   addBins[1,3]<-j
   bins<-rbind(bins,addBins)
 }
}

bins[nextn,3]<-length(rr)

bins[bins[,3]==bins[,2],3]<-bins[bins[,3]==bins[,2],3]+2
bins[(bins[,3]-bins[,2])<2,3]<-bins[(bins[,3]-bins[,2])<2,3]+1
#bins[bins[,3]==bins[,2],1]<-1
#bins[(bins[,3]-bins[,2])<2,1]<-1


  lengthBinsOld<-0
      while(length(bins)!=lengthBinsOld)
        {lengthBinsOld<-length(bins)
         bins<-joinAdjacent(bins)
        }



## actual peak detection for each bin featuring a relevant signal

for(n in 1:length(bins[,1]))
{
if(bins[n,1]==0)
{
kong<-ams.pp.buildModelMatrix(bins[n,2],bins[n,3],50,spectrum,EFA=FALSE,charge=c(1,2,3,4,5,6),width=0.00008,maxRangeSize=10)
fooResult<-predict(kong,EFA=FALSE,spectrum,includingZero=FALSE,gdf=TRUE,bic.steps=0)
fight<-ams.pp.listPeaks(fooResult,kong$iCollection)
resultList<-rbind(resultList,fight)
}
}
#save result and run postprocessing
save(resultList,file="0resultList_width_0.00008.rda")
runPostProc(0.00008,0,30,"0resultList_width_0.00008.rda",pathOut="",mz.bin=spectrum[,1], chargestates=c(1,2,3,4,5,6))
load("pp_resultList__30_0.08.RDA")
pp_resultList_BS<-cbind(spectrum[pp_resultList[,1],1],pp_resultList[,2:3])
save(pp_resultList_BS,file="0pp_resultList_g_30width0.08.RDA")




##simulation


load("mz.bins.500-700.equispaced11.rda")


runSim<-function(SNR,Nu,thresh,width,path,charge, limit,EFA=TRUE,hierarchical=FALSE,gdf=TRUE,bic.steps=0,maxRangeSize=10){
	#make result vectors
	resultList<-matrix(nrow=0,ncol=4)
	lengthResultList<-vector(length=Nu)
	inCorrectList<-list()
	correct<-vector(length=Nu)
	tooMany<-vector(length=Nu)
	errsRss.tp<-vector(length=Nu)
	sumFrzl<-vector(length=Nu)
	sumBeta<-vector(length=Nu)               
	numberSimple<-vector(length=Nu)


       #load data
	load("relevantRegions-uniprot_sprot-HUMAN-tryptic-equispaced11-500-700-charge-1-N-500-k-20.rda")

	load("uniprot_sprot-HUMAN-tryptic-equispaced11-500-700-charge-1-N-500-k-20.rda")

	load(paste("uniprot_sprot-HUMAN-tryptic-equispaced11-500-700-charge-1-N-500-k-20-spectra-SNR-",SNR,".rda",sep=""))
	
	load("mz.bins.500-700.equispaced11.rda") 

	#define bins
	for(i in 1:Nu)
	{
		resultList2<-matrix(nrow=0,ncol=4)

		bins<-matrix(nrow=1,ncol=3)
		bins[1,1]<-rr[i,1]
		bins[1,2]<-1
		nextn<-1
		addBins<-matrix(nrow=1,ncol=3)

		for (j in 2:length(rr[i,]))
		{  if(rr[i,j]!=rr[i,j-1])
	         {bins[nextn,3]<-j-1
	          nextn<-nextn+1
	          addBins[1,1]<-rr[i,j]
	          addBins[1,2]<-j
	          addBins[1,3]<-j
		  bins<-rbind(bins,addBins)
	         }
		}

		bins[nextn,3]<-length(rr[i,])
		dat<-cbind(mz.bins,N$S[i,])

		#for each bin, run strap
		for(n in 1:length(bins[,1]))
			{
				if(bins[n,1]==0)
				{
					kong<-ams.pp.buildModelMatrix(bins[n,2],bins[n,3],100,dat,EFA=EFA,charge=charge,width=width,maxRangeSize=maxRangeSize)
					fooResult<-predict(kong,EFA=EFA,dat,includingZero=hierarchical,gdf=gdf,bic.steps=bic.steps)
					fight<-ams.pp.listPeaks(fooResult,kong$iCollection)
					resultList2<-rbind(resultList2,fight)
				}
			}	

		result<-resultList2
		

		print(result)

		#obtaining ground Truth
		breaks<-ams.pp.bins2breaks(mz.bins)
		m0<-findInterval(samples$info[5][((i-1)*20+1):(i*20),1],breaks)
		zInfo<-samples$info[6][((i-1)*20+1):(i*20),1]
		intensity<-samples$info[7][((i-1)*20+1):(i*20),1]
		groundTruth<-matrix(0,nrow=20,ncol=3)
		#Calculating tpr and fpr
		for(j in 1:20)
			{groundTruth[j,]<-cbind(m0[j],zInfo[j],intensity[j])
			}
		print(groundTruth)
		est.TPs.idx<-matrix(0,nrow=0,ncol=3)
		idfyerSet<-NULL
		for (k in 1:20)
		{
		  if (sum(cbind(groundTruth[k,1],groundTruth[k,2]) %in% result[,1:2])==2)
		    {correct[i]<-correct[i]+1
		     idfyer<-which(   (result[,1]==groundTruth[k,1]) & (result[,2]==groundTruth[k,2]) )
		     est.TPs.idx<-rbind(est.TPs.idx,cbind(result[idfyer,1],result[idfyer,2],result[idfyer,3]))
		     idfyerSet<-cbind(idfyerSet,idfyer)
		    }
		    else
		    {inCorrectList<-c(inCorrectList, cbind(groundTruth[k,], i))
		    }
		}
		if(!is.null(idfyerSet))
		{
			falseResult<-result[-idfyerSet,,drop=F]
		}
		else
		{
			falseResult<-result
		}
		print(falseResult)
		tooMany[i]<-length(falseResult[,1])

		#store results
		#resultList<-list(resultList,result)
		resultList<-rbind(resultList,result)
		lengthResultList[i]<-length(result)/4
		#beta gewurschtel quantitativ
		beta <- NULL;
		for (n in 1:length(est.TPs.idx[,1])) {
			# find corresponding peak in sample set
			 c.idx <- which(groundTruth[,1] == est.TPs.idx[n,1]);
			 true.beta <- sum(groundTruth[c.idx,3]); # guard against multiple peaks
			 beta <- rbind(beta, c(est.TPs.idx[n,3], true.beta))
		}
		betasum <- apply(beta,2,sum);
		errsRss.tp [i]<- sum((beta[,1]/betasum[1] - beta[,2]/betasum[2])^2);
		# linus
		frzl <- rbind(beta, cbind(result[!(result[,1] %in% est.TPs.idx[,1]),3] ,rep(0,sum(!(result[,1] %in% est.TPs.idx[,1])))))#tooMany[i])));
		frzl <- t(t(frzl)/apply(frzl,2,sum));
		#browser();
		errsRss.tp[i]<-(sum((frzl[,1]-frzl[,2])^2) / (1/correct[i]*errsRss.tp[i]*(length(result[,1]))))-1#ncol(pq.truth))) -1
		sumFrzl[i]<-sum((frzl[,1]-frzl[,2])^2)
		sumBeta[i]<-sum((beta[,1]/betasum[1] - beta[,2]/betasum[2])^2);
	
	
	
	}
	save(correct,file=paste(path,"correct_",thresh,"_",SNR,"_",width*1000,".RDA",sep=""))
	save(tooMany,file=paste(path,"tooMany_",thresh,"_",SNR,"_",width*1000,".RDA",sep=""))
	save(inCorrectList,file=paste(path,"inCorrectList_",thresh,"_",SNR,"_",width*1000,".RDA",sep=""))
	save(sumFrzl,file=paste(path,"sumFrzl_",thresh,"_",SNR,"_",width*1000,".RDA",sep=""))
	save(sumBeta,file=paste(path,"sumBeta_",thresh,"_",SNR,"_",width*1000,".RDA",sep=""))
	save(resultList,file=paste(path,"resultList_",thresh,"_",SNR,"_",width*1000,".RDA",sep=""))
	save(lengthResultList,file=paste(path,"lengthResultList_",thresh,"_",SNR,"_",width*1000,".RDA",sep=""))
	save(errsRss.tp,file=paste(path,"errs_",thresh,"_",SNR,"_",width*1000,".RDA",sep=""))
	save(numberSimple,file=paste(path,"numberSimple_",thresh,"_",SNR,"_",width*1000,".RDA",sep=""))
	
	
	#calculate big N for PPV calculation
	Num<-vector(length=Nu)
	for (i in 1:Nu)
		{Num[i]<-5500-sum(rr[i,]==1)
		}
	Num<-2*Num
	TP<-correct
	FP<-tooMany
	FN<-50-TP
	TN<-Num - (FP+TP+FN)

	ACC<-(TP+TN)/Num
	NPV<-TN/(TN+FN)
	PPV<-TP/(TP+FP)
	quantError<-sumFrzl/(1/correct*sumBeta*(correct+tooMany))


}


## adjusted postprocessing for ground truth comparison 
runPostProc<-function(Nu,SNR,thresh,width,quant,g,path,threshCounts=0,mz.bin=mz.bins, chargestates=c(2,3,4,5),recalc=FALSE){
load(paste(path,"resultList_",thresh,"_",SNR,"_",width*1000,".RDA",sep=""))
load(paste(path,"lengthResultList_",thresh,"_",SNR,"_",width*1000,".RDA",sep=""))
inCorrectList<-list()
correct<-vector(length=Nu)
tooMany<-vector(length=Nu)
errsRss.tp<-vector(length=Nu)
sumFrzl<-vector(length=Nu)
sumBeta<-vector(length=Nu)
pp_resultList<-matrix(nrow=0,ncol=4)
pp_lengthResultList<-vector(length=Nu)
charges<-chargestates

##recalc the breakpoints
if(recalc)
{breakCountMatrix<-matrix(ncol=max(chargestates),nrow=length(mz.bin)+1)
	 for (i in charges)
		{breakCountMatrix[,i]<-compareSwissprot(i,mz.bin)
		}
	save(breakCountMatrix,file=paste(path,"/breakCountMatrix.RDA",sep=""))
}
else
{	load(paste(path,"/breakCountMatrix.RDA",sep=""))
}

if(chargestates[1]==2)
{
load("uniprot_sprot-HUMAN-tryptic-equispaced11-500-700-charges-2-5-N-500-k-20.rda")
load(paste("uniprot_sprot-HUMAN-tryptic-equispaced11-500-700-charges-2-5-N-500-k-20-spectra-SNR-",SNR,".rda",sep=""))
}
else
{
load("uniprot_sprot-HUMAN-tryptic-equispaced11-500-700-charge-1-N-500-k-20.rda")
load(paste("uniprot_sprot-HUMAN-tryptic-equispaced11-500-700-charge-1-N-500-k-20-spectra-SNR-",SNR,".rda",sep=""))
}

for(i in 1:Nu)
{
      resultListSample<-resultList[max((sum(lengthResultList[min(i-1,1):(i-1)])+1),1):sum(lengthResultList[1:i]),]
      resultListSample<-resultListSample[resultListSample[,1]>0,]
      for (k in chargestates)
	{ resultListSampleCharge<-resultListSample[resultListSample[,2]==k,,drop=F]
	  otherResults<-resultListSample[resultListSample[,2]!=k,]
         resultListSample<- rbind(resultListSampleCharge[breakCountMatrix[resultListSampleCharge[,1],k]>threshCounts, ],otherResults)
       }

       result<-resultListSample
       result<-ams.pp.pl.slidingMaxAdjusted(result,FALSE,g)
       result<-result[result[,1]>0,]
       xCombined<-NULL
	for(j in 1:length(result[,1]))
          {
    		if(result[j,4]>.5)
			{model<-buildAveragineModelShort(mz.bin[result[j,1]],mz.bin, result[j,2],S="Y")
		         xNew<-Re(fft(model[[1]],inverse=TRUE))
	     		  xNew<-c(xNew[max((model[[2]]-result[j,1]),0):length(xNew)],rep(0,max(0,(model[[2]]-result[j,1]))))
			}
		else	{model<-buildAveragineModelShort(mz.bin[result[j,1]],mz.bin, result[j,2],S="N")
		        xNew<-Re(fft(model[[1]],inverse=TRUE))
	   		 xNew<-c(xNew[max((model[[2]]-result[j,1]),0):length(xNew)],rep(0,max(0,(model[[2]]-result[j,1]))))
			}

             xCombined<-cbind(xCombined,xNew)
	  }
	print(dim(xCombined))
	print(length(N$S[i,]))
      
       
       xCombined<-xCombined[,result[,3]>.Machine$double.eps]



      #obtaining ground Truth
       breaks<-ams.pp.bins2breaks(mz.bin)
	m0<-findInterval(samples$info[5][((i-1)*20+1):(i*20),1],breaks)
	zInfo<-samples$info[6][((i-1)*20+1):(i*20),1]
	intensity<-samples$info[7][((i-1)*20+1):(i*20),1]
	groundTruth<-matrix(0,nrow=20,ncol=3)
	#Calculating tpr and fpr
	for(j in 1:20)
	{groundTruth[j,]<-cbind(m0[j],zInfo[j],intensity[j])
	}
	print(groundTruth)
	est.TPs.idx<-matrix(0,nrow=0,ncol=3)
	idfyerSet<-NULL
	for (k in 1:20)

      {
          if (sum(cbind(groundTruth[k,1],groundTruth[k,2]) %in% result[,1:2])==2)
            {correct[i]<-correct[i]+1
            idfyer<-which(   (result[,1]==groundTruth[k,1]) & (result[,2]==groundTruth[k,2]) )
            est.TPs.idx<-rbind(est.TPs.idx,cbind(result[idfyer,1],result[idfyer,2],result[idfyer,3]))
            idfyerSet<-cbind(idfyerSet,idfyer)
            }
          else
           {inCorrectList<-list(inCorrectList, groundTruth[k,], i)
           }
     }

     falseResult<-result[-idfyerSet,]
     
     print(falseResult)
     if(length(falseResult)==4)
	{tooMany[i]<-1
	}
     else{
         tooMany[i]<-length(falseResult[,1])
	 }
     #beta gewurschtel quantitativ
     beta <- NULL;
     for (n in 1:length(est.TPs.idx[,1])) {
     # find corresponding peak in sample set
        c.idx <- which(groundTruth[,1] == est.TPs.idx[n,1]);
        true.beta <- sum(groundTruth[c.idx,3]); # guard against multiple peaks
        beta <- rbind(beta, c(est.TPs.idx[n,3], true.beta))
      }
      betasum <- apply(beta,2,sum);
      errsRss.tp [i]<- sum((beta[,1]/betasum[1] - beta[,2]/betasum[2])^2);
      # linus
      frzl <- rbind(beta, cbind(result[!(result[,1] %in% est.TPs.idx[,1]),3] ,rep(0,sum(!(result[,1] %in% est.TPs.idx[,1])))))#tooMany[i])));
      #print(frzl)
      frzl <- t(t(frzl)/apply(frzl,2,sum));
      errsRss.tp[i]<-(sum((frzl[,1]-frzl[,2])^2) / (1/correct[i]*errsRss.tp[i]*(length(result[,1]))))-1#ncol(pq.truth))) -1
      sumFrzl[i]<-sum((frzl[,1]-frzl[,2])^2)
      sumBeta[i]<-sum((beta[,1]/betasum[1] - beta[,2]/betasum[2])^2);
      
      pp_resultList<-rbind(pp_resultList,result)
      pp_lengthResultList[i]<-length(result)/4	
    }

save(correct,file=paste(path,"pp_correct_",thresh,"_",g,"_",quant,"_",SNR,"_",width*1000,".RDA",sep=""))
save(tooMany,file=paste(path,"pp_tooMany_",thresh,"_",g,"_",quant,"_",SNR,"_",width*1000,".RDA",sep=""))
save(inCorrectList,file=paste(path,"pp_inCorrectList_",thresh,"_",g,"_",quant,"_",SNR,"_",width*1000,".RDA",sep=""))
save(pp_resultList,file=paste(path,"pp_resultList_",thresh,"_",g,"_",quant,"_",SNR,"_",width*1000,".RDA",sep=""))
save(pp_lengthResultList,file=paste(path,"pp_lengthResultList_",thresh,"_",g,"_",quant,"_",SNR,"_",width*1000,".RDA",sep=""))
save(errsRss.tp,file=paste(path,"pp_errs_",thresh,"_",g,"_",quant,"_",SNR,"_",width*1000,".RDA",sep=""))

}




# now run simulation

runSim(100,500,0,.0001,"/results/",c(1),1e-26,EFA=FALSE,hierarchical=FALSE,gdf=TRUE,bic.steps=0,maxRangeSize=10)
runSim(50,500,0,.0001,"/results/",c(1),1e-26,EFA=FALSE,hierarchical=FALSE,gdf=TRUE,bic.steps=0,maxRangeSize=10)
runSim(25,500,0,.0001,"/results/",c(1),1e-26,EFA=FALSE,hierarchical=FALSE,gdf=TRUE,bic.steps=0,maxRangeSize=10)
runSim(10,500,0,.0001,"/results/",c(1),1e-26,EFA=FALSE,hierarchical=FALSE,gdf=TRUE,bic.steps=0,maxRangeSize=10)
runSim(5,500,0,.0001,"/results/",c(1),1e-26,EFA=FALSE,hierarchical=FALSE,gdf=TRUE,bic.steps=0,maxRangeSize=10)


runPostProc(500,100,0,.0001,0,3,"/results/",0,mz.bin=mz.bins, chargestates=c(1),recalc=TRUE)
runPostProc(500,50,0,.0001,0,3,"/results/",0,mz.bin=mz.bins, chargestates=c(1),recalc=FALSE)
runPostProc(500,25,0,.0001,0,3,"/results/",0,mz.bin=mz.bins, chargestates=c(1),recalc=FALSE)
runPostProc(500,10,0,.0001,0,3,"/results/",0,mz.bin=mz.bins, chargestates=c(1),recalc=FALSE)
runPostProc(500,5,0,.0001,0,3,"/results/",0,mz.bin=mz.bins, chargestates=c(1),recalc=FALSE)
