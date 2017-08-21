
exactAveragine<-function(mass)
{
  numberAvs<-mass/111.125
  numberHs<-7.7583*numberAvs
  numberCs<-4.9384*numberAvs
  numberNs<-1.3577*numberAvs
  numberOs<-1.4773*numberAvs
  numberSs<-0.0417*numberAvs
  return(c(numberHs,numberCs, numberNs, numberOs,numberSs))
}
lightAveragine<-function(mass)
{
  numberAvs<-mass/109.7878
  numberHs<-7.7583*numberAvs
  numberCs<-4.9384*numberAvs
  numberNs<-1.3577*numberAvs
  numberOs<-1.4773*numberAvs
  return(c(numberHs,numberCs, numberNs, numberOs,0))
}
sulveragine<-function(mass)
{
  numberAvs<-mass/417.5719
  numberHs<-27*numberAvs
  numberCs<-13*numberAvs
  numberNs<-3*numberAvs
  numberOs<-6*numberAvs
  numberSs<-3*numberAvs
  nonSul<-lightAveragine(mass)
  result<-(6*nonSul+c(numberHs,numberCs, numberNs, numberOs,numberSs))/7
  return(result)
}





# calculates pdf for remaining fractional part of atom abundance (so abundance-floor(abundance))
# molecule is either S, C ,H, N, or O
# exactness is given as the inverse of the required exactness
# correction indicates whether a pdf should be produced (so corrected for area  under curve=1


# builds fractional binomial distribution by building a hierarchical model of 2 binomial distribution with amount being
# the probability to decide for either one
# only for amount in (0,1)
getDistributionHBin<-function(molecule,amount, charge)
{  
  if((amount<0)|(amount>1)){
    print("fractional binomial is only defined for amount in (0,1)")
    return()
  }
  if (molecule=="H"){
    S1<-matrix(c(1.0078246,	2.0141021,0.99985,	0.00015),nrow=2,ncol=2)
    } # H
  if (molecule=="C")
  {S1<-matrix(c(12.0000000, 	13.0033554,0.988930,	0.011070 ),nrow=2,ncol=2)} # C
  if (molecule=="N")
  {S1<-matrix(c(14.0030732, 	15.0001088,0.996337,	0.003663 ),nrow=2,ncol=2)} # N
  if (molecule=="O")
  {S1<-matrix(c(15.9949141, 	16.9991322, 	17.9991616,0.997590,	0.000374,	0.002036 ),nrow=3,ncol=2)} # O
  if (molecule=="S")
  {S1<-matrix(c(31.972070, 	32.971456, 	33.967866, 	34,	35.967080,0.9502,	0.0075,		0.0421, 0,	0.0002 ),nrow=5,ncol=2)} # S
  S1[,1]<-S1[,1]/charge         
  S1new<-S1
  S1new[1,1]<-S1[1,1]*amount           #correcting the mass
  S1new[,1]<-S1[,1]-S1[1,1]+S1new[1,1]
  S0<-S1new                            #distribution for no atom
  S0[,2]<-0
  S0[1,2]<-1                           #all mass at zero
  Snew<-S1new                          #return object
  Snew[,2]<-amount*S1new[,2]+(1-amount)*S0[,2]
  if (molecule=="S")                   #exclude p()=0
  {Snew<-Snew[c(1,2,3,5),]
  }
  return(Snew)
}

#adapts scale of a data set to match requirements of exactntess
#Funktion funktioniert noch nicht fuer den Fall hdistr, also wenn
#nicht 0er eingefuegt werden muessen, sondern vielmehr Daten ausgesucht


fftDiscrete<-function(distr1, distr2,limit=1e-26)
{kmax<-length(distr1[,1])+length(distr2[,1])
distrNew<-matrix(0,nrow=kmax,ncol=2)
for (k in 2:kmax)
{for (i in 1:min((k-1),length(distr1[,1])))
{if (  ((k-i)>=0)  &&  ((k-i)<=length(distr2[,1])))
{distrNew[k,2]=distrNew[k,2]+distr1[i,2]*distr2[(k-i),2]
distrNew[k,1]=distrNew[k,1]+distr1[i,2]*distr2[(k-i),2]*(distr1[i,1]+distr2[(k-i),1])
}
}
}
distrNew<-distrNew[distrNew[,2]>0,]   
distrNew[,1]<-distrNew[,1]/distrNew[,2]
distrNew<-distrNew[distrNew[,2]>limit,]
return(distrNew)
}

#function that builds an exact averagine model by convoluting
#mercury result for next lowest integer model and the continuous binomials
#for each element

buildAveragineModelShort<-function(mass, masses, charge, limit=1e-26,width=.00005,S="Ave") #width needs to be set visually
{stopifnot(require(amsmercury))
  
  if (S=="Ave")
  {
    stoich<-exactAveragine(mass*1.000641*charge)
  }
  if (S=="Y")
  {
    stoich<-sulveragine(mass*1.000641*charge)
  }
  if (S=="N")
  {
    stoich<-lightAveragine(mass*1.000641*charge)
  }
  basic<-ams.mercury.mercury(stoich, charge, limit)
  hdistr<-getDistributionHBin("H",(stoich[1]-floor(stoich[1])),charge)
  cdistr<-getDistributionHBin("C",(stoich[2]-floor(stoich[2])),charge)
  ndistr<-getDistributionHBin("N",(stoich[3]-floor(stoich[3])),charge)
  odistr<-getDistributionHBin("O",(stoich[4]-floor(stoich[4])),charge)
  sdistr<-getDistributionHBin("S",(stoich[5]-floor(stoich[5])),charge)
  
  
  
  result<-fftDiscrete(hdistr,cdistr,limit)
  result<-fftDiscrete(result,ndistr,limit)
  result<-fftDiscrete(result,odistr,limit)
  result<-fftDiscrete(result,sdistr,limit)
  result<-fftDiscrete(basic,result,limit)
  
  masses<-ams.pp.bins2breaks(masses)
  sig <- rep(0, length(masses))
  position<-vector(length=length(result[,1]))
  position<-findInterval(result[,1],masses)
  for (i in 1:length(result[,1]))
  {sig[position[i]]<-sig[position[i]]+result[i,2]
  }
  minPosition<-min(position)
  
  lpf<-vector(length=length(masses))
  lastResult<-1
  i<-floor(length(masses)/2)
  lastResult<-1
  lpf[1]<-lastResult
  change<-1
  lpfi<-dnorm(masses[i],masses[i],width*masses[i])
  while((lastResult>limit)&(change<i))
  {lpf[length(lpf)+1-change]<-dnorm(masses[i-change],masses[i],width*masses[i])/lpfi
  lpf[1+change]<-dnorm(masses[i+change],masses[i],width*masses[i])/lpfi
  lastResult<-min(lpf[1+change],lpf[length(lpf)+1-change])
  change<-change+1
  }
  res<-fft(sig)*fft(lpf)/length(lpf)
  return(list(res,minPosition))
}




findMaxPeakAndPosition<-function(data)
{       maxValue<-max(data[,2])
position<-max(c(1:length(data[,2]))[data[,2]==maxValue])
return(cbind(maxValue,position))
}

findMaxPeakAndPositionUns<-function(data, unsolved)
{
  maxValue<-max(data[unsolved,2])
  position<-max(unsolved[data[unsolved,2]==maxValue])
  return(cbind(maxValue,position))
}   

joinAdjacent<-function(bins)
{binsNew<-bins   #copy current binning structure
nextN<-1
for (i in 1:(length(bins[,1])-1))  #for all bins
{if ((binsNew[nextN,1]==0&bins[i+1,1]==0)&(binsNew[nextN,3]>=(bins[i+1,2])-1))   #if 2 adjacent are both 0
{    
  binsNew[nextN,3]<-bins[i+1,3]
  if(binsNew[nextN,2]>bins[i+1,2])
  {binsNew[nextN,2]<-bins[i+1,2]
  }
}
  else{ #otherwise no change taking place
    nextN<-nextN+1
    binsNew[nextN,]<-bins[i+1,]
  }
} #return new Bins
return(binsNew[1:nextN,])
}
