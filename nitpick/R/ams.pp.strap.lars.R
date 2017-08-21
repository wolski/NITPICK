"ams.pp.strap.lars" <-
function(x, y, type = c("plasso"), trace = FALSE, Gram, 
		eps = .Machine$double.eps,  max.steps, use.Gram = TRUE, EFA=FALSE,
		preselection, bic.steps=5, gdf=FALSE,ignores=NULL,includingZero=FALSE)
{	###### gdf=TRUE uses generalized degrees of freedom, gdf=FALSE uses number of active components for df estimation
	###### bic.steps=0 uses comparison with minimum possible BIC
	stopifnot(require(iwrlars))


	### program automatically centers and standardizes predictors.
	###
	### Original program by Brad Efron September 2001
	### Recoded by Trevor Hastie November 2001
	### Computational efficiency December 22, 2001
	### Bug fixes and singularities February 2003
	### Conversion to R April 2003
	### Copyright Brad Efron and Trevor Hastie
	###
	### Extension for "plasso" by Bernhard Renard, Marc Kirchner, January 2007
	
	call <- match.call()
	type <- match.arg(type)
	TYPE <- switch(type, plasso = "PLASSO")
	if(trace)
		cat(paste(TYPE, "sequence\n"))
	nm <- dim(x)
	n <- nm[1]
	m <- nm[2]
	im <- inactive <- seq(m)
	one <- rep(1, n)
	vn <- dimnames(x)[[2]]	
	### Center x and y, and scale x, and save the means and sds
	meanx <- drop(one %*% x)/n
	x <- scale(x, meanx, FALSE)	# centers x
	normx <- sqrt(drop(one %*% (x^2)))
	nosignal<-normx/sqrt(n) < eps
	ignores<-NULL
	if(any(nosignal)) {
		# ignore variables with too small a variance
		ignores<-c(ignores,im[nosignal])
		inactive<-im[-ignores]
		normx[nosignal]<-eps*sqrt(n)
		if(trace)
			cat("LARS Step 0 :\t", sum(nosignal), "Variables with Variance < \eps; dropped for good\n")	#
	}
	#else ignores <- NULL #singularities; augmented later as well
	names(normx) <- NULL
	x <- scale(x, FALSE, normx)	# scales x
	if(use.Gram & missing(Gram)) {
		if(m > 500 && n < m)
			cat("There are more than 500 variables and n<m;\nYou may wish to restart and set use.Gram=FALSE\n")
		if(trace)
			cat("Computing X'X .....\n")
		Gram <- t(x) %*% x	#Time saving
	}
	mu <- mean(y)
	y <- drop(y - mu)
	Cvec <- drop(t(y) %*% x)
	ssy <- sum(y^2)
	residuals <- y
	if(missing(max.steps))
		max.steps <- 8*min(m, n-1)
	beta <- matrix(0, max.steps + 1, m)	# beta starts at 0
	Gamrat <- NULL
	arc.length <- NULL
	R2 <- 1
	RSS <- ssy
	first.in <- integer(m)
	active <- NULL	# maintains active set
	actions <- as.list(seq(max.steps))	
						# a signed index list to show what comes in and out
	drops <- FALSE	# to do with type=="lasso" or "forward.stagewise"
	Sign <- NULL	# Keeps the sign of the terms in the model
	R <- NULL	###
	k <- 0

	#NN# added Indicator Plasso in while and setting Indicator accordingly
	maxDf <- min(m-length(ignores), n-1)
	if (type=="plasso"){
		maxDf <- maxDf - 1; # skip last step for positive LASSO
	}

f<-try(ssq.beta <- nnls.fit(x[,preselection,drop=F], y),silent=TRUE)

if(class(f)=="try-error") 
{noisefactor=0
 while(class(f)=="try-error"&noisefactor<100) 
    {noise <- 1e8*noisefactor*.Machine$double.eps * matrix(rnorm(nrow(x)*length(preselection), 0, 1), nrow(x), length(preselection))
     f<-try(ssq.beta <- nnls.fit(x[,preselection,drop=F]+noise, y),silent=TRUE)
     noisefactor <- noisefactor+1
     print(noisefactor)
    }
 sigmasq <- 1e99
 if(noisefactor==100)
  {beta<-rep(0,dim(x)[2])
   return(beta)
  }
}
else
{
 	sigmasq <- mean((y-x[,preselection,drop=F] %*% ssq.beta)^2)
	noisefactor <- 0
} 
message("sigmasq: ", sigmasq)
	# regularize, if necessary:
	# add a small amount of noise

while (!is.finite(sigmasq)||(sigmasq >  var(y))) {
	message(noisefactor)
		sigmasq <- 0
		for (frzl in 1:10) {
		message(" -- ", frzl)
		noise <- 1e8*noisefactor*.Machine$double.eps * matrix(rnorm(nrow(x)*length(preselection), 0, 1), nrow(x), length(preselection))
		ssq.beta <- nnls.fit(x[,preselection,drop=F]+noise, y)
		sigmasq <- sigmasq + mean((y-x[,preselection,drop=F] %*% ssq.beta)^2)
		}
		sigmasq <- sigmasq / 10
		noisefactor <- noisefactor +1
	}
	if(sigmasq >  var(y))
	 {
	 	 sigmasq<-.5*var(y)
	 }	 
message("adjusted sigmasq: ", sigmasq, " with noise factor ", noisefactor)
	##### Added necessary precalculations:
sigmasq.Old<-sigmasq      


if(max(ssq.beta)<=.Machine$double.eps)
{return()
}

#	ssq.beta <-as.matrix( nnls.fit(as.matrix(x[,preselection2]), y))
#	sigmasq <- min(mean((y-as.matrix(x[,preselection2]) %*% (ssq.beta))^2),sigmasq)


if(gdf)
{
######(**)	
	##### Added necessary precalculations:
        # length of NNLS solution (full model) and the estimated df for that solutions in order to standardize solutions
	lengthNNLS<-sum(ssq.beta>.Machine$double.eps)  #number of positive coefficients in full models
	covEnd<-(x[,preselection,drop=F]	 %*% ssq.beta)*y
	dfEnd<-sum(covEnd)/sigmasq
}	
	predErrorFull<-sigmasq




	bic <- NULL;

if(includingZero)
{
	bic[1]<-n/predErrorFull*var(y)
}
	###### Added
	nnls.err <-NULL;
	###### Added
	df <-NULL;
	minbic.idx <- 0;
	biccount <- 0;

	### MAIN LOOP ###################################
	while((k < max.steps) & (length(active) < maxDf)) {
		action <- NULL
		k <- k + 1
		C <- Cvec[inactive]
		### identify the largest nonactive gradient
		if (type=="plasso") {
			Cmax <- max(C)
			if (Cmax <= 0)
				break; # nowhere left where we can move without violating the constraint
		} else {
			Cmax <- max(abs(C))
		}
		### Check if we are in a DROP situation
		if(!any(drops)) {
			if (type=="plasso") {
				new <- C >= Cmax - eps # see Efron et. al (2003), section 3.
			} else {
				new <- abs(C) >= Cmax - eps
			}
			C <- C[!new]	# for later
			new <- inactive[new]	# Get index numbers

			### We keep the choleski R  of X[,active] (in the order they enter)
			for(inew in new) {
				if(use.Gram) {
					R <- updateR(
						Gram[inew, inew], 
						R, 
						drop(Gram[inew, active]), 
						Gram=TRUE,
						eps=eps
					)
				} else {
					R <- updateR(
						x[, inew],
						R,
						x[, active],
						Gram=FALSE,
						eps=eps
					)
				}
				if(attr(R, "rank") == length(active)) {
					##singularity; back out
					nR <- seq(length(active))
					R <- R[nR, nR, drop = FALSE]
					attr(R, "rank") <- length(active)
					ignores <- c(ignores, inew)
					action <- c(action,  - inew)
					if(trace)
						cat("LARS Step", k, ":\t Variable", inew, 
							"\tcollinear; dropped for good\n")
				} else {
					if(first.in[inew] == 0)
						first.in[inew] <- k
					active <- c(active, inew)
					#NN# if added
					if (type=="plasso") {
						Sign <- c(Sign, 1)
					} else {
						Sign <- c(Sign, sign(Cvec[inew]))
					}
					action <- c(action, inew)
					if(trace)
						cat("LARS Step", k, ":\t Variable", inew, 
							"\tadded\n")
				}
			}
		} else {
			action <-  - dropid
		}
		Gi1 <- backsolve(R, backsolvet(R, Sign))	
		### Now we have to do the forward.stagewise dance
		### This is equivalent to NNLS
		dropouts<-NULL
		A <- 1/sqrt(sum(Gi1 * Sign))
		w <- A * Gi1	# note that w has the right signs
		if(!use.Gram) u <- drop(x[, active, drop = FALSE] %*% w)	###
		### Now we see how far we go along this direction before the
		### next competitor arrives. There are several cases
		###
		### If the active set is all of x, go all the way
		if(length(active) >=  min(n-1, m - length(ignores) ) ) {
			gamhat <- Cmax/A
		} else {
			if(use.Gram) {
				a <- drop(w %*% Gram[active,  - c(active,ignores), drop = FALSE])
			}
			else {
				a <- drop(u %*% x[,  - c(active, ignores), drop=FALSE])
			}
			#NN# if plasso
			if (type=="plasso"){
				gam <- c((Cmax - C)/(A - a))
			} else {
				gam <- c((Cmax - C)/(A - a), (Cmax + C)/(A + a))
			}	
			### Any dropouts will have gam=0, which are ignored here
			#NN# if plasso
			if (type=="plasso"){
				gamhat <- min(gam[gam > eps])
			} else {
				gamhat <- min(gam[gam > eps], Cmax/A)
			}	
		}
		if(type == "lasso"|type=="plasso") {
			dropid <- NULL
			b1 <- beta[k, active]	# beta starts at 0
			z1 <-  - b1/w
			zmin <- min(z1[z1 > eps], gamhat)
			if(zmin < gamhat) {
				gamhat <- zmin
				drops <- z1 == zmin
			}
			else drops <- FALSE
		}
		beta[k + 1,  ] <- beta[k,  ]
		beta[k + 1, active] <- beta[k + 1, active] + gamhat * w
		if(use.Gram) {
			Cvec <- Cvec - gamhat * Gram[, active, drop = FALSE] %*% w
		}
		else {
			residuals <- residuals - gamhat * u
			Cvec <- drop(t(residuals) %*% x)
		}
		Gamrat <- c(Gamrat, gamhat/(Cmax/A))
		arc.length <- c(arc.length, gamhat)	
		### Check if we have to drop any guys
		if((type == "lasso"|type=="plasso") && any(drops)) {
			dropid <- seq(drops)[drops] # convert logical -> numbers
			for(id in rev(dropid)) {
				if(trace)
				cat("Lasso Step", k+1, ":\t Variable",
					active[id], "\tdropped\n")
				R <- downdateR(R, id)
			}
			dropid <- active[drops]	# indices from 1:m
			beta[k+1,dropid]<-0  # added to make sure dropped coef is zero
			active <- active[!drops]
			Sign <- Sign[!drops]
		}
		if(!is.null(vn))
			names(action) <- vn[abs(action)]
		actions[[k]] <- action
	####Added for Extended Fractional Averagine

if(EFA)
{
  if(is.finite(active[k]))
  { if(((active[k]%%2)==0))
        ignores<-c(ignores,(active[k]-1))
    if(((active[k]%%2)==1))
        ignores<-c(ignores,(active[k]+1))
  }      
}
	#### End Change for Extended Fractional Averagine

		inactive <- im[ - c(active, ignores)]

		##### Change in df estimation (now using Ye) and in stopping criterion
                ##### Added: if and gdf estimation
		if(gdf){
			if(is.finite((max(beta))))
				{
         xbetahat<-x%*%beta[k+1,] #calculate X Betahat for given k
				}
			else
        {
         print("non-finite max(beta)")
         nnls.beta<-nnls.fit(x[,which(beta[k+1,]>0),drop=F], y);
			   xbetahat<-x[,which(beta[k+1,]>0),drop=F]%*%nnls.beta
			  }
				
			cov<-xbetahat*y          #calculated cov(X Betahat(i), y(i))



			df[k]<-sum(cov)/sigmasq   
			df[k]<-df[k]*(lengthNNLS/dfEnd) #standardize such that full model has df of full model
					  # to account for innaccuracy in cov estimation
		}
		else {
			df[k]<-sum(beta[k+1,]!=0)
		}
		##### Addition end
		# get beta from nnls.fit and calcualte rss
		beta.nnls <- nnls.fit(x[,which(beta[k+1,]>0),drop=F], y);
		nnls.estimate <- x[,which(beta[k+1,]>0),drop=F] %*% beta.nnls;
		nnls.residual <- y - nnls.estimate;
		nnls.err[k] <- mean(nnls.residual^2);

		######try if this helps:
		if (nnls.err[k] < sigmasq)
			{nnls.err[k] <- sigmasq
			}

    if(includingZero)
    {
                ###### change bic to include df estimate from further up
		     bic[k+1] <- (n/sigmasq)*(nnls.err[k]+((log(n)*df[k]*sigmasq)/n))
    }
    else
    {
     bic[k] <- (n/sigmasq)*(nnls.err[k]+((log(n)*df[k]*sigmasq)/n))
    }
                ###### Added if and Exact Criterion
		if (bic.steps==0){
			#calculate lowest possible bic at that point (so current df, but min possible pred. error)
			#so pred Error of full model
			possibleBicMin<-(n/sigmasq)*(predErrorFull+((log(n)*df[k]*sigmasq)/n))

			
			if ( ( k>1 ) && ( bic[minbic.idx]<possibleBicMin ) && ((sigmasq.Old-sigmasq)/sigmasq <10  ) ) {
				break
			}
			else{#minBic is searched and stored
				minbic.idx <- max(which(bic == min(bic)));
			}
				

		}
		##### Addition ended
		else{	sigmasq.Old<-sigmasq
			minbic.idx.old <- minbic.idx;
		#browser();
			minbic.idx <- max(which(bic == min(bic)));
			if (minbic.idx == minbic.idx.old) {
				biccount <- biccount + 1;
				if (biccount == bic.steps)
					break
			} else {
				biccount <- 1;
			}
		}
	}
#	beta <- beta[seq(k + 1),  ]	#
if(includingZero)
{
beta <- beta[minbic.idx,,drop=F] 
}
else
{
beta <- beta[minbic.idx+1,,drop=F] 
}

#	dimnames(beta) <- list(paste(0:k), vn)	### Now compute RSS and R2
#	if(trace)
#		cat("Computing residuals, RSS etc .....\n")
#	residuals <- y - x %*% t(beta)
	beta <- scale(beta, FALSE, normx)
#	RSS <- apply(residuals^2, 2, sum)
	#Problem with RSS for plasso by the fact that the last estimates do give some weird results
#	R2 <- 1 - RSS/RSS[1]
#	Cp <- ((n - k - 1) * RSS)/rev(RSS)[1] - n + 2 * seq(k + 1)
# 	object <- list(call = call, type = TYPE, R2 = R2, RSS = RSS, Cp = Cp, 
# 			actions = actions[seq(k)], entry = first.in, Gamrat = Gamrat, 
# 			arc.length = arc.length, Gram = if(use.Gram) Gram else NULL, 
# 			beta = beta, mu = mu, normx = normx, meanx = meanx, bic=bic)
# 	class(object) <- "lars"
        ###### added return information of k
	object <- list(call=call, type="plasso", beta=beta, steps=minbic.idx,k=k,bic=bic,df=df,nnls.err=nnls.err,vary=var(y),sigma=sigmasq,n=n);
	object
}
