"nnls.fit" <-
function(x,y,wsqrt=1,eps=0,rank.tol=1e-07) {
  ## Purpose: Nonnegative Least Squares (similar to the S-Plus function
  ## with the same name) with the help of the R-library quadprog

 ## ------------------------------------------------------------------------
  ## Attention:
  ## - weights are square roots of usual weights
  ## - the constraint is coefficient>=eps

 ## ------------------------------------------------------------------------
  ## Author: Marcel Wolbers, July 99
  ##
#========================================================================
  require ("quadprog")
  m <- NCOL(x)
  if (length(eps)==1) eps <- rep(eps,m)
  x <- x * wsqrt
  y <- y * wsqrt
  #sometimes a rescaling of x and y helps (if solve.QP.compact fails otherwise)
  xscale <- apply(abs(x),2,mean)
  yscale <- mean(abs(y))
  x <- t(t(x)/xscale)
  y <- y/yscale
  Rinv <- backsolve(qr.R(qr(x)),diag(m))
  cf <- solve.QP.compact(Dmat=Rinv,dvec=t(x)%*%y,Amat=rbind(rep(1,m)),
                   Aind=rbind(rep(1,m),1:m),bvec=eps*xscale/yscale,
                         factorized=TRUE)$sol
  cf <- cf*yscale/xscale  #scale back
  cf
}

