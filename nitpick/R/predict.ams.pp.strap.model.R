predict.ams.pp.strap.model<-function(model,EFA=TRUE,data,bic.steps=0,gdf=TRUE,includingZero=FALSE)
{   xCombined<-model$xCombined
    if(!is.null(xCombined))
    {
     beginWithOffset<-model$beginWithOffset
     endWithOffset<-model$endWithOffset
     preselection<-model$preselection
    
     larsResult<-NULL
     larsResult$beta<-NULL
     if(length(preselection)==0)
       { larsResult$beta<-0
       }
     else
       {
         larsResult<-ams.pp.strap.lars(xCombined,data[beginWithOffset:endWithOffset,2], EFA=EFA,preselection=preselection,bic.steps=bic.steps,gdf=gdf,ignores=ignores,includingZero=includingZero)
       }
     if(is.null(larsResult)||is.null(larsResult$beta)|| !is.finite(larsResult$beta) ||(sum(larsResult$beta)==0))
       { xlassoTotal<-NULL
          nnls<-NULL
          beta<-NULL
       }
     else
       { xlassoTotal<-c(1:dim(xCombined)[2])[larsResult$beta>0]
         nnls<-nnls.fit(as.matrix(xCombined[,xlassoTotal]),data[beginWithOffset:endWithOffset,2])
         beta<-rep(0,dim(xCombined)[2])
         beta[xlassoTotal]<-nnls
       }
     }
     else
     {
      beta<-vector(length=0)
     }   
     return(beta)
}