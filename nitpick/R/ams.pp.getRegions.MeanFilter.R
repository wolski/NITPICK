ams.pp.getRegions.MeanFilter<-function (spectrum, mz.nIntervals, mz.neighborhood, noise.threshold.factor)
{
    embedding <- embed(c(rep(0, (mz.neighborhood - 1)/2), spectrum[,
        2], rep(0, (mz.neighborhood - 1)/2)), mz.neighborhood)
    if(dim(embedding)[1]<dim(spectrum)[1])
      {embedding<-rbind(embedding,matrix(0,nrow=(dim(spectrum)[1]-dim(embedding)[1]),ncol=dim(embedding)[2]))
      }
    means <- apply(embedding, 1, mean,na.rm=TRUE)
    region.labels <- cut(spectrum[, 1], mz.nIntervals, labels = F)
    relevantRegions <- NULL
    for (i in 1:mz.nIntervals) {
        tmp.means <- means[region.labels == i]
        noise.mean <- min(tmp.means,na.rm=TRUE)
        relevantRegions <- c(relevantRegions, !(tmp.means > (noise.mean *
            noise.threshold.factor)))
    }
    return(relevantRegions)
}
