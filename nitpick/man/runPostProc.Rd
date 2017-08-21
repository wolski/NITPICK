\name{runPostProc}
\alias{runPostProc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{postprocessing based on local maxima detection in neighborhood}
\description{
postprocessing by applying local maxima detection}
\usage{
runPostProc(width,quant=0,g=3,pathIn,pathOut,mz.bin=mz.bins, chargestates=c(1,2,3,4,5,6),database=FALSE,recalc=FALSE,databasePath=NULL,threshCounts=0)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{width}{peak width}}
  \item{quant}{lower quantile of the intensity of peaks which should be disregarded}
  \item{g}{size of the neighborhood (also called G)}
  \item{pathIn}{Path where results from peak detection are to be found}
  \item{pathOut}{swhere results should be stored}
  \item{mz.bin}{binning scheme of the spectrum}
  \item{chargestates}{charges occuring in the spectrum}
  \item{database}{if mz-values should be compared to possible values from database}
  \item{recalc}{if computation of possible values shoud be retrieved}
  \item{databasePath}{where database is to be found}
  \item{threshCount}{how often peptide at a mz-bin needs to be counted so that mass is accepted in postprocessing}
}
\details{
  FILL ME OUT.
}
\value{
	A peak list, consisting of all local maxima, according to the 
	given parameterization.
}
%\references{ ~put references to the literature/web site here ~ }
\author{Marc Kirchner \email{marc.kirchner@iwr.uni-heidleberg.de}}
\examples{}
\keyword{multivariate}
