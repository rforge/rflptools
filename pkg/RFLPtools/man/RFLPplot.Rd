\name{RFLPplot}
\alias{RFLPplot}
\title{
  Function to plot RFLP data.
}
\description{
  Given RFLP data is plotted where the samples are sorted according
  to the corresponding dendrogram.
}
\usage{
RFLPplot(x, nrBands, nrMissing, distfun = dist, 
         hclust.method = "complete", mar.bottom = 5, 
         cex.axis = 0.5, colBands)
}
\arguments{
  \item{x}{ data.frame with RFLP data; see \code{\link{RFLPdata}}.}
  \item{nrBands}{ if not missing, then only samples with the specified number
  of bands are considered. }
  \item{nrMissing}{ if not missing, then it is assumed that some bands may be missing.
  That is, all samples with number of bands in nrBands, nrBands+1, ..., nrBands+nrMissing
  are considered.}
  \item{distfun}{ function computing the distance with default \code{dist};
  see \code{\link[stats]{dist}}.}
  \item{hclust.method}{ method used for hierarchical clustering;
  see \code{\link[stats]{hclust}}.}
  \item{mar.bottom}{ bottom margin of the plot; see \code{\link[graphics]{par}}. }
  \item{cex.axis}{ size of the x-axis annotation. }
  \item{colBands}{ color for the bands. Has to be of length 1 or number of samples. 
  If missing, \code{"Set1"} of \pkg{RColorBrewer} is used; see
  \code{\link[RColorBrewer]{brewer.pal}}. }
}
\details{
  RFLP data is plotted. The samples are sorted according to the corresponding
  dendrogram which is computed via function \code{hclust}.
  
  The option to specify \code{nrMissing} may be useful, if gel image quality is low, 
  and the detection of bands is doubtful.
}
\value{invisible}
%\references{ ??? }
\author{ 
  Fabienne Flessa \email{Fabienne.Flessa@uni-bayreuth.de},\cr
  Alexandra Kehl \email{Alexandra.Kehl@botgarten.uni-tuebingen.de},\cr
  Matthias Kohl \email{Matthias.Kohl@stamats.de} 
}
%\note{}
\seealso{ \code{\link{RFLPdata}}, \code{\link[stats]{dist}} }
\examples{
data(RFLPdata)
par(mfrow = c(1,2))
plot(hclust(RFLPdist(RFLPdata, nrBands = 3)), cex = 0.7)
RFLPplot(RFLPdata, nrBands = 3, mar.bottom = 6, cex.axis = 0.8)

par(mfrow = c(1,2))
plot(hclust(RFLPdist2(RFLPdata, nrBands = 9, nrMissing = 1)), cex = 0.7)
RFLPplot(RFLPdata, nrBands = 9, nrMissing = 1, mar.bottom = 6, cex.axis = 0.8)


distfun <- function(x) dist(x, method = "maximum")
par(mfrow = c(1,2))
plot(hclust(RFLPdist(RFLPdata, nrBands = 3, distfun = distfun), 
            method = "average"), cex = 0.7, cex.lab = 0.7)
RFLPplot(RFLPdata, nrBands = 3, distfun = distfun, hclust.method = "average", 
         mar.bottom = 6, cex.axis = 0.8)
}
\keyword{multivariate}
\keyword{cluster}
