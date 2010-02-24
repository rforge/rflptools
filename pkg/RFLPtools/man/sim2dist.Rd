\name{sim2dist}
\alias{sim2dist}
\title{
  Convert similarity matrix to dist object.
}
\description{
  Function to convert similarity matrix to object of S3 class \code{"dist"}.
}
\usage{
sim2dist(x, maxSim = 1)
}
\arguments{
  \item{x}{ symmetric matrix: similarity matrix. }
  \item{maxSim}{ maximum similarity possible. }
}
\details{
  Similarity is converted to distance by \code{maxSim - x}. 
  The resulting matrix is converted to an object of S3 class \code{"dist"}
  by \code{as.dist}
}
\value{
  Object of S3 class \code{"dist"} is returned; see \code{\link[stats]{dist}}.
}
%\references{ ??? }
\author{ 
  Fabienne Flessa \email{Fabienne.Flessa@uni-bayreuth.de},\cr
  Alexandra Kehl \email{Alexandra.Kehl@uni-bayreuth.de},\cr
  Matthias Kohl \email{Matthias.Kohl@stamats.de} 
}
%\note{}
\seealso{ \code{\link{BLASTdata}}, \code{\link{simMatrix}} }
\examples{
data(BLASTdata)

## without sequence range
\dontrun{
res <- simMatrix(BLASTdata)
}

## with sequence range
range(BLASTdata$alignment.length)
res1 <- simMatrix(BLASTdata, sequence.range = TRUE, Min = 100, Max = 450)
res2 <- simMatrix(BLASTdata, sequence.range = TRUE, Min = 500)

## visualize similarity matrix
library(MKmisc)
corPlot(res2, minCor = 0, 
        labels = colnames(res2), title = "(Dis-)Similarity Plot")


## or
library(lattice)
myCol <- colorRampPalette(brewer.pal(8, "RdYlGn"))(128)
levelplot(res2, col.regions = myCol,
          at = do.breaks(c(0, max(res2)), 128),
          xlab = "", ylab = "",
          ## Rotate label of x axis
          scales = list(x = list(rot = 90)),
          main = "(Dis-)Similarity Plot")

## convert to distance
res.d <- sim2dist(res2)

## hierarchical clustering
plot(hclust(res.d))
}
\keyword{multivariate}
\keyword{cluster}
