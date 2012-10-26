\name{RFLPdist}
\alias{RFLPdist}
\title{
  Compute distances for RFLP data.
}
\description{
  Within each group containing RFLP-samples exhibiting a equal number of bands, 
  the distance between the molecular weights is computed.
}
\usage{
RFLPdist(x, distfun = dist, nrBands)
}
\arguments{
  \item{x}{ data.frame with RFLP data; see \code{\link{RFLPdata}}.}
  \item{distfun}{ function computing the distance with default \code{dist};
  cf. \code{\link[stats]{dist}}.}
  \item{nrBands}{ if not missing, then only samples with the specified number
  of bands are considered. }
}
\details{
  For each number of bands the given distance between the molecular weights is computed. 
  The result is a named list of distances where the names correspond to the number of 
  bands which occur in each group.
  
  If \code{nrBands} is specified only samples with this number of bands are considered.
}
\value{
  A named list with the distances; see \code{\link[stats]{dist}}.
  
  In case \code{nrBands} is not missing, an object of S3 class \code{dist}.
}
\references{
  Poussier, Stephane; Trigalet-Demery, Danielle; Vandewalle, Peggy; Goffinet, Bruno; Luisetti, Jacques; Trigalet, Andre.
  Genetic diversity of Ralstonia solanacearum as assessed by PCR-RFLP of the hrp gene region, AFLP and 16S rRNA sequence 
  analysis, and identification of an African subdivision.
  Microbiology 2000 146:1679-1692 
  
  Matsumoto, Masaru; Furuya, Naruto; Takanami, Yoichi; Matsuyama, Nobuaki.  
  RFLP analysis of the PCR-amplified 28S rDNA in Rhizoctonia solani. 
  Mycoscience 1996 37:351 - 356
}
\author{ 
  Fabienne Flessa \email{Fabienne.Flessa@uni-bayreuth.de},\cr
  Alexandra Kehl \email{Alexandra.Kehl@botgarten.uni-tuebingen.de},\cr
  Matthias Kohl \email{Matthias.Kohl@stamats.de} 
}
%\note{}
\seealso{ \code{\link{RFLPdata}}, \code{\link[stats]{dist}} }
\examples{
## Euclidean distance
data(RFLPdata)
res <- RFLPdist(RFLPdata)
names(res) ## number of bands
res$"6"

RFLPdist(RFLPdata, nrBands = 6)

## Other distances
res1 <- RFLPdist(RFLPdata, distfun = function(x) dist(x, method = "manhattan"))
res2 <- RFLPdist(RFLPdata, distfun = function(x) dist(x, method = "maximum"))
res[[1]]
res1[[1]]
res2[[1]]

## cut dendrogram at height 50
clust4bd <- hclust(res[[2]])
cgroups50 <- cutree(clust4bd, h=50)
cgroups50

## or
library(MKmisc)
res3 <- RFLPdist(RFLPdata, distfun = corDist)
res3$"9"

## hierarchical clustering
par(mfrow = c(2,2))
plot(hclust(res[[1]]), main = "Euclidean distance")
plot(hclust(res1[[1]]), main = "Manhattan distance")
plot(hclust(res2[[1]]), main = "Maximum distance")
plot(hclust(res3[[1]]), main = "Pearson correlation distance")


## Similarity matrix
library(MKmisc)
myCol <- colorRampPalette(brewer.pal(8, "RdYlGn"))(128)
ord <- order.dendrogram(as.dendrogram(hclust(res[[1]])))
temp <- as.matrix(res[[1]])
simPlot(temp[ord,ord], col = rev(myCol), minVal = 0, 
        labels = colnames(temp), title = "(Dis-)Similarity Plot")


## or
library(lattice)
levelplot(temp[ord,ord], col.regions = rev(myCol),
          at = do.breaks(c(0, max(temp)), 128),
          xlab = "", ylab = "",
          ## Rotate label of x axis
          scales = list(x = list(rot = 90)),
          main = "(Dis-)Similarity Plot")

## multidimensional scaling
loc <- cmdscale(res[[5]])
x <- loc[,1]
y <- -loc[,2]
plot(x, y, type="n", xlab="", ylab="", xlim = 1.05*range(x), main="Multidemsional scaling")
text(x, y, rownames(loc), cex=0.8)
}
\keyword{multivariate}
\keyword{cluster}
