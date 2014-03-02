\name{linCombDist}
\alias{linCombDist}
\title{ Linear Combination of Distances }
\description{
  This function computes linear combinations of distances.
}
\usage{
linCombDist(x, distfun1, w1, distfun2, w2, diag = FALSE, upper = FALSE)
}
\arguments{
  \item{x}{object which is passed to \code{distfun1} and \code{distfun2}.}
  \item{distfun1}{function used to compute an object of class \code{"dist"}.}
  \item{w1}{weight for result of \code{distfun1}.}
  \item{distfun2}{function used to compute an object of class \code{"dist"}.}
  \item{w2}{weight for result of \code{distfun2}.}
  \item{diag}{ see \code{\link[stats]{dist}} }
  \item{upper}{ see \code{\link[stats]{dist}} }
}
\details{
  This function computes and returns the distance matrix computed by
  a linear combination of two distance matrices.
}
\value{
  \code{linCombDist} returns an object of class \code{"dist"}; cf. \code{\link[stats]{dist}}.
}
\references{
  Flessa, F., Kehl, A., Kohl, M.
  Analysing diversity and community structures using PCR-RFLP: a new software application.
  Molecular Ecology Resources 2013 Jul; 13(4):726-33. 
}
\author{ Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\examples{
## assume a shift in the measured bands
M <- rbind(c(550, 500, 300, 250), c(510, 460, 260, 210),
           c(700, 650, 450, 400), c(550, 490, 310, 250))
dist(M)
diffDist(M)

## convex combination of dist and diffDist
linCombDist(M, distfun1 = dist, w1 = 0.5, distfun2 = diffDist, w2 = 0.5)

## linear combination
linCombDist(M, distfun1 = dist, w1 = 2, distfun2 = diffDist, w2 = 5)

## maximum distance
linCombDist(M, distfun1 = function(x) dist(x, method = "maximum"), w1 = 0.5, 
            distfun2 = function(x) diffDist(x, method = "maximum"), w2 = 0.5)
            
data(RFLPdata)
distfun <- function(x) linCombDist(x, distfun1 = dist, w1 = 0.1, distfun2 = diffDist, w2 = 0.9)
par(mfrow = c(2, 2))
plot(hclust(RFLPdist(RFLPdata, nrBands = 3, distfun = distfun)), cex = 0.7, cex.lab = 0.7)
RFLPplot(RFLPdata, nrBands = 3, distfun = distfun, mar.bottom = 6, cex.axis = 0.8)
plot(hclust(RFLPdist(RFLPdata, nrBands = 3)), cex = 0.7, cex.lab = 0.7)
RFLPplot(RFLPdata, nrBands = 3, mar.bottom = 6, cex.axis = 0.8)
}
\keyword{multivariate}
