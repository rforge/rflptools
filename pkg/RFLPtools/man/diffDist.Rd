\name{diffDist}
\alias{diffDist}
\title{ Distance Matrix Computation }
\description{
  This function computes and returns the distance matrix computed by
  using the specified distance measure to compute the distances between
  the rows of a data matrix. Instead of the row values as in the case of 
  \code{\link[stats]{dist}}, the successive differences of the row values
  are used.
}
\usage{
diffDist(x, method = "euclidean", diag = FALSE, upper = FALSE, p = 2)
}
\arguments{
  \item{x}{a numeric matrix, data frame or \code{"dist"} object.}
  \item{method}{the distance measure to be used. This must be one of
    \code{"euclidean"}, \code{"maximum"}, \code{"manhattan"},
    \code{"canberra"}, \code{"binary"} or \code{"minkowski"}.
    Any unambiguous substring can be given.}
  \item{diag}{logical value indicating whether the diagonal of the
    distance matrix should be printed by \code{print.dist}.}
  \item{upper}{logical value indicating whether the upper triangle of the
    distance matrix should be printed by \code{print.dist}.}
  \item{p}{The power of the Minkowski distance.}
}
\details{
  This function computes and returns the distance matrix computed by
  using the specified distance measure to compute the distances between
  the rows of a data matrix. Instead of the row values as in the case of 
  \code{\link[stats]{dist}}, the successive differences of the row values
  are used.
  
  It's a simple wrapper function arround \code{\link[stats]{dist}}. For
  more details about the distances we refer to \code{\link[stats]{dist}}.
  
  The function may be helpful, if there is a shift w.r.t.\ the measured
  bands; e.g.\ \code{c(550, 500, 300, 250)} vs.\ \code{c(510, 460, 260, 210)}.
}
\value{
  \code{diffDist} returns an object of class \code{"dist"}; cf. \code{\link[stats]{dist}}.
}
%\references{}
\author{ Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\examples{
## assume a shift in the measured bands
M <- rbind(c(550, 500, 300, 250), c(510, 460, 260, 210),
           c(550, 500, 300, 200))
dist(M)
diffDist(M)
}
\keyword{multivariate}
