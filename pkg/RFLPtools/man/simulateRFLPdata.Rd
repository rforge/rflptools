\name{simulateRFLPdata}
\alias{simulateRFLPdata}
\title{
  Simulate RFLP data.
}
\description{
  Simulates RFLP data for comparions of algorithms.
}
\usage{
simulateRFLPdata(N = 10, nrBands = 3:12, bandCenters = seq(100, 800, by = 100),
                 delta = 50, refData = FALSE)
}
\arguments{
  \item{N}{ integer: number samples which shall be simulated per number of bands.}
  \item{nrBands}{ integer: vector of number of bands.}
  \item{bandCenters}{ numeric: vector of band centers.}
  \item{delta}{ numeric: uniform distribution with \code{min = bandCenter - delta} and \code{max = bandCenter + delta}
                is used.}
  \item{refData}{ logical: if TRUE, additonal columns \code{Taxonname} and \code{Accesion} are generated.}
}
\details{
  The function can be used to simulate RFLP data. For every number of band specified in 
  \code{nrBands} a total number of \code{N} samples are generated. 
  
  First the band centers are randomly selected (with replacement) from \code{bandCenter} 
  which form the centers of intervals of length \code{2*delta}. From these intervals 
  uniform random numbers are drawn leading to randomly generated RFLP data.
}
\value{
  A data frame with \code{N*length(nrBands)} observations on the following four variables
  \describe{
    \item{\code{Sample}}{character: sample identifier. }
    \item{\code{Band}}{integer: band number. }
    \item{\code{MW}}{integer: molecular weight. }
    \item{\code{Enzyme}}{character: enzyme name. }
  }
  is generated. If \code{refData = TRUE} then the following two additional variables
  are added.
  \describe{
    \item{\code{Taxonname}}{character: taxon name. }
    \item{\code{Accession}}{character: accession number. }
  }
}
%\references{}
\author{ 
  Matthias Kohl \email{Matthias.Kohl@stamats.de} 
}
%\note{}
\seealso{ \code{\link{RFLPdata}}, \code{\link{RFLPref}} }
\examples{
  simData <- simulateRFLPdata()
}
\keyword{multivariate}
