\name{RFLPlod}
\alias{RFLPlod}
\title{
  Remove bands below LOD
}
\description{
  Function to exclude bands below a given LOD.
}
\usage{
RFLPlod(x, LOD)
}
\arguments{
  \item{x}{ data.frame with RFLP data.}
  \item{LOD}{ threshold for low-bp bands.}
}
\details{
  Low-bp bands may be regarded as unreliable. Function 
  \code{RFLPlod} can be used to exclude such bands, which
  are likely to be absent in some other samples, before
  further analyses.
}
\value{
  A \code{data.frame} with variables
  \describe{
    \item{\code{Sample}}{character: sample identifier. }
    \item{\code{Band}}{integer: band number. }
    \item{\code{MW}}{integer: molecular weight. }
    \item{\code{Gel}}{character: gel identifier. }
  }
}
\references{
  Flessa, F., Kehl, A., Kohl, M.
  Analysing diversity and community structures using PCR-RFLP: a new software application.
  Molecular Ecology Resources 2013 Jul; 13(4):726-33. 
}
\author{ 
  Fabienne Flessa \email{Fabienne.Flessa@uni-bayreuth.de},\cr
  Alexandra Kehl \email{Alexandra.Kehl@botgarten.uni-tuebingen.de},\cr
  Matthias Kohl \email{Matthias.Kohl@stamats.de} 
}
%\note{}
\seealso{ \code{\link{RFLPdata}} }
\examples{
data(RFLPdata)
## remove bands with MW smaller than 60
RFLPdata.lod <- RFLPlod(RFLPdata, LOD = 60)
par(mfrow = c(1, 2))
RFLPplot(RFLPdata, nrBands = 4, ylim = c(40, 670))
RFLPplot(RFLPdata.lod, nrBands = 4, ylim = c(40, 670))
title(sub = "After applying RFLPlod")
}
\keyword{manip}
