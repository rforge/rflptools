\name{RFLPrefplot}
\alias{RFLPrefplot}
\title{
  Function for a visual comparison of RFLP samples with reference samples.
}
\description{
  Given RFLP samples are plotted together with reference samples and
  sorted by their distance to the reference sample.
}
\usage{
RFLPrefplot(x, ref, distfun = dist, nrBands, mar.bottom = 5, 
            cex.main = 1.2, cex.axis = 0.5, devNew = FALSE, 
            colBands)
}
\arguments{
  \item{x}{ data.frame with RFLP data; e.g. \code{\link{RFLPdata}}.}
  \item{ref}{ data.frame with RFLP reference data; e.g. \code{\link{RFLPref}}.}
  \item{distfun}{ function computing the distance with default \code{dist};
  see \code{\link[stats]{dist}}.}
  \item{nrBands}{ if not missing, then only samples with the specified number
  of bands are considered. }
  \item{mar.bottom}{ bottom margin of the plot; see \code{\link[graphics]{par}}. }
  \item{cex.main}{ size of the plot title. }
  \item{cex.axis}{ size of the x-axis annotation. }
  \item{devNew}{ logical. Open new graphics device for each plot. }
  \item{colBands}{ color for the bands. Has to be of length 1 or number of samples. 
  If missing, \code{"Set1"} of \pkg{RColorBrewer} is used; see
  \code{\link[RColorBrewer]{brewer.pal}}. }
}
\details{
  Given RFLP samples are plotted together with reference samples and
  sorted by their distance to the reference sample.
}
\value{invisible}
%\references{ ??? }
\author{ 
  Fabienne Flessa \email{Fabienne.Flessa@uni-bayreuth.de},\cr
  Alexandra Kehl \email{Alexandra.Kehl@botgarten.uni-tuebingen.de},\cr
  Matthias Kohl \email{Matthias.Kohl@stamats.de} 
}
%\note{}
\seealso{ \code{\link{RFLPplot}} }
\examples{
data(RFLPdata)
data(RFLPref)
dev.new(width = 12)
RFLPrefplot(RFLPdata, RFLPref, nrBands = 4, cex.axis = 0.5)
dev.off()

RFLPrefplot(RFLPdata, RFLPref, nrBands = 6, cex.axis = 0.8)
RFLPrefplot(RFLPdata, RFLPref, nrBands = 9, cex.axis = 0.8)

RFLPrefplot(RFLPdata, RFLPref[RFLPref$Sample == "Ni_29_A3",], nrBands = 4, cex.axis = 0.7)

Dir <- system.file("extdata", package = "RFLPtools") # input directory 
filename <- file.path(Dir, "AZ091016_report.txt")
RFLP1 <- read.rflp(file = filename)
RFLP2 <- RFLPqc(RFLP1)

dev.new(width = 12)
RFLPrefplot(RFLP1, RFLPref, nrBands = 4, cex.axis = 0.8)
dev.off()

RFLPrefplot(RFLP1, RFLPref, nrBands = 5, cex.axis = 0.8)
}
\keyword{multivariate}
