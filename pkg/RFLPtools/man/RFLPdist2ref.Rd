\name{RFLPdist2ref}
\alias{RFLPdist2ref}
\title{
  Compute distance between RFLP data and RFLP reference data.
}
\description{
  Function to compute distance between RFLP data and RFLP reference data.
}
\usage{
RFLPdist2ref(x, ref, distfun = dist, nrBands, LOD = 0)
}
\arguments{
  \item{x}{ data.frame with RFLP data; e.g. \code{\link{RFLPdata}}.}
  \item{ref}{ data.frame with RFLP reference data; e.g. \code{\link{RFLPref}}.}
  \item{distfun}{ function computing the distance with default \code{dist};
  cf. \code{\link[stats]{dist}}.}
  \item{nrBands}{ only samples and reference samples with this number
  of bands are considered. }
  \item{LOD}{ threshold for low-bp bands. }
}
\details{
  For each sample with \code{nrBands} bands the distance to each reference 
  sample with \code{nrBands} bands is computed. The result is a matrix with 
  the corresponding distances where rows represent the samples and columns
  the reference samples.
  
  If \code{LOD > 0} is specified, all values below \code{LOD} are removed before the 
  distances are calculated. This applies to \code{x} and \code{ref}.
}
\value{
  A matrix with distances.
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
\seealso{ \code{\link{RFLPdata}}, \code{\link[stats]{dist}} }
\examples{
## Euclidean distance
data(RFLPdata)
data(RFLPref)
nrBands(RFLPref)
RFLPdist2ref(RFLPdata, RFLPref, nrBands = 4)
RFLPdist2ref(RFLPdata, RFLPref, nrBands = 6)

Dir <- system.file("extdata", package = "RFLPtools") # input directory 
filename <- file.path(Dir, "AZ091016_report.txt")
RFLP1 <- read.rflp(file = filename)
RFLP2 <- RFLPqc(RFLP1)
nrBands(RFLP2)
RFLPdist2ref(RFLP1, RFLPref, nrBands = 4)
RFLPdist2ref(RFLP1, RFLPref, nrBands = 5)
}
\keyword{multivariate}
\keyword{cluster}
