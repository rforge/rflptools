\name{RFLPqc}
\alias{RFLPqc}
\title{
  Quality control for RFLP data
}
\description{
  Function to perform quality control for RFLP data based on a
  comparison between the total length of the digested PCR
  amplification product and the sum of the fragment lengths. If the
  sum is smaller or larger than the PCR amplification product
  (within a certain range to define), the samples can be excluded
  from further analyses. This function is helpful for data sets
  containig faint or uncertain bands. It is necessary to include 
  the total length of the PCR amplification product for each sample 
  as largest fragment in the data set, see \code{\link{RFLPdata}}.
}
\usage{
RFLPqc(x, rm.band1 = TRUE, QC.lo = 0.8, QC.up = 1.07, QC.rm = FALSE)
}
\arguments{
  \item{x}{ data.frame with RFLP data.}
  \item{rm.band1}{ logical: remove first band.}
  \item{QC.lo}{ numeric: a real number in (0,1).}
  \item{QC.up}{ numeric: a real number larger than 1.}
  \item{QC.rm}{ logical: remove samples with unsufficient quality.}
}
\details{
  In case the first band corresponds to the total length of the fragment one can perform
  a quality control comparing the length of the first band with the sum of the lengths
  of the remaining bands for each sample. If the sum is smaller than \code{QC.lo} times 
  the length of the first band or larger than \code{QC.up} times the length of the first 
  band, respectively, a text message is printed. 
  
  If \code{rm.band1 = TRUE} band 1 of all samples is removed and the remaining band 
  numbers are reduced by 1. 
  
  If \code{QC.rm = TRUE} samples of insufficient quality are entirely removed from the 
  given data and the resulting \code{data.frame} is returned.
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
\seealso{ \code{\link{RFLPdata}}, \code{\link{RFLPdist}} }
\examples{
Dir <- system.file("extdata", package = "RFLPtools") # input directory 
filename <- file.path(Dir, "AZ091016_report.txt")
RFLP1 <- read.rflp(file = filename)
str(RFLP1)

RFLP2 <- RFLPqc(RFLP1, rm.band1 = FALSE) # identical to RFLP1
identical(RFLP1, RFLP2)

RFLP3 <- RFLPqc(RFLP1)
str(RFLP3)

RFLP4 <- RFLPqc(RFLP1, rm.band1 = TRUE, QC.rm = TRUE)
str(RFLP4)
}
\keyword{manip}
