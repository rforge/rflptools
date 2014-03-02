\name{read.rflp}
\alias{read.rflp}
\title{
  Read RFLP data
}
\description{
  Function to read RFLP data (e.g. generated with software package Gene
  Profiler 4.05 (Scanalytics Inc.)) for DNA fragment analysis and
  genotyping, and exported to a text file.
}
\usage{
read.rflp(file)
}
\arguments{
  \item{file}{ character: RFLP file to read in.}
}
\details{
  The function reads data from a text file which was generated e.g. with the 
  software package Gene Profiler 4.05 (Scanalytics Inc.) for DNA fragment
  analysis and genotyping. The data file contains sample identifier (\code{Sample}), 
  band number (\code{Band}), molecular weight (\code{MW}) and gel identifier (\code{Gel}) 
  (see \code{\link{RFLPdata}}).
  
  If gel identifier \code{Gel} is missing it is extracted from the sample identifier 
  \code{Sample}.
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
filename <- file.path(Dir, "RFLPexample.txt")
RFLP1 <- read.rflp(file = filename)
str(RFLP1)

filename <- file.path(Dir, "AZ091016_report.txt")
RFLP2 <- read.rflp(file = filename)
str(RFLP2)
}
\keyword{file}
\keyword{connection}
