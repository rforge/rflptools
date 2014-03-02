\name{RFLPcombine}
\alias{RFLPcombine}
\title{
  Combine RFLP data sets
}
\description{
  Function to combine an arbitrary number of RFLP data sets.
}
\usage{
RFLPcombine(\dots)
}
\arguments{
  \item{\dots}{ two or more data.frames with RFLP data.}
}
\details{
  The data sets are combined using \code{\link{rbind}}.
  
  If data sets with identical sample identifiers are given, the
  identifiers are made unique using \code{\link{make.unique}}.
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
res <- RFLPcombine(RFLPdata, RFLPdata, RFLPdata)
RFLPplot(res, nrBands = 4)
}
\keyword{manip}
