\name{refDataGerm}
\alias{refDataGerm}
\docType{data}
\title{ Example data set from GERM software }
\description{
  This is the reference data taken from the GERM software. 
}
\usage{data(refDataGerm)}
\format{
  A data frame with 250 observations on the following six variables
  \describe{
    \item{\code{Sample}}{character: sample identifier. }
    \item{\code{Enzyme}}{character: enzyme used.}
    \item{\code{Band}}{integer: band number. }
    \item{\code{MW}}{integer: molecular weight. }
    \item{\code{Genus}}{character: genus of sample. }
    \item{\code{Species}}{character: species of sample. }
  }
}
\details{
  See GERM software.
}
\source{
  The data set was taken from the GERM software (table 'Example Data').
}
\references{
  Ian A. Dickie, Peter G. Avis, David J. McLaughlin, Peter B. Reich. 
  Good-Enough RFLP Matcher (GERM) program.
  Mycorrhiza 2003, 13:171-172.
}
\examples{
data(refDataGerm)
str(refDataGerm)
}
\keyword{datasets}
