\name{FragMatch}
\alias{FragMatch}
\title{
  Compute matches for RFLP data via FragMatch.
}
\description{
  Compute matches for RFLP data using FragMatch - a program for the analysis of DNA 
  fragment data.
}
\usage{
FragMatch(newData, refData, maxValue = 1000, errorBound = 25,
          weight = 1, na.rm = TRUE)
}
\arguments{
  \item{newData}{ data.frame with new RFLP data; see \code{\link{newDataGerm}}.}
  \item{refData}{ data.frame with reference RFLP data; see \code{\link{refDataGerm}}.}
  \item{maxValue}{ numeric: maximum value for which the error bound is applied. 
                   Can be a vector of length larger than 1.}
  \item{errorBound}{ numeric: error bound corresponding to \code{maxValue}. 
                     Can be a vector of length larger than 1.}
  \item{weight}{ numeric: weight for weighting partial matches; see details section.}
  \item{na.rm}{ logical: indicating whether NA values should be stripped before the 
                computation proceeds.}
}
\details{
  A rather simple algorithm which consists of counting the number of matches where
  it is considered a match if the value is inside a range of +/- \code{errorBound}.
  
  If there is more than one enzyme, one can use weights to give the partial perfect 
  matches for a certain enzyme a higher (or also smaller) weight.
}
\value{
  A character matrix with entries of the form \code{"a_b"} which means that there
  were \code{a} out of \code{b} possible matches.
}
\references{
  T. A. Saari,  S. K. Saari, C. D. Campbell, I. J Alexander, I. C. Anderson.
  FragMatch - a program for the analysis of DNA fragment data.
  Mycorrhiza 2007, 17:133-136
}
\author{ 
  Matthias Kohl \email{Matthias.Kohl@stamats.de} 
}
%\note{}
\seealso{ \code{\link{newDataGerm}}, \code{\link{refDataGerm}} }
\examples{
  data(refDataGerm)
  data(newDataGerm)
  
  res <- FragMatch(newDataGerm, refDataGerm)
}
\keyword{multivariate}
