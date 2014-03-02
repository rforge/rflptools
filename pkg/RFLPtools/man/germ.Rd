\name{germ}
\alias{germ}
\title{
  Compute matches for RFLP data via GERM.
}
\description{
  Compute matches for RFLP data using the Good-Enough RFLP Matcher (GERM) program.
}
\usage{
germ(newData, refData, parameters = list("Max forward error" = 25,
                                         "Max backward error" = 25,
                                         "Max sum error" = 100,
                                         "Lower measurement limit" = 100), 
     method = "joint", na.rm = TRUE)
}
\arguments{
  \item{newData}{ data.frame with new RFLP data; see \code{\link{newDataGerm}}.}
  \item{refData}{ data.frame with reference RFLP data; see \code{\link{refDataGerm}}.}
  \item{parameters}{ list of the four program parameters of GERM; see details section.}
  \item{method}{ matching and ranking method used for computation; see details section.}
  \item{na.rm}{ logical: indicating whether NA values should be stripped before the computation proceeds.}
}
\details{
  There are four matching and ranking methods which are \code{"joint"}, \code{"forward"}, 
  \code{"backward"}, and \code{"sum"}. For more details see Dickie et al. (2003).
  
  The parameters of the GERM software are:
  \code{"Max forward error"}: Used if "matching and ranking method" is set to \code{"forward"} or \code{"joint"}.
  \code{"Max backward error"}: Used if "matching and ranking method" is set to \code{"backward"} or \code{"joint"}.
  \code{"Max sum error"}: Used for matching if "matching and ranking method" is set to \code{"sum"}.
  \code{"Lower measurement limit"}: The lower bound of measurements (often 100 or 50, depending on ladder used).
}
\value{
  A named list with the results.
}
\references{
  Ian A. Dickie, Peter G. Avis, David J. McLaughlin, Peter B. Reich. 
  Good-Enough RFLP Matcher (GERM) program.
  Mycorrhiza 2003, 13:171-172.
}
\author{ 
  Matthias Kohl \email{Matthias.Kohl@stamats.de} 
}
%\note{}
\seealso{ \code{\link{newDataGerm}}, \code{\link{refDataGerm}} }
\examples{
  data(refDataGerm)
  data(newDataGerm)
  
  ## Example 1
  res1 <- germ(newDataGerm[1:7,], refDataGerm)
  
  ## Example 2
  res2 <- germ(newDataGerm[8:15,], refDataGerm)
  
  ## Example 3
  res3 <- germ(newDataGerm[16:20,], refDataGerm)
  
  ## all three examples in one step
  res.all <- germ(newDataGerm, refDataGerm)
}
\keyword{multivariate}
