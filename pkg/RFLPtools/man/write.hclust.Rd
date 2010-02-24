\name{write.hclust}
\alias{write.hclust}
\title{
  Cut a hierarchical cluster tree and write cluster identifiers to a 
  text file.
}
\description{
  The tree obtained by a hierarchical cluster analysis is cut into
  groups by using \code{\link[stats]{cutree}} and the results are
  exported to a text file.
}
\usage{
write.hclust(x, file, prefix, h = NULL, k = NULL, append = FALSE, dec = ",")
}
\arguments{
  \item{x}{ object of class \code{hclust}: result of hierarchical cluster analysis
    computed via function \code{\link[stats]{hclust}}.}
  \item{file}{ either a character string naming a file or a connection open
    for writing.  \code{""} indicates output to the console.}
  \item{prefix}{ character. Information about the cluster analysis.}
  \item{h}{numeric scalar or vector with heights where the tree should
    be cut.}
  \item{k}{ an integer scalar or vector with the desired number of groups.}
  \item{append}{ logical. Only relevant if \code{file} is a character string.  If
    \code{TRUE}, the output is appended to the file.  If \code{FALSE}, any
    existing file of the name is destroyed.}
  \item{dec}{ the string to use for decimal points in numeric or complex
    columns: must be a single character. }
}
\details{
  The results are written to file by a call to \code{\link[utils]{write.table}}
  where the columns in the resulting file are seperated by tabulators 
  (i.e. \code{sep="\t"}) and no row names are exported (i.e. \code{row.names = FALSE}).
}
%\references{ ??? }
\author{ 
  Fabienne Flessa \email{Fabienne.Flessa@uni-bayreuth.de},\cr
  Alexandra Kehl \email{Alexandra.Kehl@uni-bayreuth.de},\cr
  Matthias Kohl \email{Matthias.Kohl@stamats.de} 
}
%\note{}
\seealso{ \code{\link[utils]{write.table}}, \code{\link[stats]{cutree}} }
\examples{
data(RFLPdata)
res <- RFLPdist(RFLPdata, nrBands = 4)
cl <- hclust(res)
write.hclust(cl, file = "Test.txt", prefix = "Bd4", h = 50)

res <- RFLPdist2(RFLPdata, nrBands = 4, nrMissing = 1)
cl <- hclust(res)
write.hclust(cl, file = "Test.txt", append = TRUE, prefix = "Bd4_Mis1", h = 60)
}
\keyword{file}
\keyword{connection}
