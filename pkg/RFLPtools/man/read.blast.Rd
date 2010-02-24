\name{read.blast}
\alias{read.blast}
\title{
  Read BLAST data
}
\description{
  Function to read BLAST data generated with standalone BLAST from NCBI.
}
\usage{
read.blast(file)
}
\arguments{
  \item{file}{ character: BLAST file to read in.}
}
\details{
  The function reads data which was generated with standalone BLAST from NCBI.
  Pairwise similarities of DNA sequences are calculated among all sequences 
  to analyse applying Standalone Blast
  (cf. ftp://ftp.ncbi.nih.gov/blast/executables/release/) 
  with the parameters -m 8 -r 2 -G 5 -E 2.
}
\value{
  A \code{data.frame} with variables
  \describe{
    \item{\code{query.id}}{character: sequence identifier. }
    \item{\code{subject.id}}{character: subject identifier. }
    \item{\code{identity}}{numeric: identity between sequences (in percent). }
    \item{\code{alignment.length}}{integer: number of nucleotides. }
    \item{\code{mismatches}}{integer: number of mismatches. }
    \item{\code{gap.opens}}{integer: number of gaps. }
    \item{\code{q.start}}{integer: query sequence start. }
    \item{\code{q.end}}{integer: query sequence end. }
    \item{\code{s.start}}{integer: subject sequence start. }
    \item{\code{s.end}}{integer: subject sequence end. }
    \item{\code{evalue}}{numeric: evalue. }
    \item{\code{bit.score}}{numeric: score value. }
  }
}
\references{
  Standalone Blast download: ftp://ftp.ncbi.nih.gov/blast/executables/release/
  
  Blast News: http://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastNews
}
\author{ 
  Fabienne Flessa \email{Fabienne.Flessa@uni-bayreuth.de},\cr
  Alexandra Kehl \email{Alexandra.Kehl@uni-bayreuth.de},\cr
  Matthias Kohl \email{Matthias.Kohl@stamats.de} 
}
%\note{}
\seealso{ \code{\link{BLASTdata}}, \code{\link{simMatrix}} }
\examples{
Dir <- system.file("extdata", package = "RFLPtools") # input directory 
filename <- file.path(Dir, "BLASTexample.txt")
BLAST1 <- read.blast(file = filename)
str(BLAST1)
}
\keyword{file}
\keyword{connection}
