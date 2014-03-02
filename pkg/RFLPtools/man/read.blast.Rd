\name{read.blast}
\alias{read.blast}
\title{
  Read BLAST data
}
\description{
  Function to read BLAST data generated with standalone BLAST from NCBI.
}
\usage{
read.blast(file, sep = "\t")
}
\arguments{
  \item{file}{ character: BLAST file to read in.}
  \item{sep}{ the field separator character. Values on each line of the file are 
  separated by this character. Default \code{"\t"}.}
}
\details{
  The function reads data which was generated with standalone BLAST from NCBI;
  see ftp://ftp.ncbi.nih.gov/blast/executables/release/.
  
  Possible steps:\cr
  1) Install NCBI BLAST\cr
  2) Generate and import database(s)\cr
  3) Apply BLAST with options \code{outfmt} and \code{out}; e.g.\cr
  \code{blastn -query Testquery -db Testdatabase -outfmt 6 -out out.txt}\cr
  or\cr
  \code{blastn -query Testquery -db Testdatabase -outfmt 10 -out out.csv}\cr
  One can also call BLAST from inside R by using function \code{system}\cr
  system("blastn -query Testquery -db Testdatabase -outfmt 6 -out out.txt")\cr
  4) Read in the results\cr
  \code{test.res <- read.blast(file = "out.txt")}\cr
  or\cr
  \code{test.res <- read.blast(file = "out.csv", sep = ",")}
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
\seealso{ \code{\link{BLASTdata}}, \code{\link{simMatrix}} }
\examples{
Dir <- system.file("extdata", package = "RFLPtools") # input directory 
filename <- file.path(Dir, "BLASTexample.txt")
BLAST1 <- read.blast(file = filename)
str(BLAST1)
}
\keyword{file}
\keyword{connection}
