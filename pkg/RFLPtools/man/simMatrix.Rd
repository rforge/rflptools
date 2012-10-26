\name{simMatrix}
\alias{simMatrix}
\title{
  Similarity matrix for BLAST data.
}
\description{
  Function to compute similarity matrix for all-vs-all BLAST results of rDNA sequences 
  generated with standalone BLAST from NCBI or local BLAST implemented in BioEdit.
}
\usage{
simMatrix(x, sequence.range = FALSE, Min, Max)
}
\arguments{
  \item{x}{ data.frame with BLAST data; see \code{\link{BLASTdata}}.}
  \item{sequence.range}{ logical: use sequence range.}
  \item{Min}{ minimum sequence length. }
  \item{Max}{ maximum sequence length. }
}
\details{
    The given BLAST data is used to compute a similarity matrix using the following algorithm:
    First, the length of each sequence (LS) comprised in the input data file is extracted. 
    If there is more than one comparison for one sequence including different parts of the 
    respective sequence, that one with maximum base length is chosen. 
    Subsequently, the number of matching bases (mB) is calculated by multiplying 
    two variables comprised in the BLAST output: the identity between sequences (\%) and 
    the number of nucleotides divided by 100. The, resulting value is rounded to integer. 
    Furthermore, the similarity is calculated by dividing mB by LS. 
    Finally, the similarity matrix including all sequences is built. 
    If the similarity of a combination is not shown in the BLAST report file (because 
    the similarity was lower than 70\%), this comparison is included in the similarity 
    matrix with the result zero.
}
\value{
  Similarity matrix.
}
\references{
  Standalone Blast download: ftp://ftp.ncbi.nih.gov/blast/executables/release/
  
  Blast News: http://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastNews
  
  BioEdit v7.0.9: Tom Hall, Ibis Biosciences; http://www.mbio.ncsu.edu/BioEdit/bioedit.html

  Persoh, D., Melcher, M., Flessa, F., Rambold, G.: 
  First fungal community analyses of endophytic ascomycetes associated with Viscum album ssp. austriacum 
  and itshost Pinus sylvestris. Fungal Biology 2010, in press.
}
\author{ 
  Fabienne Flessa \email{Fabienne.Flessa@uni-bayreuth.de},\cr
  Alexandra Kehl \email{Alexandra.Kehl@botgarten.uni-tuebingen.de},\cr
  Matthias Kohl \email{Matthias.Kohl@stamats.de} 
}
%\note{}
\seealso{ \code{\link{BLASTdata}}, \code{\link{sim2dist}} }
\examples{
data(BLASTdata)

## without sequence range
## code takes some time
\dontrun{
res <- simMatrix(BLASTdata)
}

## with sequence range
range(BLASTdata$alignment.length)
res1 <- simMatrix(BLASTdata, sequence.range = TRUE, Min = 100, Max = 450)
res2 <- simMatrix(BLASTdata, sequence.range = TRUE, Min = 500)
}
\keyword{multivariate}
\keyword{cluster}

