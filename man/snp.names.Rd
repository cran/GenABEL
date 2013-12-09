\name{snp.names}
\alias{snp.names}
\title{extracts names of SNPs in a region}
\description{
Based on boundary conditions specified and (or) chromosome 
selects SNP names in the region
}
\usage{
snp.names(data, begin, end, chromosome)
}
\arguments{
  \item{data}{object of class \code{\link{gwaa.data-class}}, \code{\link{snp.data-class}}, 
		\code{\link{scan.gwaa-class}} or \code{\link{check.marker-class}}}
  \item{begin}{Start position (or name of the first SNP)}
  \item{end}{End-position or name of last SNP}
  \item{chromosome}{Chromosome code}
}
\details{
	Any of the arguments, except the \code{data} can be missing
}
\value{
	A vector of names of SNPs located in the region
}
%\references{}
\author{Yurii Aulchenko}
%\note{}
\seealso{\code{\link{snp.data-class}}}
\examples{
require(GenABEL.data)
data(srdta)
snp.names(srdta, begin = 50000, end = 100000)
snp.names(srdta, begin = 50000, end = 100000, chromosome = "1")

# does not make sense with these data:
snp.names(srdta, begin = 50000, end = 100000, chromosome = "X") 

# again makes sense: 
snp.names(srdta, end = 100000)
snp.names(srdta, begin = 2200000)

# show summary for SNPs in region between 50,000 and 100,000
a <- snp.names(srdta, begin = 50000, end = 100000)
summary(srdta@gtdata[,a])
}
\keyword{misc}
