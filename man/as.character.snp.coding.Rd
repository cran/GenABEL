\name{as.character.snp.coding}
\alias{as.character.snp.coding}
\title{Attempts to convert internal snp.coding-class to character}
\description{
	A function to convert an object of \code{\link{snp.coding-class}}
	to "character"
}
\usage{
	\method{as.character}{snp.coding}(x, ...)
}
\arguments{
  \item{x}{An object of \code{\link{snp.coding-class}}}
  \item{\dots}{...}
}
%\details{
%}
\value{
	A vector containing actual (nucleotide) coding, for corresponding 
	SNPs, in character format
}
%\references{}
\author{Yurii Aulchenko}
%\note{
%}
\seealso{
\code{\link{as.character.snp.strand}},
\code{\link{as.character.snp.data}},
\code{\link{as.double.snp.data}},
\code{\link{as.hsgeno}},
\code{\link{as.genotype.snp.data}}
}
\examples{
require(GenABEL.data)
data(srdta)
as.character(srdta@gtdata@coding[1:5])
}
\keyword{misc}% at least one, from doc/KEYWORDS
