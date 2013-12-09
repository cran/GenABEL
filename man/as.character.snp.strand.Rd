\name{as.character.snp.strand}
\alias{as.character.snp.strand}
\title{Attempts to convert internal strand-class to character}
\description{
	A function to convert an object of \code{\link{snp.strand-class}}
	to "character"
}
\usage{
	\method{as.character}{snp.strand}(x, ...)
}
\arguments{
  \item{x}{An object of \code{\link{snp.strand-class}}}
  \item{\dots}{...}
}
%\details{
%}
\value{
	A vector containing strand ("+", "-" or "u"), for corresponding 
	SNPs, in character format
}
%\references{}
\author{Yurii Aulchenko}
%\note{
%}
\seealso{
\code{\link{as.character.snp.coding}},
\code{\link{as.character.snp.data}},
\code{\link{as.double.snp.data}},
\code{\link{as.hsgeno}},
\code{\link{as.genotype.snp.data}}
}
\examples{
require(GenABEL.data)
data(srdta)
as.character(srdta@gtdata@strand[1:5])
}
\keyword{misc}% at least one, from doc/KEYWORDS
