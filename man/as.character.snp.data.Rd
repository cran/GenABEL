\name{as.character.snp.data}
\alias{as.character.snp.data}
\title{Attempts to convert snp.data to character}
\description{
	A function to convert an object of \code{\link{snp.data-class}}
	to "character"
}
\usage{
	\method{as.character}{snp.data}(x, ...)
%	as.character.snp.data(x, ...)
}
\arguments{
  \item{x}{An object of \code{\link{snp.data-class}}}
  \item{\dots}{...}
}
%\details{
%}
\value{
	A matrix containing genotypes in character format
}
%\references{}
\author{Yurii Aulchenko}
%\note{
%}
\seealso{
\code{\link{as.double.snp.data}},
\code{\link{as.hsgeno}},
\code{\link{as.genotype.snp.data}}
}
\examples{
require(GenABEL.data)
data(srdta)
as.character(srdta@gtdata[1:5,1:10])
}
\keyword{misc}% at least one, from doc/KEYWORDS
