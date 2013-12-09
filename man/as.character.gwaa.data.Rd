\name{as.character.gwaa.data}
\alias{as.character.gwaa.data}
\title{Attempts to convert genotypic part of gwaa.data to character}
\description{
	A function to convert @gtdata slot of an object of 
	\code{\link{gwaa.data-class}} to "character"
}
\usage{
	\method{as.character}{gwaa.data}(x, ...)
}
\arguments{
  \item{x}{An object of \code{\link{gwaa.data-class}}}
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
\code{\link{as.character.snp.data}},
\code{\link{as.double.gwaa.data}},
\code{\link{as.double.snp.data}},
\code{\link{as.hsgeno}},
\code{\link{as.genotype.gwaa.data}},
\code{\link{as.genotype.snp.data}}
}
\examples{
require(GenABEL.data)
data(srdta)
as.character(srdta[1:5,1:10])
}
\keyword{misc}% at least one, from doc/KEYWORDS
