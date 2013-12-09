\name{as.genotype.snp.data}
\alias{as.genotype.snp.data}
\title{Attempts to convert snp.data to "genotype"}
\description{
	A function to convert an object of \code{\link{snp.data-class}}
	to "genotype" data frame
}
\usage{
	\method{as.genotype}{snp.data}(x, ...)
}
\arguments{
  \item{x}{An object of \code{\link{snp.data-class}}}
  \item{\dots}{...}
}
%\details{
%}
\value{
	A data-frame containing genotypes consumable by "genetics" library
}
%\references{}
\author{Yurii Aulchenko}
%\note{
%}
\seealso{
\code{\link{as.character.snp.data}},
\code{\link{as.double.snp.data}},
\code{\link{as.hsgeno}}
}
\examples{
require(GenABEL.data)
data(srdta)
as.genotype(srdta@gtdata[1:5,1:10])
}
\keyword{misc}% at least one, from doc/KEYWORDS
