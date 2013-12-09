\name{as.genotype.gwaa.data}
\alias{as.genotype.gwaa.data}
\title{Attempts to convert gwaa.data to "genotype"}
\description{
	A function to convert @gtdata slot of an object of \code{\link{gwaa.data-class}}
	to "genotype" data frame
}
\usage{
	\method{as.genotype}{gwaa.data}(x, ...)
}
\arguments{
  \item{x}{An object of \code{\link{gwaa.data-class}}}
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
\code{\link{as.character.gwaa.data}},
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
as.genotype(srdta[1:5,1:10])
}
\keyword{misc}% at least one, from doc/KEYWORDS
