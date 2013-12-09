\name{as.genotype}
\alias{as.genotype}
\title{Attempts to convert object to "genotype"}
\description{
	A function to convert an object 
	to "genotype" data frame
}
\usage{
	as.genotype(x, ...)
}
\arguments{
  \item{x}{An object of \code{\link{snp.data-class}}}
  \item{\dots}{...}
}
%\details{
%}
\value{
	A data-frame containing "genotype" data class, consumable by "genetics" library
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
as.genotype(srdta@gtdata[1:5,1:10])
}
\keyword{misc}% at least one, from doc/KEYWORDS
