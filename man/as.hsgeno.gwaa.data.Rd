\name{as.hsgeno.gwaa.data}
\alias{as.hsgeno.gwaa.data}
\title{Attempts to convert gwaa.data to "hsgeno"}
\description{
	A function to convert @gtdata slot of an object of \code{\link{gwaa.data-class}}
	to "hsgeno" data frame
}
\usage{
	\method{as.hsgeno}{gwaa.data}(x, ...)
}
\arguments{
  \item{x}{An object of \code{\link{gwaa.data-class}}}
  \item{\dots}{...}
}
%\details{
%}
\value{
	A data-frame containing alleles, consumable by "haplo.stats" library
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
as.hsgeno(srdta[1:5,1:10])
}
\keyword{misc}% at least one, from doc/KEYWORDS
