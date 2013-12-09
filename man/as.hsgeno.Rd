\name{as.hsgeno}
\alias{as.hsgeno}
\title{Attempts to convert object to "hsgeno"}
\description{
	A function to convert an object 
	to "hsgeno" data frame, to be used by 
	"haplo.stats" library
}
\usage{
	as.hsgeno(x, ...)
}
\arguments{
  \item{x}{An object of \code{\link{snp.data-class}} or \code{\link{gwaa.data-class}}}
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
\code{\link{as.character.snp.data}},
\code{\link{as.double.snp.data}},
\code{\link{as.genotype.snp.data}}
}
\examples{
require(GenABEL.data)
data(srdta)
as.hsgeno(srdta[1:5,1:3])
as.hsgeno(srdta@gtdata[1:5,1:3])
}
\keyword{misc}% at least one, from doc/KEYWORDS
