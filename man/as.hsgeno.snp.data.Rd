\name{as.hsgeno.snp.data}
\alias{as.hsgeno.snp.data}
\title{Attempts to convert snp.data to "hsgeno"}
\description{
	A function to convert an object of \code{\link{snp.data-class}}
	to "hsgeno" data frame
}
\usage{
	\method{as.hsgeno}{snp.data}(x, ...)
}
\arguments{
  \item{x}{An object of \code{\link{snp.data-class}}}
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
as.hsgeno(srdta@gtdata[1:5,1:10])
}
\keyword{misc}% at least one, from doc/KEYWORDS
