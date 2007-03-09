\name{as.hsgeno}
\alias{as.hsgeno}
\title{Attempts to convert object to "hsgeno"}
\description{
	A function to convert an object 
	to "hsgeno" data frame
}
\usage{
	as.hsgeno(x, ...)
}
\arguments{
  \item{x}{An object of \code{\link{snp.data-class}}}
  \item{\dots}{...}
}
\details{
}
\value{
	A data-frame containing alleles, consumable by "haplo.stats" library
}
\references{}
\author{Yurii Aulchenko}
\note{
}
\seealso{
\code{\link{as.character.snp.data}},
\code{\link{as.double.snp.data}},
\code{\link{as.genotype.snp.data}}
}
\examples{
data(srdta)
as.hsgeno(srdta@gtdata[1:5,1:10])
}
\keyword{misc}% at least one, from doc/KEYWORDS
