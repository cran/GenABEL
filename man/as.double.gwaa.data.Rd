\name{as.double.gwaa.data}
\alias{as.double.gwaa.data}
\title{Attempts to convert gwaa.data to double}
\description{
	A function to convert an object of \code{\link{gwaa.data-class}}
	to "double"
}
\usage{
	\method{as.double}{gwaa.data}(x, ...)
}
\arguments{
  \item{x}{An object of \code{\link{gwaa.data-class}}}
  \item{\dots}{...}
}
%\details{
%}
\value{
	A matrix containing genotypes in double (numeric) format
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
as.double(srdta[1:5,1:10])
}
\keyword{misc}% at least one, from doc/KEYWORDS
