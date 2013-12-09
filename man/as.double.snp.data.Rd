\name{as.double.snp.data}
\alias{as.double.snp.data}
\title{Attempts to convert snp.data to double}
\description{
	A function to convert an object of \code{\link{snp.data-class}}
	to "double"
}
\usage{
	\method{as.double}{snp.data}(x, ...)
}
\arguments{
  \item{x}{An object of \code{\link{snp.data-class}}}
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
\code{\link{as.character.snp.data}},
\code{\link{as.hsgeno}},
\code{\link{as.genotype.snp.data}}
}
\examples{
require(GenABEL.data)
data(srdta)
as.double(srdta@gtdata[1:5,1:10])
}
\keyword{misc}% at least one, from doc/KEYWORDS
