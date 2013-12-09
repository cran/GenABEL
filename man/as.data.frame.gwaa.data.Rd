\name{as.data.frame.gwaa.data}
\alias{as.data.frame.gwaa.data}
\title{Attempts to convert snp.data to "hsgeno"}
\description{
	A function taking @phdata part (data.frame) of the object of 
	\code{\link{gwaa.data-class}} 
}
\usage{
	\method{as.data.frame}{gwaa.data}(x, ...)
}
\arguments{
  \item{x}{An object of \code{\link{data.frame-class}}}
  \item{\dots}{...}
}
\details{
	Use is mainly internal
}
\value{
	A data-frame containing phenotypic data
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
as.data.frame(srdta[1:5,])
}
\keyword{misc}% at least one, from doc/KEYWORDS
