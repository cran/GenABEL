\name{summary.check.marker}
\alias{summary.check.marker}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Summary of check.marker object}
\description{
Provides cross-tabulation summarising number 
of marker which did not pass this or that criteria
}
\usage{
	\method{summary}{check.marker}(object, ...)
}
\arguments{
  \item{object}{object of class \code{\link{check.marker-class}}}
  \item{\dots}{additional arguments (not used)}
}
%\details{}
\value{
	A list containing 2 tables: per-marker and per-person 
	inconsistencies
}
%\references{ ~put references to the literature/web site here ~ }
\author{Yurii Aulchenko}
%\note{ ~~further notes~~ }
\seealso{
\code{\link{check.marker}},
\code{\link{check.marker-class}}
}
\examples{
require(GenABEL.data)
data(srdta)
mc <- check.marker(srdta,ids=c(1:500))
summary(mc)
}
\keyword{htest}% at least one, from doc/KEYWORDS
