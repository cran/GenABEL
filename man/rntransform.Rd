\name{rntransform}
\alias{rntransform}
\title{Rank-transformation to normality}
\description{
Rank-transformation to normality.
}
\usage{
rntransform(var)
}
\arguments{
  \item{var}{Vector containing a variable to be transformed}
}
\details{
Rank-transformation to normality generates perfectly normal distribution 
from ANY distribution, unless many/heavy ties are present.
}
\value{
Vector containing transformed variable
}
%\references{}
\author{Yurii Aulchenko}
%\note{ ~~further notes~~ 
%
%}
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
	# uniformly distributed variable
	x <- round(runif(200)*100)
	# get 7 missing values
	x[round(runif(7,min=1,max=100))] <- NA
	# transform
	y <- rntransform(x)
	# test normality of the original and transformed var
	shapiro.test(x)
	shapiro.test(y)
	# plot histogram
	par(mfcol=c(2,1))
	hist(x)
	hist(y)
}
\keyword{utilities}
