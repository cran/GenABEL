\name{rntransform}
\alias{rntransform}
\title{Rank-transformation to normality}
\description{
Rank-transformation to normality of a variable or residuals from GLM analysis.
}
\usage{
rntransform(formula,data,family=gaussian)
}
\arguments{
  \item{formula}{GLM formula for the variable to be transformed, or just the variable}
  \item{data}{data.frame or gwaa.data object containing the data}
  \item{family}{GLM family}
}
\details{
Rank-transformation to normality generates perfectly normal distribution 
from ANY distribution, unless many/heavy ties are present in variable 
(or residuals, if formula is used).

When formula is supplied, this procedure first calls \code{\link{ztransform}}, 
and then applies rank transformation to residuals.
}
\value{
Vector containing transformed variable, distributed as standard normal.
}
%\references{}
\author{Yurii Aulchenko}
%\note{ ~~further notes~~ 
%
%}
\seealso{\code{\link{ztransform}}}
\examples{
	# uniformly distributed variable
	x <- round(runif(200)*100)
	# get 7 missing values
	x[round(runif(7,min=1,max=100))] <- NA
	# Z-transform
	y0 <- ztransform(x)
	# Rank-transform to normality
	y1 <- rntransform(x)
	# test normality of the original and transformed var
	shapiro.test(x)
	shapiro.test(y0)
	shapiro.test(y1)
	# plot histogram
	par(mfcol=c(3,1))
	hist(x)
	hist(y0)
	hist(y1)
}
\keyword{utilities}
