\name{ztransform}
\alias{ztransform}
\title{Transformation to standard Normal}
\description{
Transformation of a variable or residuals from GLM analysis to standard Normal.
}
\usage{
ztransform(formula,data,family=gaussian)
}
\arguments{
  \item{formula}{GLM formula for the variable to be transformed, or just the variable}
  \item{data}{data.frame or gwaa.data object containing the data}
  \item{family}{GLM family}
}
\details{
Transformation to normality generates a variable which has mean zero and variance 
of one. If formula used, residuals from regression model are scaled to standard Normal.
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
	# tests with genetic data
	require(GenABEL.data)
	data(srdta)
	Zqt1 <- ztransform(qt1,srdta)
	Zqt1sexA <- ztransform(qt1~sex,srdta)
}
\keyword{utilities}
