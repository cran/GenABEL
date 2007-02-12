\name{catable}
\alias{catable}
\title{function to generate summary table for quantitative data}
\description{
	This function makes a table with number of observations 
	which fall between user-defined categories
}
\usage{
catable(data, categories = c(quantile(data,c(0.01,0.1,0.5,0.9,0.99))), cumulative = F)
}
\arguments{
  \item{data}{A vector of numerics}
  \item{categories}{vector containing desired cut-off levels}
  \item{cumulative}{whether cumulative distribution should be shown}
}
\details{
}
\value{
	table with number and proportion of observations falling between 
	categories
}
\references{}
\author{Yurii Aulchenko}
\note{
}
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
	data(srdta)
	callr <- summary(srdta@gtdata)[,"CallRate"]
	catable(callr,c(0.93,0.95,0.99))
	catable(callr)
	catable(callr,cum=TRUE)
}
\keyword{distribution}
