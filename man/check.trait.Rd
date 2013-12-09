\name{check.trait}
\alias{check.trait}
\title{function to do primitive trait quality control}
\description{
This function check for outliers (using FDR framework) 
and plots the raw data.
}
\usage{
check.trait(trait, data, fdrate = 0.05, graph = TRUE, binshow = FALSE, 
		qoption = "bh95")
}
\arguments{
  \item{trait}{name (or list of names) of trait(s) to be checked}
  \item{data}{gwaa.data object or data frame containing the trait}
  \item{fdrate}{false discovery rate to apply for QC}
  \item{graph}{if graphical output should be produced}
  \item{binshow}{if binary traits should be plotted}
  \item{qoption}{how to compute q-values (not implemented, currently using only BH95)}
}
\details{
The P-value that a particulat measurment is an outlier is compted as folowing.
Consider trait vector Y with particulat \eqn{i^{th}} measurment denodet as \eqn{y_i}. 
Let Y(-i) is vector, which is the same as Y, except that \eqn{i^{th}} measurment
is dropped. Then Chi-square for measurment i is computed as 

\deqn{Chi_{i} = (mean(Y(-i)) - y_i)^2/var(Y(-i))}

P-value is computed using 1 d.f., and the vector of P-values enters FDR 
computation procedure (BH95 by default).
}
\value{
	No value returned, output is made to the screen and graphical device.
}
%\references{}
\author{Yurii Aulchenko}
%\note{
%}
\seealso{
	\code{\link{check.marker}}
}
\examples{
require(GenABEL.data)
data(srdta)
check.trait("qt3",data=srdta)
n <- names(srdta@phdata)
check.trait(n,data=srdta)
}
\keyword{misc}% at least one, from doc/KEYWORDS
