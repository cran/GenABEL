\name{estlambda}
\alias{estlambda}
\title{Estimate the inflation factor for a distribution of P-values}
\description{
Estimate the inflation factor for a distribution of P-values or 1df chi-square test.
The major use of this procedure is the Genomic Control, but can also be used to 
visualise the distribution of P-values coming from other tests.
}
\usage{
estlambda(data, plot = TRUE, proportion = 1.0, ... )
}
\arguments{
  \item{data}{A vector of reals. If all are <=1, it is assumed that this 
		is a vector of P-values, else it is treated as a vector 
		of chi-squares with 1 d.f.}
  \item{plot}{Wether the prot should be presented}
  \item{proportion}{The proportion of lowest P (Chi2) to be used when 
		estimating the inflation factor Lambda}
 \item{...}{arguments passed to \code{\link{plot}} function}
}
%\details{
%}
\value{
	A list with elements
	\item{estimate}{Estimate of Lambda}
	\item{se}{Standard error of the estimate}
}
%\references{}
\author{Yurii Aulchenko}
%\note{
%}
\seealso{
\code{\link{ccfast}},
\code{\link{qtscore}}
}
\examples{
data(srdta)
pex <- summary(srdta@gtdata)[,"Pexact"]
estlambda(pex)
a <- ccfast("bt",srdta)
lambda(a)
}
\keyword{htest}% at least one, from doc/KEYWORDS
