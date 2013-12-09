\name{qvaluebh95}
\alias{qvaluebh95}
\title{Computes Benjamini-Hochberg (95) q-value}
\description{
Computes Benjamini-Hochberg (95) q-value
}
\usage{
qvaluebh95(p, fdrate=0.1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{p}{vector containing p-values}
  \item{fdrate}{desired FRD}
}
%\details{
%}
\value{
	A list containing
	\item{pass}{Is true if this P-value passed specified FDR}
	\item{qvalue}{qvalue}
}
%\references{}
\author{Yurii Aulchenko}
%\note{}
%\seealso{
%}
\examples{
require(GenABEL.data)
data(srdta)
a<-qtscore(qt2,data=srdta)
qv <- qvaluebh95(a[,"P1df"])
plot(a[,"Position"],-log10(qv$qvalue))
}
\keyword{htest}% at least one, from doc/KEYWORDS
