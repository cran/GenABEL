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
\details{
}
\value{
	List
	\section{Names}{
	\item{pass}{Is true if this P-value passed specified FDR}
	\item{qvalue}{qvalue}
	} 
}
%\references{}
\author{Yurii Aulchenko}
%\note{}
\seealso{
}
\examples{
data(srdta)
a<-ccfast("bt",data=srdta)
qv <- qvaluebh95(a$P1df)
plot(a$map,-log10(qv$qvalue))
}
\keyword{htest}% at least one, from doc/KEYWORDS
