\name{add.plot.scan.gwaa}
\alias{add.plot.scan.gwaa}
\title{function to plot additional GWAA results}
\description{
Add plot of results of GWA analysis
}
\usage{
add.plot.scan.gwaa(x, ..., df = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{object of type scan.gwaa-class, as returned by \code{\link{scan.glm}}, \code{\link{qtscore}},
		\code{\link{ccfast}}, \code{\link{emp.ccfast}}, \code{\link{emp.qtscore}}, or \code{\link{scan.haplo}}}
  \item{\dots}{additional arguments to be passed to plot}
  \item{df}{P-value at which df to add (1 or 2)}
}
%\details{}
\value{
	No value returned.
}
%\references{}
\author{Yurii Aulchenko}
%\note{}
\seealso{
	\code{\link{plot.scan.gwaa}}, 
	\code{\link{snp.subset}}, 
	\code{\link{scan.glm}}, 
	\code{\link{qtscore}}, 
	\code{\link{ccfast}},
	\code{\link{emp.qtscore}}, 
	\code{\link{emp.ccfast}},
	\code{\link{scan.haplo}}
}
\examples{
data(srdta)
a <- ccfast("bt",srdta,snps=c(1:100))
plot(a)
a1 <- qtscore("bt~CRSNP",srdta,snps=c(1:100))
add.plot.scan.gwaa(a1,col="red",type="l")
}
\keyword{hplot}% at least one, from doc/KEYWORDS
