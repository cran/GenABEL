\name{plot.scan.gwaa}
\alias{plot.scan.gwaa}
\title{function to plot GWAA results}
\description{
Plots results of GWA analysis
}
\usage{
plot.scan.gwaa(x, y, ..., df=1,ylim)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{object of type \code{\link{scan.gwaa-class}}, as returned by \code{\link{scan.glm}}, \code{\link{qtscore}},
		\code{\link{ccfast}}, \code{\link{emp.ccfast}}, \code{\link{emp.qtscore}}, or \code{\link{scan.haplo}}}
  \item{y}{this argument is not used}
  \item{\dots}{additional arguments to be passed to plot}
  \item{df}{Plot results of 1 or 2-df test (1, 2). Could be also "Pc1df" (for GC corrected P-values), "Pc2df" 
		(for robust genomic control of the 2 df test) or "all" (all results)}
  \item{ylim}{ylim, same as in the standard plot function}
}
%\details{}
\value{
	No value returned.
}
%\references{}
\author{Yurii Aulchenko}
%\note{}
\seealso{
	\code{\link{scan.gwaa-class}}, 
	\code{\link{add.plot}}, 
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
a <- qtscore(bt,srdta,snps=c(1:250))
plot(a)
add.plot(a,df="Pc1df",col="green")
}
\keyword{hplot}% at least one, from doc/KEYWORDS
